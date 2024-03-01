# -------------------  FUNÇÃO ----------------- #


## lendo arquivo 
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]
dados <- read.table(LIB, header = FALSE, sep = "")


## função calculo de somia dos cromossomos onde o valor será correspondente a
#profundiade mediana do cromossomo dividido pela profundidade mediana da cobertura
somia <- function(df, grupo, depth){
  library(dplyr)
  library(tidyverse)
  library(stringr)
  library(conflicted)
  conflict_prefer("filter", "dplyr")
  res <- df %>%
    group_by({{grupo}}) %>%
    dplyr::filter(str_detect({{grupo}}, "NC_")) %>%
    summarise(prof_mediana = median({{depth}})) %>%
    mutate(mediana_cob = median(prof_mediana)) %>% 
    mutate(somia = prof_mediana/mediana_cob) %>%
    mutate(somia2 = somia*2)
  
  
  return(res)
}

resultados <- as.data.frame(somia(dados, V1, V3))

#arredondar para uma casa após a virgula 
resultados$somia2 <- round(resultados$somia2, 1)


## substituindo o código pelo numero do cromossomo no genoma
#length(resultados$V1)
resultados$V1[resultados$V1 == c("NC_009277.2", "NC_009386.2", "NC_009387.2", "NC_009388.2", "NC_009389.2", "NC_009390.2", "NC_009391.2", "NC_009392.2", "NC_009393.2", "NC_009394.2", "NC_009395.2", "NC_009396.2", "NC_009397.2", "NC_009398.2", "NC_009399.2", "NC_009400.2", "NC_009401.2", "NC_009402.2", "NC_009403.2", "NC_009404.2", "NC_009405.2", "NC_009406.2", "NC_009407.2", "NC_009408.2", "NC_009409.2", "NC_009410.2", "NC_009411.2", "NC_009412.2", 
                                 "NC_009413.2", "NC_009414.2", "NC_009415.2", "NC_009416.2", "NC_009417.2", "NC_009418.2", "NC_009419.2", "NC_009420.2")] <- c("chr_6", "chr_1", "chr_2", "chr_3", "chr_4", "chr_5", "chr_7", "chr_8", "chr_9", "chr_10", "chr_11", "chr_12", "chr_13", "chr_14", "chr_15", 
                                                                                                                                                               "chr_16", "chr_17", "chr_18", "chr_19", "chr_20", "chr_21", "chr_22", "chr_23", "chr_24", "chr_25", "chr_26", "chr_27", "chr_28", "chr_29", "chr_30", "chr_31", "chr_32", "chr_33", 
                                                                                                                                                               "chr_34", "chr_35", "chr_36")
print(resultados)
#salvar
#write.table(resultados, file = "~/Área de Trabalho/leishmania/depth/resultados.txt", sep = "")


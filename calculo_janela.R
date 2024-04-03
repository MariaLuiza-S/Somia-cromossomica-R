######### cobertura janela de mil pb por cromossomo ############# 
#dados <- read.table("~/Desktop/leishmania/depth/SRR25744243.txt", header = FALSE, sep = "") 

# direcionado o terminal a leitura do arquivo localizado a frente do script
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]
dados <- read.table(LIB, header = FALSE, sep = "")

library(dplyr)
library(stringr)
##### cálculo da janela ############


# vetor com os códigos dos cromossomos para serem utilizados no loop
cromossomos = c("NC_009386.2", "NC_009387.2", "NC_009388.2", "NC_009389.2", "NC_009390.2", "NC_009277.2",
                "NC_009391.2", "NC_009392.2", "NC_009393.2", "NC_009394.2", "NC_009395.2", "NC_009396.2",
                "NC_009397.2", "NC_009398.2", "NC_009399.2", "NC_009400.2", "NC_009401.2", "NC_009402.2",
                "NC_009403.2", "NC_009404.2", "NC_009405.2", "NC_009406.2", "NC_009407.2", "NC_009408.2",
                "NC_009409.2", "NC_009410.2", "NC_009411.2", "NC_009412.2", "NC_009413.2", "NC_009414.2",
                "NC_009415.2", "NC_009416.2", "NC_009417.2", "NC_009418.2", "NC_009419.2", "NC_009420.2"
                )

## data frame vazio que será utilizado para armazenar o resultado final
chr_all <- NULL

## incio do loop - itera sobre cada elemento do vetor cromossomo as medidas estatísticas
for (chr in 1:length(cromossomos)){
  
  ## data frame temporário para armazenamento dos cromossomos
  chr_tmp <- (dados) %>% 
    group_by(V1) %>%
    filter(str_detect(V1, cromossomos[chr])) %>%
    mutate(cromossomo = cromossomos[chr]) %>%
    as.data.frame()
  
  ## data frame temporário para armazentamento das análises estatísticas (Média, mediana e desvio padrao) a cada mil pb por cromossomos
  # calculo do numero total de grupos (número total de pares de bases / 1000) e criacao de um vetor de 1 ate o numero total de grupos para cada cromossomo
  # cada numero do vetor sera repetido 1000 vezes e os primeiros elementos de cada vetor distribuidos a seus respectivos grupos
  Jchr_tmp <- chr_tmp %>%
    mutate(intervalo = rep(1:ceiling(n() / 1000), each = 1000)[1:n()]) %>%
    group_by(intervalo, cromossomo) %>%
    summarise(
      media = mean(V3),
      mediana = median(V3),
      desvio_padrao = sd(V3),
      .groups = "drop"
    ) %>%
  as.data.frame()
  
  
  ## armazenamento dos resultados 
  chr_all <- rbind(chr_all, Jchr_tmp)  
  
}


## exlusao de datas frames que não serao mais utilizados
rm(chr_tmp, Jchr_tmp, dados)

### substituindo o código dos cromossomos por números
chr_all$cromossomo <- recode(chr_all$cromossomo, "NC_009386.2" = "chr_1", "NC_009387.2" = "chr_2", "NC_009388.2" = "chr_3", "NC_009389.2" = "chr_4", "NC_009390.2" = "chr_5", "NC_009277.2" = "chr_6",
                             "NC_009391.2" = "chr_7", "NC_009392.2" = "chr_8", "NC_009393.2" = "chr_9", "NC_009394.2" = "chr_10", "NC_009395.2" = "chr_11", "NC_009396.2" = "chr_12",
                             "NC_009397.2" = "chr_13", "NC_009398.2" = "chr_14", "NC_009399.2" = "chr_15", "NC_009400.2" = "chr_16", "NC_009401.2" = "chr_17", "NC_009402.2" = "chr_18",
                             "NC_009403.2" = "chr_19", "NC_009404.2" = "chr_20", "NC_009405.2" = "chr_21", "NC_009406.2" = "chr_22", "NC_009407.2" = "chr_23", "NC_009408.2" = "chr_24",
                             "NC_009409.2" = "chr_25", "NC_009410.2" = "chr_26", "NC_009411.2" = "chr_27", "NC_009412.2" = "chr_28", "NC_009413.2" = "chr_29", "NC_009414.2" = "chr_30",
                             "NC_009415.2" = "chr_31", "NC_009416.2" = "chr_32", "NC_009417.2" = "chr_33", "NC_009418.2" = "chr_34", "NC_009419.2" = "chr_35", "NC_009420.2" = "chr_36")

#write.table(chr_all, file = "~/Área de Trabalho/leishmania/depth/LIB.txt", sep = "")
options(max.print = 1000000)
print(chr_all)


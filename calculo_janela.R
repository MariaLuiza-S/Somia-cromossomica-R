######### cobertura janela de mil pb por cromossomo ############# 
dados <- read.table("~/Desktop/leishmania/depth/SRR25744243.txt", header = FALSE, sep = "") 


library(dplyr)
library(stringr)


## filtrar cada cromossomo de forma individual 
chr_1 <- dados %>% 
  group_by(V1) %>%
  filter(str_detect(V1, "NC_009386.2"))


## medidas estatísticas (media, mediana e desvio padrao) para cada mil pb do cromossomo
# calculo do numero total de grupos e criacao de um vetor de 1 ate o numero total de grupos;
# cada numero do vetor sera repetido 1000 vezes e os primeiros elementos de cada vetor distribuidos a seus respectivos grupos

resultados <- chr_1 %>%
  mutate(chr1 = rep(1:ceiling(n() / 1000), each = 1000)[1:n()]) %>%
  group_by(chr1) %>%
  summarise(
    media = mean(V3),
    mediana = median(V3),
    desvio_padrao = sd(V3)
  )


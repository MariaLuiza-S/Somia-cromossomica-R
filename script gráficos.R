######################## carregando os pacotes ###############
library(pheatmap)
library(tidyverse)
library(stringr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gridExtra)


# leituras dos arquivos da pasta com final txt
lista  <- list.files(pattern="txt")

# reunindo os arquivos em uma lista
arquivos <- lapply(lista, function(x) read.table(x, header = TRUE, sep = "", dec = ".")) 

# criando um data frame com a uni?o dos arquivos da lista
dados    <- do.call("cbind", arquivos) 



## renomeando coluna V1 para n?mero dos cromossomos
rename_duplicates <- function(dados) {
  col_names <- colnames(dados)
  somia_indices <- which(col_names == "V1")
  
  for (i in seq_along(somia_indices)) {
    col_names[somia_indices[i]] <- paste0("cromossomos_", i)
  }
  
  colnames(dados) <- col_names
  return(dados)
}

# novo dataframe com a coluna cromossomos renomeada
DADOS <- rename_duplicates(dados)


### renomear amostras
## trecho a ser deletado do nomes das bibliotecas
delet <- "_somia.txt"

# renomear colunas duplicadas com base nos nomes dos arquivos
rename_duplicates <- function(DADOS, lista, delet) {
  col_names <- colnames(DADOS)
  somia_indices <- which(col_names == "somia2")
  
  lista_modificada <- gsub(delet, "", lista)
  
  for (i in seq_along(somia_indices)) {
    col_names[somia_indices[i]] <- lista_modificada[i]
  }
  
  
  colnames(DADOS) <- col_names
  return(DADOS)
}


# aplicando a mudan?a de nomes no data frame
somia <- rename_duplicates(DADOS, lista, delet)

# selecionar apenas a coluna do cromossomo e de somia - excluindo as colunas de profundidades
ploidia <- select(somia, cromossomos_1, starts_with("ERR"), starts_with("SRR"))


############# spearman ##############
#################### transpostar a matriz ##########################
tploidia <- as.data.frame(t(ploidia))
teste <- as.data.frame(tploidia[-1,])
str(teste)

# mudar para num?rico
teste <- teste %>% mutate_if(is.character, as.numeric)
library(GGally)

# metodo de spearman
correlation_matrix <- cor(teste, method = "spearman")

melted_corr <- melt(correlation_matrix)


melted_corr_lower <- melted_corr[melted_corr$Var1 > melted_corr$Var2, ]

melted_corr_diagonal <- melted_corr[melted_corr$Var1 == melted_corr$Var2, ]


ggplot(data = melted_corr_lower, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name="") +
  geom_text(aes(label = round(value, 1)), color = "black", size = 2) +
  
  # Adicionar rótulos ao longo da diagonal
  geom_text(data = melted_corr_diagonal, aes(x = Var1, y = Var2, label = Var1), 
            color = "black", size = 3, angle = 45, hjust = 1, vjust = -0.5) +
  
  # Manter os rótulos do eixo x e personalizar a aparência
  scale_x_discrete(limits = unique(melted_corr$Var1)) +
  scale_y_discrete(limits = unique(melted_corr$Var2)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),  # Remover rótulos padrão do eixo y
    axis.text.x = element_text(angle = 0, hjust = 1, size = 10),  # Manter rótulos no eixo x
    axis.title.x = element_blank(), 
    axis.ticks.y = element_blank(),  # Remover ticks do eixo y
    panel.grid = element_blank()
  ) +
  coord_fixed()

x11()


# tornado a coluna de cromossomos n?merica
ploidia$cromossomos_1[ploidia$cromossomos_1 == c("chr_6", "chr_1", "chr_2", "chr_3", "chr_4", "chr_5", "chr_7", "chr_8", "chr_9", "chr_10", 
                                                 "chr_11", "chr_12", "chr_13", "chr_14", "chr_15", "chr_16", "chr_17", "chr_18", "chr_19", "chr_20", 
                                                 "chr_21", "chr_22", "chr_23", "chr_24", "chr_25", "chr_26", "chr_27", "chr_28", "chr_29", "chr_30", 
                                                 "chr_31", "chr_32", "chr_33", "chr_34", "chr_35", "chr_36")] <- c(6, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                                                                                                   21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36)

str(ploidia)
ploidia$cromossomos_1 <- as.numeric(ploidia$cromossomos_1)

### exluir a coluna de cromossomos pra plotar heatmap
ploidiaR <- ploidia[,-1]



################## HEATMAP ##################
pheatmap(ploidiaR, cluster_rows = FALSE, show_colnames = F, 
         main = "Estimativa de ploidia")


# mudar a formata??o dos dados para construir boxplot - adicionar coluna 'variable'
box <- reshape2::melt(ploidia, id.var = 'cromossomos_1')

# transformar os cromossomos em factor e atribuir a ordem
box$cromossomos_1 <- factor(box$cromossomos_1, levels = c("1", "2", "3", "4", "5", "6", 
                                                          "7", "8", "9", "10", "11", "12", 
                                                          "13", "14", "15", "16", "17", "18", 
                                                          "19", "20", "21", "22", "23", "24", 
                                                          "25", "26", "27", "28", "29", "30", 
                                                          "31", "32", "33", "34", "35", "36"))





### arquivo fenotipos
grupos <- read.csv("D:/PESQUISA/Grupos.csv", header = T, sep = ",", dec = ".", fileEncoding = "UTF-8")

colnames(grupos)[1] <- "variable"
colnames(grupos)[2] <- "Grupos"
colnames(grupos)[3] <- "Fármacos"

dados_box <- inner_join(box, grupos)




# ajustando associações ao data frame
amostras <- data.frame(grupos[,1])
grupos <- data.frame(grupos[,-1])

### associando a amotra ao fenotipo
row.names(grupos) <- make.unique(amostras[, 1])
colnames(grupos)[1] <- "Grupos"


# cor
gp_colors <- list(Grupos = c("Controle" = "deepskyblue",
                             "Resistência" = "coral1"), Fármacos = c("Terbafina" = "sandybrown", 
                                                                  "S/F" = "beige",
                                                                  "Anfoterecina B" = "gold1",
                                                                  "Miltefosina" = "yellowgreen",
                                                                  "Tiofeno" = "lightpink", 
                                                                  "Antimônio" = "violet", 
                                                                  "Paromicina" = "sienna1", 
                                                                  "Alopurinol" = "orange",
                                                                  "Fluorouracil" = "skyblue"))

heat_plot <- pheatmap(ploidiaR,
                      cluster_rows = F, cluster_cols = T, 
                     # clustering_distance_cols = 'euclidean',
                    #  clustering_distance_rows = 'euclidean',
                      annotation_colors = gp_colors,
                     clustering_method = 'ward.D', 
                      annotation_col = grupos,
                      annotation_names_row = T, 
                      annotation_names_col = F,
                      fontsize_row = 10,         
                      fontsize_col = 7,         
                      angle_col = 45, 
                      show_colnames = F, show_rownames = T,
                      main = "Estimativa do número de cópias cromossômicas") 
x11()




library(pheatmap)


str(grupos)
unique(grupos$Fármacos)

 ############## BOXPLOT ###############
boxplot <- dados_box[,-5]

x11()

ggplot(boxplot, aes(x=as.factor(cromossomos_1), y=value, fill=Grupos)) +
  geom_boxplot(position=position_dodge(0.8), colour="black", outlier.colour = "black",
               outlier.shape = 1, outlier.size = 1.5) +
  ylab("Ploidia") +
  ylim(0, 8) +
  xlab("Cromossomos") +
  labs(title = "Estimativa de Ploidia - Controle e Resistência") +
  theme_bw() +
  theme_classic() +
  scale_fill_manual(values=c("blue", "red"), 
                    name="Grupo", 
                    labels=c("Controle", "Resistência"))

############### VARIA??O (SD) ###############

# verificando mediana e desvio padr?o por cromossomos e amostra
# filtrando amostras em que a ploidia ? 2x mais ou menos que o desvio padr?o
sd <- boxplot %>%
  group_by(cromossomos_1) %>%
  mutate(
    desvio_padrao = sd(value), 
    mediana = median(value), 
    media = mean(value),
    IIQ = IQR(value),
    Q1 = quantile(value, 0.25), 
    Q3 = quantile(value, 0.75)) %>%
  mutate(OutM = ifelse(value < (mediana - 2* desvio_padrao) | value  > (mediana + 2* desvio_padrao), 1, 0), 
         outliers = ifelse(value < Q1 - 1.5*IIQ | value > Q3 + 1.5*IIQ, 1, 0))



# filtrando os outliers
outliers <- sd %>% filter(outliers == 1)


amostras_outliers <- outliers %>%
  group_by(variable) %>%
  summarise(media_ploidia = mean(value, na.rm = TRUE))










########################## FILTRAR APENAS PELO CORTE
#lista <- list.files(pattern = "txt")

#arquivos <- lapply(lista, function(x) read.table(x, header = TRUE, sep = "", dec = "."))

#arquivos_filtrados <- Filter(function(df) all(df$media >= 20), arquivos)

#arquivos_filtrados <- lista[sapply(lista, function(x) {
 # df <- read.table(x, header = TRUE, sep = "", dec = ".")
  #any(df$media < 10)
})]

# Exibir os arquivos filtrados
#print(arquivos_filtrados)


# Combinar os data frames filtrados em um ?nico data frame
#dados <- do.call("cbind", arquivos_filtrados)
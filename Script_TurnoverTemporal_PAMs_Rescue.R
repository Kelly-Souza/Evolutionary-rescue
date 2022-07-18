setwd("C:/Users/lets/Desktop/RESCUE/Test_final_2")
rm(list=ls())

list.files(pattern = ".txt")
library(betapart)

P_1 <- read.table("PAM_final_RescueBio8eBio16_Trailing_Stable_Leading.txt", h = T)
P_1 <- P_1[,-c(1:2)]

for(i in 1:ncol(P_1)){
 P_1[,i] <- ifelse(test = P_1[,i]>= 1, yes = 1, no = 0)
 print(i)
}

P_0 <- read.table("PAM_filtered_amphibian.txt", h = T)
P_0 <- P_0[,-c(1:2)]
dim(P_0)

list.files()

colnames(P_0)
colnames(P_1) <- colnames(P_0)

head(P_0[1:6,1:6])
head(P_1[1:6,1:6])

?beta.temp

Turnover <- numeric(length = nrow(P_0))
Nestedness <- numeric(length = nrow(P_0))
Dissimilarity <- numeric(length = nrow(P_0))

for(i in 1:nrow(P_0)){
  Temp <- beta.temp(x = P_0[i,], y = P_1[i,], index.family="sorensen") 
  #Turnover[i]<- beta.temp(x = P_0[i,-c(1,2)], y = P_1[i,], index.family="sorensen")$beta.sim 
  Turnover[i]<- Temp$beta.sim 
  Nestedness[i]<- Temp$beta.sne 
  Dissimilarity[i]<- Temp$beta.sor 
  print(i)
}

  
is.na(Dissimilarity) == is.na(Nestedness)
is.na(Dissimilarity) == is.na(Turnover)
is.na(Nestedness) == is.na(Turnover)
# nao sao iguais as posicoes com NA
# XYZTurnover <- cbind(P_0[,1:2],ifelse(is.nan(Turnover), yes = NA, no = Turnover))
# write.table(XYZTurnover, "XYZTurnover.txt")


PAM <- read.table(file = "PAM_filtered.txt", header = T)
coords <- PAM[,1:2]
head(coords)
XYZTurnover_Nestedness_Dissimilarity <- cbind(coords,Turnover,  Nestedness, Dissimilarity)
write.table(x = XYZTurnover_Nestedness_Dissimilarity, file = "XYZTurnover_Nestedness_DissimilarityBio8Bio16.txt")
head(XYZTurnover_Nestedness_Dissimilarity)

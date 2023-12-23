################################################################################
rm(list = ls())
gc()
################################################################################
# Process image features
LUAD.image <- read.csv("LUAD_type1_features.csv")
row.names(LUAD.image) <- LUAD.image[,1]
# Process gene features
LUAD.gene <- read.table("LUAD_type2_features.txt")
LUAD.gene <- t(LUAD.gene)
LUAD.gene <- as.data.frame(LUAD.gene)
names(LUAD.gene) <- LUAD.gene[1,]
LUAD.gene <- LUAD.gene[-c(1,2),]
na.split1 <- strsplit(LUAD.gene[,1], split = '-')
for (j in 1:length(na.split1)) {
  LUAD.gene[j,1] <- paste(na.split1[[j]][1],na.split1[[j]][2],na.split1[[j]][3], sep = ".")
}
LUAD.gene <- LUAD.gene[!duplicated(LUAD.gene[,1]),]
row.names(LUAD.gene) <- LUAD.gene[,1]

ID <- intersect(LUAD.image[,1],LUAD.gene[,1])
data.x <- LUAD.image[match(ID,LUAD.image[,1]),-1]
data.z <- LUAD.gene[match(ID,LUAD.gene[,1]),-1]
rm(LUAD.image,LUAD.gene,na.split1,j)
gc()
################################################################################
# Data cleaning on gene features
n <- length(ID)
# remove missing genes
pp <- dim(data.z)[2]
NA.num <- rep(0,pp)
for (j in 1:pp) {
  NA.num[j] <- sum(is.na(data.z[,j]))
}
data.z <- data.z[,which(NA.num == 0)]
# remove genes who have high repetition rate 
pp <- dim(data.z)[2]
dup.num <- rep(0,pp)
for (j in 1:pp) {
  dup.num[j] <- sum(duplicated(data.z[,j]))
}
data.z <- data.z[,which(dup.num < 0.1*n)]
# replace outliers by max and min values. 
pp <- dim(data.z)[2]
var.num <- rep(0,pp)
for (j in 1:pp) {
  a <- as.numeric(data.z[,j])
  a.out <- which(a < mean(a) - 6*sd(a) | a > mean(a) + 6*sd(a) | a > 10)
  if(length(a.out) > 0){
    a[a.out] <- as.numeric(a[a.out] > 0) * max(a[-a.out]) + as.numeric(a[a.out] < 0) * min(a[-a.out])
  }
  data.z[,j] <- a
}
# if(!requireNamespace("BiocManager", quietly = TRUE))      
#   install.packages("BiocManager") 
# BiocManager::install("KEGGREST") 
library(KEGGREST)
library(psych)
# Select genes on Lung cancer pathway
gs <- keggGet('hsa05223')
genes <- unlist(lapply(gs[[1]]$GENE, function(x) strsplit(x,';'))) 
genelist <- data.frame(genes[1:length(genes)%%3 ==2])[,1]
n.ID <- match(genelist,names(data.z))
n.ID <- n.ID[!is.na(n.ID)]
data.z <- data.z[,n.ID]

rm(list=setdiff(ls(),c("data.x","data.z","ID")))
gc()

data.patient <- read.csv("LUAD_clinical.csv")
na.split2 <- strsplit(data.patient[,1], split = '-')
for (j in 1:length(na.split2)) {
  data.patient[j,1] <- paste(na.split2[[j]][1],na.split2[[j]][2],na.split2[[j]][3], sep = ".")
}
row.names(data.patient) <- data.patient[,1]

id <- intersect(ID, rownames(data.patient)[which(data.patient$FEV1_PERCENT_REF_PREBRONCHOLIATOR != "[Not Available]")])
data.x <- data.x[id,]
data.z <- data.z[id,]
data.y <- data.frame(FEV1 = as.numeric(data.patient[id,]$FEV1_PERCENT_REF_PREBRONCHOLIATOR))
data.x.scale <- as.data.frame(scale(data.x))
data.z.scale <- as.data.frame(scale(data.z))
data.y.scale <- as.data.frame(scale(data.y))
id <- row.names(data.x)
load('freq.RData')
data <- list(data.x = data.x, data.z = data.z, data.y = data.y, data.x.scale = data.x.scale, data.z.scale = data.z.scale,
             data.y.scale = data.y.scale, freq = freq, id = id)
save(data, file = "LUADdata116.RData")






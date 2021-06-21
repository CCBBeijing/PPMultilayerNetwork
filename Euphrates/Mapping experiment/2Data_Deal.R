marker <- read.csv("HuyangMap-Genotype.csv",check.names = F)
marker_name <- read.csv("HuyangMap-Genotype.csv",check.names = F)
marker <- as.matrix(marker[,-1:-3])
marker_name <- as.matrix(marker_name[,1:3])
marker[which(marker=="ll")] <- 0
marker[which(marker=="lm")] <- 1
marker[which(marker=="nn")] <- 0
marker[which(marker=="np")] <- 1
marker[which(marker=="hh")] <- 0
marker[which(marker=="hk")] <- 1
marker[which(marker=="kk")] <- 2
marker[which(marker=="--")] <- 9

CK <- as.matrix(read.csv("0CK_data.csv",check.names = F))
Salt <- as.matrix(read.csv("0Salt_data.csv",check.names = F))
ck_PL_root <- CK[which(CK=="primaryroot",arr.ind = T)[,1],c(1,17:31)]  
salt_PL_root <- Salt[which(Salt=="primaryroot",arr.ind = T)[,1],c(1,17:31)]   
ck_PL_root <-apply(ck_PL_root , 2, as.numeric) 
salt_PL_root <-apply(salt_PL_root , 2, as.numeric)

ck_PLr_root <- c()
for (n in unique(ck_PL_root[,1])) {
  
  ck_PLr_root <- rbind(ck_PLr_root, apply(matrix(ck_PL_root[which(ck_PL_root[,1]==n),],ncol = dim(ck_PL_root)[2]), 2, mean,na.rm=T) ) 
} 

colnames(ck_PLr_root) <- c("line",(0:(dim(ck_PLr_root)[2]-2)))   

salt_PLr_root <- c()
for (n in unique(salt_PL_root[,1])) {
  
  salt_PLr_root <- rbind(salt_PLr_root, apply(matrix(salt_PL_root[which(salt_PL_root[,1]==n),],ncol = dim(salt_PL_root)[2]), 2, mean,na.rm=T) ) 
} 
colnames(salt_PLr_root) <- c("line",(0:(dim(salt_PLr_root)[2]-2)))  

unique(c(ck_PLr_root[,1],salt_PLr_root[,1]))
T_marker <- marker[,colnames(marker)%in%unique(c(ck_PLr_root[,1],salt_PLr_root[,1]))]


ck_PLr_root  <- subset(ck_PLr_root,ck_PLr_root[,1]%in%salt_PLr_root[,1])
salt_PLr_root <- subset(salt_PLr_root,salt_PLr_root[,1]%in%ck_PLr_root[,1])
T_marker <- T_marker[,colnames(T_marker)%in%c(ck_PLr_root[,1])]

ck_PLr_root <- ck_PLr_root[order(ck_PLr_root[,1]),]
salt_PLr_root <- salt_PLr_root[order(salt_PLr_root[,1]),]

diff_PLr_root <- cbind(ck_PLr_root[,1],ck_PLr_root[,-1] - salt_PLr_root[,-1])

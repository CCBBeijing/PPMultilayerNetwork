diameter_ck <- read.csv("diameter_ck.csv",header=T,check.names=F)
diameter_salt <- read.csv("diameter_salt.csv",header=T,check.names=F)
hysnp <- read.csv("hysnp.csv",header=T,check.names=F)#snp data


##################################data processing ###########################################################################
an.load <- function(diameter_salt=diameter_salt ,
                    diameter_ck=diameter_ck,
                    hysnp=hysnp){
  r.snp <- colnames(hysnp)
  salt.s <- c()
  for(i in 1:114){
    salt.s <- c(salt.s,which(diameter_ck[i,1]==diameter_salt[,1]))
  }
  
  disalt <- diameter_salt[salt.s,]
  
  
  ck.s <- c()
  for(i in 1:82){
    ck.s <- c(ck.s,which(diameter_salt[i,1]==diameter_ck[,1]))
  }
  
  dick <- diameter_ck[ck.s,]
  
  
  
  snp.s2 <- c()
  for(i in 1:81){
    snp.s2 <- c(snp.s2,which(disalt[i,1]==r.snp))
  }
  
  hysnp <-hysnp[,snp.s2]
  
  rr <- colnames(hysnp)
  snp.s3 <- c()
  for(i in 1:73){
    snp.s3 <- c(snp.s3,which(rr[i]==dick[,1]))
  }
  
  
  diam_ck <- dick[snp.s3,]
  
  snp.s4 <- c()
  for(i in 1:73){
    snp.s4 <- c(snp.s4,which(rr[i]==disalt[,1]))
  }
  diam_salt <- disalt[snp.s4,]
  t <- seq(20,120,20)
  res <- list(diam_salt=diam_salt,diam_ck=diam_ck,hysnp=hysnp,t=t)
  res
  
}

dat <- an.load(diameter_salt=diameter_salt ,
               diameter_ck=diameter_ck,
               hysnp=hysnp)

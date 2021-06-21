library(splines)
library(orthogonalsplinebasis)
library(grplasso)
library(parallel)

cluster_Lpar <- c() #对每个类进行平均拟合
for (n in 1:length(all_cluster_value)) {
  cat(n,"\n") 
  cluster_Lpar <-  rbind(cluster_Lpar,smooth.optim_cluster(times=seq(13,78,length.out = 30),para=rep(0.1,6),effect_value=all_cluster_value[[n]]))
  cat("\n")
}

cluster_meanValue_450 <- c() #将每个聚类的均值放在一起,450个时间点
for (n in 1:dim(cluster_Lpar)[1]) {
  
  cluster_meanValue_450 <- cbind(cluster_meanValue_450,Legendre.model(seq(13,78,length.out = 450),cluster_Lpar[n,]))
  
}
colnames(cluster_meanValue_450) <- c(1:403)

cluster_meanValue_D_450 <- c() #将每个聚类的均值的勒让德导数放在一起,450个时间点
for (n in 1:dim(cluster_Lpar)[1]) {
  
  cluster_meanValue_D_450 <- cbind(cluster_meanValue_D_450,dLegendre.model(seq(13,78,length.out = 450),cluster_Lpar[n,]))
  
}
colnames(cluster_meanValue_D_450) <- c(1:403)

cluster_varsel <- varsel1(X=cluster_meanValue_450,Y=cluster_meanValue_D_450,tt=seq(13,78,length.out = 450))
colnames(cluster_varsel$connect) <- c(1:403)


cluster_list <- list()##整合变量选择结果
for (n in 1:403) {
  
  ttemp <- list()
  ttemp[[1]] <- n
  ttemp[[2]] <- names(which(cluster_varsel$connect[n,]==1))
  cluster_list[[n]] <- ttemp
}

cluster_mean <- c() #对每个类进行平均
for (n in 1:length(all_cluster_value)) {
  cat(n,"\n") 
  cluster_mean <-  rbind(cluster_mean,apply(all_cluster_value[[n]], 2,mean))
  cat("\n")
}

cluster_optim <- optim.parallel(connect=cluster_varsel$connect,effect=t(cluster_mean),
                                n.cores=1,proc=ode.optim,order=6,times=seq(13,78,length=30),nstep=29)

inter_effect <- regasso(connect1=cluster_varsel$connect,gene=cluster_mean,interaction=cluster_optim)
out_ModuleInform(cluster_list,cluster_optim,inter_effect,ClusterNum)



library(MASS)
#功能作图
PL_FunMap <- get_VG_FunMap(ck_PLr_root,salt_PLr_root,diff_PLr_root,T_marker,c(13:78))
PL_VG_30 <- get_VG(PL_FunMap[[2]],T_marker,seq(13,78,length.out = 30))
rownames(PL_VG_30) <- c(1:8305)

manhat_LR <- read.csv("snp_inform.csv",header = T)
manhat_LR[,4] <- PL_FunMap[[1]]

##功能聚类
FunC_order5 <- GeneClusterBIC(mExpression=(PL_VG_30), times=seq(13,78,length.out = 30),G=c(1:500),orderLOP=5)
ClusterNum <- length(table(FunC_order5$cluster[1,]))
Toplayer_cluster <- FunC_order5$cluster[1,]
all_cluster_value <- get_cluster(FunC_order5$cluster[1,],ClusterNum,PL_VG_30)



Logistic_diff = function(A,t){get_miu3(A[1:3],t)-get_miu3(A[4:6],t)}
rrloss <- function(yy,par,t){
  sum((yy-Logistic_diff(par,t))^2)
}
get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno[,-1],2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))
}
get_par <- function(hydata,marker_data,t){
  control <- dat$diam_ck
  stress <- dat$diam_salt
  diff_vg <- matrix(NA,ncol = 18,nrow =dim(marker_data)[1] ) 
  for (a in 1:dim(marker_data)[1]) {
    
    AA <-  as.numeric(which(marker_data[a,]==1))
    aa <- as.numeric(which(marker_data[a,]==0))
    Aa <- as.numeric(which(marker_data[a,]==2))
    
    par1 <-     colMeans(hydata[AA,])
    par0 <-     colMeans(hydata[aa,])
    par2 <-     colMeans(hydata[Aa,])
    if(length(Aa)==0){ 
      op1 <- optim(par = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t)),rrloss,t=t,yy=par1,method = "L-BFGS-B",lower = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t))*0.1, upper = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t))*3)
      op2 <- optim(par = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t)),rrloss,t=t,yy=par0,method = "L-BFGS-B",lower = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t))*0.1, upper = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t))*5)
      Vg <- c(op1$par,op2$par) 
    } else{
      op1 <- optim(par = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t)),rrloss,t=t,yy=par1,method = "L-BFGS-B",lower = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t))*0.1, upper = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t))*3)
      op2 <- optim(par = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t)),rrloss,t=t,yy=par0,method = "L-BFGS-B",lower = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t))*0.1, upper = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t))*5)
      op3 <- optim(par = c(get_initial_par(control[Aa,],t),get_initial_par(stress[Aa,],t)),rrloss,t=t,yy=par2,method = "L-BFGS-B",lower = c(get_initial_par(control[Aa,],t),get_initial_par(stress[Aa,],t))*0.1, upper = c(get_initial_par(control[Aa,],t),get_initial_par(stress[Aa,],t))*5)
      
      
      
      Vg <- c(op1$par,op2$par,op3$par) 
      
    }
    
    
    diff_vg[a,1:length(Vg)] <- Vg
    cat(a,"finished","\n")
    
  }
  
  
  return(diff_vg) 
}

get_par <- function(hydata,marker_data,t){
  control <- dat$diam_ck
  stress <- dat$diam_salt
  diff_vg <- matrix(NA,ncol = 18,nrow =dim(marker_data)[1] ) 
  for (a in 1:dim(marker_data)[1]) {
    AA <-  as.numeric(which(marker_data[a,]==1))
    aa <- as.numeric(which(marker_data[a,]==0))
    Aa <- as.numeric(which(marker_data[a,]==2))
    
    par1 <-     colMeans(hydata[AA,])
    par0 <-     colMeans(hydata[aa,])
    par2 <-     colMeans(hydata[Aa,])
    if(length(Aa)==0){ 
      op1 <- optim(par = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t)),rrloss,t=t,yy=par1,method = "BFGS")
      op2 <- optim(par = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t)),rrloss,t=t,yy=par0,method = "BFGS")
      Vg <- c(op1$par,op2$par) 
    } else{
      op1 <- optim(par = c(get_initial_par(control[AA,],t),get_initial_par(stress[AA,],t)),rrloss,t=t,yy=par1,method = "BFGS")
      op2 <- optim(par = c(get_initial_par(control[aa,],t),get_initial_par(stress[aa,],t)),rrloss,t=t,yy=par0,method = "BFGS")
      op3 <- optim(par = c(get_initial_par(control[Aa,],t),get_initial_par(stress[Aa,],t)),rrloss,t=t,yy=par2,method = "BFGS")
      
      
      
      Vg <- c(op1$par,op2$par,op3$par) 
      
    }
    
    
    diff_vg[a,1:length(Vg)] <- Vg
    cat(a,"finished","\n")
    
  }
  
  
  return(diff_vg) 
}

hypar <- get_par(hydata=(dat$diam_ck[,-1]-dat$diam_salt[,-1]),marker_data=dat$hysnp,t=dat$t)

get_VG <- function(FunMap_par,marker_data,t){
  
  T_marker <-marker_data#marker_data[,colnames(marker_data)%in%c(pheno_diff[,1])]
  
  diff_vg <- c() 
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(which(T_marker[a,]==1))
    aa <- as.numeric(which(T_marker[a,]==0))
    Aa <- as.numeric(which(T_marker[a,]==2))
    all <- as.numeric(which(T_marker[a,]!=9))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    p1 <- (NAA*2+NAa)/((NAA+NAa+Naa)*2) #A鍩哄洜棰戠巼
    p0 <- (Naa*2+NAa)/((NAA+NAa+Naa)*2) #a鍩哄洜棰戠巼
    
    mean_AA <- Logistic_diff(hypar[a,1:6],t)
    mean_aa <- Logistic_diff(hypar[a,7:12],t)
    
    #mean_AA <- Logistic_diff(FunMap_par[,a][8:13],t)
    #mean_aa <- Logistic_diff(FunMap_par[,a][2:7],t)
    AE <- (mean_AA - mean_aa)
    
    if(NAa==0){ Vg <- 2*p1*p0*(AE^2)  } else{
      mean_Aa <- Logistic_diff(hypar[a,13:18],t)
      #mean_Aa <- Logistic_diff(FunMap_par[,a][14:19],t)
      AE <- (mean_AA - mean_aa)/2
      DE <- mean_Aa - (mean_AA + mean_aa)/2
      Vg <- 2*p1*p0*((AE + (p1 - p0)*DE)^2) + 4*p1*p1*p0*p0*DE*DE
      
    }
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  colnames(diff_vg) <- seq(min(t),max(t),length.out = length(t))
  rownames(diff_vg) <- c(1:dim(marker_data)[1])
  return(sqrt(diff_vg)) 
}



hy_genetic_effect <- get_VG(hypar,dat$hysnp,t=seq(20,120,4))
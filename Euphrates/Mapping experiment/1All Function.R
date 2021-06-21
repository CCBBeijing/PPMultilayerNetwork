###FUNMAP
Logistic <- function(par,t){
  
  y <- (par[1]/(1+par[2]*exp(-par[3]*t))-(par[4]*exp(-par[5]*t)))
  return(y)
  
} 
SAD1 <- function(par, times = t) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi<- par[1]   
  v2 <- par[2]
  
  tmp <- (1-phi^2)   
  sigma <- array(1, dim=c(n,n))   
  
  for(i in 1:n)   
  {     
    sigma[i,i:n] <- phi^( c(i:n) - i ) * (1-phi^(2*i))/tmp     
    sigma[i:n,i] <- sigma[i,i:n]   
  }   
  sigma <- sigma * abs(v2)  
  
  return(sigma)
}
Likelihood = function(pheno,t,par){
  miu = Logistic(par[1:5],t)
  sigma = SAD1(par = c(par[6],par[7]),times = t)
  L0 = c()
  L0 = sum(dmvnorm(pheno,miu,sigma,log = T))
  return(-L0)
}
optim_deal <- function(pheno_deal,t){
  library(mvtnorm)
  itime <- 100
  itimes <- 1
  #par0 <- c(get_initial_par(pheno_ck,t),get_initial_par(pheno_salt,t),0.1,0.1,0.1,0.1)
  
  #par0 <- c( 19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,1.0660809,1.8476307,1.0659458,3.7009838 )
  par0 <- c(24.52088643,11.64758911,0.06179968,8.74285591,0.09998942,0.1,0.1)
  #par0 <- c(19.3898481,2.9533945,0.1775538,1.86229345  ,0.01507654,0.1,0.1)
  repeat{
    a <- optim(par=par0,Likelihood,pheno=pheno_deal,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood,pheno=pheno_deal,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H0 Finished",'\n')
  b
} 
Likelihood_loss = function(pheno,t,par){
  miu = Logistic(par,t)
  L0 = c()
  L0 = sum((apply(pheno, 2, mean)-miu)^2)
  return(L0)
}
optim_loss <- function(pheno,t){
  library(mvtnorm)
  itime <- 100
  itimes <- 1
  par0 <- c(16.19012753, -28.48257768,  41.94282944,  22.09364587 ,  0.02529833 ) ##Stress
  repeat{
    a <- optim(par=par0,Likelihood_loss,pheno=pheno,t=t,control=list(maxit=20000),method = "BFGS")
    
    b <- optim(a$par,Likelihood_loss,pheno=pheno,t=t,control=list(maxit=20000),method = "BFGS")
    
    cat("Logistic",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) == 0 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  b
  
}

###FUNMAP_diff
Logistic_diff <- function(par,t){
  
  y <- (par[1]/(1+par[2]*exp(-par[3]*t))-(par[4]*exp(-par[5]*t))) - (par[6]/(1+par[7]*exp(-par[8]*t))-(par[9]*exp(-par[10]*t)))
  return(y)
  
} 

SAD1_diff <- function(par, times = t) {   
  n <- ifelse (is.vector(times), length(times), NCOL(times) )   
  phi_ck<- par[1]   
  v2_ck <- par[2]
  phi_salt<- par[3]   
  v2_salt <- par[4] 
  
  tmp_ck <- (1-phi_ck^2)   
  sigma_ck <- array(1, dim=c(n,n))   
  
  tmp_salt <- (1-phi_salt^2)   
  sigma_salt <- array(1, dim=c(n,n)) 
  
  
  for(i in 1:n)   
  {     
    sigma_ck[i,i:n] <- phi_ck^( c(i:n) - i ) * (1-phi_ck^(2*i))/tmp_ck     
    sigma_ck[i:n,i] <- sigma_ck[i,i:n]   
    
    sigma_salt[i,i:n] <- phi_salt^( c(i:n) - i ) * (1-phi_salt^(2*i))/tmp_salt     
    sigma_salt[i:n,i] <- sigma_salt[i,i:n] 
    
  }   
  sigma_ck <- sigma_ck * abs(v2_ck)  
  sigma_salt <- sigma_salt * abs(v2_salt)
  sigma <- sigma_ck + sigma_salt
  
  return(sigma); 
}
get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno,2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),1,1)
}  
get_initial_par_unnorm <- function(pheno,t){
  mean0 <- apply(pheno,2,mean)  
  c(max(mean0),
    max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])),1,0.1)
}  


Likelihood_diff = function(pheno,t,par){
  miu = Logistic_diff(par,t)
  sigma = SAD1_diff(par = c(par[11],par[12],par[13],par[14]),times = t)
  L0 = c()
  L0 = sum(dmvnorm(pheno,miu,sigma,log = T))
  return(-L0)
}##pheno每一行为一个系号的表型数据
Likelihood_diff_H1_2 = function(pheno_AA,pheno_aa,t,par){
  
  miu1 = Logistic_diff(par[1:10],t)
  miu1.1 = Logistic_diff(par[11:20],t)
  
  sigma = SAD1_diff(par = c(par[21],par[22],par[23],par[24]),times = t)
  L1 = c()
  L1 = sum(dmvnorm(pheno_AA,miu1,sigma,log = T))
  L1.1 = c()
  L1.1 = sum(dmvnorm(pheno_aa,miu1.1,sigma,log = T))
  
  L0 <- L1 + L1.1 
  
  return(-L0)
}
Likelihood_diff_H1_3 = function(pheno_AA,pheno_aa,pheno_Aa,t,par){
  
  miu1 = Logistic_diff(par[1:10],t)
  miu1.1 = Logistic_diff(par[11:20],t)
  miu1.2 = Logistic_diff(par[21:30],t)
  
  sigma = SAD1_diff(par = c(par[31],par[32],par[33],par[34]),times = t)
  L1 = c()
  L1 = sum(dmvnorm(pheno_AA,miu1,sigma,log = T))
  L1.1 = c()
  L1.1 = sum(dmvnorm(pheno_aa,miu1.1,sigma,log = T))
  L1.2 = c()
  L1.2 = sum(dmvnorm(pheno_Aa,miu1.2,sigma,log = T))
  
  L0 <- L1 + L1.1 + L1.2
  return(-L0)
}

optim_diff <- function(pheno_ck,pheno_salt,pheno_diff,t){
  library(mvtnorm)
  itime <- 20
  itimes <- 1
  #par0 <- c(get_initial_par(pheno_ck,t),get_initial_par(pheno_salt,t),0.1,0.1,0.1,0.1)
  
  #par0 <- c( 19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,1.0660809,1.8476307,1.0659458,3.7009838 )
  par0 <- c(1.579665e+01, 7.406909e+00, 4.694319e-02, 1.812807e+02, 1.575449e+02, 3.618391e+00, 2.159658e+00, 7.963673e-02,
            1.011235e+02, 4.254119e+01, 1.028639e+00, 6.418383e-06, 1.065976e+00, 5.551713e+00)
  repeat{
    a <- optim(par=par0,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood_diff,pheno=pheno_diff,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H0 Finished",'\n')
  b
}  ##H0情况
optim_diff_H1_2 <- function(pheno_ck,pheno_salt,pheno_AA_diff,pheno_aa_diff,t){
  
  itime <- 20
  itimes <- 1
  #par0 <- c(19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,1.0660809,1.8476307,1.0659458,3.7009838)*1.1
  par0 <- c(1.579665e+01, 7.406909e+00, 4.694319e-02, 1.812807e+02, 1.575449e+02, 3.618391e+00, 2.159658e+00, 7.963673e-02,
            1.011235e+02, 4.254119e+01, 
            1.579665e+01, 7.406909e+00, 4.694319e-02, 1.812807e+02, 1.575449e+02, 3.618391e+00, 2.159658e+00, 7.963673e-02,
            1.011235e+02, 4.254119e+01, 
            1.028639e+00, 6.418383e-06, 1.065976e+00, 5.551713e+00)
  repeat{
    a <- optim(par=par0,Likelihood_diff_H1_2,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood_diff_H1_2,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H1_1 Finished",'\n')
  b
  
}  ##两基因型
optim_diff_H1_3 <- function(pheno_ck,pheno_salt,pheno_AA_diff,pheno_aa_diff,pheno_Aa_diff,t){
  
  itime <- 10
  itimes <- 1
  #par0 <- c(19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,19.3898481,2.9533945,0.1775538,280.3381369,172.3222069,5.8190876,1.3704297,0.9021405,-9.3612689, 2.5313579,1.0660809,1.8476307,1.0659458,3.7009838)
  par0 <- c(1.579665e+01, 7.406909e+00, 4.694319e-02, 1.812807e+02, 1.575449e+02, 3.618391e+00, 2.159658e+00, 7.963673e-02,
            1.011235e+02, 4.254119e+01, 
            1.579665e+01, 7.406909e+00, 4.694319e-02, 1.812807e+02, 1.575449e+02, 3.618391e+00, 2.159658e+00, 7.963673e-02,
            1.011235e+02, 4.254119e+01, 
            1.579665e+01, 7.406909e+00, 4.694319e-02, 1.812807e+02, 1.575449e+02, 3.618391e+00, 2.159658e+00, 7.963673e-02,
            1.011235e+02, 4.254119e+01,
            1.028639e+00, 6.418383e-06, 1.065976e+00, 5.551713e+00)
  
  repeat{
    a <- optim(par=par0,Likelihood_diff_H1_3,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,pheno_Aa=pheno_Aa_diff,t=t,control=list(maxit=1000))
    
    b <- optim(a$par,Likelihood_diff_H1_3,pheno_AA=pheno_AA_diff,pheno_aa=pheno_aa_diff,pheno_Aa=pheno_Aa_diff,t=t,control=list(maxit=1000))
    
    #cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$par-b$par) < 0.005 )||itimes == itime){ #itimes越高精度越高
      break
    }else{
      par0 <- b$par
    }
  }
  cat("H1_3 Finished",'\n')
  b
  
}  ##三基因型

get_VG_FunMap <- function(pheno_ck,pheno_salt,pheno_diff,marker,t){
  
  library(mvtnorm)
  diff_LR <- c() 
  diff_par <- c()
  for (a in 1:dim(marker)[1]) {
    
    AA <- as.numeric(names(which(marker[a,]==1)))
    aa <- as.numeric(names(which(marker[a,]==0)))
    Aa <- as.numeric(names(which(marker[a,]==2)))
    all <- as.numeric(names(which(marker[a,]!=9)))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    all_diff <- pheno_diff[pheno_diff[,1]%in%all,]
    
    AA_diff <- pheno_diff[pheno_diff[,1]%in%AA,]
    
    aa_diff <- pheno_diff[pheno_diff[,1]%in%aa,]
    
    optim_all_diff <- optim_diff(pheno_ck[,-1:-2],pheno_salt[,-1:-2],pheno_diff=all_diff[,-1:-2],t)
    
    
    if(NAa==0){ 
      
      optim_AA_aa_diff <- optim_diff_H1_2(pheno_ck[,-1:-2],pheno_salt[,-1:-2],AA_diff[,-1:-2],aa_diff[,-1:-2],t)
      
      LR_diff <- -2*(-optim_all_diff$value + optim_AA_aa_diff$value)  
      par_diff <- c(optim_all_diff$par,optim_AA_aa_diff$par,rep(NA,10))
      
      
    } else{
      
      Aa_diff <- pheno_diff[pheno_diff[,1]%in%Aa,]
      
      optim_AA_aa_Aa_diff <- optim_diff_H1_3(pheno_ck[,-1:-2],pheno_salt[,-1:-2],AA_diff[,-1:-2],aa_diff[,-1:-2],Aa_diff[,-1:-2],t)
      
      LR_diff <- -2*(-optim_all_diff$value + optim_AA_aa_Aa_diff$value)
      par_diff <- c(optim_all_diff$par,optim_AA_aa_Aa_diff$par)
    }
    diff_LR <- c(diff_LR,LR_diff)
    diff_par <- rbind(diff_par,par_diff)
    cat("SNP",a,"finished","\n")
    
  }
  
  return(list(diff_LR,diff_par) ) 
}

get_VG <- function(FunMap_par,marker_data,t){
  
  diff_vg <- c() 
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(names(which(marker_data[a,]==1)))
    aa <- as.numeric(names(which(marker_data[a,]==0)))
    Aa <- as.numeric(names(which(marker_data[a,]==2)))
    all <- as.numeric(names(which(marker_data[a,]!=9)))
    
    NAA <- length(AA)
    Naa <- length(aa)
    NAa <- length(Aa)
    
    p1 <- (NAA*2+NAa)/((NAA+NAa+Naa)*2) #A基因频率
    p0 <- (Naa*2+NAa)/((NAA+NAa+Naa)*2) #a基因频率
    
    mean_AA <- Logistic_diff(FunMap_par[a,][15:24],t)
    mean_aa <- Logistic_diff(FunMap_par[a,][25:34],t)
    AE <- (mean_AA - mean_aa)/2 
    
    if(NAa==0){ Vg <- 2*p1*p0*(AE^2)  } else{
      
      mean_Aa <- Logistic_diff(FunMap_par[a,][35:44],t)
      
      DE <- mean_Aa - (mean_AA + mean_aa)/2
      Vg <- 2*p1*p0*((AE + (p1 - p0)*DE)^2) + 4*p1*p1*p0*p0*DE*DE
      
    }
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  colnames(diff_vg) <- seq(min(t),max(t),length.out = length(t))
  
  return(sqrt(diff_vg)) 
}###根据求得的参数计算遗传标准差

get_cluster <- function(cutree,cluster_num,plast_Value){
  
  allcluster <- list()    
  for (a in 1:cluster_num) {
    
    allcluster[[a]] <- as.data.frame(plast_Value)[which(cutree==c(1:cluster_num)[a]),]
    
  }
  allcluster
} #将每一类的样本表型数据整合在一起

##变量选择与微分方程用
Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL ){
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}
dLegendre.model <-function( t, mu, tmin=NULL, tmax=NULL ){
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1]*0 + 1*mu[2];
  if (np.order>=2)
    L <- L + 0.5 * (6 * ti)* mu[3] ;
  if (np.order>=3)
    L <- L +0.5 * (15 * ti ^ 2 - 3)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)* mu[5];
  if (np.order>=5)
    L <- L + 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)*mu[6];
  if (np.order>=6)
    L <- L + (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *ti)* mu[7];
  if (np.order>=7)
    L <- L + (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)* mu[8];
  return(L);
}

Likehood_Legendre_ind <- function(times,para,E){
  
  sum( (E -Legendre.model(t=times,mu=para))^2)
  
}
Likehood_Legendre_cluster <- function(times,para,E){
  
  sum((apply(E, 2, mean)-Legendre.model(t=times,mu=para))^2)
  
} #类，先平均再最小二乘

smooth.optim_ind <- function(times,para,effect_value){
  
  mean_par <- c()
  
  L <- optim(para,Likehood_Legendre_ind,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_ind,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_ind,E=effect_value,times=times,method="BFGS")
  mean_par <- L$par
  mean_par
}
smooth.optim_cluster <- function(times,para,effect_value){
  
  mean_par <- c()
  
  L <- optim(para,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  L <- optim(L$par,Likehood_Legendre_cluster,E=effect_value,times=times,method="BFGS")
  mean_par <- L$par
  mean_par
}

get_cluster <- function(cutree,cluster_num,plast_Value){
  
  allcluster <- list()    
  for (a in 1:cluster_num) {
    
    allcluster[[a]] <- as.data.frame(plast_Value)[which(cutree==c(1:cluster_num)[a]),]
    
  }
  allcluster
} #将每一类的样本表型数据整合在一起

varsel1 <- function(X,Y,tt,order=6){
  
  nobs = nrow(X)
  ndim = ncol(X)
  dfo = rep(order-1,ndim)
  index = rep(1:ndim,times=dfo)
  aa2 <- c()
  for(k in 1:ndim){
    aa1 <- c()
    for(j in 1:(order-1)){
      aa <- c()
      for(i in 1:length(tt)){
        aa <- c(aa,Legendre.modelve(tt[i],np.order = j,tmin = min(tt), tmax = max(tt)))
      }
      aa1 <- cbind(aa1,aa*X[,k])
    }
    aa2 <- cbind(aa2,aa1)
  }
  
  
  Xc = scale(aa2,center=T,scale=T)
  n = nrow(Xc)
  
  connect = matrix(0,nrow=ndim,ncol=ndim)
  coefest = matrix(0,nrow=sum(dfo),ncol=ndim)
  regfun = vector("list",length=ndim)
  for(i in 1:ndim)
  {
    inde <- index[-which(index==i)]
    XXc <- Xc[,-which(index==i)]
    
    yc <- Y[,i]-mean(Y[,i])
    
    out1 <- GrpLasso(X=XXc,y=yc,index=inde,lambda=30,crit="BIC")
    var.grp <- out1$var.select  # genes selected
    coef.grp <- out1$coef
    
    ### Adaptive Group Lasso
    index.adp <- index[is.element(index,var.grp)]
    W.adp = sapply(1:length(var.grp),function(j) sqrt(sum(coef.grp[index.adp==var.grp[j]]^2)))
    Xc.adp = Xc[,is.element(index,var.grp)]
    Xcs.adp = scale(Xc.adp,center=F,scale=rep(1/W.adp,times=dfo[var.grp]))
    init.adp = coef.grp/rep(W.adp,times=dfo[var.grp])
    lambda = lambdamax(Xcs.adp,yc,index=index.adp,coef.ini=init.adp,
                       penscale=sqrt,center=F,standardize=F,model=LinReg())*0.7^(1:18)
    out2 = GrpLasso(X=Xc.adp,y=yc,W=W.adp,index=index.adp,ini=coef.grp,
                    lambda=lambda,crit="BIC")
    var.adp = out2$var.select
    coef.adp = out2$coef
    connect[i,var.adp] <-  1
    coefest[is.element(index,var.adp),i] <-  coef.adp
    regfun[[i]] <-  sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    connect[i,var.adp] = 1
    coefest[is.element(index,var.adp),i] = coef.adp
    regfun[[i]] = sapply(var.adp,function(ii) Xc[,is.element(index,ii)]%*%coefest[is.element(index,ii),i])
    cat("var=",i,var.adp,"\n")
  }
  return(list(connect=connect,regfun=regfun,coefest=coefest))
}

Legendre.modelve <- function(t, np.order, tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- NA
  L <- 1;
  if (np.order >= 2)
    L <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L <- 0.5 * (15 * ti ^ 2 - 3)
  if (np.order >= 4)
    L <- 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)
  if (np.order >= 5)
    L <- 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L <- (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                       ti)
  if (np.order >= 7)
    L <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                       ti ^ 2 - 35) 
  return(L);
}

GrpLasso <- function(X,y,W=NULL,index,ini=rep(0,ncol(X)),lambda=NULL,crit=c("BIC","EBIC"),center=F){
  if(center==T){
    y = y-mean(y)
    X = scale(X,center=T,scale=F)
  }
  n = nrow(X)
  ind = unique(index)
  p = length(ind)
  dfo = sapply(1:p,function(j) sum(index==ind[j]))
  
  # fit model for a sequence of penalty parameters
  if(!is.null(W)){
    X = scale(X,center=F,scale=rep(1/W,times=dfo))
    ini = ini/rep(W,times=dfo)
  }
  
  # set up the candidates for penalty parameter
  if(is.null(lambda)){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:20)
  }else if(length(lambda)==1){
    lambda = lambdamax(X,y,index=index,coef.ini=ini,penscale=sqrt,center=F,
                       standardize=F,model=LinReg())*0.9^(1:lambda)
  }
  
  fit = grplasso(X,y,index=index,lambda=lambda,model=LinReg(),center=F,
                 standardize=F,coef.ini=ini,penscale=sqrt,
                 control=grpl.control(update.hess="lambda",tol=10^-8,trace=0))
  
  # calculate BIC/EBIC
  nlambda = length(lambda)
  rss = sapply(1:nlambda,function(j) sum((y-fit$fitted[,j])^2))
  var.select = sapply(1:nlambda,function(j) unique(index[abs(fit$coef[,j])>0]))
  dfo.lambda = sapply(1:nlambda,function(j) sum(dfo[is.element(ind,var.select[[j]])]))
  if(crit!="BIC" & crit!="EBIC"){
    cat("Error: Criterion not implemented. Reset to BIC!\n")
    crit = "BIC"
  }
  if(crit=="BIC"){
    bic = log(rss)+dfo.lambda*log(n)/n
  }else if(crit == "EBIC"){
    bic = log(rss)+dfo.lambda*log(n)/n+0.5*dfo.lambda*log(p)/n
  }
  
  # select model with smallest value of selection criterion
  id.ss = which.min(bic)
  var.ss = var.select[[id.ss]]
  fit.ss = fit$fitted[,id.ss]
  coef.ss = fit$coef[,id.ss]
  if(!is.null(W)){
    coef.ss = coef.ss*rep(W,times=dfo)
  }
  coef.ss = coef.ss[is.element(index,var.ss)]
  
  return(list(var.select=var.ss,coefficients=coef.ss,fitted=fit.ss,BIC=bic,
              lambda=lambda,fit.org=fit))
}


optim.parallel <- function(connect,effect,n.cores,proc,order,times,nstep){
  
  diag(connect) <- 1
  
  LL <- L_derivative(nt=times,nstep=nstep,order=order)
  
  nx <- dim(effect)[2]
  
  grp <- floor(nx/n.cores)
  grp.i <- c()
  if(n.cores==1){
    grp.i <- c(grp.i,rep(1,nx))
  }else{
    for(ii in 1:n.cores){
      if(ii==n.cores){
        grp.i <- c(grp.i,rep(ii,nx-grp*(ii-1)))
      }else{
        grp.i <- c(grp.i,rep(ii,grp))
      }
    }
  }
  
  grp.ii <- unique(grp.i)
  
  res.list <- mclapply(grp.ii, function(i)
  {
    y.c <- 	which(grp.i==i)
    A <- sapply(y.c, proc, connect=connect,effect=effect,LL=LL,nstep=nstep,order=order,times=times);
    return (unlist(A));
  }, mc.cores=n.cores )
  
  res1 <- do.call("c", res.list)
  res2 <- parallel.data.optim(res1,connect,times)
  return(res2)
}
L_derivative <- function(nt,nstep,order){
  
  stp <- (max(nt)-min(nt))/nstep
  res <- c()
  for(j in 1:nstep){
    
    tg1 <- Legendre.model1((j-1)*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg2 <- Legendre.model1(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg3 <- Legendre.model1(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg4 <- Legendre.model1(j*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tmp1 <- rbind(tg1,tg2,tg3,tg4)
    res <- rbind(res,tmp1)
  }
  res
}
Legendre.model11 <- function(t, np.order,tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- 1;
  if (np.order >= 2)
    L[2] <- 0.5 * (6 * ti) 
  if (np.order >= 3)
    L[3] <- 0.5 * (15 * ti ^ 2 - 3) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * 4 * ti ^ 3 - 60 * ti) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *
                         ti) 
  if (np.order >= 7)
    L[7] <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *
                          ti ^ 2 - 35)
  return(L);
}
Legendre.model1 <- function(t, np.order,tmin = NULL, tmax = NULL){
  u <- -1;
  v <- 1;
  if (is.null(tmin))
    tmin <- min(t);
  if (is.null(tmax))
    tmax <- max(t);
  ti    <- u + ((v - u) * (t - tmin)) / (tmax - tmin);
  L <- rep(NA,np.order)
  L[1] <- ti;
  if (np.order >= 2)
    L[2] <- 0.5 * (3 * ti^2 - 1) 
  if (np.order >= 3)
    L[3] <- 0.5 * (5 * ti ^ 3 - 3 * ti) 
  if (np.order >= 4)
    L[4] <-  0.125 * (35 * ti ^ 4 - 30 * ti^2 + 3) 
  if (np.order >= 5)
    L[5] <-  0.125 * (63 * ti ^ 5 - 70 * ti ^ 3 + 15*ti)
  if (np.order >= 6)
    L[6] <-(1 / 16) * (231 * ti ^ 6 - 315 * ti ^ 4 + 105 * ti^2 -5) 
  
  return(L);
}
ode.optim <- function(y.c,connect,effect,LL,nstep,order,times){
  self <- y.c 
  indexx <- which(connect[y.c,]==1)
  #para <- rep(0.00001,length(indexx)*(order-1))
  para <- rep(0.001,length(indexx)*(order-1))
  #res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
  #            LL=LL,method="BFGS",control=list(maxit=2000,trace=T))
  res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
               LL=LL,method = "BFGS",control=list(maxit=4000,trace=F))
  cat("Gene=",y.c," ",res$value,"\n")
  A <- ode.sovle.ind(NG=(effect),res$par,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,LL=LL,self=self)
  return(A)
}

fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){
  
  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL,self = self)
  sum((NG[,self]-rowSums(odes))^2) ##最小二乘

  }

ode.sovle.ind <- function(NG,fitpar,nconnect,nt,order,nstep,LL,self){
  stp <- (max(nt)-min(nt))/nstep
  index <- which(nconnect==1)
  
  ind.par <- matrix(fitpar[1:(length(index)*(order-1))],ncol=order-1,byrow=T)
  allrep <- matrix(rep(0,length(index)),nrow=1)
  allrep[which(index==self)] <- NG[1,self]
  nn <- 1
  for(j in 1:nstep){
    tg1 <- (rowSums(t(apply(ind.par,1,"*",LL[nn,])))*NG[j,index])
    tg2 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+1,])))*NG[j,index])
    tg3 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+2,])))*NG[j,index])
    tg4 <- (rowSums(t(apply(ind.par,1,"*",LL[nn+3,])))*NG[j,index])
    tmp <- allrep[j,] +stp*(tg1+2*tg2+2*tg3+tg4)/6
    allrep <- rbind(allrep,tmp)
    nn <- nn + 4
  }
   
  self_name <- colnames(NG)[self]
  no_name <- which(colnames(allrep)==self_name)
  allrep[,no_name] <- allrep[,no_name]
  allrep
}


parallel.data.optim <- function(rd,nm,ntt){
  
  nrd <- matrix(rd,nrow=length(ntt))
  nn <- dim(nm)[1]
  ki <- 0
  allist <- list()
  for(i in 1:nn){
    iii <- (which(nm[i,]==1))
    iiil <- length(iii)
    tmp.d <- nrd[,(ki+1):(ki+iiil)]
    if(is.matrix(tmp.d)){
      colnames(tmp.d) <- iii
    }else{
      names(tmp.d) <- iii
    }
    
    allist[[i]] <- tmp.d
    ki <- ki + iiil
  }
  
  return(allist)
}

regasso <- function(connect1,gene,interaction){
  
  diag(connect1) <- 0
  ng <- dim(gene)[1]
  allcor <- list()
  for(i in 1:ng){
    a1 <- gene[i,]
    nng1 <- (interaction[[i]])
    if(!is.matrix(nng1)){
      next
    }else{
      nng <- as.matrix(nng1[,-which(colnames(nng1)==i)])
      corr <- c()
      for(j in 1:dim(nng)[2]){
        corr <- c(corr,cor(a1,nng[,j]))
      }
      allcor[[i]] <- corr
    }
  }
  for (q in 1:length(allcor)) {
    connect1[q,which(connect1[q,]==1)] <- allcor[[q]]
  }
  
  return(connect1)
}
regasso_moudule <- function(connect1,gene,interaction){
  aaa <- colnames(connect1)
  diag(connect1) <- 0
  ng <- dim(gene)[1]
  allcor <- list()
  for(i in 1:ng){
    a1 <- gene[i,]
    nng1 <- (interaction[[i]])
    if(!is.matrix(nng1)){
      next
    }else{
      nng <- as.matrix(nng1[,-which(colnames(nng1)==aaa[i])])
      corr <- c()
      for(j in 1:dim(nng)[2]){
        corr <- c(corr,cor(as.numeric(a1),as.numeric(nng[,j])))
      }
      allcor[[i]] <- corr
    }
  }
  
  for (q in 1:length(allcor)) {
    connect1[q,which(connect1[q,]==1)] <- allcor[[q]]
  }
  
  return(connect1)
}


##解微分方程，求得个SNP互作数据
optim_inter <- function(all_cluster_value,ind_Lpar11,norder){
  
  library(splines)
  library(orthogonalsplinebasis)
  library(MASS)
  library(grplasso)
  library(parallel)
  
  ind_Lpar1 <- ind_Lpar11
  clusterAA <- all_cluster_value
  
  
  ind_LparAA <- ind_Lpar1[as.numeric(rownames(clusterAA)),]
  rownames(ind_LparAA) <-  rownames(clusterAA) #调出先前求出的该类个体的勒让德拟合参数
  
  ttt <- seq(13,78,length.out = (dim(ind_LparAA)[1]+5))
  cAA_mean <- apply(ind_LparAA,1, Legendre.model,t=ttt)
  cAA_d <- apply(ind_LparAA,1, dLegendre.model,t=ttt)
  
  clusterAA_50 <- varsel1(X=cAA_mean,Y=cAA_d,tt=ttt)
  colnames(clusterAA_50$connect) <- as.numeric(rownames(clusterAA))
  rownames(clusterAA_50$connect) <- as.numeric(rownames(clusterAA))
  
  
  clusterAA_list <- list()
  for (n in 1:dim(ind_LparAA)[1]) {
    
    ttemp <- list()
    ttemp[[1]] <-  rownames(clusterAA)[n]
    ttemp[[2]] <- names(which(clusterAA_50$connect[n,]==1))
    clusterAA_list[[n]] <- ttemp
  }
  
  
  tryyAA <- optim.parallel(connect=clusterAA_50$connect,effect=t(all_cluster_value),
                           n.cores=1,proc=ode.optim,order=norder,times=seq(13,78,length=30),nstep=29)
  
  for (e in 1:length(tryyAA)) {
    
    colnames(tryyAA[[e]]) <- rownames(clusterAA)[as.numeric(colnames(tryyAA[[e]]))]
    
  } #重命名
  
  ##计算相关系数
  inter_effect <- regasso_moudule(connect1=clusterAA_50$connect,gene=all_cluster_value,interaction=tryyAA)
  
  
  
  for (n in 1:length(clusterAA_list)) {
    
    clusterAA_list[[n]][[3]] <- c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0])[order(as.numeric(names(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0]))))]
    
  }
  
  
  ##计算相关系数
  inter_effect <- regasso_moudule(connect1=clusterAA_50$connect,gene=all_cluster_value,interaction=tryyAA)
  
  
  
  for (n in 1:length(clusterAA_list)) {
    
    
    nono <- which(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0])[order(as.numeric(names(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0]))))]==0)
    
    clusterAA_list[[n]][[3]] <- colSums(get_integrate(out=tryyAA[[n]],t=seq(13,78,length=30)))
    
    clusterAA_list[[n]][[3]][nono] <- 0
  } 
  
  
  
  ##assign(paste("optim_cluster_",which_cluster,"list",sep=""),clusterAA_list, envir = .GlobalEnv)
  optim_cluster <- list()
  optim_cluster[[1]] <- clusterAA_list
  optim_cluster[[2]] <- tryyAA
  optim_cluster
}


get_integrate <- function(out,t){
  s <- t[2] - t[1]
  c <- c()
  for (n in 1:length(out[1,])) {
    a <- out[,n]
    b <- c()
    for (i in 1:(length(out[,1])-1)) {
      b <- c(b,a[i] * s)
    }
    c <- cbind(c,b)
  }
  colnames(c) <- colnames(out)
  c/100
}
##根据上述求得的数据输出画图用的数据
out_Netdata <- function(optim_cluster,which_cluster,clusterAA){
  
  clusterAA_list <- optim_cluster[[1]]
  afterAA <- c()
  effect_selfAA <- c()
  for (i in 1:length(clusterAA_list)){
    dep <- clusterAA_list[[i]][[1]]
    ind <- clusterAA_list[[i]][[2]]
    effectAA <- clusterAA_list[[i]][[3]]
    effect_selfAA <- rbind(effect_selfAA,effectAA[which(names(effectAA)==dep)])
    effectAA <-effectAA[-which(names(effectAA)==dep)]
    one <- c()
    for (j in 1:length(ind)) {
      if(effectAA[j] >= 0){
        type <- 1
      }else{
        type <- 0
      }
      one <- rbind(one,c(ind[j],dep,abs(effectAA[j]),type))
    }
    afterAA <- rbind(afterAA,one)
    
  }
  eageAA <- cbind(afterAA[,-3],c(0:(length(afterAA[,1])-1)),afterAA[,3])
  colnames(eageAA) <- c("source","target","colour","Id","Weight")
  
  ttype <- c()
  for (n in 1:length(effect_selfAA)) {
    
    if(effect_selfAA[n] >= 0){
      ttype[n] <- 1
    }else{
      ttype[n] <- 0
    }
    
  }
  effect_selfAA <- cbind(rownames(clusterAA),paste("SNP",rownames(clusterAA),sep = ""),abs(effect_selfAA),ttype)
  colnames(effect_selfAA) <- c("ID","Label","weight","type")
  
  
  write.csv(effect_selfAA,paste("node",which_cluster,".csv",sep = ""), row.names = F,fileEncoding = "UTF-8")
  write.csv(eageAA,paste("eage",which_cluster,".csv",sep = ""), row.names = F,fileEncoding = "UTF-8")
  
  get_eageInform(eageAA,clusterAA) 
  
}
out_Netdata_cytoscape <- function(optim_cluster,which_cluster,clusterAA){
  
  clusterAA_list <- optim_cluster[[1]]
  afterAA <- c()
  effect_selfAA <- c()
  for (i in 1:length(clusterAA_list)){
    dep <- clusterAA_list[[i]][[1]]
    ind <- clusterAA_list[[i]][[2]]
    effectAA <- clusterAA_list[[i]][[3]]
    effect_selfAA <- rbind(effect_selfAA,effectAA[which(names(effectAA)==dep)])
    effectAA <-effectAA[-which(names(effectAA)==dep)]
    one <- c()
    for (j in 1:length(ind)) {
      if(effectAA[j] >= 0){
        type <- 1
      }else{
        type <- 0
      }
      one <- rbind(one,c(ind[j],dep,abs(effectAA[j]),type))
    }
    afterAA <- rbind(afterAA,one)
    
  }
  eageAA <- cbind(afterAA[,-3],c(0:(length(afterAA[,1])-1)),afterAA[,3])
  colnames(eageAA) <- c("source","target","colour","Id","Weight")
  
  new <- eageAA
  new[,1] <- paste("S",new[,1],sep = "")
  new[,2] <- paste("S",new[,2],sep = "")
  write.csv(new,paste("M",which_cluster,".csv",sep = ""), row.names = F,fileEncoding = "UTF-8")
  
  get_eageInform(eageAA,clusterAA) 
  
}

get_eageInform <- function(eage,cluster){
  
  eage_sati <- c()
  for (n in rownames(cluster)) {
    
    
    if(any(which(eage[,1]==n))){
      a1 <- eage[which(eage[,1]==n),]
      a1.0 <- 0
      a1.1 <- 0
      if(!is.matrix(a1)){
        if(a1[3]==0){
          a1.0 <- 1
          a1.1 <- 0}
        else{
          a1.0 <- 0
          a1.1 <- 1}}
      else{
        for (q in 1:dim(a1)[1]) {
          a2 <- as.matrix(a1)[q,]
          if(a2[3]==0){
            a1.0 <- a1.0+1
            a1.1 <- a1.1+0
          }
          else{
            a1.0 <- a1.0+0
            a1.1 <- a1.1+1
          }
        }
      }
      a3 <- c(a1.0,a1.1)
      
    }
    else{a3 <- c(0,0)}
    
    
    b1.0 <- 0
    b1.1 <- 0
    b1 <- eage[which(eage[,2]==n),]
    if(!is.matrix(b1)){
      if(b1[3]==0){
        b1.0 <- 1
        b1.1 <- 0
      }
      else{
        b1.0 <- 0
        b1.1 <- 1
      }
    }
    else{
      for (w in 1:dim(b1)[1]) {
        b2 <- b1[w,]
        if(b2[3]==0){
          b1.0 <- b1.0+1
          b1.1 <- b1.1+0
        }
        else{
          b1.0 <- b1.0+0
          b1.1 <- b1.1+1}
      }
    }
    b3 <- c(b1.0,b1.1)
    cc <- c(a3,b3)
    eage_sati <- rbind(eage_sati,cc)
    
  }
  eage_sati <- cbind(rownames(cluster),eage_sati)
  colnames(eage_sati) <- c("no","sour_con","sour_pro","tar_con","tar_pro")
  rownames(eage_sati) <- c(1:dim(cluster)[1])
  eage_sati
  
  
}

out_ModuleInform <- function(cluster_list,cluster_optim,inter_effect,ClusterNum){
  
  for (n in 1:length(cluster_list)) {
    
    cluster_list[[n]][[3]] <- c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0])[order(as.numeric(names(c(inter_effect[n,][n],inter_effect[n,][inter_effect[n,]!=0]))))]
    
  }
  
  
  after <- c()
  effect_self <- c()
  for (i in 1:length(cluster_list)){
    dep <- cluster_list[[i]][[1]]
    ind <- cluster_list[[i]][[2]]
    effect <- cluster_list[[i]][[3]]
    effect_self <- rbind(effect_self,effect[which(names(effect)==dep)])
    effect <-effect[-which(names(effect)==dep)]
    one <- c()
    for (j in 1:length(ind)) {
      if(effect[j] >= 0){
        type <- 1
      }else{
        type <- 0
      }
      one <- rbind(one,c(ind[j],dep,abs(effect[j]),type))
    }
    after <- rbind(after,one)
    
  }
  eage <- cbind(after[,-3],c(0:(length(after[,1])-1)),after[,3])
  colnames(eage) <- c("source","target","colour","Id","Weight")
  
  ttype <- c()
  for (n in 1:length(effect_self)) {
    
    if(effect_self[n] >= 0){
      ttype[n] <- 1
    }else{
      ttype[n] <- 0
    }
    
  }
  effect_self <- cbind(c(1:length(cluster_list)),paste("M",c(1:403),sep = ""),abs(effect_self),ttype)
  colnames(effect_self) <- c("ID","Label","weight","type")
  
  rownames(cluster_mean) <- c(1:ClusterNum)
  get_eageInform(eage,cluster_mean)
  
  write.csv(effect_self,"Module-node.csv", row.names = F,fileEncoding = "UTF-8")
  new <- eage
  new[,1] <- paste("M",new[,1],sep = "")
  new[,2] <- paste("M",new[,2],sep = "")
  write.csv(new,"Module-eage.csv", row.names = F,fileEncoding = "UTF-8")
  
  Link_Module <- get_eageInform(eage,cluster_mean)
  Link_Module[,1] <- paste("M",Link_Module[,1],sep = "")
  write.csv(Link_Module,"Link_Module.csv", row.names = F,fileEncoding = "UTF-8")
  
  
}

# FunCluster

LgdP <- expression( tt,
                    ( 3* tt^2 - 1 )/2 , 
                    ( 5 * tt^3 - 3* tt )/2, 
                    ( 35 * tt^4 - 30 * tt^2 + 3)/8,
                    ( 63 * tt^5 - 70 * tt^3 + 15 * tt )/8,
                    ( 231 * tt^6 - 315 * tt^4 + 105 * tt^2 - 5)/16,
                    ( 429 * tt^7 - 693 * tt^5 + 315 * tt^3 - 35 * tt)/16,
                    ( 6435 * tt^8 - 12012 * tt^6 + 6930 * tt^4 - 1260 * tt^2 + 35)/128,
                    ( 12155 * tt^9 - 25740 * tt^7 + 18018 * tt^5 - 4620 * tt^3 + 315 * tt)/128,
                    ( 46189 * tt^10 - 109395 * tt^8 + 90090 * tt^6 - 30030 * tt^4 + 3465 * tt^2 - 63)/256 )

GetMR <- function(rho,times)
{
  MR <- matrix(1,length(times),length(times))
  for ( i in 1:length(times)){
    for(j in 1:length(times)){
      MR[i,j]= rho^(abs(times[j] - times[i]))
    }
  }
  return (MR)
}
GetMX <- function(times,r)
{
  tnum = length(times) 
  X <- matrix(1,tnum,r+1)
  
  for(t in 1:tnum ){
    tt <- -1 + 2*(times[t] - times[1])/(times[tnum] - times[1])
    for(i in 1:r){
      X[t,i+1] <- eval(LgdP[i])
    }
  }
  return (X)
}
GetInitPij <- function(N,J)
{
  P <- matrix(1/J,N,J)
  for (i in 1:N){
    P[i,] <- rnorm(J, mean=1/J, sd= 0.2 * 1/J )
    P[i,] <- P[i,]/sum(P[i,])
  }
  
  return (P)
}
GetMeanMatrix <- function(J,times,P,X,Asdata,InvMSigema)
{
  m <- matrix(NA,J,length(times))
  N <- length(Asdata[,1])
  r <- length(X[1,])
  
  xInvSigema <- t(X) %*% InvMSigema
  xInvSigemax <- xInvSigema%*% X
  
  mU <- matrix(NA, J, ncol(X))
  for( j in 1:J){
    ud <- matrix(0, r, r)
    for( i in 1: N){
      ud <- ud + P[i,j]*xInvSigemax
    }
    ubd <- matrix(0, r, 1)
    for( i in 1: N){
      ubd <- ubd + P[i,j]*( xInvSigema %*% (Asdata[i,]) )
    }
    uj <- ginv(ud) %*% ubd
    m[j,] <- X %*% uj
    mU[j,] <- uj
  }
  return(list(M = m, U = mU))
}
GetNewSsquare <- function(Asdata,m,MR,times,P,J)
{
  N <- length(Asdata[,1])
  
  InvMR <- ginv(MR)
  newSsquare <- 0
  for(i in 1:N){
    SumJ <- 0
    for(j in 1:J){
      yi_mj <- Asdata[i,]-m[j,]
      SumJ <- SumJ + P[i,j] * ((yi_mj) %*% InvMR %*% (yi_mj) )
    }
    newSsquare <- newSsquare + SumJ
  }
  newSsquare <- as.numeric(newSsquare/(length(times)*N))
  
  return(newSsquare)
}
GetNewRho.b <- function(rho,rhoDir)
{
  
  newrho <- as.numeric(rho + 0.005*rhoDir)##0.004
  if (newrho > 0.99) newrho <- 0.99
  if (newrho < 0) newrho <- 0
  
  return (newrho)
}
GetNewRho<- function(Asdata,m,MR,times,P,J,rho,Ssquare)
{
  N <- length(Asdata[,1])
  newrho <- 0
  for(i in 1:N){
    SumJ <- 0
    for(j in 1:J){
      yi_mj <- Asdata[i,]-m[j,]
      Item1 <- (1/(1 - rho*rho))*((yi_mj) %*% MR %*% (yi_mj) )
      Item2 <- 0
      for(k in 2:(length(times)-1) )
        Item2 <- Item2 + (yi_mj[k]^2)
      Item2 <- Item2 * rho
      Item3 <- 0
      for(k in 1:(length(times)-1) )
        Item2 <- Item3 + yi_mj[k] * yi_mj[k+1]
      SumJ <- SumJ + P[i,j] * (Item1 + Item2 - Item3)
    }
    newrho <- newrho + SumJ
  }
  newrho <- as.numeric(newrho/( (length(times)-1)* N * Ssquare))
  
  if(abs(newrho) >= 1) return( sign(newrho)*.5)
  else return(newrho)
}
GetLikelihood <- function(Asdata,m,omiga,InvMSigema,DetMSigema,P,J,times) 
{
  N <- length(Asdata[,1])
  
  LogDetMSigema <- log(DetMSigema)/2
  LogM2Pi <- length(times)*log(2*pi)/2
  
  oneterm <- function(i, j) {
    f <- function(i,j)
      P[i,j]*(log(omiga[j]) - LogM2Pi - LogDetMSigema
              - ( ((Asdata[i,]-m[j,])) %*% InvMSigema %*% (Asdata[i,]-m[j,])) /2)
    mapply(f, i, j)
  }
  tmp <- outer(1:N, 1:J, oneterm)
  tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
  return(sum(tmp))
}
StepE <- function(Asdata,m,omiga,InvMSigema,DetMSigema,P,J,times)
{
  TwoPiExp <- (2*pi) ^ ( length(times)/2 )
  TwoPiExp <- TwoPiExp * DetMSigema
  
  N <- length(Asdata[,1])
  
  tmp <- rep(0,N)
  for( i in 1:N){
    Fi <- rep(0,J)
    for( j in 1:J){
      yi_mj <- Asdata[i,]-m[j,]
      Fi[j] = exp( ( (yi_mj) %*% InvMSigema %*% (yi_mj) ) / -2) / TwoPiExp 
    }
    OmigaF <- omiga %*% Fi
    P[i,] <- (omiga * Fi) / OmigaF
    tmp[i] <- log(OmigaF)
  }
  tmp[!is.finite(tmp)] <- min(tmp[is.finite(tmp)])
  P[!is.finite(P)] <- max(P[is.finite(P)]) ######
  Likelihood <- sum(tmp)
  
  return(list(P = P, Likelihood = Likelihood))
}

StepM.b <- function(Asdata,m,MR,times,Ssquare,P,rho,rhoDir,J,rpt)
{
  newSsquare <- GetNewSsquare(Asdata,m,MR,times,P,J)
  if (rpt > 0)
    newrho <- GetNewRho.b(rho,rhoDir)
  else
    newrho <- rho
  
  return( c(newSsquare, newrho))
}


RunEM.C <- function(Asdata,times,rho,Ssquare,X,P,MR,MSigema,InvMSigema,DetMSigema,omiga,m,J,r)
{
  IncLimit <- 3
  REPEAT_LIMIT <- 30  #200
  LIKELIHOOD_DIFF <- 0.5 #0.01
  
  rpt <- 1
  Likelihood <- -Inf
  
  rhoDir <- 1
  rhoIncCount <- 0
  
  while(TRUE){
    cat("rpt= ",rpt)
    OldLikelihood <- Likelihood
    EResult <- StepE(Asdata,m,omiga,InvMSigema,DetMSigema,P,J,times)
    P <- EResult$P
    Likelihood <- EResult$Likelihood
    
    if( (abs(OldLikelihood - Likelihood) < LIKELIHOOD_DIFF) ){
      #cat("quit due to likelihood\n")
      #cat("LIKELIHOOD_DIFF:",LIKELIHOOD_DIFF,"\n")
      break
    }
    
    if( rpt >= REPEAT_LIMIT ){
      #cat("quit due to rpt\n")
      break
    }
    
    if ( Likelihood >= OldLikelihood){
      #rhoIncCount <- 0
    }else{
      #cat("decrease, limit:", IncLimit, "\n")
      rhoIncCount <- rhoIncCount + 1
      if (rhoIncCount >= IncLimit){
        rhoIncCount <- 0
        rhoDir <- rhoDir * -1
      }
    }
    
    newpars <- StepM.b(Asdata,m,MR,times,Ssquare,P,rho,rhoDir,J,rpt)
    #     cat("newpars:\n")
    #     print(newpars)
    Ssquare <- newpars[1]
    rho <- newpars[2]
    
    MR <- GetMR(rho,times)
    MSigema <- Ssquare * MR
    # print(MSigema)
    InvMSigema <- ginv(MSigema)
    DetMSigema <- (det(MSigema))^0.5
    
    N <- length(Asdata[,1])
    omiga <- colSums(P)/ N 
    
    rMeans <- GetMeanMatrix(J,times,P,X,Asdata,InvMSigema)
    m <- rMeans$M
    mU <- rMeans$U
    rpt <- rpt + 1
  }
  
  return(list(rho=rho,Ssquare=Ssquare,Likelihood=Likelihood,m=m,mU=mU,P=P))
}


InitAndRunEM <- function(mExpression, times, J=5, r=4)  # J: Genes; r :order of legendre polynominal
{
  rho <- 0.8
  Ssquare <- 20
  
  X <- GetMX(times,r)
  N <- length(mExpression[,1])
  P <- GetInitPij(N,J)
  
  MR <- GetMR(rho,times)
  MSigema <- Ssquare * MR
  #print(MSigema)
  InvMSigema <- ginv(MSigema)
  DetMSigema <- (det(MSigema))^0.5
  
  omiga <- colSums(P)/ N
  m <- GetMeanMatrix(J,times,P,X,mExpression,InvMSigema)$M
  
  EMResults <- RunEM.C(mExpression,times,rho,Ssquare,X,P,MR,MSigema,InvMSigema,DetMSigema,omiga,m,J,r) 
  EMResults
}

GeneCluster <- function(mExpression, times, NumberOfCluster,orderLOP){ 
  r=orderLOP
  J=NumberOfCluster
  
  if(!requireNamespace("MASS",quietly = TRUE)){
    stop("MASS needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  
  if (r>10) {warning("range of r is 1 to 10")}
  if (r<0) {warning("r must be greater than 0")}
  if (length(times)!= ncol(mExpression)) {warning(" length of time vector is different to column length ")}
  if (J>nrow(mExpression)) {warning("number of cluster larger than number of genes")}
  if (!is.numeric(mExpression)) {warning("data is not numeric ")}
  
  
  lgd=r
  
  
  a1 <- times
  LineColors <- colors()
  
  opar <- par(no.readonly = TRUE)
  par(col.axis="blue", mar=c(4, 4, 2.5, 0.25),pty = "m",lwd=1)
  
  Asdata <- mExpression
  
  tmp<-InitAndRunEM(mExpression, times, J,r)
  
  Pij <-tmp$P
  m <- tmp$m
  vu <-tmp$mU
  
  indx <- max.col(Pij)
  grp <- unique(indx)
  grp <- order(grp) # added by Yaqun Wang, 1/10/2017
  
  TransColors <- rainbow(length(grp),alpha = 0.12)
  SolidColors <- rainbow(length(grp),alpha = 0.8)
  
  vY <- c(min(mExpression), max(mExpression), rep(0,length(a1) - 2))
  plot(a1, vY , type = "n",xlab="Time", ylab="Expression")
  TitleStr <- sprintf("Summary of Gene Expression( Clusters: %d )",length(grp))
  title(TitleStr, font.main=3)
  
  for( i in 1:length(grp) ){
    y18 <- as.numeric(m[grp[i],])
    par(lwd = 2)
    lines(a1, y18, col=SolidColors[i])
    par(lwd = 1)
  }
  
  for( i in 1:length(grp) ){
    indxTemp <- which(indx == grp[i])
    
    plot(a1, vY , type = "n",xlab="Time", ylab="Expression")
    TitleStr <- sprintf("Gene Expression of Cluster%d, Genes:%d",grp[i],length(indxTemp)) # Modified by Yaqun wang, i => grp[i]
    title(TitleStr, font.main=1,cex.main = 1)
    
    
    for(j in 1: length(indxTemp) ){
      lines(a1,Asdata[indxTemp[j], ],col=TransColors[i])
    }
    par(lwd = 3)
    y18 <- as.numeric(m[grp[i],])
    lines(a1, y18, col=SolidColors[i])
    par(lwd = 1)
  }
  par(opar)
  
  return(list(MeanExpression=m,LOPCoefficient =vu,Classifications=indx))
}
GeneClusterBIC <- function( mExpression, times, G , orderLOP){ 
  if(!requireNamespace("MASS",quietly = TRUE)){
    stop("MASS needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  r=orderLOP
  lgd=r
  if (r>10) {warning("range of r is 1 to 10")}
  if (r<0) {warning("r must be greater than 0")}
  if (length(times)!=ncol(mExpression)) {warning(" length of time vector is different to column length ")}
  if (!is.numeric(mExpression)) {warning("data is not numeric ")}
  
  EMResults <- matrix(0, length(G) , 4)
  i <- 1
  N_cluster <- c()
  for(j in G){
    cat("G: ", j, "\n")
    EMResults[i , 1] <- j
    tmp<-InitAndRunEM(mExpression, times, J=j,r=lgd)
    EMResults[i , 2:4] <- c(tmp$rho,tmp$Ssquare,tmp$Likelihood)
    Pij <-tmp$P####
    indx <- max.col(Pij)####
    N_cluster <- rbind(N_cluster,indx)####
    i <- i + 1
  }
  m <- ncol(mExpression)
  N <- nrow(mExpression)
  
  BIC <- -2 * EMResults[ , 4] + log(N) * (EMResults[ , 1] * (m + 1) + 1)
  BIC <- rbind(EMResults[ , 1], BIC)
  
  markSmallest=FALSE
  startN=1
  par(mar=(c(4,4,2,2)+.1),cex=1.3)
  Nk <- ncol(BIC)
  plot(BIC[1, startN : Nk], BIC[2, startN : Nk] / 1e+3 , type="b",pch = 5 , col= "blue",
       lty = "solid",main="",xlab = "# of cluster", ylab = "BIC", cex=0.9)
  ind <- which( BIC[2, ] == min(BIC[2, ]) )
  points(BIC[1, ind], BIC[2, ind] / 1e+3, 
         col = "deeppink", pch = 8 )
  mtext(expression("x"*10^3),line = 0.3,adj = 0,cex=1.2)
  
  optimal=which(BIC[2,]==min(BIC[2,])) 
  
  return(list(BIC=BIC,optimal=optimal,cluster=N_cluster))
  
}
GetPt <- function(t, times, r){
  tnum <- length(times)
  Pt<- rep(1,r+1)
  tt <- -1 + 2 * (t - times[1])/(times[tnum] - times[1])
  for(i in 1:r){
    Pt[i+1] <- eval(LgdP[i])
  }
  return(Pt)
}

GetYHat <- function(xx, LOPCoefficient, times, r){  
  yhat <- matrix(0,nrow(LOPCoefficient),length(xx))
  for(i in 1:length(xx)){
    Pt <- GetPt(xx[i],times, r)
    yhat[,i] <- Pt %*% t(LOPCoefficient)
  }
  return(yhat)
}

GeneClusterInterp <- function(LOPCoefficient, OriginalTime, outLen = 20){ 
  times=OriginalTime
  
  r <- length(LOPCoefficient[1,]) - 1
  newtime  <- seq(min(times),max(times), len=outLen)
  yHat <- GetYHat(newtime, LOPCoefficient, times, r)
  return(rbind(newtime,yHat))
}
ModifyMatrix <-function(LOPCoefficient,times,lowCut,upCut){
  
  cut=((max(times)-min(times))+1)*10
  output2=GeneClusterInterp(LOPCoefficient,times,cut)
  newoutput<-matrix(NA,nrow(LOPCoefficient),cut)
  for(i in 2:nrow(output2)){
    if((output2[i,1]> upCut )|(output2[i,1]< lowCut)){newoutput[i-1,]=output2[i,]}
    else{
      for(j in 1:(ncol(output2)-1)){
        if (((output2[i,j] < upCut )&(output2[i,j] > lowCut)) & ((output2[i,j+1] > upCut)|(output2[i,j+1] < lowCut))){
          newoutput[i-1,1:(cut-j)]=output2[i,(j+1):cut]
          break}
      }
    }
  }
  
  rownames(newoutput)=seq(1,nrow(LOPCoefficient))
  
  na_count<- apply(newoutput, 1, function(z) sum(is.na(z)))
  
  Modify=newoutput[!(na_count==cut),]
  
  Modify=na.omit(t(Modify))
  
  Modify=t(Modify)
  
  na_count=data.frame(na_count)
  
  return(list(na_count=na_count,Modify=Modify))
}

GeneClusterNet <- function(mExpression,times, orderLOP,alpha1=0.5,alpha2=0.05, realign=F, cutoff=c(lowCut=-0.35,upCut=0.2),NumberOfCluster=0,sLabels=NULL)
{ 
  if(!requireNamespace("G1DBN",quietly = TRUE)){
    stop("G1DBN needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  if(!requireNamespace("igraph",quietly = TRUE)){
    stop("igraph needed for this function to work. Please install it.",
         call.=FALSE)
  }
  
  
  
  data=mExpression
  
  J=NumberOfCluster
  
  if(NumberOfCluster==0){J=GeneClusterBIC(data, times, G = c(1:15), orderLOP)$optimal}
  
  vu=GeneCluster(data, times, J, orderLOP)$LOPCoefficient
  
  if(realign){Xn2=ModifyMatrix(vu,times,cutoff[1],cutoff[2])$Modify}
  
  else {Xn2=GeneClusterInterp(vu, times, outLen = (max(times)-min(times)+1))[-1,]}
  
  kk=(Xn2[,1]>0)
  
  color=rep("yellow",length(kk))
  for(i in 1:length(kk)){
    if(kk[i]=="FALSE"){color[i]="red"} 
  }
  
  Xn=t(Xn2)
  S1 <- DBNScoreStep1(Xn, method='ls')
  S2 <- DBNScoreStep2(S1$S1ls, data=Xn, method='ls', alpha1=alpha1)
  G2 <- BuildEdges(S2,threshold=alpha2,dec=FALSE)
  Step2InferredNet<- BuildNetwork(G2,1:nrow(Xn2))
  
  
  #   if(NumberOfCluster==0) {sLabel=seq(1:nrow(Xn2))
  #   }else {if(nrow(Xn2)==J){sLabel=sLabels}
  #    else {sLabel=seq(1:nrow(Xn2))
  #    }
  #   }
  
  if(NumberOfCluster==0) {sLabel=seq(1:nrow(Xn2))
  }else {if( (nrow(Xn2)==J) & (length(sLabels) !=0 )){sLabel=sLabels} # Nodified by Yaqun 1/9/2017
    else {sLabel=seq(1:nrow(Xn2))
    }
  }
  
  
  if(realign){sLabel=row.names(Xn2)
  
  rownames(Step2InferredNet$AdjMatrix)=row.names(Xn2)
  colnames(Step2InferredNet$AdjMatrix)=row.names(Xn2)
  
  rownames(Step2InferredNet$Score)=row.names(Xn2)
  colnames(Step2InferredNet$Score)=row.names(Xn2)
  
  }
  
  g <- graph.adjacency(t(Step2InferredNet$AdjMatrix), mode ='directed')
  g <- simplify(g, remove.multiple = T, remove.loops = T)
  #V(g)$color <- color
  V(g)$color <- 'deepskyblue3'
  E(g)$color <- 'deeppink'
  E(g)$arrow.width <- 0.5
  plot(g, layout=layout.circle, vertex.label=sLabel, vertex.label.dist=0,vertex.size=15)
  title(main="DBN Inferred network") 
  
  
  return(list(Score=Step2InferredNet$Score ,AdjMatrix=Step2InferredNet$AdjMatrix))
}


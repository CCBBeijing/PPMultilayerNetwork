

snp=read.csv( "S_SNP.csv",header = T,row.names = 1)#snp data
rmos=read.csv( "log10S_revise.csv",header = T,row.names = 1)#Data of Staphylococcus aureus monoculture
rmoe=read.csv( "log10ES_revise.csv",header = T,row.names = 1)#Data of Staphylococcus aureus co-culture
t <- c(0,0.5,1,1.5,2,4,6,8,10,12,16,20,24,30,36)#Time point of phenotypic data measurement

#################co_functional mapping
library(mvtnorm)

library(pbapply)
get_miu3 =function(B,t){B[1]/(1+exp((4*B[2]*(B[3]-t)/B[1])+2))}#Logistic  model
get_r2=function(B,t,y){
  sum((y-get_miu3(B,t))^2)
}
SAD1_get_matrix = function(par, times = t, options=list()) {   
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
  return(sigma); 
}#SAD1 model
get_u = function(A,t){get_miu3(A[1:3],t)-get_miu3(A[4:6],t)}#Mean vector model of phenotypic plasticity
get_sig= function(par,times = t, options=list()){SAD1_get_matrix(par = c(par[1],par[2]),times = t) + SAD1_get_matrix(par = c(par[3],par[4]),times = t) 
}#Variance matrix model

get_initial_par <- function(pheno,t){
  mean0 <- apply(pheno[,],2,mean)  
  c(max(mean0),
    (max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)]))),
    t[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]-mean0[which.max(((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))]/max((mean0[-1]-mean0[-length(mean0)])/(t[-1]-t[-length(t)])))
}#The initial parameters were calculated according to biological significance
get_initial_par2 <- function(pheno,t){
  
  mcurve =colMeans(pheno)
  dd <- optim(par=get_initial_par(pheno,t),get_r2,t=t,y=mcurve)
  return(dd$par)
  
}#The initial parameters were calculated according to biological significance

H0 = function(yt,t,par){
  miu=get_u(par[1:6],t)
  sigma=get_sig(par = par[7:10],times=t)
  L0 = c()
  L0 = sum(dmvnorm(t(yt),miu,sigma,log = TRUE))
  return(-L0)
}
H1 = function(yt,t,m0,m1,par){
  
  
  
  
  
  p1.0 <- yt[,c(m0)]
  p1.1 <- yt[,c(m1)]
  
  miu0 = get_u(par[1:6],t)
  miu1 = get_u(par[7:12],t)
  sigma = get_sig(par = par[13:16],times=t)
  
  L1.1 = sum(dmvnorm(t(p1.0),miu0,sigma,log = TRUE))
  L1.0=sum(dmvnorm(t(p1.1),miu1,sigma,log = TRUE))
  L1= L1.1+L1.0
  return(-L1)
}
inpar <- c(get_initial_par2(rmos,t),get_initial_par2(rmoe,t),0.5,0.02,0.5,0.02)

parl0 <- optim(par =inpar ,H0,yt = t(rmos)-t(rmoe),t = t )
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)
parl0 <- optim(par = parl0$par,H0,yt = t(rmos-rmoe),t = t)



optim_diff <- function(pheno_ck0,pheno_salt0,pheno_ck1,pheno_salt1,pheno_diff,t,m0,m1,xx){
  
  itime <- 100
  itimes <- 1
  
  par0 <-as.numeric(c(parl0$par[1:6]*xx,parl0$par[1:6]*xx,parl0$par[7:10]))
  repeat{
    a <- optim(par = par0 ,H1,yt = pheno_diff, t = t,m0=m0,m1=m1)
    
    b <- optim(a$par,H1,yt = pheno_diff, t = t,m0=m0,m1=m1)
    
    
    cat("Logistic_diff",itimes,b$value,'\n')
    
    itimes <- itimes + 1
    
    if(all( abs(a$value-b$value) < 1 )||itimes == itime){ #itimesԽ?߾???Խ??
      break
    }else{
      par0 <- b$par
    }
  }
  b
  
}

get_lr <- function(SNP,control,stress,xx){
  vnp=SNP
  m0 <- which( vnp == 0)
  m1 <- which( vnp == 1)
  pheno_ck <- control[,-1]
  pheno_salt <- stress[,-1]
  pheno_ck0 <- pheno_ck[,m0]
  pheno_ck1 <- pheno_ck[,m1]
  pheno_salt0 <- pheno_salt[,m0]
  pheno_salt1 <- pheno_salt[,m0]
  
  parl1 <- optim_diff(pheno_ck0,pheno_salt0,pheno_ck1,pheno_salt1,(control[,-1]-stress[,-1]),t,m0,m1,xx)
  LR_plastic3 <- c(-2*(-parl0$value+parl1$value),parl1$par)
}

LR_plastic<- apply(snp,1,get_lr,control=cbind(t,t(rmos)),stress=cbind(t,t(rmoe)),xx=1)






#########################################The genetic standard deviation was calculated according to the result parameters of functional mapping
get_VG <- function(jin2,jin5,marker_data,t,LR){
  
  diff_vg <- c() 
  
  
  
  
  
  for (a in 1:dim(marker_data)[1]) {
    
    AA <- as.numeric(which(marker_data[a,]==1))
    aa <- as.numeric(which(marker_data[a,]==0))
    
    NAA <- length(AA)
    Naa <- length(aa)
    
    p1 <- (NAA*2)/((NAA+Naa)*2) #A????Ƶ??
    p0 <- (Naa*2)/((NAA+Naa)*2) #a????Ƶ??
    
    mean_AA <- get_u(LR[2:7,a],t)
    mean_aa <- get_u(LR[8:13,a],t)
    AE <- (mean_AA - mean_aa)/2 
    
    Vg <- 2*p1*p0*(AE^2)  
    diff_vg <- rbind(diff_vg,Vg)
    cat(a,"finished","\n")
    
  }
  
  
  
  return(sqrt(diff_vg)) 
}

varh2tx <- get_VG(cbind(t,t(rmos)),cbind(t,t(rmoe)),snp,c(1:max(t)),LR_plastic)
row.names(varh2tx) <- c(1:nrow(snp))



##############################################functional cluster

ykmeans <- function(data,k,iter.max=10){
  
  data=rbind(colMeans(data),data)
  l <- nrow(data)
  mdist <- dist(data,method = "euclidean")
  mdist <- as.numeric(mdist[1:l-1])
  smdist <- order(mdist,decreasing = T)[1:k]
  
  
  center <- data[smdist,]
  ycl <- kmeans(data,centers = center,iter.max = iter.max)
  return(ycl)
}#The optimized kemeans method is used to get the initial parameters

requiredPackages = c("mvtnorm","reshape2","pbapply",'parallel','orthopolynom')
for(packages in requiredPackages){
  if(!require(packages,character.only = TRUE)) install.packages(packages)
  require(packages,character.only = TRUE)
}

age <- c(1:max(t))#times
get_init_par <- function(data,k){
  
  f3  <- function(y) {
    
    lop <- legendre.polynomials(n=13,normalized = F)
    lop_matrix <- as.matrix(as.data.frame(polynomial.values(polynomials = lop,x=scaleX(1:36,u=-1,v=1))))
    colnames(lop_matrix) <- paste0("L",0:13)
    lop_fit <- lm(y~lop_matrix[,2:14]) 
    return(as.numeric(lop_fit$coefficients))
  }
  
  
  
  
  
  #????k-means???ó?ʼ????--------------------------------------
  
  init_cluster <- ykmeans(data,k,iter.max = 100)
  pro <- table(init_cluster$cluster)/nrow(data) #????ģ?͸??ʸ???
  
  cuM <- init_cluster$centers
  cusd <- diag(cov(df))
  
  init_curve_para <-  t(apply(cuM, 1, f3))
  init_sd_para <- c(mean(cusd),0.4) #SAD1
  init_pro <- pro
  
  #????mclust
  #init_cluster <- mclust::densityMclust(df,k)
  #cuM <- t(init_cluster$parameters$mean)
  #init_pro <- init_cluster$parameters$pro
  #init_sd_para <- c(mean(cusd),0.5)
  
  return_object <- list(init_sd_para,init_curve_para,init_pro)
  names(return_object)<-c("init_sd_par","init_curve","init_pro")
  return(return_object)
}#get initial parameters

get_cluster <- function(data,k,input){
  
  requiredPackages = c("mvtnorm","reshape2")
  for(packages in requiredPackages){
    if(!require(packages,character.only = TRUE)) install.packages(packages)
    require(packages,character.only = TRUE)
  }
  Delta <- 100; iter <- 0; itermax <- 20;
  SAD1_get_matrix <- function(par,data){
    p <-  ncol(data)
    v2 <- par[1]
    phi <- par[2]
    tmp <- (1-phi^2)
    sigma <- array(dim=c(p,p))
    for(i in 1:p){
      sigma[i,i:p] <- phi^( c(i:p) - i ) * (1-phi^(2*i ))/tmp
      sigma[i:p,i] <- sigma[i,i:p]}
    sigma <- sigma*abs(v2)
    return(sigma)
  } 
  AR1_get_matrix <- function(par,data){
    n <- ncol(data)
    rho <- par[2]
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                      (1:n - 1))
    return(par[1]*rho^exponent)
  }
  Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
  {
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
  mle <- function(par,data,prob){
    par1 <- par[1:2]
    par2 <- matrix(par[-c(1:2)],nrow = k,ncol = 2)
    temp_S <- sapply(1:k, function(c) dmvnorm(data,
                                              Legendre.model(t=c(1:36),mu=par2[c,]),
                                              SAD1_get_matrix(par1,data))*prob[c] )
    LL <- sum(-log(rowSums(temp_S)))
    return(LL)
  }
  while ( Delta > 10 && iter <= itermax ) {
    # initiation
    if(iter == 0){
      init_sd_para <- input[[1]]
      init_curve_para <- input[[2]]
      pro <- input[[3]]
    }
    #E step, calculate the posterior probability
    old_par <- c(init_sd_para,init_curve_para)
    LL_mem <- mle(old_par,data,pro)
    mvn.c <- sapply(1:k, function(c) dmvnorm(data,
                                             Legendre.model(t=c(1:36),mu=init_curve_para[c,]),
                                             SAD1_get_matrix(init_sd_para,data))*pro[c] )
    omega <- mvn.c/rowSums(mvn.c)
    #M step, calculate parameters
    pro <- colSums(omega)/sum(omega)
    
    new_par <- try(optim(old_par, mle, data=data, prob=pro, method = "Nelder-Mead"))
    if ('try-error' %in% class(new_par))
      break
    L_Value <- new_par$value
    init_sd_para <- new_par$par[1:2]
    init_curve_para <- matrix(new_par$par[-c(1:2)],nrow = k)
    Delta <- abs(L_Value-LL_mem)
    if (Delta > 1200) break
    
    cat('\n',"iter=",iter,"LL=",L_Value,"Delta=",Delta,'\n')
    iter <- iter+1; LL_mem <- L_Value
    
  } 
  
  BIC <- 2*(L_Value)+log(nrow(data))*length(old_par)
  
  cluster <- apply(omega,1,which.max)
  #clustered_df <- data.frame(cbind(row.names(data),df_fitted,cluster))
  clustered_df <- data.frame(cbind(row.names(data),df,cluster))
  colnames(clustered_df) <- c("row.names(data)",age,"cluster")
  long_df <- melt(clustered_df,id.vars=c("row.names(data)","cluster"))
  long_df[,4] <- as.numeric(long_df[,4])
  long_df[,3] <- as.numeric(long_df[,3])
  colnames(long_df) <- c("gene","cluster","time","fpkm")
  
  clustered_df <- clustered_df[,-1]
  return_object <- list(init_sd_para,init_curve_para,pro,LL_mem,BIC,clustered_df)
  names(return_object)<-c("sd_par", "curve_par", "pro", "LL", "BIC", "clustered_data")
  
  return(return_object)
  
}#get cluster results

get_kcl <- function(data,clust){
  
  k=length( table(clust))
  kcl <- list()
  for (i in 1:k) { posi<- which(clust==i)
  kcl[[i]] <- data[posi,]
  }
  
  return(kcl)
  
  
  
}#Get the cluster data results 
df <- varh2tx
input=get_init_par(df,9)

a <- get_cluster(df,9,input) 
table(as.numeric(a$clustered_data[,37]))
fitall <-as.numeric( a$clustered_data[,37])


fclall <- get_kcl(varh2tx,fitall)#get the modules



#########################################The variables are selected by Lasso, and the fourth-order Runge Kutta is used to solve the differential equation to construct the network


library(splines)
library(orthogonalsplinebasis)
library(MASS)
library(grplasso)
library(parallel)

LMall <- function(NX,nt,nstep=30,order){
  
  stp <- (max(nt)-min(nt))/nstep
  res <- c()
  for(j in 1:nstep){
    
    tg1 <- Legendre.model11((j-1)*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg2 <- Legendre.model11(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg3 <- Legendre.model11(j*stp/2+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tg4 <- Legendre.model11(j*stp+1,np.order=order-1,tmin=min(nt), tmax=max(nt))
    tmp1 <- rbind(tg1,tg2,tg3,tg4)
    res <- rbind(res,tmp1)
  }
  res
}

fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){
  
  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL,self=self)
  sumode <-rowSums(odes)
  sumode[1] <- sumode[1]+NG[1,self]
  sum((NG[,self]-(sumode))^2
      +0*sum(para^2)
  )
}
ode.optim <- function(y.c,connect,effect,LL,nstep,order,times){
  self <- y.c 
  indexx <- which(connect[y.c,]==1)
  para <- rep(0.0001,length(indexx)*(order-1))
  res <- optim(para,fitPKM,NG=(effect),self=y.c,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,
               LL=LL,method="BFGS",control=list(maxit=200,trace=T)
  )
  cat("Gene=",y.c," ",res$value,"\n")
  A <- ode.sovle.ind(NG=(effect),res$par,nconnect=connect[y.c,],nt=times,order=order,nstep=nstep,LL=LL,self=self)
  return(A)
}


ode.sovle.ind <- function(NG,fitpar,nconnect,nt,order,nstep,LL,self){
  stp <- (max(nt)-min(nt))/nstep
  index <- which(nconnect==1)
  
  ind.par <- matrix(fitpar[1:(length(index)*(order-1))],ncol=order-1,byrow=T)
  allrep <- matrix(rep(0,length(index)),nrow=1)
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
  if(length(which(allrep[,no_name]<0))==0){
    return(allrep)
  } else{return(allrep*500)}
  
}
smL <- function(times,para,DS1){
  
  sum((DS1-Legendre.model(t=times,mu=para))^2)
}

Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
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

Legendre.model11 <- function(t, np.order,tmin = NULL, tmax = NULL)
{
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
    L[6] <-(1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 * ti) 
  if (np.order >= 7)
    L[7] <- (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)
  if (np.order>=8)
    L[8] <-  (1/128) * (6435 * 8 * ti ^ 7 - 12012 * 6 * ti ^ 5 + 6930 * 4 * ti ^ 3 - 1260 * 2 * ti)
  if (np.order>=9)
    L[9] <-  (1/128) * (12155 * 9 * ti ^ 8 - 25740 * 7 * ti ^ 6 + 18018 * 5 * ti ^ 4 - 4620 * 3 * ti ^ 2 + 315)
  if (np.order>=10)
    L[10] <-  (1/256) * (46189 * 10 * ti ^ 9 - 109395 * 8 * ti ^ 7 + 90090 * 6 * ti ^ 5 - 30030 * 4 * ti ^ 3 + 3465 * 2 * ti)
  if (np.order>=11)
    L[11] <-  (1/512) * (176358 * 11 * ti ^ 10 - 461890 * 9 * ti ^ 8 + 437580 * 7 * ti ^ 6 - 180180 * 5 * ti ^ 4 + 30030 * 3 * ti^ 2 - 1386)
  if (np.order>=12)
    L[12] <-  660.1943 * 12 * ti ^ 11 -1894.471 * 10 * ti ^ 9 + 2029.79 * 8 * ti ^ 7  - 997.0898 * 6 * ti ^ 5 + 219.9463 * 4 * ti^ 3 - 17.5957* 2 * ti
  if (np.order>=13)
    L[13] <-  1269.604 * 13 * ti ^ 12 -3961.166 * 11 * ti ^ 10 + 4736.177 * 9 * ti ^ 8  - 2706.387 * 7 * ti ^ 6 + 747.8174 * 5 * ti^ 4 - 87.97852* 3 * ti^ 2 + 2.932617
  if (np.order>=14)
    L[14] <-  2448.523 * 14 * ti ^ 13 -8252.429 * 12 * ti ^ 11 + 10893.21 * 10 * ti ^ 9  - 7104.265 * 8 * ti ^ 7 + 2368.088 * 6 * ti^ 5 - 373.9087* 4 * ti^ 3 + 21.99463* 2 * ti 
  if (np.order>=15)
    L[15] <-  4733.811 * 15 * ti ^ 14 -17139.66 * 13 * ti ^ 12 + 24757.29 * 11 * ti ^ 10  - 18155.34 * 9 * ti ^ 8 + 7104.265 * 7 * ti^ 6 - 1420.853* 5 * ti^ 4 + 124.6362* 3 * ti ^ 2 - 3.14209
  if (np.order>=16)
    L[16] <-  9171.759 * 16 * ti ^ 15 -35503.58 * 14 * ti ^ 13 + 55703.9 * 12 * ti ^ 11  - 45388.36 * 10 * ti ^ 9 + 20424.76 * 8 * ti^ 7 - 4972.986* 6 * ti^ 5 + 592.0221* 4 * ti ^ 3 - 26.70776* 2 * ti 
  if (np.order>=17)
    L[17] <-  17804 * 17 * ti ^ 16 -73374.07 * 15 * ti ^ 14 + 124262.5 * 13 * ti ^ 12  - 111407.8 * 11 * ti ^ 10 + 56735.45 * 9 * ti^ 8 - 16339.81* 7 * ti^ 6 + 2486.493* 5 * ti ^ 4 - 169.1492* 3 * ti ^ 2 +3.33847
  if (np.order>=18)
    L[18] <-  34618.89 * 18 * ti ^ 17 -151334 * 16 * ti ^ 15 + 275152.8 * 14 * ti ^ 13  - 269235.5 * 12 * ti ^ 11 + 153185.7 * 10 * ti^ 9 - 51061.91* 8 * ti^ 7 + 9531.556* 6 * ti ^ 5 - 888.0331* 4 * ti ^ 3 +31.71547* 2 * ti 
  if (np.order>=19)
    L[19] <-  67415.74 * 19 * ti ^ 18 -311570 * 17 * ti ^ 16 + 605336.1 * 15 * ti ^ 14  - 642023.1 * 13 * ti ^ 12 + 403853.3 * 11 * ti^ 10 - 153185.7* 9 * ti^ 8 + 34041.27* 7 * ti ^ 6 - 4084.952* 5 * ti ^ 4 +222.0083* 3 * ti ^ 2 -3.523941
  if (np.order>=20)
    L[20] <-  131460.7 * 20 * ti ^ 19 -640449.5 * 18 * ti ^ 17 + 1324173 * 16 * ti ^ 15  - 1513340 * 14 * ti ^ 13 + 1043288 * 12 * ti^ 11 - 444238.6* 10 * ti^ 9 + 114889.3* 8 * ti ^ 7 - 17020.64* 6 * ti ^ 5 +1276.548* 4 * ti ^ 3 -37.00138* 2 * ti 
  if (np.order>=21)
    L[21] <-  256661.4 * 21 * ti ^ 20 -1314607 * 19 * ti ^ 18 + 2882023 * 17 * ti ^ 16  - 3531127 * 15 * ti ^ 14 + 2648345 * 13 * ti^ 12 - 1251945* 11 * ti^ 10 + 370198.8* 9 * ti ^ 8 - 65651.02* 7 * ti ^ 6 +6382.738* 5 * ti ^ 4 -283.6773* 3 * ti^ 2 + 3.700138 
  if (np.order>=22)
    L[22] <-  501656.3 * 22 * ti ^ 21 -2694944 * 20 * ti ^ 19 + 6244383 * 18 * ti ^ 17  - 8165732 * 16 * ti ^ 15 + 6620863 * 14 * ti^ 13 - 3442849* 12 * ti^ 11 + 1147616* 10 * ti ^ 9 - 237985* 8 * ti ^ 7 +28722.32* 6 * ti ^ 5 -1772.983* 4 * ti^ 3 + 42.55159* 2 * ti
  if (np.order>=23)
    L[23] <-  981501.4 * 23 * ti ^ 22 -5518219 * 21 * ti ^ 20 + 13474721 * 19 * ti ^ 18  - 18733149 * 17 * ti ^ 16 + 16331463 * 15 * ti^ 14 - 9269209* 13 * ti^ 12 + 3442849* 11 * ti ^ 10 - 819725.9* 9 * ti ^ 8 +118992.5* 7 * ti ^ 6 -9574.107* 5 * ti^ 4 + 354.5966* 3 * ti^ 2 -3.868326
  if (np.order>=24)
    L[24] <-  1922107 * 24 * ti ^ 23 -11287266 * 22 * ti ^ 21 + 28970650 * 20 * ti ^ 19  - 42669950 * 18 * ti ^ 17 + 39807941 * 16 * ti^ 15 - 24497195* 14 * ti^ 13 + 10041643* 12 * ti ^ 11 - 2705096* 10 * ti ^ 9 +461095.8* 8 * ti ^ 7 -46274.85* 6 * ti^ 5 + 2393.527* 4 * ti^ 3 -48.35408* 2 * ti
  if (np.order>=25)
    L[25] <-  3767330 * 25 * ti ^ 24 -23065284 * 23 * ti ^ 22 + 62079965 * 21 * ti ^ 20  - 96568835 * 19 * ti ^ 18 + 96007388 * 17 * ti^ 16 - 63692706* 15 * ti^ 14 + 28580061* 13 * ti ^ 12 - 8607122* 11 * ti ^ 10 +1690685* 9 * ti ^ 8 -204931.5* 7 * ti^ 6 + 13882.46* 5 * ti^ 4 -435.1867* 3 * ti ^ 2 + 4.029506
  if (np.order>=26)
    L[26] <-  7389762 * 26 * ti ^ 25 -47091621 * 24 * ti ^ 23 + 132625380 * 22 * ti ^ 21  - 217279879 * 20 * ti ^ 19 + 229350983 * 18 * ti^ 17 - 163212560* 16 * ti^ 15 + 79615883* 14 * ti ^ 13 - 26538628* 12 * ti ^ 11 +5917397* 10 * ti ^ 9 -845342.4* 8 * ti^ 7 + 71726.02* 6 * ti^ 5 -3155.104* 4 * ti ^ 3 + 54.39834* 2 * ti
  
  return(L);
}

optim.parallel <- function(connect,effect,n.cores,proc,order,times,nstep){
  
  diag(connect) <- 1
  nt1 <- min(times)
  nt2 <- max(times)
  
  LL <- LMall(NX=1,nt=seq(nt1,nt2,(nt2-nt1)/nstep),nstep=nstep,order=order)
  
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

get_cl_meancurve <- function(data,clust){
  
  
  library(ggplot2)
  library(reshape2)
  library(scales)
  k=length( table(clust))
  fclmean <- c()
  yansee <- hue_pal()(k)
  for (i in 1:k) {
    
    fclm <- colMeans(data[[i]])
    
    
    fclmean <- rbind(fclmean,fclm)
  }
  
  rownames(fclmean) <- c(1:k)
  
 
  
  
  
  
  marker <- t(fclmean)
  n <- dim(marker)[1]
  marker_list <- list()
  name <- colnames(marker)
  
  
  library(glmnet)
  
  for (col in 1:length(marker[1,])) {
    #tim <- proc.time()
    
    
    m <- marker[,col]
    M <- marker[,-col]
    
    x_matrix <- M
    ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 10,alpha = 0)
    best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
    fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 100,alpha = 1,penalty.factor = 1/ abs(best_ridge_coef),keep = TRUE)
    best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
    marker_list_one <- list()
    marker_list_one[[1]] <- name[col]#??һ???б???ֱ??qtl??????
    marker_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]#?ڶ????б??Ǽ???qtl??????
    marker_list_one[[3]] <- best_alasso_coef1@x[-1]#???????б??Ǳ?��ѡ??ϵ??
    
    marker_list[[col]] <- marker_list_one}
  
  for (i in 1:length(marker_list)) {    
    if(length(marker_list[[i]][[3]])>1){
      
      
      lsr <- quantile(abs(marker_list[[i]][[3]]))[2]
      marker_list[[i]][[2]] <-  marker_list[[i]][[2]][which(abs(marker_list[[i]][[3]])>lsr)]
      marker_list[[i]][[3]] <-   marker_list[[i]][[3]][which(abs(marker_list[[i]][[3]])>lsr)]}
    
  } 
  
  marker_list
  
  
  lasso_1 <- matrix(0,nrow=dim(marker)[2],ncol=dim(marker)[2])
  colnames(lasso_1) <- colnames(marker)
  rownames(lasso_1) <- colnames(marker)
  for(c in 1:k){
    
    for(d in 1:length(marker_list[[c]][[2]])){
      lasso_1[c,which(name== marker_list[[c]][[2]][d])]=0.5
      
    }
  }
  
  lasso=lasso_1
  for (i in 1:k) { lasso[i,i] <- 1
  
  }
  library(corrplot)
  
  cor_plot <-  corrplot( lasso,shade.col = NA,tl.col = "black",tl.srt = 45)
  for(c in 1:k){
    
    for(d in 1:length(marker_list[[c]][[2]])){
      lasso_1[c,which(name== marker_list[[c]][[2]][d])]=1
      
    }
  }
  F2014_H.odee <- optim.parallel(connect=lasso_1,effect=t(fclmean),
                                 n.cores=1,proc=ode.optim,order=13,times=c(1:36),nstep=c(35))
  fc_res <- list(
    
  )
  fc_res[[1]] <-  F2014_H.odee
  fc_res[[2]] <-   cor_plot
  fc_res[[3]] <-  marker_list
 
  
  names(fc_res)<-c("ode_res","lasso_res_plot","lasso_res"
  )
  return(fc_res) 
  
  
  
  
  
  
  
  
}

fclall_res <- get_cl_meancurve(fclall,fitall)#get modules network




#################get submodules and snp network

get_final_res <- function(data){
  
  marker <- t(data)
  n <- dim(marker)[1]
  marker_list <- list()
  name <- colnames(marker)
  
  
  library(glmnet)
  
  for (col in 1:length(marker[1,])) {
    #tim <- proc.time()
    
    
    m <- marker[,col]
    M <- marker[,-col]
    
    x_matrix <- M
    ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 10,alpha = 0)
    best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
    fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 100,alpha = 1,penalty.factor = 1 / abs(best_ridge_coef),keep = TRUE)
    best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
    marker_list_one <- list()
    marker_list_one[[1]] <- name[col]#??һ???б???ֱ??qtl??????
    marker_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]#?ڶ????б??Ǽ???qtl??????
    marker_list_one[[3]] <- best_alasso_coef1@x[-1]#???????б??Ǳ?��ѡ??ϵ??
    marker_list[[col]] <- marker_list_one
    
  }
  for (i in 1:length(marker_list)) {    
    if(length(marker_list[[i]][[3]])>1){
      
      
      lsr <- quantile(abs(marker_list[[i]][[3]]))[2]
      marker_list[[i]][[2]] <-  marker_list[[i]][[2]][which(abs(marker_list[[i]][[3]])>lsr)]
      marker_list[[i]][[3]] <-   marker_list[[i]][[3]][which(abs(marker_list[[i]][[3]])>lsr)]}
    
  } 
  marker_list
  lasso_1 <- matrix(0,nrow=dim(marker)[2],ncol=dim(marker)[2])
  colnames(lasso_1) <- colnames(marker)
  rownames(lasso_1) <- colnames(marker)
  k=length(marker_list)
  for(c in 1:length(marker_list)){
    
    for(d in 1:length(marker_list[[c]][[2]])){
      lasso_1[c,which(name== marker_list[[c]][[2]][d])]=0.5
      
    }
  }
  
  lasso=lasso_1
  for (i in 1:k) { lasso[i,i] <- 1
  
  }
  library(corrplot)
  
  cor_plot <-  corrplot( lasso,shade.col = NA,tl.col = "black",tl.srt = 45)
  for(c in 1:k){
    
    for(d in 1:length(marker_list[[c]][[2]])){
      lasso_1[c,which(name== marker_list[[c]][[2]][d])]=1
      
    }
  }
  #row.names(data)=as.numeric(1:nrow(data))
  F2014_H.odee <- optim.parallel(connect=lasso_1,effect=t(data),
                                 n.cores=1,proc=ode.optim,order=13,times=c(1:36),nstep=c(35))
  fc_res <- list(
    
  )
  fc_res[[1]] <-    cor_plot
  fc_res[[2]] <-  marker_list


  fc_res[[3]] <- F2014_H.odee
  names(fc_res)<-c(
    "lasso_res_plot","lasso_res","ode_res")
  return(fc_res) 
  
}

df <- fclall[[6]]
input=get_init_par(df,57)
f6 <- get_cluster(df,57,input) 
fit6 <- as.numeric(f6$clustered_data[,37])
fcl6 <- get_kcl(fclall[[6]],fit6)
fcl6_36res <- get_final_res(fcl6[[36]])



df <- fclall[[1]]
input=get_init_par(df,14)
f1 <- get_cluster(df,14,input) 
fit1 <- as.numeric(f1$clustered_data[,37])
fcl1 <- get_kcl(fclall[[1]],fit1)
fcl1_4res <- get_final_res(fcl1[[4]])

fitPKM <- function(para,NG,self,nconnect,nt,order,nstep,LL){
  
  odes <- ode.sovle.ind(NG,para,nconnect,nt,order,nstep,LL,self=self)
  sumode <-rowSums(odes)
  sumode[1] <- sumode[1]+NG[1,self]
  sum((NG[,self]-(sumode))^2
      +0.01*sum(para^2)
  )
}
df <- fclall[[3]]
input=get_init_par(df,11)
f3 <- get_cluster(df,11,input) 
fit3 <- as.numeric(f3$clustered_data[,37])
fit3[which(fit3==11)]=5
fcl3 <- get_kcl(fclall[[3]],fit3)

fcl3_5res <- get_final_res(fcl3[[5]])


fcl3_3res <- get_final_res(fcl3[[3]])

fcl3_7res <- get_final_res(fcl3[[7]])




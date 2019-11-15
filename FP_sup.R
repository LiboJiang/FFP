
mle_curve <- function(pheno,times){
  
  pl <- log(pheno)
  p <- as.numeric(colMeans(pl,na.rm=T))
  pp1 <- p[1:23];pp2 <- p[24:31];
  cpara1 <- AIC_L(pheno=pp1,times=1:23)
  cpara2 <- AIC_L(pheno=pp2,times=24:31)
  spara <- c(0.5,sd(pl,na.rm=T))
  para <- c(spara,cpara1,cpara2)
  
  mle_fun <- function(para,y,times){
    
    sp <- para[1:2]
    vp <- para[-c(1:2)]
    vp1 <- vp[1:4]
    mu1 <- Legendre.model(times[1:23],vp1)
    vp2 <- vp[5:7]
    mu2 <- Legendre.model(times[24:31],vp2)
    cov <- AR1.get_mat(sp,times)
    pv <- sum(-dmvnorm(y,c(mu1,mu2),cov,log=T),na.rm=T)
    return(pv)
  }
  
  L0 <- optim(para,mle_fun,y=pl,times=times,method="BFGS",control=list(trace=F,maxit=20000))
  return(L0)
  
}





mle_H1 <- function(dat,times=1:31){
  
  
  marker <- dat$genotype
  pheno <- dat$phenotype
  y <- log(pheno)
  
  r0 <- mle_curve(pheno,times)
  
  nm <- dim(marker)[1]
  
  npara <- c(r0$par,r0$par[-c(1:2)])
  res <-c()
  for(i in 1:nm){
    snp <- as.numeric(marker[i,])
    snpi <- as.numeric(names(table(snp)))
    LH1 <- optim(npara,mle_g,y=y,times=times,snp=snp,snpi=snpi,method="BFGS",control=list(trace=F,maxit=20000))
    LR <- 2*(r0$value-LH1$value)
    pvalue <- pchisq(LR,7,lower.tail = F)
    rr <- c(LR,pvalue,LH1$value,LH1$par)
    cat("snp=",i,"rr=",rr,"\n")
    res <- rbind(res,rr)
  }
  return(res)
  
}


mle_g <- function(npara,y,times,snp,snpi){
  
  snpara <- npara[1:2]
  vnpara <- npara[-c(1:2)]
  cov1 <- AR1.get_mat(snpara,times)
  
  L1 <- 0
  for(i in 1:length(snpi)){
    index <- which(snp==snpi[i])
    yy <- y[index,]
    vp1 <- vnpara[(7*(i-1)+1):(7*i)]
    vp11 <- vp1[1:4];vp12 <- vp1[5:7]
    mu11 <- Legendre.model(times[1:23],vp11)
    mu12 <- Legendre.model(times[24:31],vp12)
    pv1 <- sum(-dmvnorm(yy,c(mu1,mu2),cov1,log=T),na.rm=T)
    L1 <- L1 + pv1
  }
  return(L1)
}













AR1.get_mat <- function(par0, times, traits=1, options=list())
{
  par<-par0;
  if (class(par0)=="list")
    par <- unlist(par0);
  
  t_len <- length( times );
  
  Ar.1 <- array(0, dim=c(t_len*traits,t_len*traits));
  for (i0 in 1:traits)
    for (i1 in 1:traits)
    {
      if (i0==i1)
        for (k0 in 1:t_len)
          for (k1 in 1:t_len)
          {
            Ar.1[(i0-1)*t_len+k0,(i1-1)*t_len+k1] <- par[i0*2]^2 * par[i0*2-1]^abs( k0 - k1 );
          }
    }
  
  return(Ar.1);
}

AIC_L <- function(pheno,times){
  
  od <- 12
  val <- c()
  for(nn in 2:od){
    
    para1 <- rep(0.5,nn)
    LL <- optim(para1,smL,DS1=pheno,times=times,method="BFGS")
    
    val <- c(val,2*nn+LL$value)
  }
  index <- c(2:od)[(which(val==min(val)))]
  para1 <- rep(0.5,index+1)
  LL <- optim(para1,smL,DS1=pheno,times=times,method="BFGS")
  return(LL$par)
}



smL <- function(times,para,DS1){
  
  sum((DS1-Legendre.model(t=times,mu=para))^2)
}


Legendre.model <- function( t, mu, tmin=NULL, tmax=NULL )
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



Genetic_eff <- function(dat,ret){
  
  ug <- unique(ret[,1])
  dt <- seq(1,23,0.1)
  dt1 <- seq(1,13,0.1)
  dt2 <- seq(14,23,0.1)
  dt3 <- seq(24,31,0.1)
  
  g_e1 <- c();g_e2 <- c();g_e3 <- c()
  for(i in ug){
    ni <- which(ret[,1]==i)[1]
    gpar <- ret[ni,-c(1:5)]
    g1par <- gpar[1:7];g2par <- gpar[8:14]
    g_e1 <- cbind(g_e1,(Legendre.model(dt,g1par[1:4])-Legendre.model(dt,g2par[1:4]))[1:121])
    g_e2 <- cbind(g_e2,(Legendre.model(dt,g1par[1:4])-Legendre.model(dt,g2par[1:4]))[122:221])
    g_e3 <- cbind(g_e3,Legendre.model(dt3,g1par[5:7])-Legendre.model(dt3,g2par[5:7]))
  }
  
  QQ <- c()
  for(i in 1:length(ug)){
    QQ <- c(QQ,rep(paste0("Q",i),length(which(ret[,1]==ug[i]))))
  }
  ndat <- cbind(QL=QQ,dat$info)
  
  
  return(list(dt1=dt1,dt2=dt2,dt3=dt3,ge1=g_e1,ge2=g_e2,ge3=g_e3,info=ndat))
}

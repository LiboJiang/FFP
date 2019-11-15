











load("dat.RData")
source("FP_sup.R")


H0 <- mle_curve(pheno=dat$phenotype,times=1:31)




ret <- mle_H1(dat,times=1:31)


GE <- Genetic_eff(dat,ret)


source("../tmp-file/WQ/WQ_sup.R")



require(orthogonalsplinebasis)
require(grplasso)
library(parallel)

#s1
nt1 <- seq(min(GE$dt1),max(GE$dt1),length=30)
stage2_dt1 <- smooth.optim(times=GE$dt1,para=rep(.1,5),y=t(GE$ge1),nt=nt1)

stage3_dt1 <- varsel(X=t(stage2_dt1$smooth.d),Y=t(stage2_dt1$dsmooth.d),tt=nt1)


d1.odee <- optim.parallel(connect=stage3_dt1$connect,effect=t(stage2_dt1$smooth.d),
                           n.cores=4,proc=ode.optim,order=6,times=nt1,nstep=29)

d1_res <- interType(con=stage3_dt1$connect,alle=d1.odee,sme=stage2_dt1$smooth.d)

aaad1 <- regasso(connect1=stage3_dt1$connect,gene=stage2_dt1$smooth.d,interaction=d1.odee)

#s2


nt2 <- seq(min(GE$dt2),max(GE$dt2),length=30)
stage2_dt2 <- smooth.optim(times=GE$dt2,para=rep(.1,5),y=t(GE$ge2),nt=nt2)

stage3_dt2 <- varsel(X=t(stage2_dt2$smooth.d),Y=t(stage2_dt2$dsmooth.d),tt=nt2)


d2.odee <- optim.parallel(connect=stage3_dt2$connect,effect=t(stage2_dt2$smooth.d),
                          n.cores=4,proc=ode.optim,order=6,times=nt2,nstep=29)

d2_res <- interType(con=stage3_dt2$connect,alle=d2.odee,sme=stage2_dt2$smooth.d)

aaad2 <- regasso(connect1=stage3_dt2$connect,gene=stage2_dt2$smooth.d,interaction=d2.odee)



#s3


nt3 <- seq(min(GE$dt3),max(GE$dt3),length=30)
stage2_dt3 <- smooth.optim(times=GE$dt3,para=rep(.1,5),y=t(GE$ge3),nt=nt3)

stage3_dt3 <- varsel(X=t(stage2_dt3$smooth.d),Y=t(stage2_dt3$dsmooth.d),tt=nt3)


d3.odee <- optim.parallel(connect=stage3_dt3$connect,effect=t(stage2_dt3$smooth.d),
                          n.cores=4,proc=ode.optim,order=6,times=nt3,nstep=29)

d3_res <- interType(con=stage3_dt3$connect,alle=d3.odee,sme=stage2_dt3$smooth.d)

aaad3 <- regasso(connect1=stage3_dt3$connect,gene=stage2_dt3$smooth.d,interaction=d3.odee)


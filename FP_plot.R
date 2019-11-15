


p1 <- as.numeric(colMeans(log(dat$phenotype),na.rm=T))

plot(1:31,p1)


para1 <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5)

LL <- optim(para1,smL,DS1=p1[24:31],times=24:31,method="BFGS")

LL


plot(24:31,p1[24:31],ylim=c(3,7))
lines(seq(24,31,0.1),Legendre.model(t=seq(24,31,0.1),mu=LL$par))


lines(seq(1,31,0.1),Legendre.model(t=seq(1,31,0.1),mu=L0$par[-c(1:2)]),col="red")


H0_plot <- function(dat,par1){
  
  np <- dim(dat$phenotype)[1]
  dt1 <- seq(1,23,0.1)
  dt2 <- seq(24,31,0.1)
  pdf("Figure1.pdf",width=6,height=4)
  par(mar=c(4,4,1,1),oma=c(0.5,0.5,0.5,0.5))
  plot(NA,NA,xlim=c(1,31),ylim=c(0,8),xlab="Time (day)",ylab="Transpiration rate",cex.lab=1.5,mgp=c(2.1,0.6,0))
  for(i in 1:np){
    lines(1:23,log(dat$phenotype[i,1:23]),col="#CFCFCF50")
    lines(24:31,log(dat$phenotype[i,24:31]),col="#CFCFCF50")
  }
  rect(1,-10,13,13,100,border=NA,col="#EE636350")
  rect(13,-10,23,100,border=NA,col="#3CB37150")
  rect(24,-10,31,100,border=NA,col="#63B8FF50")
  #points(1:31,colMeans(log(dat$phenotype),na.rm=T),col="red")
  lines(dt1,Legendre.model(dt1,par1$par[3:6]),col="#EE6363",lwd=2)
  lines(dt2,Legendre.model(dt2,par1$par[7:9]),col="#EE6363",lwd=2)
  text(7,7.5,"Normal",cex=1.2,col="#EE6363",font=2)
  text(18,7.5,"Drought",cex=1.2,col="#3CB371",font=2)
  text(27.5,7.5,"Recover",cex=1.2,col="#63B8FF",font=2)
  Arrows(1.8,0.5,12.2,0.5,lwd=2,col="#EE6363",arr.type="triangle",code=3)
  text(7,1,"12 days")
  Arrows(13.8,0.5,22.2,0.5,lwd=2,col="#EE6363",arr.type="triangle",code=3)
  text(18,1,"10 days")
  Arrows(24.8,0.5,30.2,0.5,lwd=2,col="#EE6363",arr.type="triangle",code=3)
  text(27.5,1,"7 days")
  dev.off()
  
}


which(ret==max(ret[,1]))




H1_plot <- function(dat,ret,index){
  
  gpar <- ret[index,-c(1:5)]
  g1par <- gpar[1:7];g2par <- gpar[8:14]
  np <- dim(dat$phenotype)[1]
  snp <- as.numeric(dat$genotype[index,])
  snp1 <- which(snp==0);snp2 <- which(snp==1);
  dt1 <- seq(1,23,0.1)
  dt2 <- seq(24,31,0.1)
  pdf("Figure3.pdf",width=6,height=4)
  par(mar=c(4,4,1,1),oma=c(0.5,0.5,0.5,0.5))
  plot(NA,NA,xlim=c(1,33),ylim=c(0,8.2),xlab="Time (Day)",ylab="Transpiration rate",cex.lab=1.5,mgp=c(2.1,0.6,0))
  for(i in 1:length(snp1)){
    lines(1:23,log(dat$phenotype[snp1[i],1:23]),col="#C1FFC130")
    lines(24:31,log(dat$phenotype[snp1[i],24:31]),col="#C1FFC130")
  }
  for(i in 1:length(snp2)){
    lines(1:23,log(dat$phenotype[snp2[i],1:23]),col="#1E90FF30")
    lines(24:31,log(dat$phenotype[snp2[i],24:31]),col="#1E90FF30")
  }
  #points(1:31,colMeans(log(dat$phenotype),na.rm=T),col="red")
  lines(dt1,Legendre.model(dt1,g1par[1:4]),col="green",lwd=2)
  lines(dt2,Legendre.model(dt2,g1par[5:7]),col="green",lwd=2)
  lines(dt1,Legendre.model(dt1,g2par[1:4]),col="blue",lwd=2)
  lines(dt2,Legendre.model(dt2,g2par[5:7]),col="blue",lwd=2)
  text(32.5,6,expression(italic(QQ)),cex=1.0,col="green")
  text(32.5,5,expression(italic(qq)),cex=1.0,col="blue")
  rect(1,-10,13,13,100,border=NA,col="#EE636350")
  rect(13,-10,23,100,border=NA,col="#3CB37150")
  rect(24,-1,31,100,border=NA,col="#63B8FF50")
  #text(7,8.1,"Chr11:97.3-103.4 cM",cex=1.0,col="black")
  mtext("Chr11:97.3-103.4 cM",3,cex=1.0,col="black")
  text(7,7.5,"Normal",cex=1.2,col="#EE6363",font=2)
  text(18,7.5,"Drought",cex=1.2,col="#3CB371",font=2)
  text(27.5,7.5,"Recover",cex=1.2,col="#63B8FF",font=2)
  Arrows(1.8,0.5,12.2,0.5,lwd=2,col="#EE6363",arr.type="triangle",code=3)
  text(7,1,"12 days")
  Arrows(13.8,0.5,22.2,0.5,lwd=2,col="#EE6363",arr.type="triangle",code=3)
  text(18,1,"10 days")
  Arrows(24.8,0.5,30.2,0.5,lwd=2,col="#EE6363",arr.type="triangle",code=3)
  text(27.5,1,"7 days")
  dev.off()
  
}




man_plot <- function(dat,ret){
  
  pvs <- -log10(ret[,2])
  allpos <- as.numeric(dat$info[,3])
  chr_len <- c()
  for(chr in 1:12){
    ids <- which(dat$info[,2] == chr)
    chr_len <- c( chr_len, max(allpos[ids]) )
  }

  leiji_len <- cumsum( c(0,chr_len))
  position <- rep(0, length(pvs))
  for(chr in 1:12){
    ids <- which(dat$info[,2] == chr);
    position[ids] <- allpos[ids]+ leiji_len[chr];
  }
  
  pdf("Figure2.pdf",width=6,height=4)
  par(mar=c(3,3,1,1),oma=c(0.5,0.5,0.5,0.5))
  
  threshold <- -log10(0.01/length(pvs))
  cat("threshold:", threshold, "\n");
  zuigao <- NULL;
  if( threshold >= max(pvs,na.rm=TRUE) ){ zuigao <- threshold; }else{ zuigao <- max(pvs,na.rm=TRUE); }
  si <- which(pvs>threshold)
  sigsnp <- cbind(dat$info[si,],pvalue=ret[si,2])
  write.csv(sigsnp,file="sigsnp.csv")
  plot(NA, NA, type="l", xlim=c(0, max(position)), ylim=c(0, zuigao*1.15 ),
       xlab="Chromosome", ylab=expression(paste( "-log"[10], "(", italic(p), ")", sep="" )),
       xaxs="i", yaxs="i",col="#EE9572",lwd=2,
       axes=FALSE, mgp=c(1.8,0.5,0), cex.lab=1.3)
  for(i in 2:(length(leiji_len)-1)){
    segments(leiji_len[i],-10,leiji_len[i],100)
  }
  lines(position,pvs,lwd=2,col="#EE9572")
  eachs <-  zuigao * 1.1 / 3;
  ats <- c(0, eachs, 2*eachs, 3*eachs);
  mid_pos <- rep(0, 12);
  for(i in 1:12){  mid_pos[i] <- (leiji_len[i] + leiji_len[i+1])/2;  }
  axis(side=1, at=mid_pos[1:12], labels=c(1:12), tick=FALSE, cex.axis=0.8, mgp=c(0.3,0.3,0) );#, "scaffold"
  axis(side=1, at=position, labels=rep("",length(position)), tick=0.2, cex.axis=0.8,col="#E5E5E5",lwd=0.2);#, "scaffold"
  axis(side=2, at=ats, labels=round(ats,1), lwd=0, lwd.tick=1.2, tick=TRUE, mgp=c(0.5,0.5,0), cex.axis=1.1);
  abline(h=0,lty=1,col="black",lwd=1.2);
  if(threshold>max(pvs, na.rm=TRUE)){
    abline(h=threshold*1.15,lty=1,col="black",lwd=1.2);
  }else{
    abline(h=max(pvs, na.rm=TRUE)*1.15,lty=1,col="black",lwd=1.2);
  }
  abline(v=0,lty=1,col="black",lwd=1.2);
  abline(h=threshold,col="black",lwd=1.2,lty=2);
  abline(v=leiji_len[13],lty=1,col="black",lwd=1.2);
  dev.off()
  
}

#!/usr/bin/env Rscript
library("argparse")
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--prefix", action="store")
parser$add_argument("--maxp", default=100, action="store", type="integer")
parser$add_argument("--recessive", default=FALSE, action="store_true")
parser$add_argument("--removenull", default=FALSE, action="store_true")
args <- parser$parse_args()

dat<-as.data.frame( read.delim(paste(args$prefix, ".burden.txt", sep=""), header=T, stringsAsFactors = F, sep="\t") )

if(args$recessive==TRUE){
  dat$log10p_obs<-(-log10(dat$P_REC))
  plotname<-paste(args$prefix, ".REC.QQ.png", sep="")
}else{
  dat$log10p_obs<-(-log10(dat$P_DOM))
  plotname<-paste(args$prefix, ".QQ.png", sep="")
}

dat<-subset(dat, log10p_obs>0)

dat = dat[!is.nan(dat$log10p_obs),]

dat<-dat[order(dat$log10p_obs),]
set.seed(123)
dat$log10p_exp <- sort(-log10(runif(nrow(dat))))

#Calculate slope of line 
p0<-0
u0=suppressWarnings(max(dat[dat$log10p_obs==0,]$log10p_exp))
tmp.u0=suppressWarnings(max(dat[dat$log10p_obs==0,]$log10p_exp))
u0=ifelse( is.infinite( tmp.u0 ) == T , 0 , tmp.u0 )
tmp<-subset(dat, log10p_obs!=0)
p95<-tmp[round(nrow(tmp[tmp$log10p_obs!=0,]))*0.95,]$log10p_obs
u95<-tmp[round(nrow(tmp[tmp$log10p_obs!=0,]))*0.95,]$log10p_exp
slope<-(p95-p0)/(u95-u0)
yint<-p95-slope*u95

maxp<-ceiling(max(dat$log10p_obs))
if(maxp>args$maxp){maxp<-args$maxp}

#Generate the actual plot
options(bitmapType='cairo')
png(plotname, width=500, height=500)
par(mar=c(4,5,2,2))
plot(x=dat$log10p_exp, y=dat$log10p_obs, 
xlim=c(0,maxp), ylim=c(0, maxp), 
cex=0.5, pch=19,cex.lab=1.5, cex.axis=1.3,xlab=expression(Expected ~ -log[10](p)), ylab=expression(Observed ~ -log[10](p)), main="")
abline(yint, slope, col="blue", lty=2)
abline(0,1)
mtext(paste("lambda =", round(slope, digits=5)), side=3)
invisible(dev.off())

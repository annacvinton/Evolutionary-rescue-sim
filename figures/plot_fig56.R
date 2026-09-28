#!/usr/bin/env Rscript
## Fig 5: without a gradient there is no rescue. Fig 6: robustness composite.
f5 <- read.csv("fig5_data.csv")
pdf("fig5_homogeneous_v3.pdf", width=8.4, height=3.9)
par(mfrow=c(1,2), mar=c(4.4,4.4,2.6,1), mgp=c(2.7,.7,0), tcl=-.3, las=1)
cols <- c("#888888","#2a7f8f","#2a7f8f","#2a7f8f")
bp <- barplot(f5$baseline, names.arg=f5$slope, col=cols, ylab="pre-perturbation abundance",
        xlab="gradient slope", main="no gradient: largest populations...", cex.main=.95)
bp <- barplot(f5$persistence, names.arg=f5$slope, col=cols, ylim=c(0,0.7),
        ylab="fraction of runs persisting", xlab="gradient slope",
        main="...and almost no rescue", cex.main=.95)
text(bp[1], f5$persistence[1]+0.04, "2 of 100", cex=.8)
dev.off()
cat("wrote fig5\n")

dr <- read.csv("fig6_draws.csv"); kn <- read.csv("fig6_kernel.csv")
vb <- read.csv("viability_final.csv")
pdf("fig6_robustness_v3.pdf", width=10.5, height=3.6)
par(mfrow=c(1,3), mar=c(4.2,4.2,2.6,1), mgp=c(2.5,.7,0), tcl=-.3, las=1)
plot(NA, xlim=c(-.3,4.3), ylim=c(.45,.75), xaxt="n", xlab="autocorrelation (patch SD 2)",
     ylab="recovery (survivors)", main="five landscape replicates agree", cex.main=1)
axis(1,c(0,2,4))
for (w in unique(dr$draw)) { q<-dr[dr$draw==w,]; lines(q$ac,q$recover,type="b",pch=19,col=adjustcolor("#d1762f",.65)) }
plot(NA, xlim=c(-.3,4.3), ylim=c(.3,.8), xaxt="n", xlab="autocorrelation (patch SD 2)",
     ylab="persistence", main="kernel shape is irrelevant", cex.main=1)
axis(1,c(0,2,4))
for (k in 0:1) { q<-kn[kn$kernel==k,]; q<-q[order(q$cond),]
  lines(c(0,2,4),q$persist,type="b",pch=c(19,17)[k+1],lwd=2,col=c("#555555","#2a7f8f")[k+1]) }
legend("bottomright",bty="n",pch=c(19,17),col=c("#555555","#2a7f8f"),
       legend=c("exponential (original)","Gaussian (v3)"),cex=.9)
v <- vb[vb$slope %in% c(.6,.8,1,1.2),]
tab <- table(factor(v$slope), factor(v$klass, levels=c("viable","declining","non-viable")))
barplot(t(tab), beside=FALSE, col=c("#1d5a73","#7fb2c9","#f4f0e6"),
        xlab="gradient slope", ylab="treatment cells", main="every analysed cell certified by controls", cex.main=1)
legend("topright",bty="n",fill=c("#1d5a73","#7fb2c9","#f4f0e6"),
       legend=c("viable","slowly declining","non-viable"),cex=.85)
dev.off()
cat("wrote fig6\n")

#!/usr/bin/env Rscript
# Fig 3: autocorrelation raises abundance, and persistence follows.
g <- read.csv("fig3_panels.csv")
sdcol <- c("1"="#8f6f2a","2"="#2a5a8f")
pdf("fig3_abundance_mediation_v3.pdf", width=8.8, height=4.0)
par(mfrow=c(1,2), mar=c(4,4.4,2.6,1), mgp=c(2.6,.7,0), tcl=-.3, las=1)
for (panel in 1:2) {
  ylab <- if(panel==1) "pre-perturbation abundance (cell median)" else "fraction of runs persisting"
  main <- if(panel==1) "autocorrelation sets abundance..." else "...persistence shifts modestly in the viable interior"
  yl <- if(panel==1) c(350,1350) else c(.35,.75)
  plot(NA, xlim=c(-.3,4.3), ylim=yl, xaxt="n", xlab="spatial autocorrelation", ylab=ylab, main=main, cex.main=1)
  axis(1, c(0,2,4))
  for (sd in c(1,2)) {
    q <- g[g$patch_sd==sd,]; q<-q[order(q$ac),]; cc <- sdcol[as.character(sd)]
    if (panel==1) { arrows(q$ac,q$base_lo,y1=q$base_hi,angle=90,code=3,length=.03,col=cc)
                    lines(q$ac,q$base,type="b",pch=19,lwd=2.2,col=cc) }
    else          { arrows(q$ac,q$persist-q$pse,y1=q$persist+q$pse,angle=90,code=3,length=.03,col=cc)
                    lines(q$ac,q$persist,type="b",pch=19,lwd=2.2,col=cc) }
  }
  legend("bottomright", bty="n", lwd=2.2, pch=19, col=sdcol, legend=c("patch SD 1","patch SD 2"), cex=.9)
}
dev.off()
cat("wrote fig3\n")

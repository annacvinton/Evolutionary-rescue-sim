#!/usr/bin/env Rscript
# Figure 2: the three currencies of rescue. Rows = outcomes, columns = factors.
# Reads phase_profiles_v3.csv (cluster-mean +/- SE). Base R only.
p <- read.csv("phase_profiles_v3.csv")
facs <- c("dispersal","slope","mutation","autocorrelation")
outs <- c("resistance","recovery","persistence")
xlabs <- c(dispersal="dispersal SD", slope="gradient slope",
           mutation="mutational SD", autocorrelation="autocorrelation (patch SD 2)")
cols <- c(resistance="#2a7f8f", recovery="#d1762f", persistence="#555555")
pdf("fig2_phase_profiles_v3.pdf", width=10.5, height=7.6)
par(mfrow=c(3,4), mar=c(3.6,3.9,1.6,.8), mgp=c(2.3,.7,0), tcl=-.3, las=1, oma=c(0,1.4,1.8,0))
for (o in outs) {
 qq <- p[p$outcome==o,]
 yl <- range(c(qq$mean-2*qq$se, qq$mean+2*qq$se)); yl <- yl + c(-.08,.08)*diff(yl)
 for (f in facs) {
  q <- p[p$factor==f & p$outcome==o,]
  q <- q[order(q$level),]
  ylim <- yl
  plot(seq_len(nrow(q)), q$mean, type="b", pch=19, lwd=2.2, col=cols[o],
       xlim=c(.7,nrow(q)+.3), ylim=ylim, xaxt="n",
       xlab=if(o=="persistence") xlabs[f] else "", ylab="", main="")
  axis(1, seq_len(nrow(q)), q$level)
  arrows(seq_len(nrow(q)), q$mean-q$se, y1=q$mean+q$se, angle=90, code=3, length=.03, col=cols[o])
  if (f=="dispersal") mtext(o, side=2, line=2.9, cex=.8, las=0, font=2)
  if (o=="resistance") mtext(xlabs[f], side=3, line=.4, cex=.75, font=2)
 }
}
mtext("Each factor acts on a different outcome", outer=TRUE, line=.3, cex=1, font=2)
dev.off()
cat("wrote fig2_phase_profiles_v3.pdf\n")

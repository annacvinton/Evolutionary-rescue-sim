#!/usr/bin/env Rscript
# Fig 3: the landscape acts through abundance. Fig 4: where populations persist at all.
ce <- read.csv("cells_persist_base_v3.csv")
vb <- read.csv("viability_labeled.csv")
accol <- c("0"="#c15a2a","2"="#b0a135","4"="#2a8f8f")

pdf("fig3_abundance_mediation_v3.pdf", width=9.2, height=4.2)
par(mfrow=c(1,2), mar=c(4,4.2,2.4,1), mgp=c(2.5,.7,0), tcl=-.3, las=1)
q <- ce[ce$patch_sd>0 & ce$disp==3 & ce$mutsd==0.75 & ce$slope==0.8,]   # centre cell: abundance varies only via the landscape
plot(q$base, q$persist, log="x", pch=19, cex=.55, col=adjustcolor(accol[as.character(q$ac)],.55),
     xlab="pre-perturbation abundance (cell median)", ylab="fraction of runs persisting",
     main="persistence follows abundance;\nautocorrelation adds nothing beyond it (slope 0.8, disp 3, mut on)", cex.main=.95)
lo <- loess(persist~log(base), q); xs <- exp(seq(log(min(q$base)),log(max(q$base)),len=80))
lines(xs, predict(lo, data.frame(base=xs)), lwd=2.4, col="grey25")
legend("topleft", bty="n", pch=19, col=accol, legend=paste("AC",names(accol)), cex=.85)
b <- aggregate(base~ac+patch_sd, ce[ce$patch_sd>0,], median)
bp <- matrix(b$base[order(b$patch_sd,b$ac)], nrow=3)
barplot(bp, beside=TRUE, col=accol, names.arg=c("patch SD 1","patch SD 2"),
        ylab="pre-perturbation abundance", main="autocorrelation sets abundance,\nmost strongly on rough terrain", cex.main=.95)
legend("topright", bty="n", fill=accol, legend=paste("AC",names(accol)), cex=.85)
dev.off()
cat("wrote fig3\n")

pdf("fig4_viability_map_v3.pdf", width=10, height=3.4)
par(mfrow=c(1,4), mar=c(5.2,4.6,2.4,.8), las=1)
conds <- c("flat","SD1/AC0","SD1/AC2","SD1/AC4","SD2/AC0","SD2/AC2","SD2/AC4")
pal <- colorRampPalette(c("#f4f0e6","#7fb2c9","#1d5a73"))(21)
for (sl in c(0.6,0.8,1.0,1.2)) {
  q <- vb[vb$slope==sl & vb$mutsd==0.75,]
  M <- matrix(NA,3,7,dimnames=list(c(1.5,3,6),conds))
  for (i in 1:nrow(q)) M[as.character(q$disp[i]), q$label[i]] <- q$surv[i]
  image(1:7,1:3,t(M),col=pal,zlim=c(0,1),axes=FALSE,xlab="",ylab="dispersal SD",
        main=paste("slope",sl), cex.main=1)
  axis(1,1:7,conds,las=2,cex.axis=.75); axis(2,1:3,c(1.5,3,6))
  for(i in 1:3) for(j in 1:7) if(!is.na(M[i,j])) text(j,i,sprintf("%.1f",M[i,j]),cex=.65,col=ifelse(M[i,j]>.6,"white","grey20"))
}
dev.off()
cat("wrote fig4\n")

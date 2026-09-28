#!/usr/bin/env Rscript
# Supplementary figures S1-S4 from run_summary_v3.csv + viability_final.csv.
# S1 slope x dispersal (the dominant interaction) on all three outcomes
# S2 variance-infusion interactions on recovery (with dispersal; with magnitude)
# S3 magnitude-of-change profiles on all three outcomes
# S4 interaction effect-size heatmap (15 interactions x 3 outcomes)
d  <- read.csv("run_summary_v3.csv"); via <- read.csv("viability_final.csv")
K5 <- c("slope","disp","mutsd","patch_sd","ac")
d  <- merge(d[d$slope > 0, ], via[, c(K5,"surv")], by = K5)
core <- d[d$slope < 1.2 & d$surv >= 0.8, ]
core$resist <- core$trough_n/core$base_n; core$recover <- (core$peak_n-core$trough_n)/core$base_n
s <- core[core$survived == 1 & core$base_n > 0 & core$peak_t > core$trough_t, ]
cols3 <- c("#2a7f8f","#d1762f","#555555"); names(cols3) <- c("resistance","recovery","persistence")
pal   <- c("#7fb2c9","#3d7f9a","#1d5a73")

## S1
pdf("figS1_slope_x_dispersal_v3.pdf", width=10, height=3.6)
par(mfrow=c(1,3), mar=c(4.2,4.4,2.4,1), mgp=c(2.6,.7,0), tcl=-.3, las=1)
for (o in c("resistance","recovery","persistence")) {
  dat <- if (o=="persistence") core else s
  yv  <- c(resistance="resist",recovery="recover",persistence="survived")[o]
  ag  <- aggregate(dat[[yv]], dat[c("slope","disp")], mean)
  plot(NA, xlim=c(1,6.5), ylim=range(ag$x)*c(.9,1.08), xlab="dispersal SD", ylab=o, main=paste(o,"by slope"), log="x", xaxt="n")
  axis(1, c(1.5,3,6))
  for (i in 1:3) { sl <- c(.6,.8,1)[i]; q <- ag[ag$slope==sl,]; q<-q[order(q$disp),]
    lines(q$disp, q$x, type="b", pch=19, lwd=2.2, col=pal[i]) }
  if (o=="resistance") legend("topleft", bty="n", lwd=2.2, pch=19, col=pal, legend=paste("slope",c(.6,.8,1)), cex=.9)
}
dev.off(); cat("wrote S1\n")

## S2
pdf("figS2_variance_infusion_interactions_v3.pdf", width=8.4, height=3.8)
par(mfrow=c(1,2), mar=c(4.2,4.4,2.4,1), mgp=c(2.6,.7,0), tcl=-.3, las=1)
mc <- c("0"="#999999","0.75"="#d1762f")
ag <- aggregate(s$recover, s[c("mutsd","disp")], mean)
plot(NA, xlim=c(1,6.5), ylim=c(0,1.1), log="x", xaxt="n", xlab="dispersal SD", ylab="recovery", main="variance infusion x dispersal")
axis(1,c(1.5,3,6))
for (m in c(0,.75)) { q<-ag[ag$mutsd==m,]; q<-q[order(q$disp),]; lines(q$disp,q$x,type="b",pch=19,lwd=2.2,col=mc[as.character(m)]) }
legend("topleft",bty="n",lwd=2.2,pch=19,col=mc,legend=c("no variance input","variance input"),cex=.9)
ag <- aggregate(s$recover, s[c("mutsd","pert_treat")], mean)
plot(NA, xlim=c(8.7,13.3), ylim=c(0,1.1), xlab="magnitude of environmental change", ylab="recovery", main="variance infusion x magnitude")
for (m in c(0,.75)) { q<-ag[ag$mutsd==m,]; q<-q[order(q$pert_treat),]; lines(q$pert_treat,q$x,type="b",pch=19,lwd=2.2,col=mc[as.character(m)]) }
dev.off(); cat("wrote S2\n")

## S3
pdf("figS3_magnitude_profiles_v3.pdf", width=10, height=3.4)
par(mfrow=c(1,3), mar=c(4.2,4.4,2.4,1), mgp=c(2.6,.7,0), tcl=-.3, las=1)
for (o in c("resistance","recovery","persistence")) {
  dat <- if (o=="persistence") core else s
  yv  <- c(resistance="resist",recovery="recover",persistence="survived")[o]
  ag  <- aggregate(dat[[yv]], list(P = dat$pert_treat), mean)
  plot(ag$P, ag$x, type="b", pch=19, lwd=2.2, col=cols3[o], xlab="magnitude of environmental change", ylab=o, main=paste(o,"by magnitude"))
}
dev.off(); cat("wrote S3\n")

## S4
ie <- read.csv("interactions_final.csv", row.names=1)
ie <- ie[order(-ie$recovery), ]
pdf("figS4_interaction_heatmap_v3.pdf", width=6.2, height=5.6)
par(mar=c(4,12,3,1), las=1)
M <- as.matrix(ie[, c("resistance","recovery","persistence_dR2")]); colnames(M) <- c("resistance","recovery","persistence")
M <- M[nrow(M):1, ]
image(1:3, 1:nrow(M), t(M), col=colorRampPalette(c("#f7f4ec","#7fb2c9","#1d5a73"))(30), axes=FALSE, xlab="", ylab="",
      main="variance explained by each pairwise interaction", cex.main=.95)
axis(1, 1:3, colnames(M)); axis(2, 1:nrow(M), rownames(M), cex.axis=.8)
for (i in 1:nrow(M)) for (j in 1:3) text(j, i, sprintf("%.1f", M[i,j]), cex=.7, col=ifelse(M[i,j]>8,"white","grey20"))
mtext("partial eta-squared (%) for the phases; delta pseudo-R-squared (pp) for persistence", side=1, line=2.6, cex=.75)
dev.off(); cat("wrote S4\n")

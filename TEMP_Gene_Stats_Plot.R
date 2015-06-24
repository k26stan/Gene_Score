
# TAB <- read.table( "Downloads/20150604_Phased/Gene_Stats.txt",sep="\t",header=T )
# TAB <- read.table( "Downloads/20150604_Phased/Chr22/Gene_Stats.1.txt",sep="\t",header=T )
TAB.1 <- read.table( "Downloads/20150604_Phased/Chr22/Gene_Stats.1.txt",sep="\t",header=T )
TAB.2 <- read.table( "Downloads/20150604_Phased/Chr22/Gene_Stats.txt",sep="\t",header=T )
TAB <- rbind( TAB.1, TAB.2[401:nrow(TAB.2),] )
LIM <- c(0,7)
COLS <- c("dodgerblue2","chartreuse2","firebrick2","gold2")
plot( 0,0,type="n",xlim=LIM,ylim=LIM,xlab="-log10(p.exp)",ylab="-log10(p.obs)")
abline(0,1) ; abline(h=seq(0,10,1),lty=2,col="grey50") ; abline(v=seq(0,10,1),lty=2,col="grey50")
WHICH <- which(!is.na(TAB$P_CH))
points( -log10( 1:length(WHICH)/length(WHICH) ), -log10(sort(TAB$P_CH[WHICH])), col=COLS[1],pch="+" )
WHICH <- which(!is.na(TAB$P_HAP))
points( -log10( 1:length(WHICH)/length(WHICH) ), -log10(sort(TAB$P_HAP[WHICH])),col=COLS[2],pch="+" )
WHICH <- which(!is.na(TAB$P_BURD))
points( -log10( 1:length(WHICH)/length(WHICH) ), -log10(sort(TAB$P_BURD[WHICH])), col=COLS[3],pch="+" )
WHICH <- which(!is.na(TAB$P_BURD.ex))
points( -log10( 1:length(WHICH)/length(WHICH) ), -log10(sort(TAB$P_BURD.ex[WHICH])), col=COLS[4],pch="+" )
legend( "bottomright", fill=COLS,legend=c("Comp_Het","Haplo","Burd","Burd_Ex"),bg="white" )
abline(h=-log10(.01),lty=2,col="chocolate2",lwd=2) ; abline(v=-log10(.01),lty=2,col="chocolate2",lwd=2)

WHICH <- which( TAB$P_CH<.01 | TAB$P_HAP<.01 | TAB$P_BURD<.01 | TAB$P_BURD.ex<.01 )
TAB[ WHICH, c("GTX","n.VAR","P_CH","P_HAP","P_BURD","P_BURD.ex") ]

quartz()
par(mfrow=c(1,2))
hist( TAB$AREA, breaks=seq(-20,20,.5),xlim=extendrange(TAB$AREA) )
hist( TAB$AREA.ex, breaks=seq(-20,20,.5),xlim=extendrange(TAB$AREA.ex) )

load( "Downloads/20150604_Phased/Gene_Stats.Rdata" )
EXON <- COMPILE$Exon
FULL <- COMPILE$Full


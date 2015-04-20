## Rscript to plot Phased Variants for any given Gene ##
## October 14, 2014 ##
## Kristopher Standish ##

# Pull out Command Line
LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/Phased/20141014_Testing/TEMP_GENES.list","1")
# LINE <- c("/projects/janssen/Phased/20141014_Testing/TEMP_GENES.2.list","6")
# LINE <- c("/Users/kstandis/Data/Burn/Phased/20141014_Testing/TEMP_GENES.22.list","22")

# Parse Thru Command Line
INPUT_1.split <- strsplit(LINE[1],".",fixed=T)[[1]]
if ( INPUT_1.split[length(INPUT_1.split)]=="list" ) {
	PathToGeneList <- LINE[1]
	INPUT_1.type <- "List"
}else{ 
	GENE <- LINE[1] 
	INPUT_1.type <- "Gene" }
Which_Chrom <- LINE[2]
Which_Chr <- paste("chr",Which_Chrom,sep="")

###########################################################
## LOAD DATA ##############################################
###########################################################
print("Loading Data") 

Which_Samp <- 1

## Get Current Directory
PathToHome <- getwd()
PathToTab <- paste(PathToHome,"/Rout_",Which_Chrom,"_",Which_Samp,"_Table.txt",sep="")
PathToLoc <- paste(PathToHome,"/Rout_",Which_Chrom,"_",Which_Samp,"_Location.Rdata",sep="")
PathToMrg <- paste(PathToHome,"/Rout_",Which_Chrom,"_",Which_Samp,"_Merged.txt",sep="")
# PathToGeneTable <- "/Users/kstandis/Data/Tools/GG-Gene_Names_DB.txt"
PathToGeneTable <- "/home/kstandis/RV/Tools/GG-Gene_Names_DB.txt"

## Load Gene List
if ( INPUT_1.type=="List" ) { GENE_LIST <- as.character( read.table(PathToGeneList, header=F)[,1] ) }

## Load Table
TAB <- read.table(PathToTab, sep="\t",header=T)

## Load Locations
load(PathToLoc)
LOC <- LOCATIONS.save

## Load Merged Tables
MRG <- read.table(PathToMrg, sep="\t",header=T)

## Load Gene Table
GG_colClasses <- c( rep("character",2),rep("numeric",5),rep("character",11) )
GG <- read.table(PathToGeneTable, sep="\t",header=T, comment.char="",fill=T,quote="",colClasses=GG_colClasses)

###########################################################
## PULL OUT GENE DATA #####################################
###########################################################
print("Data Loaded; Creating Function")

PLOT_GENE <- function(GENE) {
## Set Path To Save To
PathToSave1 <- paste(PathToHome,"/Plot_Gene_Map_",GENE,".jpeg",sep="")
PathToSave2 <- paste(PathToHome,"/Plot_Gene_VarCount_",GENE,".jpeg",sep="")

## Get Gene Data
# print("Filtering Data by Gene")
TAB.g <- TAB[ which( rownames(TAB)==GENE & TAB$CHR==Which_Chr ), ]
LOC.g <- LOC[[ GENE ]]
MRG.g <- MRG[ which( MRG$GENE_ID==GENE & MRG$CHR==Which_Chrom ), ]
GG.g <- GG[ which( GG$hg19.kgXref.geneSymbol==GENE & GG$hg19.knownGene.chrom==Which_Chr ), ]

###########################################################
## PLOT GENE ##############################################
###########################################################

## Compile Parameters
NUM_TRANS <- nrow(GG.g)
NUM_VARS <- nrow(MRG.g)
COLS <- c("cadetblue1","dodgerblue1","chartreuse1","mediumpurple2","sienna1","firebrick1") # Orig
COLS <- c("slateblue3","cadetblue1","steelblue2","firebrick2","goldenrod2","springgreen3")
X_LIM <- extendrange( c(min(GG.g$hg19.knownGene.txStart),max(GG.g$hg19.knownGene.txEnd)), f=.2 )
Y_LIM <- c( 0, NUM_TRANS+1 )

## Plot Gene Map
# print("Plotting Gene Map")
jpeg(PathToSave1, height=500+50*NUM_TRANS,width=1400+NUM_VARS, pointsize=28+.1*NUM_TRANS)
plot(0,0,type="n", xlim=X_LIM,ylim=Y_LIM, xlab=paste("Chromosome",Which_Chrom,"Position"),ylab="Transcript", main=GENE, yaxt="n" )
## Loop through all Transcripts
for ( t in 1:NUM_TRANS ) {
	trans <- strsplit( as.character(GG.g$X.hg19.knownGene.name[t]),".",fixed=T )[[1]][1]
	## Specify Y-Values
	Y <- t

	## Get Variant Positions
	HAP_1 <- MRG.g$POS[ which(MRG.g$Hap_1==1) ]
	HAP_2 <- MRG.g$POS[ which(MRG.g$Hap_2==1) ]

	## Get Outer Coordinates
	TX_RNG <- c( GG.g$hg19.knownGene.txStart[t], GG.g$hg19.knownGene.txEnd[t] )
	CD_RNG <- c( GG.g$hg19.knownGene.cdsStart[t], GG.g$hg19.knownGene.cdsEnd[t] )

	## Get Exon Coordinates
	EX_B <- as.numeric( strsplit( as.character(GG.g$hg19.knownGene.exonStarts), "," )[[t]] )
	EX_E <- as.numeric( strsplit( as.character(GG.g$hg19.knownGene.exonEnds), "," )[[t]] )

	## Plot this Shiz
	arrows( TX_RNG[1],Y,TX_RNG[2],Y, code=3,angle=90, lwd=5,col=COLS[1] )
	abline( h=Y, col="grey50",lwd=1 )
	polygon( CD_RNG[c(1,2,2,1)],c(-.1,-.1,.1,.1)+Y, col=COLS[2],lwd=1 )
	for ( e in 1:length(EX_B) ) {
		X_COORDS <- c( EX_B[e],EX_E[e],EX_E[e],EX_B[e] )
		polygon( X_COORDS,c(-.2,-.2,.2,.2)+Y, col=COLS[3],lwd=1 )	
	}
	# arrows( TX_RNG[1],Y,TX_RNG[2],Y, code=2,angle=35, lwd=5,col=COLS[1] )
	arrows( TX_RNG[1],Y,TX_RNG[2],Y, code=3,angle=90,length=0, lwd=5,col=COLS[1] )

	if ( length(HAP_1)>0 ) { arrows( HAP_1,Y-.3,HAP_1,Y, lwd=2,col=COLS[5],length=.1 ) }
	if ( length(HAP_2)>0 ) { arrows( HAP_2,Y+.3,HAP_2,Y, lwd=2,col=COLS[4],length=.1 ) }
	if ( length(intersect(HAP_1,HAP_2))>0 ) { arrows( intersect(HAP_1,HAP_2),Y-.3,intersect(HAP_1,HAP_2),Y, lwd=2,col=COLS[6],length=.1 )
	arrows( intersect(HAP_1,HAP_2),Y+.3,intersect(HAP_1,HAP_2),Y, lwd=2,col=COLS[6],length=.1 ) }
	text( X_LIM[1], Y-.4, labels=trans, pos=4, cex=.7)
} # Close Transcripts Loop
dev.off()

## Plot Number of Variants of each Type on each Strand
if ( nrow(MRG.g)>0 ) {
	r <- floor(sqrt(NUM_TRANS)) ; c <- ceiling(sqrt(NUM_TRANS)) ; if (r*c < NUM_TRANS) { r <- r+1 }
	jpeg(PathToSave2, height=400+400*r,width=1200+200*c, pointsize=32)
	par(mfrow=c(r,c))
	for ( t in 1:NUM_TRANS ) {
		trans <- strsplit( as.character(GG.g$X.hg19.knownGene.name[t]),".",fixed=T )[[1]][1]
		WHICH_SLASH <- sapply( lapply( strsplit( as.character(MRG.g$Gene),"///" ), function(x) grep(trans,x) ), "[", 1)
		# print("Compiling Variant Counts by Strand")
		BAR_COLS <- COLS[4:6]
		STRAND_COUNTS <- array(0 , c( 3,ncol(LOC.g) ) )
		colnames(STRAND_COUNTS) <- c("Downstream","UTR","Exon","Intron","Upstream","Nons","ncRNA") # colnames(LOC.g)
		rownames(STRAND_COUNTS) <- c("Hap_1","Hap_2","Both")
		for ( l in 1:nrow(LOC.g) ) {
			if ( !is.na( WHICH_SLASH[l] ) ) {
				LOC.ID <- paste( strsplit( rownames(LOC.g)[l], ":" )[[1]][1:2], collapse=":" )
				# Which row of "MRG.g" to refer to
				MRG.IND <- which( MRG.g$ID == LOC.ID )
				# In what element is the variant located?
					# LOC.COL <- which( LOC.g[l,]!=0 )
				LOC.trans <- strsplit( as.character(MRG.g$Location)[ MRG.IND ],"///" )[[1]][ WHICH_SLASH[l] ]
				if ( grepl("noncoding_rna",LOC.trans,ignore.case=T) ) { LOC.COL <- "ncRNA" }
				if ( LOC.trans=="Downstream" ) { LOC.COL <- "Downstream" }
				if ( LOC.trans=="Upstream" ) { LOC.COL <- "Upstream" }
				if ( grepl("intron",LOC.trans,ignore.case=T) ) { LOC.COL <- "Intron" }
				if ( grepl("utr",LOC.trans,ignore.case=T) ) { LOC.COL <- "UTR" }
				if ( grepl("exon",LOC.trans,ignore.case=T) ) { LOC.COL <- "Exon"
					COD.trans <- strsplit( as.character(MRG.g$Coding_Impact)[ MRG.IND ],"///" )[[1]][ WHICH_SLASH[l] ]
					if ( grepl("nons",COD.trans,ignore.case=T) ) { LOC.COL <- "Nons" }
				}
				# Which strand is the variant on?
				WHICH_HAP <- which( MRG.g[ MRG.IND, c("Hap_1","Hap_2") ]==1 )
				if ( length(WHICH_HAP)==2 ) { HAP.IND <- 3 }else{ HAP.IND <- WHICH_HAP }
				# Add one to whichever strand it's on
				STRAND_COUNTS[ HAP.IND,LOC.COL ] <- STRAND_COUNTS[ HAP.IND,LOC.COL ] + 1
			}
		}
		# print("Making Barplot of Variant Counts")

		barplot(STRAND_COUNTS, main=paste(GENE,"-",trans,), xlab="",ylab="Variant Count", col=BAR_COLS, beside=T, las=2)
		# if ( colSums(STRAND_COUNTS)["ncRNA"] > 0 ) {
			# LEGEND.X <- 21 }else{ LEGEND.X <- 25 }
		# legend( LEGEND.X, max(STRAND_COUNTS), fill=BAR_COLS, legend=rownames(STRAND_COUNTS), cex=.7 )
		if ( t==1 ) { legend( 4, max(STRAND_COUNTS), fill=BAR_COLS, legend=rownames(STRAND_COUNTS), cex=.7 ) }
	} # Close "NUM_TRANS" loop
	dev.off()
} # Close "if ( nrow(LOC.g)>0 )"

} # Close "PLOT_GENE" Function

## Loop thru Genes to actually plot these
if ( INPUT_1.type=="List" ) { 
	print("Creating Plots")
	loop <- 0
	for (GENE in GENE_LIST) {
		loop <- loop + 1
		PLOT_GENE(GENE)
		if ( loop %% 5 == 0 ) { print( paste("Finished",loop,"of",length(GENE_LIST)) ) }
	} # Close GENE loop
}else{
	print(paste("Creating Plot for",GENE))
	PLOT_GENE(GENE)
} # Close "IF"

###########################################################
## END OF DOC #############################################
###########################################################




# PLOT_GENE("FAM83F")
# PLOT_GENE("FBXW4P1")
# PLOT_GENE("GAB4")
# PLOT_GENE("GAL3ST1")
# PLOT_GENE("GAS2L1")
# PLOT_GENE("FLJ32756")
# PLOT_GENE("FBXO7")
# PLOT_GENE("GAS2L1")
# PLOT_GENE("GAS2L1")



# FAM83F
# FBLN1
# FBXO7
# FBXW4P1
# FLJ23584
# FLJ32756
# FLJ41941
# FOXRED2
# GAB4
# GAL3ST1
# GALR3
# GAS2L1
# GATSL3


#  "A4GALT"  "abParts" "ACO2"    "ACR"     "ADM2"    "ADORA2A"

# > head(rownames(TAB),40)
#  [1] "AK022914" "AK056135" "DQ573684" "POTEH"    "OR11H1"   "DQ571479"
#  [7] "CCT8L2"   "BC038197" "XKR3"     "GAB4"     "CECR7"    "IL17RA"  
# [13] "CECR6"    "BC021738" "CECR5"    "CECR1"    "CECR2"    "SLC25A18"
# [19] "ATP6V1E1" "BCL2L13"  "BID"      "MICAL3"   "FLJ41941" "PEX26"   
# [25] "TUBA8"    "USP18"    "GGT3P"    "BC112340" "BC051721" "DGCR6"   
# [31] "PRODH"    "DGCR5"    "DGCR9"    "DGCR10"   "DGCR2"    "DGCR14"  
# [37] "TSSK2"    "GSC2"     "SLC25A1"  "CLTCL1" 



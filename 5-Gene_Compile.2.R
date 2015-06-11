## Script to Loop Through Genes and Compile/Plots Things ##
## Called by Gene_Score Pipeline ##
## February 16, 2015 ##
## Kristopher Standish ##


## Usage ##
# Rscript 5-Gene_Compile.R <Path/To/Out_Dir> <Path/To/Pheno> <Covs_Command>

# 1774 - SNAR-A3_uc010ybt.1 

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( "/projects/janssen/Phased/20150218_Run","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
# LINE <- c( "/projects/janssen/Phased/20150603_Test_Genes","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
# LINE <- c( "/projects/janssen/Phased/20150604_Chr22","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
PathToOut <- LINE[1]
PathToPheno <- LINE[2]
Cov_List <- LINE[3]

# Specify other Paths
PathToCovFile <- paste( PathToOut, "Assoc/Cov_w_PCs.txt", sep="/" )
PathToGeneCoords <- paste( PathToOut, "Gene_Info.txt", sep="/" )
PathToGenes <- paste( PathToOut, "Genes/", sep="/" )
PathToAssoc <- paste( PathToOut, "Assoc/", sep="/" )

# Check for proper parsing
print(paste( "Output:", PathToOut ))
print(paste( "Pheno:", PathToPheno ))
print(paste( "Cov List:", Cov_List ))
print(paste( "Cov File:", PathToCovFile ))
print(paste( "Gene Coords:", PathToGeneCoords ))
print(paste( "Gene Files:", PathToGenes ))
print(paste( "Association:", PathToAssoc ))

# Load Necessary Libraries
library(gplots)

###############################################################
## LOAD DATA ##################################################
###############################################################

## Load Phenotype/Covariate Files (if exists)
if ( file.exists( PathToPheno ) ) {
	print( "Loading Phenotype/Association Files" )
	# Load Phenotype/Covariate/Assoc Files
	PHENO <- read.table( PathToPheno, sep="\t",header=T )
	COVS.l <- read.table( PathToCovFile, sep="\t",header=T )
	# Reformat Covariate Table
	Cov_List.sp <- strsplit( Cov_List, "," )[[1]]
	COVS <- COVS.l[, c("IID",Cov_List.sp) ]
	# Merge Phenotype/Covariate Files
	PC <- merge( x=PHENO[,c("IID","Pheno")], y=COVS, by="IID" )
	PC <- PC
	PC[,"IID"] <- unlist(strsplit( as.character(PC$IID),"-"))[seq(1,2*nrow(PC),2)]
	## Load Association Results
	SNP.P <- read.table( paste(PathToAssoc,"SNP/SNP_Assoc.P",sep=""), sep="\t",header=T)
	SNP.HWE <- read.table( paste(PathToAssoc,"SNP/SNP_Assoc.hwe",sep=""), sep="",header=T)
	IND.P <- read.table( paste(PathToAssoc,"IND/IND_Assoc.P",sep=""), sep="\t",header=T)
	IND.HWE <- read.table( paste(PathToAssoc,"IND/IND_Assoc.hwe",sep=""), sep="",header=T)

	## Get Phenotype Sample List
	pheno.samps <- as.character( PC[,"IID"] )
}else{ print("No Phenotype Provided") }

## Load BIM files
print( "Loading BIM Files" )
SNP.bim <- read.table( paste(PathToOut,"/SNP_Vars.bim",sep=""), sep="\t",header=F)
colnames(SNP.bim) <- c("CHR","SNP","XXX","BP","REF","ALT")
IND.bim <- read.table( paste(PathToOut,"/IND_Vars.bim",sep=""), sep="\t",header=F)
colnames(IND.bim) <- c("CHR","SNP","XXX","BP","REF","ALT")

## Load Gene Coords
print( "Loading Gene Coordinate File" )
GENE_COORDS <- read.table( PathToGeneCoords, sep="\t", header=T, comment.char="" )
colnames(GENE_COORDS) <- gsub("X.","", colnames(GENE_COORDS), fixed=T )
colnames(GENE_COORDS) <- gsub("hg19.knownGene.","", colnames(GENE_COORDS), fixed=T )
colnames(GENE_COORDS) <- gsub("hg19.kgXref.","", colnames(GENE_COORDS), fixed=T )
GENE_COORDS <- GENE_COORDS[ which(GENE_COORDS[,"chrom"]!="chrX"), ]

###############################################################
## GET ORGANIZED ##############################################
###############################################################

## Parse Gene Coordinates
 # List of Gene_Transcripts
GTX_LIST <- paste( GENE_COORDS[,"geneSymbol"], GENE_COORDS[,"name"], sep="_" )
n.gtx <- length(GTX_LIST)
 # Transcript/Coding/Exon Start/Stop
CHR <- GENE_COORDS[,"chrom"]
TX_TAB <- data.frame( TX_S=GENE_COORDS[,"txStart"], TX_E=GENE_COORDS[,"txEnd"] )
RNG_TAB <- data.frame( RNG_S=TX_TAB[,"TX_S"]-5000, RNG_E=TX_TAB[,"TX_E"]+5000 )
CD_TAB <- data.frame( CD_S=GENE_COORDS[,"cdsStart"], CD_E=GENE_COORDS[,"cdsEnd"] )
EX_B_LIST <- strsplit( as.character(GENE_COORDS$exonStarts), "," )
EX_E_LIST <- strsplit( as.character(GENE_COORDS$exonEnds), "," )
EX_COUNT <- GENE_COORDS[,"exonCount"]

## Get out Sample Names for Haplotyped File
print( "Loading Phased Sample List" )
hap.samps <- as.character( read.table( paste(PathToOut,"/Phased.sample",sep=""), sep="",header=T )[,1] )
n.hap.samps <- length(hap.samps)
hap.colnames <- c("CHR","SNP","BP","REF","ALT", paste( rep(hap.samps, rep(2,n.hap.samps)), 1:2, sep="_" ) )

###############################################################
## LOOP THROUGH GENES #########################################
###############################################################

## Set up Objects to Compile
 # Basic Gene Info
cats.COOR <- c("GTX","CHR","TX_S","TX_E","CD_S","CD_E","EX_S","EX_E","n.EX" )
cats.VARS <- c( "n.VAR","n.SNP","n.IND","n.VAR.ph","n.SNP.ph","n.IND.ph","n.VAR.ex","n.SNP.ex","n.IND.ex" )
cats.PHAS <- c( "n.VAR.ph.ex","n.SNP.ph.ex","n.IND.ph.ex" )
cats.GWAS <- c( "AREA","AREA.ex","BEST_P","BEST_P.ex" )
cats.CH <- c( "p.COMP_HET_snp","p.COMP_HET_ind","p.COMP_HET_snp.ex","p.COMP_HET_ind.ex" )
cats.P <- c( "P_HAP","P_CH","P_BURD","P_BURD.ex" )
GENE.cats <- c( cats.COOR, cats.VARS, cats.PHAS, cats.GWAS, cats.CH, cats.P )
GENE <- array( , c(n.gtx,length(GENE.cats)) )
colnames(GENE) <- GENE.cats
GENE[,"GTX"] <- GTX_LIST
GENE[,"CHR"] <- CHR
GENE[,"TX_S"] <- TX_TAB[,"TX_S"]
GENE[,"TX_E"] <- TX_TAB[,"TX_E"]
GENE[,"CD_S"] <- CD_TAB[,"CD_S"]
GENE[,"CD_E"] <- CD_TAB[,"CD_E"]
GENE[,"EX_S"] <- as.character( GENE_COORDS$exonStarts )
GENE[,"EX_E"] <- as.character( GENE_COORDS$exonEnds )
GENE[,"n.EX"] <- EX_COUNT
 # Cohort Data
FULL <- EXON <- DAMG <- list()

PLOT_FRACTION <- 1/10 # 1/50
start_time <- proc.time()
## LOOP START ######################################
# for ( gtx in 1:n.gtx ) {
for ( gtx in 1562:n.gtx ) {
# for ( gtx in 1751:n.gtx ) {
	# gtx <- grep( "GRIN2B", GTX_LIST )
	# gtx <- grep( "IL31RA_uc011cqj.2", GTX_LIST )
	## Compile Info on Gene_Transcript
	name <- GTX_LIST[gtx]
	chr <- gsub( "chr","", as.character( CHR[gtx] ) )
	rng <- as.numeric( RNG_TAB[gtx,] )
	tx_rng <- as.numeric( TX_TAB[gtx,] )
	cd_rng <- as.numeric( CD_TAB[gtx,] )
	ex_b <- as.numeric( EX_B_LIST[[gtx]] )
	ex_e <- as.numeric( EX_E_LIST[[gtx]] )
	ex_cnt <- as.numeric( EX_COUNT[gtx] )

	## Print Status each Loop
	print(paste( "####### Starting",gtx,"of",n.gtx,"-",name,"#####",round(proc.time()-start_time,2)[3] ))
	
	####################################################
	## LOAD/PULL GTX DATA/FILES ########################
	print(paste( "Loading Files",round(proc.time()-start_time,2)[3] ))

	## Check if files exist & contain >0 Lines
	 # If not, skip to next Gene_Transcript
	if ( !file.exists(paste(PathToGenes,name,"/SNP_Vars.hwe",sep="")) | !file.exists(paste(PathToGenes,name,"/IND_Vars.hwe",sep="")) ) { next }
	if ( length(readLines(paste(PathToGenes,name,"/Phased.haps",sep="")))==0 ) { next }
	
	## Load HW Files
	snp.hwe.l <- read.table( paste(PathToGenes,name,"/SNP_Vars.hwe",sep=""), sep="",header=T )
	ind.hwe.l <- read.table( paste(PathToGenes,name,"/IND_Vars.hwe",sep=""), sep="",header=T )

	## Load Genotype Files
	snp.raw.l <- read.table( paste(PathToGenes,name,"/SNP_Vars.raw",sep=""), sep="",header=T )
	ind.raw.l <- read.table( paste(PathToGenes,name,"/IND_Vars.raw",sep=""), sep="",header=T )

	## Load Phased Haplotype File (& Rename Columns)
	hap.l <- read.table( paste(PathToGenes,name,"/Phased.haps",sep=""), sep="",header=F )
	colnames(hap.l) <- hap.colnames
	snp.hap <- hap.l[ which( hap.l$REF%in%c("C","G","T","A") & hap.l$ALT%in%c("C","G","T","A") ), ]
	ind.hap <- hap.l[ -which( hap.l$REF%in%c("C","G","T","A") & hap.l$ALT%in%c("C","G","T","A") ), ]

	## Pull out Variant Positions from BIM files (& Compile to 1 file)
	snp.bim.l <- SNP.bim[ which(SNP.bim$CHR==chr & SNP.bim$BP>=rng[1] & SNP.bim$BP<=rng[2] ), ]
	ind.bim.l <- IND.bim[ which(IND.bim$CHR==chr & IND.bim$BP>=rng[1] & IND.bim$BP<=rng[2] ), ]

	## Pull Out Single-Locus Results
	if ( file.exists(PathToPheno) ) {
		snp.p.l <- SNP.P[ which(SNP.P$CHR==chr & SNP.P$BP>=tx_rng[1]-5000 & SNP.P$BP<=tx_rng[2]+5000 ), ]
		ind.p.l <- IND.P[ which(IND.P$CHR==chr & IND.P$BP>=tx_rng[1]-5000 & IND.P$BP<=tx_rng[2]+5000 ), ]
	}

	## Load Annotation File
	# annot <- read.table( paste(PathToGenes,name,"/Annots_Short.txt",sep=""), sep="\t",header=T )
	# annot <- data.frame( TAG=paste(gsub("chr","",as.character(annot$Chromosome)),as.character(annot$End),sep="_"), annot)
	# eqtls <- read.table( paste(PathToGenes,name,"/eQTLs.txt",sep=""), sep="",header=T )

	####################################################
	## FILTER/ORGANIZE SOME TABLES #####################
	print(paste( "Filtering Tables",round(proc.time()-start_time,2)[3] ))

	## Get Intersection of Samples
	snp.raw <- snp.raw.l ; ind.raw <- ind.raw.l
	snp.raw[,"IID"] <- sapply(strsplit(as.character(snp.raw[,"IID"]),"-"),"[",1)
	ind.raw[,"IID"] <- sapply(strsplit(as.character(ind.raw[,"IID"]),"-"),"[",1)
	raw.samps <- intersect( snp.raw[,"IID"], ind.raw[,"IID"] )
	if ( file.exists(PathToPheno) ) {
		shared.samps <- sort( Reduce( intersect, list(hap.samps,raw.samps,pheno.samps) ) )
	}else{
		shared.samps <- sort( intersect( hap.samps, raw.samps ) )
	}
	n.samps <- length(shared.samps)

	## Remove unused Samples from RAW, HAP, and Pheno (if.exists) tables
	 # Raw Files
	snp.raw <- snp.raw[ which(snp.raw[,"IID"]%in%shared.samps), ]
	ind.raw <- ind.raw[ which(ind.raw[,"IID"]%in%shared.samps), ]
	 # Hap File
	hap.colnames.samp <- paste( rep(shared.samps, rep(2,n.samps)), 1:2, sep="_" )
	snp.hap <- snp.hap[ ,c( 1:grep("ALT",colnames(snp.hap)),which(colnames(snp.hap)%in%hap.colnames.samp) ) ]
	ind.hap <- ind.hap[ ,c( 1:grep("ALT",colnames(ind.hap)),which(colnames(ind.hap)%in%hap.colnames.samp) ) ]

	## Remove Superfluous Columns from HWE & Change colnames
	snp.hwe <- snp.hwe.l ; ind.hwe <- ind.hwe.l
	snp.hwe <- snp.hwe[ , -which(colnames(snp.hwe)%in%c("TEST","O.HET.","E.HET.")) ]
	colnames(snp.hwe)[ncol(snp.hwe)] <- "P_HW"
	ind.hwe <- ind.hwe[ , -which(colnames(ind.hwe)%in%c("TEST","O.HET.","E.HET.")) ]
	colnames(ind.hwe)[ncol(ind.hwe)] <- "P_HW"

	## Remove Superfluous Columns from BIM table
	snp.bim <- snp.bim.l ; ind.bim <- ind.bim.l
	snp.bim <- snp.bim[ , -which(colnames(snp.bim)%in%c("XXX")) ]
	ind.bim <- ind.bim[ , -which(colnames(ind.bim)%in%c("XXX")) ]

	####################################################
	## CREATE VARIANT TAGS & MERGE TABLES ##############
	print(paste( "Merging Tables",round(proc.time()-start_time,2)[3] ))

	## Tag for Hap tables
	tag.snp.hap <- paste( snp.hap[,"CHR"],snp.hap[,"BP"],sep="_" )
	tag.ind.hap <- paste( ind.hap[,"CHR"],ind.hap[,"BP"],sep="_" )
	 # Include into HAP tables
	snp.hap.2 <- data.frame( TAG=tag.snp.hap, snp.hap[,-which(colnames(snp.hap)%in%c("CHR","SNP","BP","REF","ALT"))] )
	ind.hap.2 <- data.frame( TAG=tag.ind.hap, ind.hap[,-which(colnames(ind.hap)%in%c("CHR","SNP","BP","REF","ALT"))] )

	## Tag for Raw tables
	tag.snp.raw <- colnames(snp.raw)[(grep("PHENOTYPE",colnames(snp.raw))+1):ncol(snp.raw)]
	tag.ind.raw <- colnames(ind.raw)[(grep("PHENOTYPE",colnames(ind.raw))+1):ncol(ind.raw)]

	## Add Raw Tag to BIM tables
	snp.bim <- data.frame( snp.bim, TYPE="snp", RAW_TAG=tag.snp.raw )
	ind.bim <- data.frame( ind.bim, TYPE="ind", RAW_TAG=tag.ind.raw )

	## Merge HWE & BIM tables
	snp.mg.1 <- merge( snp.bim, snp.hwe[,-which(colnames(snp.hwe)%in%c("CHR"))], by="SNP", all=T )
	ind.mg.1 <- merge( ind.bim, ind.hwe[,-which(colnames(ind.hwe)%in%c("CHR"))], by="SNP", all=T )

	## Merge MG.1 and P tables (if Pheno exists)
	snp.p <- snp.p.l ; ind.p <- ind.p.l
	if ( file.exists(PathToPheno) ) {
		colnames(snp.p)[ncol(snp.p)] <- "P_Assoc"
		snp.mg.2 <- merge( snp.mg.1, snp.p[c("SNP","P_Assoc")], by="SNP",all=T )
		colnames(ind.p)[ncol(ind.p)] <- "P_Assoc"
		ind.mg.2 <- merge( ind.mg.1, ind.p[c("SNP","P_Assoc")], by="SNP",all=T )
	}else{
		snp.mg.2 <- snp.mg.1
		ind.mg.2 <- ind.mg.1
	}

	## Merge MG.2 and RAW tables
	snp.raw.t <- snp.raw[,tag.snp.raw] ; rownames(snp.raw.t) <- snp.raw[,"IID"]
	snp.raw.t <- t(snp.raw.t)
	snp.mg.3 <- merge( snp.mg.2, snp.raw.t, by.x="RAW_TAG",by.y="row.names",all=T )
	ind.raw.t <- ind.raw[,tag.ind.raw] ; rownames(ind.raw.t) <- ind.raw[,"IID"]
	ind.raw.t <- t(ind.raw.t)
	ind.mg.3 <- merge( ind.mg.2, ind.raw.t, by.x="RAW_TAG",by.y="row.names",all=T )
	 # Tag for MG.3
	tag.snp.mg.3 <- paste( snp.mg.3$CHR,snp.mg.3$BP,sep="_" )
	tag.ind.mg.3 <- paste( ind.mg.3$CHR,ind.mg.3$BP,sep="_" )
	snp.mg.3 <- data.frame( TAG=tag.snp.mg.3, snp.mg.3 )
	ind.mg.3 <- data.frame( TAG=tag.ind.mg.3, ind.mg.3 )

	## Merge MG.3 and HAP tables
	snp.mg.4 <- merge( snp.mg.3,snp.hap.2, by="TAG", all=T )
	ind.mg.4 <- merge( ind.mg.3,ind.hap.2, by="TAG", all=T )

	####################################################
	## CLASSIFY VARIANTS ###############################
	print(paste( "Classifying Variants",round(proc.time()-start_time,2)[3] ))

	## Specify tables for Subsequent Analyses
	VAR <- rbind( snp.mg.4, ind.mg.4 )
	VAR <- VAR[ order(VAR$BP,decreasing=F), ]

	## Remove HW Violations
	VAR <- VAR[ which(VAR$P_HW>1e-20), ]

	## Compile which variants are in Exons
	TAG.exon <- c()
	for ( e in 1:ex_cnt ) {
		WHICH <- which( VAR[,"BP"]>ex_b[e] & VAR[,"BP"]<ex_e[e] )
		if ( length(WHICH)>0) { TAG.exon <- c( TAG.exon, as.character(VAR[WHICH,"TAG"]) ) }
	}

	## Specify Location of Variant (relative to Exon/Intron/Etc)
	LOC <- rep(1,nrow(VAR))
	LOC[ which( VAR[,"BP"]>=tx_rng[1] & VAR[,"BP"]<=tx_rng[2] ) ] <- 2
	LOC[ which( VAR[,"BP"]>=cd_rng[1] & VAR[,"BP"]<=cd_rng[2] ) ] <- 3
	LOC[ which( VAR[,"TAG"] %in% TAG.exon ) ] <- 4
	VAR <- data.frame( VAR, LOC )

	## Identify Potential Strand Issues **** WORK IN PROGRESS ****
	REF <- as.character( VAR$REF )
	ALT <- as.character( VAR$ALT )
	A1 <-  as.character( VAR$A1 )
	A2 <-  as.character( VAR$A2 )
	length(which(REF==A1)) / nrow(VAR)
	length(which(REF==A2)) / nrow(VAR)
	length(which(ALT==A1)) / nrow(VAR)
	length(which(ALT==A2)) / nrow(VAR)

	## Identify Discordant Calls b/n Phased & Raw tables
	 # Validate Genotype Frequencies w/ RAW and HAP tables
	MAF.raw <- apply( VAR[,which(colnames(VAR)%in%shared.samps)], 1, sum ) / (2*n.samps)
	MAF.hap <- apply( VAR[,which(colnames(VAR)%in%hap.colnames.samp)], 1, sum,na.rm=T ) / (2*n.samps)
	MAF <- data.frame( RAW=MAF.raw, HAP=MAF.hap, SUM=MAF.raw+MAF.hap, DIFF=MAF.raw-MAF.hap )
	FLIP.which.m <- which( MAF$SUM==1 | (MAF$SUM>.9 & abs(MAF$DIFF)>.1) ) # FLIP.which <- which( MAF$SUM==1 | (MAF$SUM>.8 & abs(MAF$DIFF)>.2) )
	 # Identify Discordant Calls for each Sample
	DISC.which.1 <- DISC.which.2 <- DISC.which.r <- list()
	for ( s in 1:n.samps ) {
		samp <- shared.samps[s]
		# print(samp)
		col.1 <- paste(samp,"1",sep="_")
		col.2 <- paste(samp,"2",sep="_")
		col.r <- samp
		disc.hap.1 <- which( VAR[,col.1]==1 & VAR[,col.r]==0 )
		disc.hap.2 <- which( VAR[,col.2]==1 & VAR[,col.r]==0 )
		disc.hap.r <- which( VAR[,col.r] > rowSums(VAR[,c(col.1,col.2)]) )
		DISC.which.1[[samp]] <- disc.hap.1
		DISC.which.2[[samp]] <- disc.hap.2
		DISC.which.r[[samp]] <- disc.hap.r
	}
	TAB.disc.1 <- table(unlist(DISC.which.1))
	TAB.disc.2 <- table(unlist(DISC.which.2))
	TAB.disc.r <- table(unlist(DISC.which.r))
	FLIP.which.1 <- names(TAB.disc.1)[which(TAB.disc.1>1)]
	FLIP.which.2 <- names(TAB.disc.2)[which(TAB.disc.2>1)]
	FLIP.which.r <- names(TAB.disc.r)[which(TAB.disc.r>1)]

	FLIP.tab <- table( c(FLIP.which.m,FLIP.which.1,FLIP.which.2,FLIP.which.r) )
	FLIP <- sort( as.numeric(names(FLIP.tab)[which(FLIP.tab>1)]) )

	## Flip Variants w/ Strand Issue
	VAR[FLIP,hap.colnames.samp] <- -VAR[FLIP,hap.colnames.samp] + 1
	# VAR.2 <- VAR
	# VAR.2[FLIP,hap.colnames.samp] <- -VAR.2[FLIP,hap.colnames.samp] + 1

	## Add MAF into 
	MAF.raw <- apply( VAR[,which(colnames(VAR)%in%shared.samps)], 1, sum ) / (2*n.samps)
	MAF.hap <- apply( VAR[,which(colnames(VAR)%in%hap.colnames.samp)], 1, sum,na.rm=T ) / (2*n.samps)
	MAF <- data.frame( RAW=MAF.raw, HAP=MAF.hap, SUM=MAF.raw+MAF.hap, DIFF=MAF.raw-MAF.hap )
	VAR <- data.frame( VAR, MAF_raw=MAF.raw, MAF_hap=MAF.hap, PHAS=as.numeric(!is.na(VAR[,hap.colnames.samp[4]])) )

	## Create Exon, SNP, and Indel Tables
	which.snp <- which(VAR$TYPE=="snp")
	which.ind <- which(VAR$TYPE=="ind")
	which.exon <- which(VAR$LOC==4)
	SNP <- VAR[ which.snp, ]
	IND <- VAR[ which.ind, ]
	VAR.exon <- VAR[ which.exon, ]
	SNP.exon <- VAR.exon[ intersect(which.snp,which.exon), ]
	IND.exon <- VAR.exon[ intersect(which.ind,which.exon), ]
	 # Get SNP, Indel, Exonic Var Counts
  	n.snp <- length( which.ind )
 	n.ind <- length( which.ind )
   	n.snp.exon <- length( intersect(which.snp,which.exon) )
   	n.ind.exon <- length( intersect(which.ind,which.exon) )
 

	####################################################
	## PLOT GENE MAP & GWAS RESULTS ####################
	print(paste( "Plotting Gene Map",round(proc.time()-start_time,2)[3] ))

	PLOT_RUNIF <- runif(1,0,1)
	## Set Plotting Parameters
	NUM_VARS <- nrow(VAR)
	COLS.gene <- c("steelblue1","cadetblue1","slateblue2")
	COLS.P.list <- c("firebrick4","steelblue3","cadetblue3","slateblue3")
	XLIM <- c(rng[1],rng[2])
	if ( file.exists(PathToPheno) ) {
 		YLIM <- c( -1, max( 4, max(-log10(VAR[,"P_Assoc"]),na.rm=T) ) )
 	}else{ YLIM <- c(-1,1) }
	Y <- 0
	if ( YLIM[2]>4 | PLOT_RUNIF<PLOT_FRACTION ) {
		jpeg( paste(PathToOut,"/Plots/Gene_Map_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
		plot(0,0,type="n", xlim=XLIM,ylim=YLIM, xlab=paste("Chromosome",chr,"Position"),ylab="-log10(p)", main=paste("Map & Single-Locus Results of:",name), yaxt="n" )
		## Plot P-Values
		if ( file.exists(PathToPheno) ) {
			WHICH <- which(!is.na( VAR$P_Assoc))
			COLS.P <- COLS.P.list[ VAR$LOC ]
			abline( h=seq(0,YLIM[2],1), lty=2, col="grey50" )
			abline( h=seq(-1,0,.1), lty=2, col="grey50" )
			axis(2, at=seq(-1,YLIM[2],1) )
			points( VAR$BP[WHICH], -log10(VAR$P_Assoc[WHICH]), col=COLS.P[WHICH], pch=c(10,16)[factor(VAR$TYPE[WHICH])], lwd=3, cex=2*(VAR$MAF_raw[WHICH])^(1/6) )
		}else{ WHICH <- which(VAR$MAF_raw>.01) }
		## Plot Gene Borders
		arrows( tx_rng[1],Y,tx_rng[2],Y, code=3,angle=90, lwd=5,col=COLS.gene[1] )
		abline( h=Y, col="grey50",lwd=1 )
		polygon( cd_rng[c(1,2,2,1)],c(-.1,-.1,.1,.1)+Y, col=COLS.gene[2],lwd=1 )
		for ( e in 1:ex_cnt ) {
			X_COORDS <- c( ex_b[e],ex_e[e],ex_e[e],ex_b[e] )
			polygon( X_COORDS,c(-.2,-.2,.2,.2)+Y, col=COLS.gene[3],lwd=1 )
		}
		arrows( tx_rng[1],Y,tx_rng[2],Y, code=3,angle=90,length=0, lwd=5,col=COLS.gene[1] )
		## Plot Allele Frequency
		arrows( VAR$BP[WHICH], 0, VAR$BP[WHICH], -VAR$MAF_raw[WHICH], col=COLS.P, angle=90 )
		dev.off()
	} # Close Gene Map "IF"

	## QQ Plot
	print( "Plotting QQ Plot" )
	if ( file.exists(PathToPheno) ) {
		NUM_VARS <- length(which( !is.na(VAR$P_Assoc) ))
		NUM_EXONIC <- length(which( !is.na(VAR$P_Assoc) & VAR$LOC==4 ))
		COLS <- c("steelblue1","slateblue1")
		COLS.4 <- gsub("1","3",COLS)
		LIM <- c( 0, max( 4, max(-log10(VAR[,"P_Assoc"]),na.rm=T) ) )
		## Calculate Observed, Expected, & Area for Plot
		 # Expected & Observed
		EXP.full <- -log10( 1:NUM_VARS / NUM_VARS )
		OBS.full <- -log10( sort( VAR$P_Assoc ) )
		IDX.full <- c(1:length(EXP.full),length(EXP.full),1)
		 # Calculate Area - Trapezoid Rule
		H_VALS <- EXP.full[2:NUM_VARS-1] - EXP.full[2:NUM_VARS]
		B_SUM <- ( OBS.full[2:NUM_VARS-1]-EXP.full[2:NUM_VARS-1] ) + ( OBS.full[2:NUM_VARS]-EXP.full[2:NUM_VARS] )
		AREA.full.traps <- .5*H_VALS*B_SUM
		AREA.full <- sum( AREA.full.traps ) # - .5*max(EXP.full)^2
		GENE[gtx,"AREA"] <- round( AREA.full, 5)
		GENE[gtx,"BEST_P"] <- 10^(min(-OBS.full,na.rm=T))
		if ( NUM_EXONIC > 0 ) {
			 # Expected & Observed
			EXP.exon <- -log10( 1:NUM_EXONIC / NUM_EXONIC )
			OBS.exon <- -log10( sort( VAR$P_Assoc[which(VAR$LOC==4)] ) )
			IDX.exon <- c(1:length(EXP.exon),length(EXP.exon),1)
			 # Calculate Area - Trapezoid Rule
			H_VALS <- EXP.exon[2:NUM_EXONIC-1] - EXP.exon[2:NUM_EXONIC]
			B_SUM <- ( OBS.exon[2:NUM_EXONIC-1]-EXP.exon[2:NUM_EXONIC-1] ) + ( OBS.exon[2:NUM_EXONIC]-EXP.exon[2:NUM_EXONIC] )
			AREA.exon.traps <- .5*H_VALS*B_SUM
			AREA.exon <- sum( AREA.exon.traps ) # - .5*max(EXP.exon)^2
			GENE[gtx,"AREA.ex"] <- round( AREA.exon, 5)
			GENE[gtx,"BEST_P.ex"] <- 10^(min(-OBS.exon,na.rm=T))
		}
		## Make Plot		
		if ( LIM[2]>4 | PLOT_RUNIF<PLOT_FRACTION ) {
			jpeg( paste(PathToOut,"/Plots/Gene_QQ_",name,".jpeg",sep=""), height=1500,width=1500, pointsize=32)
			plot(0,0,type="n", xlim=LIM,ylim=LIM, xlab="Expected -log10(p)",ylab="Observed -log10(p)", main=paste("QQ-Plot for:",name) )
			## Add Lines
			abline( h=seq(0,LIM[2],1), lty=2,lwd=1,col="grey50")
			abline( v=seq(0,LIM[2],1), lty=2,lwd=1,col="grey50")
			abline( 0,1, lty=1,lwd=2,col="black" )
			## Plot Full Set
			polygon( EXP.full[IDX.full], c(OBS.full,EXP.full[c(NUM_VARS,1)]), col=COLS[1], border=COLS.4[1], density=20,angle=45 )
			points( EXP.full, OBS.full, pch="+", col=COLS.4[1], type="o" )
			text( quantile(LIM,.7),quantile(LIM,.1), label=paste("Area:",round(AREA.full,3),"-",NUM_VARS,"Vars"), col=COLS.4[1], cex=1.2 )
			## Plot Exon Set
			if ( NUM_EXONIC > 0 ) {
				polygon( EXP.exon[IDX.exon], c(OBS.exon,EXP.exon[c(NUM_EXONIC,1)]), col=COLS[2], border=COLS.4[2], density=20,angle=-45 )
				points( EXP.exon, OBS.exon, pch="+", col=COLS.4[2], type="o" )	
				text( quantile(LIM,.7),quantile(LIM,.07), label=paste("Area:",round(AREA.exon,3),"-",NUM_EXONIC,"Vars"), col=COLS.4[2], cex=1.2 )
			}
			legend("topleft",fill=COLS,density=20,legend=c("Full","Exonic") )
			dev.off()
		} # Close QQ "IF"
	}

	####################################################
	## COMPILE STATS ABOUT GENE TRANSCRIPT #############
	print(paste( "Compiling Gene Stats",round(proc.time()-start_time,2)[3] ))
	 # By Location
	   # Ovarall
 	   # Exonic
 	   # Intronic
	   # Damaging (?)

 	## Compile Info about the Gene_Transcript
 	 # Full
 	GENE[gtx,"n.VAR"] <- nrow(VAR)
 	GENE[gtx,"n.SNP"] <- length(which( VAR$TYPE=="snp" ))
 	GENE[gtx,"n.IND"] <- length(which( VAR$TYPE=="ind" ))
 	GENE[gtx,"n.VAR.ph"] <- length(which(VAR$PHAS==1))
 	GENE[gtx,"n.SNP.ph"] <- length(which(VAR$PHAS==1 & VAR$TYPE=="snp" ))
 	GENE[gtx,"n.IND.ph"] <- length(which(VAR$PHAS==1 & VAR$TYPE=="ind" ))
 	 # Exon
 	GENE[gtx,"n.VAR.ex"] <- nrow(VAR.exon)
 	GENE[gtx,"n.SNP.ex"] <- length(which( VAR.exon$TYPE=="snp" ))
 	GENE[gtx,"n.IND.ex"] <- length(which( VAR.exon$TYPE=="ind" ))
 	GENE[gtx,"n.VAR.ph.ex"] <- length(which(VAR.exon$PHAS==1))
 	GENE[gtx,"n.SNP.ph.ex"] <- length(which(VAR.exon$PHAS==1 & VAR.exon$TYPE=="snp" ))
 	GENE[gtx,"n.IND.ph.ex"] <- length(which(VAR.exon$PHAS==1 & VAR.exon$TYPE=="ind" ))

	####################################################
	## COMPILE STATS ABOUT COHORT ######################
	print(paste( "Compiling Cohort Stats",round(proc.time()-start_time,2)[3] ))
	 # What to compile
	   # % Phased
 	   # How many hets, hom_vars pp (SNPs & Indels separate)
 	   # How many comp hets pp

	# ## Cohort - Full ###################################
 # 	colnames.FL <- c("HET_snp","HOM_VAR_snp","N_PHAS_snp","PRC_PHAS_snp","HAP1_snp","HAP2_snp","UNPH_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","N_PHAS_ind","PRC_PHAS_ind","HAP1_ind","HAP2_ind","UNPH_ind","COMP_HET_ind","COMP_HET_any")
 # 	FULL[[name]] <- array( ,c(n.samps,length(colnames.FL)) )
 # 	rownames(FULL[[name]]) <- shared.samps
 # 	colnames(FULL[[name]]) <- colnames.FL
 # 	## How many Hets & Hom_Vars?
 # 	if ( n.snp>0 ) {
 # 		FULL[[name]][,"HET_snp"] <- apply( SNP[,shared.samps], 2, function(x) length(which(x==1)) )
 # 		FULL[[name]][,"HOM_VAR_snp"] <- apply( SNP[,shared.samps], 2, function(x) length(which(x==2)) )
 # 	}else{ FULL[[name]][,"HET_snp"] <- FULL[[name]][,"HOM_VAR_snp"] <- 0 }
 # 	if ( n.ind>0 ) {
 # 		FULL[[name]][,"HET_ind"] <- apply( IND[,shared.samps], 2, function(x) length(which(x==1)) )
 # 		FULL[[name]][,"HOM_VAR_ind"] <- apply( IND[,shared.samps], 2, function(x) length(which(x==2)) )
 # 	}else{ FULL[[name]][,"HET_ind"] <- FULL[[name]][,"HOM_VAR_ind"] <- 0 }

	# ## What percent of variants were phased?
	# hap.1.cols <- hap.colnames.samp[seq(1,2*n.samps,2)]
	# hap.2.cols <- hap.colnames.samp[seq(2,2*n.samps,2)]
 # 	 # SNP
 # 	hap.snp.diffs <- SNP[,hap.1.cols] - SNP[,hap.2.cols]
 # 	# colnames(hap.snp.diffs) <- gsub( "_1","", colnames(hap.snp.diffs) )
	# hap.snp.n.hets <- apply( hap.snp.diffs, 2, function(x) length(which(x!=0)) )
 # 	FULL[[name]][,"N_PHAS_snp"] <- hap.snp.n.hets
 # 	FULL[[name]][,"PRC_PHAS_snp"] <- FULL[[name]][,"N_PHAS_snp"] / FULL[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.n.hets / FULL[[name]][,"HET_snp"]
 # 	FULL[[name]][,"HAP1_snp"] <- colSums( SNP[,hap.1.cols],na.rm=T )
 # 	FULL[[name]][,"HAP2_snp"] <- colSums( SNP[,hap.2.cols],na.rm=T )
 # 	FULL[[name]][,"UNPH_snp"] <- FULL[[name]][,"HET_snp"] - FULL[[name]][,"N_PHAS_snp"]
 # 	 # Indel
 # 	hap.ind.diffs <- IND[,hap.1.cols] - IND[,hap.2.cols]
 # 	# colnames(hap.ind.diffs) <- gsub( "_1","", colnames(hap.ind.diffs) )
	# hap.ind.n.hets <- apply( hap.ind.diffs, 2, function(x) length(which(x!=0)) )
 # 	FULL[[name]][,"N_PHAS_ind"] <- hap.ind.n.hets
 # 	FULL[[name]][,"PRC_PHAS_ind"] <- FULL[[name]][,"N_PHAS_ind"] / FULL[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.n.hets / FULL[[name]][,"HET_ind"]
 # 	FULL[[name]][,"HAP1_ind"] <- colSums( SNP[,hap.1.cols],na.rm=T )
 # 	FULL[[name]][,"HAP2_ind"] <- colSums( SNP[,hap.2.cols],na.rm=T )
 # 	FULL[[name]][,"UNPH_ind"] <- FULL[[name]][,"HET_ind"] - FULL[[name]][,"N_PHAS_ind"]
 # 	 # Which Samples have a Compound Het?
 # 	FULL[[name]][,"COMP_HET_snp"] <- as.numeric( apply( hap.snp.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 # 	FULL[[name]][,"COMP_HET_ind"] <- as.numeric( apply( hap.ind.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 #  	FULL[[name]][,"COMP_HET_any"] <- as.numeric( apply( hap.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )

	## Cohort - Exon ###################################
 	colnames.FL <- c("HET_snp","HOM_VAR_snp","N_PHAS_snp","PRC_PHAS_snp","HAP1_snp","HAP2_snp","UNPH_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","N_PHAS_ind","PRC_PHAS_ind","HAP1_ind","HAP2_ind","UNPH_ind","COMP_HET_ind","COMP_HET_any")
 	FULL[[name]] <- array( ,c(n.samps,length(colnames.FL)) )
 	rownames(FULL[[name]]) <- shared.samps
 	colnames(FULL[[name]]) <- colnames.FL
 	## SNP Stats
 	if ( n.snp>0 ) {
		# How many Hets & Hom_Vars?
 		FULL[[name]][,"HET_snp"] <- apply( SNP[,shared.samps], 2, function(x) length(which(x==1)) )
 		FULL[[name]][,"HOM_VAR_snp"] <- apply( SNP[,shared.samps], 2, function(x) length(which(x==2)) )
		# What percent of variants were phased?
		hap.1.cols <- hap.colnames.samp[seq(1,2*n.samps,2)]
		hap.2.cols <- hap.colnames.samp[seq(2,2*n.samps,2)]
	 	hap.snp.diffs <- SNP[,hap.1.cols] - SNP[,hap.2.cols]
		hap.snp.n.hets <- apply( hap.snp.diffs, 2, function(x) length(which(x!=0)) )
	 	FULL[[name]][,"N_PHAS_snp"] <- hap.snp.n.hets
	 	FULL[[name]][,"PRC_PHAS_snp"] <- FULL[[name]][,"N_PHAS_snp"] / FULL[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.n.hets / FULL[[name]][,"HET_snp"]
	 	FULL[[name]][,"HAP1_snp"] <- colSums( SNP[,hap.1.cols],na.rm=T )
	 	FULL[[name]][,"HAP2_snp"] <- colSums( SNP[,hap.2.cols],na.rm=T )
	 	FULL[[name]][,"UNPH_snp"] <- FULL[[name]][,"HET_snp"] - FULL[[name]][,"N_PHAS_snp"]
 	}else{
 		FULL[[name]][,"HET_snp"] <- FULL[[name]][,"HOM_VAR_snp"] <- FULL[[name]][,"N_PHAS_snp"] <- FULL[[name]][,"HAP1_snp"] <- FULL[[name]][,"HAP2_snp"] <- FULL[[name]][,"UNPH_snp"] <- 0
 		FULL[[name]][,"PRC_PHAS_snp"] <- NA
	}
	## Indel Stats
 	if ( n.ind>0 ) {
		# How many Hets & Hom_Vars?
 		FULL[[name]][,"HET_ind"] <- apply( IND[,shared.samps], 2, function(x) length(which(x==1)) )
 		FULL[[name]][,"HOM_VAR_ind"] <- apply( IND[,shared.samps], 2, function(x) length(which(x==2)) )
		# What percent of variants were phased?
		hap.1.cols <- hap.colnames.samp[seq(1,2*n.samps,2)]
		hap.2.cols <- hap.colnames.samp[seq(2,2*n.samps,2)]
	 	hap.ind.diffs <- IND[,hap.1.cols] - IND[,hap.2.cols]
		hap.ind.n.hets <- apply( hap.ind.diffs, 2, function(x) length(which(x!=0)) )
	 	FULL[[name]][,"N_PHAS_ind"] <- hap.ind.n.hets
	 	FULL[[name]][,"PRC_PHAS_ind"] <- FULL[[name]][,"N_PHAS_ind"] / FULL[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.n.hets / FULL[[name]][,"HET_ind"]
	 	FULL[[name]][,"HAP1_ind"] <- colSums( IND[,hap.1.cols],na.rm=T )
	 	FULL[[name]][,"HAP2_ind"] <- colSums( IND[,hap.2.cols],na.rm=T )
	 	FULL[[name]][,"UNPH_ind"] <- FULL[[name]][,"HET_ind"] - FULL[[name]][,"N_PHAS_ind"]
 	}else{
 		FULL[[name]][,"HET_ind"] <- FULL[[name]][,"HOM_VAR_ind"] <- FULL[[name]][,"N_PHAS_ind"] <- FULL[[name]][,"HAP1_ind"] <- FULL[[name]][,"HAP2_ind"] <- FULL[[name]][,"UNPH_ind"] <- 0
 		FULL[[name]][,"PRC_PHAS_ind"] <- NA
	}

	## Cohort - Exon ###################################
 	colnames.EX <- c("HET_snp","HOM_VAR_snp","N_PHAS_snp","PRC_PHAS_snp","HAP1_snp","HAP2_snp","UNPH_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","N_PHAS_ind","PRC_PHAS_ind","HAP1_ind","HAP2_ind","UNPH_ind","COMP_HET_ind","COMP_HET_any")
 	EXON[[name]] <- array( ,c(n.samps,length(colnames.EX)) )
 	rownames(EXON[[name]]) <- shared.samps
 	colnames(EXON[[name]]) <- colnames.EX
 	## SNP Stats
 	if ( n.snp.exon>0 ) {
		# How many Hets & Hom_Vars?
 		EXON[[name]][,"HET_snp"] <- apply( SNP.exon[,shared.samps], 2, function(x) length(which(x==1)) )
 		EXON[[name]][,"HOM_VAR_snp"] <- apply( SNP.exon[,shared.samps], 2, function(x) length(which(x==2)) )
		# What percent of variants were phased?
		hap.1.cols <- hap.colnames.samp[seq(1,2*n.samps,2)]
		hap.2.cols <- hap.colnames.samp[seq(2,2*n.samps,2)]
	 	hap.snp.diffs <- SNP.exon[,hap.1.cols] - SNP.exon[,hap.2.cols]
		hap.snp.n.hets <- apply( hap.snp.diffs, 2, function(x) length(which(x!=0)) )
	 	EXON[[name]][,"N_PHAS_snp"] <- hap.snp.n.hets
	 	EXON[[name]][,"PRC_PHAS_snp"] <- EXON[[name]][,"N_PHAS_snp"] / EXON[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.n.hets / EXON[[name]][,"HET_snp"]
	 	EXON[[name]][,"HAP1_snp"] <- colSums( SNP.exon[,hap.1.cols],na.rm=T )
	 	EXON[[name]][,"HAP2_snp"] <- colSums( SNP.exon[,hap.2.cols],na.rm=T )
	 	EXON[[name]][,"UNPH_snp"] <- EXON[[name]][,"HET_snp"] - EXON[[name]][,"N_PHAS_snp"]
 	}else{
 		EXON[[name]][,"HET_snp"] <- EXON[[name]][,"HOM_VAR_snp"] <- EXON[[name]][,"N_PHAS_snp"] <- EXON[[name]][,"HAP1_snp"] <- EXON[[name]][,"HAP2_snp"] <- EXON[[name]][,"UNPH_snp"] <- 0
 		EXON[[name]][,"PRC_PHAS_snp"] <- NA
	}
	## Indel Stats
 	if ( n.ind.exon>0 ) {
		# How many Hets & Hom_Vars?
 		EXON[[name]][,"HET_ind"] <- apply( IND.exon[,shared.samps], 2, function(x) length(which(x==1)) )
 		EXON[[name]][,"HOM_VAR_ind"] <- apply( IND.exon[,shared.samps], 2, function(x) length(which(x==2)) )
		# What percent of variants were phased?
		hap.1.cols <- hap.colnames.samp[seq(1,2*n.samps,2)]
		hap.2.cols <- hap.colnames.samp[seq(2,2*n.samps,2)]
	 	hap.ind.diffs <- IND.exon[,hap.1.cols] - IND.exon[,hap.2.cols]
		hap.ind.n.hets <- apply( hap.ind.diffs, 2, function(x) length(which(x!=0)) )
	 	EXON[[name]][,"N_PHAS_ind"] <- hap.ind.n.hets
	 	EXON[[name]][,"PRC_PHAS_ind"] <- EXON[[name]][,"N_PHAS_ind"] / EXON[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.n.hets / EXON[[name]][,"HET_ind"]
	 	EXON[[name]][,"HAP1_ind"] <- colSums( IND.exon[,hap.1.cols],na.rm=T )
	 	EXON[[name]][,"HAP2_ind"] <- colSums( IND.exon[,hap.2.cols],na.rm=T )
	 	EXON[[name]][,"UNPH_ind"] <- EXON[[name]][,"HET_ind"] - EXON[[name]][,"N_PHAS_ind"]
 	}else{
 		EXON[[name]][,"HET_ind"] <- EXON[[name]][,"HOM_VAR_ind"] <- EXON[[name]][,"N_PHAS_ind"] <- EXON[[name]][,"HAP1_ind"] <- EXON[[name]][,"HAP2_ind"] <- EXON[[name]][,"UNPH_ind"] <- 0
 		EXON[[name]][,"PRC_PHAS_ind"] <- NA
	}

	####################################################
	## PLOT PHASING STATS ##############################
	print(paste( "Plot Phased Stats",round(proc.time()-start_time,2)[3] ))
	 # FULL: Boxplot HET_snp & HOM_VAR_snp & HET_ind & HOM_VAR_ind all in one
	 # EXON: Boxplot HET_snp & HOM_VAR_snp & HET_ind & HOM_VAR_ind all in one
	 # FULL: PRC_PHAS_snp & PRC_PHAS_ind & PRC_COMP_HET_snp & PRC_COMP_HET_ind
	 # EXON: PRC_PHAS_snp & PRC_PHAS_ind & PRC_COMP_HET_snp & PRC_COMP_HET_ind

	## Calculate % Compound Hets
	 # Full
	PRC.dat.full <- FULL[[name]][,c("PRC_PHAS_snp","PRC_PHAS_ind")]
	PRC_COMP_HET.full <- array(,c(2,2)) ; rownames(PRC_COMP_HET.full) <- c(0,1) ; colnames(PRC_COMP_HET.full) <- c("snp","ind")
	PRC_COMP_HET.full[,"snp"] <- c( length(which(FULL[[name]][,"COMP_HET_snp"]==0)),length(which(FULL[[name]][,"COMP_HET_snp"]==1)) ) / n.samps
	PRC_COMP_HET.full[,"ind"] <- c( length(which(FULL[[name]][,"COMP_HET_ind"]==0)),length(which(FULL[[name]][,"COMP_HET_ind"]==1)) ) / n.samps
	 # Exon
	PRC.dat.exon <- EXON[[name]][,c("PRC_PHAS_snp","PRC_PHAS_ind")]
	PRC_COMP_HET.exon <- array(,c(2,2)) ; rownames(PRC_COMP_HET.exon) <- c(0,1) ; colnames(PRC_COMP_HET.exon) <- c("snp","ind")
	PRC_COMP_HET.exon[,"snp"] <- c( length(which(EXON[[name]][,"COMP_HET_snp"]==0)),length(which(EXON[[name]][,"COMP_HET_snp"]==1)) ) / n.samps
	PRC_COMP_HET.exon[,"ind"] <- c( length(which(EXON[[name]][,"COMP_HET_ind"]==0)),length(which(EXON[[name]][,"COMP_HET_ind"]==1)) ) / n.samps

	## Compile Percent w/ Compound Het
	GENE[gtx,"p.COMP_HET_snp"] <- round( PRC_COMP_HET.full["1","snp"] ,5)
	GENE[gtx,"p.COMP_HET_ind"] <- round( PRC_COMP_HET.full["1","ind"] ,5)
	GENE[gtx,"p.COMP_HET_snp.ex"] <- round( PRC_COMP_HET.exon["1","snp"] ,5)
	GENE[gtx,"p.COMP_HET_ind.ex"] <- round( PRC_COMP_HET.exon["1","ind"] ,5)

	if ( LIM[2]>4 | PLOT_RUNIF<PLOT_FRACTION ) {
		jpeg( paste(PathToOut,"/Plots/Gene_Stats_",name,".jpeg",sep=""), height=1600,width=2400, pointsize=32)
		par(mfrow=c(2,2))
		 # FULL: Boxplot
		COLS.cnt <- c("springgreen3","springgreen1","gold3","gold1")
		COLS.prc <- c("chocolate2","steelblue2")
		CNT.dat <- FULL[[name]][,c("HET_snp","HOM_VAR_snp","HET_ind","HOM_VAR_ind")]
		boxplot( CNT.dat, main=paste("Number Variants:",name), xlab="Category",ylab="# Vars in Individual", col=COLS.cnt, pch="" )
		for ( i in 1:4 ) { points( jitter(rep(i,nrow(CNT.dat)),amount=.1), CNT.dat[,i], pch="+" ) }
		 # EXON: Boxplot
		CNT.dat <- EXON[[name]][,c("HET_snp","HOM_VAR_snp","HET_ind","HOM_VAR_ind")]
		boxplot( CNT.dat, main=paste("Number Exonic Variants:",name), xlab="Category",ylab="# Exonic Vars in Individual", col=COLS.cnt, pch="" )
		for ( i in 1:4 ) { points( jitter(rep(i,nrow(CNT.dat)),amount=.1), CNT.dat[,i], pch="+" ) }
		 # FULL: Perc
		plot(0,0,type="n", xlim=c(0,4.5),ylim=c(0,1), main=paste("Percent Phased & Compound Hets:",name), xlab="",ylab="", xaxt="n",yaxt="n")
		abline( h=seq(0,1,.1), lty=2,lwd=1,col="grey50" )
		axis( 2, at=seq(0,1,.2) )
		axis( 4, at=seq(0,1,.2) )
		axis( 1, at=c(.5,1.5,3,4), labels=c("%Ph-SNP","%Ph-IND","%CH-SNP","%CH-IND") )
		boxplot( PRC.dat.full, at=c(.5,1.5), xaxt="n",yaxt="n", add=T, col=COLS.cnt[c(1,3)], pch="" )
		for ( i in 1:2 ) { points( jitter(rep(i-.5,nrow(PRC.dat.full)),amount=.1), PRC.dat.full[,i], pch="+" ) }
		barplot( PRC_COMP_HET.full, width=.8, space=c(3.25,.25), add=T, xaxt="n",yaxt="n", col=COLS.prc )
		abline( v=2.25 )
		 # EXON: Perc
		plot(0,0,type="n", xlim=c(0,4.5),ylim=c(0,1), main=paste("Percent Phased & Compound Hets:",name), xlab="",ylab="", xaxt="n",yaxt="n")
		abline( h=seq(0,1,.1), lty=2,lwd=1,col="grey50" )
		axis( 2, at=seq(0,1,.2) )
		axis( 4, at=seq(0,1,.2) )
		axis( 1, at=c(.5,1.5,3,4), labels=c("%Ph-SNP","%Ph-IND","%CH-SNP","%CH-IND") )
		boxplot( PRC.dat.exon, at=c(.5,1.5), xaxt="n",yaxt="n", add=T, col=COLS.cnt[c(1,3)], pch="" )
		for ( i in 1:2 ) { points( jitter(rep(i-.5,nrow(PRC.dat.exon)),amount=.1), PRC.dat.exon[,i], pch="+" ) }
		barplot( PRC_COMP_HET.exon, width=.8, space=c(3.25,.25), add=T, xaxt="n",yaxt="n", col=COLS.prc )
		abline( v=2.25 )
		legend( "bottomright",fill=COLS.prc,legend=c("No Comp_Het","Comp_Het") )
		dev.off()
	} # Close Compile Plot "IF"

	####################################################
	## UNIQUE HAPLOTYPES ###############################
	print(paste( "Analyzing Unique Haplotypes",round(proc.time()-start_time,2)[3] ))
	## Full
	UNIQ.full.all <- apply( VAR[which(VAR$PHAS==1),hap.colnames.samp], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.full.all <- table( UNIQ.full.all )
	 # MAF .01
	UNIQ.full.maf.01 <- apply( VAR[which(VAR$PHAS==1 & VAR$MAF_raw>.01),hap.colnames.samp], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.full.maf.01 <- table( UNIQ.full.maf.01 )
	 # MAF .05
	UNIQ.full.maf.05 <- apply( VAR[which(VAR$PHAS==1 & VAR$MAF_raw>.05),hap.colnames.samp], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.full.maf.05 <- table( UNIQ.full.maf.05 )
	## Exon
	UNIQ.exon.all <- apply( VAR.exon[which(VAR.exon$PHAS==1),hap.colnames.samp], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.exon.all <- table( UNIQ.exon.all )
	 # MAF .01
	UNIQ.exon.maf.01 <- apply( VAR.exon[which(VAR.exon$PHAS==1 & VAR.exon$MAF_raw>.01),hap.colnames.samp], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.exon.maf.01 <- table( UNIQ.exon.maf.01 )
	 # MAF .05
	UNIQ.exon.maf.05 <- apply( VAR.exon[which(VAR.exon$PHAS==1 & VAR.exon$MAF_raw>.05),hap.colnames.samp], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.exon.maf.05 <- table( UNIQ.exon.maf.05 )
	## Plot it
	COLS.full <- colorRampPalette(c("white","steelblue2","black"))(5)[2:4]
	COLS.exon <- colorRampPalette(c("white","slateblue2","black"))(5)[2:4]
	if ( PLOT_RUNIF<PLOT_FRACTION ) {
		## Barplot Haplotype Frequencies
		jpeg( paste(PathToOut,"/Plots/Gene_HaploDistrib_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
		par(mfrow=c(2,3))
		 # Full
		barplot( sort(TAB.uniq.full.all,decreasing=T), main="Unique Haplotype Freq (Full/All)",xlab="Haplotype",las=2,col=COLS.full[1],border=NA,xaxt="n" )
		barplot( sort(TAB.uniq.full.maf.01,decreasing=T), main="Unique Haplotype Freq (Full/MAF>1%)",xlab="Haplotype",las=2,col=COLS.full[2],border=NA,xaxt="n" )
		barplot( sort(TAB.uniq.full.maf.05,decreasing=T), main="Unique Haplotype Freq (Full/MAF>5%)",xlab="Haplotype",las=2,col=COLS.full[3],border=NA,xaxt="n" )
		 # Exon
		barplot( sort(TAB.uniq.exon.all,decreasing=T), main="Unique Haplotype Freq (Exon/All)",xlab="Haplotype",las=2,col=COLS.exon[1],border=NA )
		barplot( sort(TAB.uniq.exon.maf.01,decreasing=T), main="Unique Haplotype Freq (Exon/1%)",xlab="Haplotype",las=2,col=COLS.exon[2],border=NA )
		barplot( sort(TAB.uniq.exon.maf.05,decreasing=T), main="Unique Haplotype Freq (Exon/5%)",xlab="Haplotype",las=2,col=COLS.exon[3],border=NA )
		dev.off()
		## Heatmap Showing Unique Haplotypes
		HAPLOS.arr <- matrix(as.numeric(unlist(sapply( names(TAB.uniq.exon.all), function(x) strsplit( x, ""), simplify="array" ))),nrow=length(TAB.uniq.exon.all),byrow=T )
		if ( all(dim(HAPLOS.arr)>=2) ) {
			colnames(HAPLOS.arr) <- VAR.exon$TAG[which(VAR.exon$PHAS==1)]
			MAFS.cols <- colorRampPalette(c("white","firebrick2","black"))(100)[ceiling(100*VAR.exon$MAF_raw[which(VAR.exon$PHAS==1)])]
			HAPLOS.cols <- colorRampPalette(c("white","chartreuse3"))(max(TAB.uniq.exon.all))[TAB.uniq.exon.all]
			jpeg( paste(PathToOut,"/Plots/Gene_HaploUniq_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
			heatmap.2( HAPLOS.arr, main="Unique Exonic Haplotypes",xlab="Variant Position",ylab="Haplotype",RowSideColors=HAPLOS.cols,ColSideColors=MAFS.cols,scale="none",trace="none",Rowv=T,Colv=F,dendrogram="row",col=c("black",COLS.exon[1]), lhei=c(1,6),lwid=c(1,6),margins=c(7,5) )
			dev.off()	
		}
	}
	
	####################################################
	## ASSOCIATION w/ PHENOTYPE ########################
	print(paste( "Analyzing Phenotype Associations",round(proc.time()-start_time,2)[3] ))
	if ( file.exists(PathToPheno) ) {
		MOD.covs <- lm( Pheno ~ . , data=PC[,-1] )
		RES.covs <- resid( MOD.covs )
		RES.covs.2 <- data.frame( IID=PC[,1],RES=RES.covs )

		## Unique Haplotypes vs Phenotype ##
		if ( length(TAB.uniq.exon.all)>1 ) {

			## Get Haplotypes for each Person
			HAPLOS.samp <- cbind( UNIQ.exon.all[seq(1,2*n.samps,2)], UNIQ.exon.all[seq(2,2*n.samps,2)] )
			HAPLOS.samp <- t(apply( HAPLOS.samp, 1, function(x) sort(x) ))
			# Put Clinical & Haplotype Data Together
			HAPLOS.samp.2 <- data.frame( IID=rep(shared.samps,2), HAP=c(HAPLOS.samp[,1],HAPLOS.samp[,2]), stringsAsFactors=F )
			PC.haps <- merge( PC, HAPLOS.samp.2, by="IID" )
			RES.covs.3 <- merge( RES.covs.2, HAPLOS.samp.2, by="IID" )
			# Specify "Rare" Haplotypes
			HAPLOS.rare <- names( which(TAB.uniq.exon.all<=10) )
			RES.covs.3[which(RES.covs.3[,"HAP"] %in% HAPLOS.rare),"HAP"] <- "Rare"
			PC.haps[which(PC.haps[,"HAP"] %in% HAPLOS.rare),"HAP"] <- "Rare"
			# Run Analyses
			MOD.haps <- lm( Pheno ~ . , data=PC.haps[,-1] )
			P.haps <- anova(MOD.haps)["HAP","Pr(>F)"]
			GENE[gtx,"P_HAP"] <- P.haps
			# Plot Residuals from Covariates vs Haplotype
			jpeg( paste(PathToOut,"/Plots/Gene_HaploAssoc_",name,".jpeg",sep=""), height=1700,width=2000, pointsize=30)
			par(mfrow=c(2,1))
			TEMP <- boxplot( RES ~ factor(HAP), data=RES.covs.3, col=COLS.exon[1],main=paste("Residuals vs Haplotype:",name),xlab="Haplotypes",ylab="Phenotype Residuals vs Covariates",xaxt="n",pch="" )
			abline( h=seq(-5,5,1),lty=2,col="grey50" )
			TEMP <- boxplot( RES ~ factor(HAP), data=RES.covs.3, col=COLS.exon[1],main=paste("Residuals vs Haplotype:",name),xlab="Haplotypes",ylab="Phenotype Residuals vs Covariates",xaxt="n",pch="",add=T )
			points( RES ~ factor(HAP), data=RES.covs.3, pch="+" )
			text( 2,-2, label=paste("P=",formatC(P.haps,format="e",digits=2)) )
			barplot( TEMP$n, names.arg=TEMP$names, col=COLS.exon[1],main=paste("Haplotype Frequency:",name),xlab="Haplotypes",ylab="Frequency",las=2)
			dev.off()
		}

		## Compound Hets vs Phenotype ##
		PC.ch <- merge( PC, EXON[[name]], by.x="IID",by.y="row.names" )
		PC.ch.2 <- PC.ch[,c(colnames(PC),"COMP_HET_any")]
		if ( length(unique( PC.ch.2[,"COMP_HET_any"] ))>1 ) {
			MOD.ch <- lm( Pheno ~ . , data=PC.ch.2[,-1] )
			P.ch <- anova(MOD.ch)["COMP_HET_any","Pr(>F)"]
			GENE[gtx,"P_CH"] <- P.ch
			RES.ch.3 <- merge( RES.covs.2, PC.ch.2[,c("IID","COMP_HET_any")], by="IID" )
			MOD.ch.res <- lm( RES ~ COMP_HET_any , data=RES.ch.3 )
			if ( PLOT_RUNIF<PLOT_FRACTION ) {
				# Plot Residuals from Covariates vs Haplotype
				jpeg( paste(PathToOut,"/Plots/Gene_CHAssoc_",name,".jpeg",sep=""), height=1200,width=1000, pointsize=30)
				TEMP <- boxplot( RES ~ factor(COMP_HET_any), data=RES.ch.3, col=COLS.exon[1],main=paste("Residuals vs CompHet Status:",name),xlab="Compound Het Status",ylab="Phenotype Residuals vs Covariates",pch="" )
				abline( h=seq(-5,5,1),lty=2,col="grey50" )
				TEMP <- boxplot( RES ~ factor(COMP_HET_any), data=RES.ch.3, col=COLS.exon[1],main=paste("Residuals vs CompHet Status:",name),xlab="Compound Het Status",ylab="Phenotype Residuals vs Covariates",pch="",add=T )
				points( RES ~ factor(COMP_HET_any), data=RES.ch.3, pch="+" )
				text( 1.5,-2, label=paste("P=",formatC(P.ch,format="e",digits=2)) )
				dev.off()
			}
		}

		## Burden vs Phenotype ##
		 # Full
		PC.burd.full <- merge( PC, FULL[[name]], by.x="IID",by.y="row.names" )
		SCORE.var <- PC.burd.full[,"HET_snp"] + PC.burd.full[,"HOM_VAR_snp"] + PC.burd.full[,"HET_ind"] + PC.burd.full[,"HOM_VAR_ind"]
		PC.burd.full.2 <- cbind( PC.burd.full[,colnames(PC)], SCORE.var )
		if ( length(unique( PC.burd.full.2[,"SCORE.var"] ))>1 ) {
			MOD.burd.full <- lm( Pheno ~ . , data=PC.burd.full.2[,-1] )
			P.burd.full <- anova(MOD.burd.full)["SCORE.var","Pr(>F)"]
			GENE[gtx,"P_BURD"] <- P.burd.full
			RES.burd.full.3 <- merge( RES.covs.2, PC.burd.full.2[,c("IID","SCORE.var")], by="IID" )
			MOD.burd.full.res <- lm( RES ~ as.numeric(SCORE.var) , data=RES.burd.full.3 )
			if ( PLOT_RUNIF<PLOT_FRACTION ) {
				# Plot Residuals from Covariates vs Burden Stat
				jpeg( paste(PathToOut,"/Plots/Gene_BurdFullAssoc_",name,".jpeg",sep=""), height=1200,width=1200, pointsize=30)
				plot( RES ~ as.numeric(SCORE.var), data=RES.burd.full.3, col=COLS.full[1],main=paste("Residuals vs Burden Score (Full):",name),xlab="Burden Score",ylab="Phenotype Residuals vs Covariates",pch="+" )
				abline( h=seq(-5,5,1),lty=2,col="grey50" )
				abline( MOD.burd.full.res, lty=1,lwd=2,col=COLS.full[3] )
				text( quantile(RES.burd.full.3$SCORE.var,.5),quantile(RES.burd.full.3$RES,.01), label=paste("P=",formatC(P.burd.full,format="e",digits=2)) )
				dev.off()
			}
		}
		 # Exon
		PC.burd.exon <- merge( PC, EXON[[name]], by.x="IID",by.y="row.names" )
		SCORE.var <- PC.burd.exon[,"HET_snp"] + PC.burd.exon[,"HOM_VAR_snp"] + PC.burd.exon[,"HET_ind"] + PC.burd.exon[,"HOM_VAR_ind"]
		PC.burd.exon.2 <- cbind( PC.burd.exon[,colnames(PC)], SCORE.var )
		if ( length(unique( PC.burd.exon.2[,"SCORE.var"] ))>1 ) {
			MOD.burd.exon <- lm( Pheno ~ . , data=PC.burd.exon.2[,-1] )
			P.burd.exon <- anova(MOD.burd.exon)["SCORE.var","Pr(>F)"]
			GENE[gtx,"P_BURD.ex"] <- P.burd.exon
			RES.burd.exon.3 <- merge( RES.covs.2, PC.burd.exon.2[,c("IID","SCORE.var")], by="IID" )
			MOD.burd.exon.res <- lm( RES ~ as.numeric(SCORE.var) , data=RES.burd.exon.3 )
			if ( PLOT_RUNIF<PLOT_FRACTION ) {
				# Plot Residuals from Covariates vs Burden Stat
				jpeg( paste(PathToOut,"/Plots/Gene_BurdExonAssoc_",name,".jpeg",sep=""), height=1200,width=1200, pointsize=30)
				TEMP <- plot( RES ~ as.numeric(SCORE.var), data=RES.burd.exon.3, col=COLS.exon[1],main=paste("Residuals vs Burden Score (Exon):",name),xlab="Burden Score",ylab="Phenotype Residuals vs Covariates",pch="+" )
				abline( h=seq(-5,5,1),lty=2,col="grey50" )
				abline( MOD.burd.exon.res, lty=1,lwd=2,col=COLS.exon[3] )
				text( quantile(RES.burd.exon.3$SCORE.var,.5),quantile(RES.burd.exon.3$RES,.01), label=paste("P=",formatC(P.burd.exon,format="e",digits=2)) )
				dev.off()
			}
		}
	}

	####################################################
	## EVERY FEW ITERATIONS, SAVE TABLES/DATA ##########
	if ( gtx%%50==0 ) {
		print(paste( "Writing Updated Output Tables",round(proc.time()-start_time,2)[3] ))
		## Save Table that Gene Data are Compiled In
		write.table( GENE[1:gtx,], paste(PathToOut,"/Gene_Stats.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

		## Save List of FULL/EXON Cohort data
		COMPILE <- list( FULL, EXON )
		names(COMPILE) <- c("Full","Exon")
		save( COMPILE, file=paste(PathToOut,"/Gene_Stats.Rdata",sep="") )	
	}
	## Moving on...

} # Close GTX Loop

###############################################################
## SAVE COMPILED DATA #########################################
###############################################################

print(paste( "Writing Final Output Tables",round(proc.time()-start_time,2)[3] ))

## Save Table that Gene Data are Compiled In
write.table( GENE, paste(PathToOut,"/Gene_Stats.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )
GENE <- read.table( paste(PathToOut,"/Gene_Stats.txt",sep=""),sep="\t",header=T, colClasses="character" )
GENE.2 <- array( , dim(GENE) )
colnames(GENE.2) <- colnames(GENE)
for ( col in 1:ncol(GENE.2) ) { GENE.2[,col] <- GENE[,col] }
GENE <- GENE.2

## Save List of FULL/EXON Cohort data
COMPILE <- list( FULL, EXON )
names(COMPILE) <- c("Full","Exon")
save( COMPILE, file=paste(PathToOut,"/Gene_Stats.Rdata",sep="") )

###############################################################
## PLOT COMPILED STATS ########################################
###############################################################
## What to plot from "GENE" table
 # Distribution amongst genes
   # Size of Gene
     # Size of Exons
     # x vs y (?)
   # Number of Exons
   # Number of Variants
     # Num SNPs & Indels
   # Area
   # % Samples that are Compound Hets

## Heatmap Amongst Gene Variables
print("Plotting Heatmap")
WHICH_COLS <- grep("n.EX",colnames(GENE)):ncol(GENE)
COR.dat <- cor( matrix(as.numeric(c(GENE[,WHICH_COLS])),ncol=length(WHICH_COLS)), use="pairwise.complete.obs",method="spearman" )
# COR.dat <- cor( GENE[,WHICH_COLS], use="pairwise.complete.obs",method="spearman" )
WHICH_RM <- which( apply( COR.dat,1,function(x)all(is.na(x)) ))
colnames(COR.dat) <- rownames(COR.dat) <- colnames(GENE)[WHICH_COLS]
if (length(WHICH_RM)>0) { COR.dat <- COR.dat[ -WHICH_RM, -WHICH_RM ] }
COLS.list <- c("gold1","chocolate2","firebrick2","black","slateblue3","steelblue2","springgreen1")
COLS <- colorRampPalette(COLS.list)(100)
jpeg( paste(PathToOut,"/Plots/Compile_Gene_Heat.jpeg",sep=""), height=2000,width=2000, pointsize=36)
heatmap.2( COR.dat, col=COLS, trace="none" )
dev.off()

## Plot Size of Gene vs Size of Exons
print("Plotting Size Scatterplot")
COLS <- c("steelblue1","slateblue1","black")
jpeg( paste(PathToOut,"/Plots/Compile_Scatter_GeneSize.jpeg",sep=""), height=1200,width=3000, pointsize=32)
par(mfrow=c(1,3))
SIZE.tx <- as.numeric(GENE[,"TX_E"]) - as.numeric(GENE[,"TX_S"])
BRKS.tx <- seq( min(SIZE.tx), max(SIZE.tx)+20000, 20000 )
hist( SIZE.tx, breaks=BRKS.tx, col=COLS[1], main="Distribution of Transcribed Region Size", xlab="Size of Transcribed Region (bp)", ylab="# Genes" )
SIZE.ex <- sapply(lapply(strsplit(GENE[,"EX_E"],","),function(x) sum(as.numeric(x)) ),"[",1) - sapply(lapply(strsplit(GENE[,"EX_S"],","),function(x) sum(as.numeric(x)) ),"[",1)
BRKS.ex <- seq( min(SIZE.ex), max(SIZE.ex)+200, 200 )
hist( SIZE.ex, breaks=BRKS.ex, col=COLS[2], main="Distribution of Exonic Size", xlab="Size of Exonic Region (bp)", ylab="# Genes" )
plot( SIZE.tx, SIZE.ex, pch="+", col=COLS[3], main="Transcribed vs Exonic Size", xlab="Size of Transcribed Region (bp)",ylab="Size of Exonic Region (bp)" )
abline( lm( SIZE.ex ~ SIZE.tx ), lwd=2)
dev.off()

## Distribution of Number of Variants
print("Plotting Num Var Distribution")
COLS <- c("steelblue1","slateblue1","springgreen3","gold3")
BIN <- 2
jpeg( paste(PathToOut,"/Plots/Compile_Num_Vars.jpeg",sep=""), height=1200,width=2400, pointsize=32)
par(mfrow=c(1,2))
 # Full
VAR.dat <- matrix( as.numeric(GENE[,c("n.VAR","n.SNP","n.IND")]), ncol=3 )
BRKS <- seq( 0,as.numeric(max(VAR.dat,na.rm=T))+250, 250)
hist( VAR.dat[,3], breaks=BRKS, col=COLS[4], density=20,angle=-30, main="Distribution of Variant Count (Full)", xlab="Variant Count", ylab="# Genes")
hist( VAR.dat[,1], breaks=BRKS, col=COLS[1], density=20,angle=30, main="Distribution of Variant Count (Full)", xlab="Variant Count", ylab="# Genes", add=T)
hist( VAR.dat[,2], breaks=BRKS, col=COLS[3], density=20,angle=60, main="Distribution of Variant Count (Full)", xlab="Variant Count", ylab="# Genes", add=T)
 # Exon
VAR.dat <- matrix( as.numeric(GENE[,c("n.VAR.ex","n.SNP.ex","n.IND.ex")]), ncol=3 )
BRKS <- seq( 0,as.numeric(max(VAR.dat,na.rm=T))+5, 5)
hist( VAR.dat[,3], breaks=BRKS, col=COLS[4], density=20,angle=-30, main="Distribution of Variant Count (Exon)", xlab="Variant Count", ylab="# Genes")
hist( VAR.dat[,1], breaks=BRKS, col=COLS[2], density=20,angle=30, main="Distribution of Variant Count (Exon)", xlab="Variant Count", ylab="# Genes", add=T)
hist( VAR.dat[,2], breaks=BRKS, col=COLS[3], density=20,angle=60, main="Distribution of Variant Count (Exon)", xlab="Variant Count", ylab="# Genes", add=T)
dev.off()

## Distribution of AREA
print("Plotting Area Distribution")
COLS <- c("steelblue1","slateblue1")
BIN <- .1
jpeg( paste(PathToOut,"/Plots/Compile_Area_Dist.jpeg",sep=""), height=1200,width=2400, pointsize=32)
par(mfrow=c(1,2))
 # Full
AREA.dat <- as.numeric( GENE[,"AREA"] )
BRKS <- seq( min(AREA.dat,na.rm=T)-BIN, max(AREA.dat,na.rm=T)+BIN, BIN )
hist( AREA.dat, breaks=BRKS, col=COLS[1], main="Distribution of AREA (Full)", xlab="AREA", ylab="# Genes")
 # Exon
AREA.dat <- as.numeric( GENE[,"AREA.ex"] )
BRKS <- seq( min(AREA.dat,na.rm=T)-BIN, max(AREA.dat,na.rm=T)+BIN, BIN )
hist( AREA.dat, breaks=BRKS, col=COLS[2], main="Distribution of AREA (Exon)", xlab="AREA", ylab="# Genes" )
dev.off()

## Distribution of Best P-Values
print("Plotting Best P-Value Distribution")
COLS <- c("steelblue1","slateblue1")
BIN <- .2
jpeg( paste(PathToOut,"/Plots/Compile_BestP_Dist.jpeg",sep=""), height=1200,width=2400, pointsize=32)
par(mfrow=c(1,2))
 # Full
BESTP_P.dat <- -log10( as.numeric( GENE[,"BEST_P"] ) )
BRKS <- seq( min(BESTP_P.dat,na.rm=T)-BIN, max(BESTP_P.dat,na.rm=T)+BIN, BIN )
hist( BESTP_P.dat, breaks=BRKS, col=COLS[1], main="Distribution of Best P-Value (Full)", xlab="Best P-Value -log10(p)", ylab="# Genes")
 # Exon
BESTP_P.dat <- -log10( as.numeric( GENE[,"BEST_P.ex"] ) )
BRKS <- seq( min(BESTP_P.dat,na.rm=T)-BIN, max(BESTP_P.dat,na.rm=T)+BIN, BIN )
hist( BESTP_P.dat, breaks=BRKS, col=COLS[2], main="Distribution of Best P-Value (Exon)", xlab="Best P-Value -log10(p)", ylab="# Genes" )
dev.off()

## Distribution of % Samples that are Compound Het
print("Plotting Perc Compound Hets")
COLS <- c("springgreen3","gold3")
BIN <- .05
XLIM <- c( 0,1 )
BRKS <- seq(0,1,BIN)
jpeg( paste(PathToOut,"/Plots/Compile_CH_Dist.jpeg",sep=""), height=1200,width=2400, pointsize=32)
par(mfrow=c(1,2))
 # Full
CH.dat.snp <- as.numeric( GENE[,"p.COMP_HET_snp"] )
CH.dat.ind <- as.numeric( GENE[,"p.COMP_HET_ind"] )
hist( CH.dat.ind, breaks=BRKS, col=COLS[2], main="Distribution of % Compound Hets (Full)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=-45 )
hist( CH.dat.snp, breaks=BRKS, col=COLS[1], main="Distribution of % Compound Hets (Full)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=45, add=T )
 # Exon
CH.dat.snp <- as.numeric( GENE[,"p.COMP_HET_snp.ex"] )
CH.dat.ind <- as.numeric( GENE[,"p.COMP_HET_ind.ex"] )
hist( CH.dat.ind, breaks=BRKS, col=COLS[2], main="Distribution of % Compound Hets (Exon)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=-45 )
hist( CH.dat.snp, breaks=BRKS, col=COLS[1], main="Distribution of % Compound Hets (Exon)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=45, add=T )
legend( "topright", fill=COLS, density=20, legend=c("SNP","IND") )
dev.off()

###############################################################
## END OF DOC #################################################
###############################################################


# GTAB.1 <- read.table("Gene_Stats.1_950.txt",sep="\t",header=T)
# GTAB.2 <- read.table("Gene_Stats.951_1750.txt",sep="\t",header=T)
# GTAB.3 <- read.table("Gene_Stats.txt",sep="\t",header=T)

# GTAB <- rbind( GTAB.1[1:950,], GTAB.2[951:1750,], GTAB.3[1751:nrow(GTAB.3),] )
# GARR <- array( ,dim(GTAB))
# colnames(GARR) <- colnames(GTAB)
# for ( col in 1:ncol(GTAB) ) { GARR[,col] <- as.character(GTAB[,col]) }
# GENE <- GARR


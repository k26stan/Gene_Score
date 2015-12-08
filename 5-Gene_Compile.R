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
# LINE <- c( "/projects/janssen/Phased/20151016_DEL_CAND_Genes","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
# LINE <- c( "/projects/janssen/Phased/20151106_GeneEnrichment_CANDS","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
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
	if ( "pheno"%in%colnames(PHENO) ) { pheno_colname <- "pheno" }
	if ( "Pheno"%in%colnames(PHENO) ) { pheno_colname <- "Pheno" }
	colnames(PHENO)[which(colnames(PHENO)==pheno_colname)] <- "Pheno"
	PC <- merge( x=PHENO[,c("IID","Pheno")], y=COVS, by="IID" )
	PC.2 <- PC
	PC.2[,"IID"] <- unlist(strsplit( as.character(PC$IID),"-"))[seq(1,2*nrow(PC),2)]
	## Load Association Results
	SNP.P <- read.table( paste(PathToAssoc,"SNP/SNP_Assoc.P",sep=""), sep="\t",header=T)
	SNP.HWE <- read.table( paste(PathToAssoc,"SNP/SNP_Assoc.hwe",sep=""), sep="",header=T)
	IND.P <- read.table( paste(PathToAssoc,"IND/IND_Assoc.P",sep=""), sep="\t",header=T)
	IND.HWE <- read.table( paste(PathToAssoc,"IND/IND_Assoc.hwe",sep=""), sep="",header=T)
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
n.samps <- length(hap.samps)
hap.colnames <- c("CHR","SNP","BP","REF","ALT", paste( rep(hap.samps, rep(2,n.samps)), 1:2, sep="_" ) )

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

PLOT_FRACTION <- 1/50
PLOT_FRACTION <- 1
start_time <- proc.time()
## LOOP START ######################################
for ( gtx in 1:n.gtx ) {
for ( gtx in grep("CYP",GTX_LIST) ) {
for ( gtx in grep("SYCE",GTX_LIST) ) {
# for ( gtx in 401:n.gtx ) {
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
	## LOAD DATA/FILES #################################
	
	## Check if files exist & contain >0 Lines
	 # If not, skip to next Gene_Transcript
	if ( !file.exists(paste(PathToGenes,name,"/SNP_Vars.hwe",sep="")) | !file.exists(paste(PathToGenes,name,"/IND_Vars.hwe",sep="")) ) { next }
	if ( length(readLines(paste(PathToGenes,name,"/Phased.haps",sep="")))==0 ) { next }
	## Load HW Files
	print( "Loading HW Files" )
	snp.hwe <- read.table( paste(PathToGenes,name,"/SNP_Vars.hwe",sep=""), sep="",header=T )
	ind.hwe <- read.table( paste(PathToGenes,name,"/IND_Vars.hwe",sep=""), sep="",header=T )
	gtx.hwe <- rbind( snp.hwe, ind.hwe )
	## Load Genotype Files
	print( "Loading Raw GT Files" )
	snp.raw <- read.table( paste(PathToGenes,name,"/SNP_Vars.raw",sep=""), sep="",header=T )
	snp.raw[,"IID"] <- sapply( strsplit(as.character(snp.raw[,"IID"]),"-"), "[", 1 )
	snp.raw.tag <- colnames(snp.raw)[7:ncol(snp.raw)]
	ind.raw <- read.table( paste(PathToGenes,name,"/IND_Vars.raw",sep=""), sep="",header=T )
	ind.raw[,"IID"] <- sapply( strsplit(as.character(ind.raw[,"IID"]),"-"), "[", 1 )
	ind.raw.tag <- colnames(ind.raw)[7:ncol(ind.raw)]
	shared.samps <- sort( Reduce( intersect, list(hap.samps,snp.raw[,"IID"],ind.raw[,"IID"]) ) )
	## Load Phased Haplotype File (& Rename Columns)
	print( "Loading Hap Files" )
	hap <- read.table( paste(PathToGenes,name,"/Phased.haps",sep=""), sep="",header=F )
	colnames(hap) <- hap.colnames
	which.hap.ind <- union( which(nchar(as.character(hap$REF))>1), which(nchar(as.character(hap$ALT))>1) )
	hap.ind.TF <- c("snp","ind")[factor(1:nrow(hap) %in% which.hap.ind)]
	hap <- data.frame( TAG=paste(hap[,"CHR"],hap[,"BP"],sep="_"), TYPE=hap.ind.TF, hap )
	## Load Annotation File
	# annot <- read.table( paste(PathToGenes,name,"/Annots_Short.txt",sep=""), sep="\t",header=T )
	# annot <- data.frame( TAG=paste(gsub("chr","",as.character(annot$Chromosome)),as.character(annot$End),sep="_"), annot)
	# eqtls <- read.table( paste(PathToGenes,name,"/eQTLs.txt",sep=""), sep="",header=T )

	####################################################
	## ORGANIZE DATA/FILES #############################
	## Pull out Variant Positions from BIM files (& Compile to 1 file)
	print( "Sorting Bim Files" )
	gtx.snp.bim <- SNP.bim[ which(SNP.bim$CHR==chr & SNP.bim$BP>=rng[1] & SNP.bim$BP<=rng[2] ), ]
	gtx.snp.bim <- data.frame( gtx.snp.bim, TYPE=rep("snp",nrow(gtx.snp.bim)), RAW_TAG=snp.raw.tag )
	gtx.ind.bim <- IND.bim[ which(IND.bim$CHR==chr & IND.bim$BP>=rng[1] & IND.bim$BP<=rng[2] ), ]
	gtx.ind.bim <- data.frame( gtx.ind.bim, TYPE=rep("ind",nrow(gtx.ind.bim)), RAW_TAG=ind.raw.tag )
	gtx.bim <- rbind( gtx.snp.bim, gtx.ind.bim )
	gtx.bim <- gtx.bim[ order(gtx.bim[,"BP"]), ]
	## Merge BIM table w/ HWE table
	print( "Merging BIM/HW Files" )
	gtx.mg <- merge( gtx.bim[,c("SNP","BP","TYPE","RAW_TAG")], gtx.hwe, by="SNP", all=T )
	gtx.mg <- data.frame( TAG=paste(gtx.mg[,"CHR"],gtx.mg[,"BP"],sep="_"), gtx.mg )
	gtx.mg <- gtx.mg[ ,c("TAG","RAW_TAG","CHR","BP","SNP","TYPE","A1","A2","GENO","P") ]
	colnames(gtx.mg)[ncol(gtx.mg)] <- "P_HW"
	GTX <- gtx.mg
	## Remove HW Violations
	RM.hwe <- which( GTX$P_HW < 1e-20 )
	if ( length(RM.hwe)>0 ) { 
		RM.hwe.tags <- as.character( GTX$TAG[ RM.hwe ] )
		GTX <- GTX[ -RM.hwe, ]
	}
	## Compile which variants are in Exons
	print( "Pulling out Exonic Variants" )
	TAG.exon <- c()
	for ( e in 1:ex_cnt ) {
		WHICH <- which( gtx.mg[,"BP"]>ex_b[e] & gtx.mg[,"BP"]<ex_e[e] )
		if ( length(WHICH)>0) { TAG.exon <- c( TAG.exon, as.character(gtx.mg[WHICH,"TAG"]) ) }
	}
	## Pull Out Single-Locus Results
	if ( file.exists(PathToPheno) ) {
		print( "Pulling out Single-Locus Results" )
		snp.p <- SNP.P[ which(SNP.P$CHR==chr & SNP.P$BP>=tx_rng[1]-5000 & SNP.P$BP<=tx_rng[2]+5000 ), ]
		ind.p <- IND.P[ which(IND.P$CHR==chr & IND.P$BP>=tx_rng[1]-5000 & IND.P$BP<=tx_rng[2]+5000 ), ]
		gtx.p <- rbind( snp.p, ind.p )
		gtx.p <- gtx.p[order(gtx.p[,"BP"]),]
		gtx.p <- data.frame( TAG=paste(gtx.p[,"CHR"],gtx.p[,"BP"],sep="_"), gtx.p )
		# Merge with HWE Data (compile all PLink data)
		gtx.mg.p <- merge( gtx.mg, gtx.p[,c("SNP","P")], by="SNP", all=T )
		colnames(gtx.mg.p)[which(colnames(gtx.mg.p)=="P")] <- "P_Assoc"
		gtx.mg.p <- gtx.mg.p[ order(gtx.mg.p[,"BP"]), ]
		gtx.mg.p <- gtx.mg.p[ ,c("TAG","RAW_TAG","CHR","BP","SNP","TYPE","A1","A2","P_Assoc","GENO","P_HW") ]
		GTX <- gtx.mg.p
	}

	## Specify Location of Variant (relative to Exon/Intron/Etc)
	print( "Specifying Variant Locations" )
	GTX <- data.frame( GTX, LOC=rep(1,nrow(GTX)) )
	GTX[ which( GTX[,"BP"]>=tx_rng[1] & GTX[,"BP"]<=tx_rng[2] ), "LOC" ] <- 2
	GTX[ which( GTX[,"BP"]>=cd_rng[1] & GTX[,"BP"]<=cd_rng[2] ), "LOC" ] <- 3
	GTX[ which( GTX[,"TAG"] %in% TAG.exon ), "LOC" ] <- 4
	## Calculate Allele Frequencies
	print( "Calculating MAF" )
	gtx.geno.arr <- t(sapply( strsplit( as.character(GTX[,"GENO"]), "/" ), "[", 1:3 ))
	gtx.af <- matrix( as.numeric(gtx.geno.arr), ncol=3 ) %*% matrix(2:0,c(3,1)) / (2*n.samps)
	GTX <- data.frame( GTX, MAF=gtx.af )
	## Merge with Haplotype File
	GTX.ph <- merge( GTX, hap[,-which(colnames(hap)%in%c("TYPE","CHR","SNP","BP"))], by="TAG", all=T )

	####################################################
	## PLOT GENE MAP & GWAS RESULTS ####################
	PLOT_RUNIF <- runif(1,0,1)

	## GENE MAPS ##
	print( "Plotting Gene Map" )
	## Set Plotting Parameters
	NUM_VARS <- nrow(GTX)
	COLS.gene <- c("steelblue1","cadetblue1","slateblue2")
	COLS.P.list <- c("firebrick4","steelblue3","cadetblue3","slateblue3")
	XLIM <- c(rng[1],rng[2])
	if ( file.exists(PathToPheno) ) {
 		YLIM <- c( -1, max( 4, max(-log10(GTX[,"P_Assoc"]),na.rm=T) ) )
 	}else{ YLIM <- c(-1,1) }
	Y <- 0
	if ( YLIM[2]>4 | PLOT_RUNIF<PLOT_FRACTION ) {
		jpeg( paste(PathToOut,"/Plots/Gene_Map_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
		plot(0,0,type="n", xlim=XLIM,ylim=YLIM, xlab=paste("Chromosome",chr,"Position"),ylab="-log10(p)", main=paste("Map & Single-Locus Results of:",name), yaxt="n" )
		## Plot P-Values
		if ( file.exists(PathToPheno) ) {
			WHICH <- which(!is.na( GTX$P_Assoc))
			COLS.P <- COLS.P.list[ GTX$LOC ]
			abline( h=seq(0,YLIM[2],1), lty=2, col="grey50" )
			abline( h=seq(-1,0,.1), lty=2, col="grey50" )
			axis(2, at=seq(-1,YLIM[2],1) )
			points( GTX$BP[WHICH], -log10(GTX$P_Assoc[WHICH]), col=COLS.P[WHICH], pch=c(10,16)[factor(GTX$TYPE[WHICH])], lwd=3, cex=2*(GTX$MAF[WHICH])^(1/6) )
		}
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
		arrows( GTX$BP, 0, GTX$BP, -GTX$MAF, col=COLS.P, angle=90 )
		dev.off()
	} # Close Gene Map "IF"

	## QQ PLOTS ##
	print( "Plotting QQ Plot" )
	if ( file.exists(PathToPheno) ) {
		NUM_VARS <- length(which( !is.na(GTX$P_Assoc) ))
		NUM_EXONIC <- length(which( !is.na(GTX$P_Assoc) & GTX$LOC==4 ))
		COLS <- c("steelblue1","slateblue1")
		COLS.4 <- gsub("1","3",COLS)
		LIM <- c( 0, max( 4, max(-log10(GTX[,"P_Assoc"]),na.rm=T) ) )
		## Calculate Observed, Expected, & Area for Plot
		 # Expected & Observed
		EXP.full <- -log10( 1:NUM_VARS / NUM_VARS )
		OBS.full <- -log10( sort( GTX$P_Assoc ) )
		IND.full <- c(1:length(EXP.full),length(EXP.full),1)
		 # Calculate Area - Trapezoid Rule
		H_VALS <- EXP.full[2:NUM_VARS-1] - EXP.full[2:NUM_VARS]
		B_SUM <- ( OBS.full[2:NUM_VARS-1]-EXP.full[2:NUM_VARS-1] ) + ( OBS.full[2:NUM_VARS]-EXP.full[2:NUM_VARS] )
		AREA.full.traps <- .5*H_VALS*B_SUM
		AREA.full <- sum( AREA.full.traps ) # - .5*max(EXP.full)^2
		GENE[gtx,"AREA"] <- AREA.full
		GENE[gtx,"BEST_P"] <- 10^(min(-OBS.full,na.rm=T))
		if ( NUM_EXONIC > 0 ) {
			 # Expected & Observed
			EXP.exon <- -log10( 1:NUM_EXONIC / NUM_EXONIC )
			OBS.exon <- -log10( sort( GTX$P_Assoc[which(GTX$LOC==4)] ) )
			IND.exon <- c(1:length(EXP.exon),length(EXP.exon),1)
			 # Calculate Area - Trapezoid Rule
			H_VALS <- EXP.exon[2:NUM_EXONIC-1] - EXP.exon[2:NUM_EXONIC]
			B_SUM <- ( OBS.exon[2:NUM_EXONIC-1]-EXP.exon[2:NUM_EXONIC-1] ) + ( OBS.exon[2:NUM_EXONIC]-EXP.exon[2:NUM_EXONIC] )
			AREA.exon.traps <- .5*H_VALS*B_SUM
			AREA.exon <- sum( AREA.exon.traps ) # - .5*max(EXP.exon)^2
			GENE[gtx,"AREA.ex"] <- AREA.exon
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
			polygon( EXP.full[IND.full], c(OBS.full,EXP.full[c(NUM_VARS,1)]), col=COLS[1], border=COLS.4[1], density=20,angle=45 )
			points( EXP.full, OBS.full, pch="+", col=COLS.4[1], type="o" )
			text( quantile(LIM,.7),quantile(LIM,.1), label=paste("Area:",round(AREA.full,3),"-",NUM_VARS,"Vars"), col=COLS.4[1], cex=1.2 )
			## Plot Exon Set
			if ( NUM_EXONIC > 0 ) {
				polygon( EXP.exon[IND.exon], c(OBS.exon,EXP.exon[c(NUM_EXONIC,1)]), col=COLS[2], border=COLS.4[2], density=20,angle=-45 )
				points( EXP.exon, OBS.exon, pch="+", col=COLS.4[2], type="o" )	
				text( quantile(LIM,.7),quantile(LIM,.07), label=paste("Area:",round(AREA.exon,3),"-",NUM_EXONIC,"Vars"), col=COLS.4[2], cex=1.2 )
			}
			legend("topleft",fill=COLS,density=20,legend=c("Full","Exonic") )
			dev.off()
		} # Close QQ "IF"
	}

	####################################################
	## COMPILE STATS ABOUT VARIANTS/COHORT #############
	 # By Location
	   # Ovarall
 	   # Exonic
 	   # Intronic
	   # Damaging (?)
	 # What to compile
	   # % Phased
 	   # How many hets, hom_vars pp (SNPs & Indels separate)
 	   # How many comp hets pp

 	## Compile Info about the Gene_Transcript
 	print( "Compiling Gene Stats" )
 	## Gene
 	 # Full
 	GENE[gtx,"n.VAR"] <- nrow(GTX)
 	GENE[gtx,"n.SNP"] <- length(which( GTX$TYPE=="snp" ))
 	GENE[gtx,"n.IND"] <- length(which( GTX$TYPE=="ind" ))
 	GENE[gtx,"n.VAR.ph"] <- nrow(hap)
 	GENE[gtx,"n.SNP.ph"] <- length(which( hap[,"TYPE"]=="snp" ))
 	GENE[gtx,"n.IND.ph"] <- length(which( hap[,"TYPE"]=="ind" ))
 	 # Exon
 	GTX.exon <- GTX[ which(GTX$LOC==4), ]
 	GTX.ph.exon <- GTX.ph[ which(GTX.ph$LOC==4), ]
 	hap.exon <- hap[ which(hap$TAG %in% TAG.exon), ]
 	GENE[gtx,"n.VAR.ex"] <- nrow(GTX.exon)
 	GENE[gtx,"n.SNP.ex"] <- length(which( GTX.exon$TYPE=="snp" ))
 	GENE[gtx,"n.IND.ex"] <- length(which( GTX.exon$TYPE=="ind" ))
 	GENE[gtx,"n.VAR.ph.ex"] <- nrow(hap.exon)
 	GENE[gtx,"n.SNP.ph.ex"] <- length(which( hap.exon[,"TYPE"]=="snp" ))
 	GENE[gtx,"n.IND.ph.ex"] <- length(which( hap.exon[,"TYPE"]=="ind" ))

	## Cohort - Full ###################################
 	print( "Compiling Cohort Stats (Full)" )
 	colnames.FL <- c("HET_snp","HOM_VAR_snp","PRC_PHAS_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","PRC_PHAS_ind","COMP_HET_ind","COMP_HET_any","HAP1_snp","HAP2_snp","UNPH_snp","HAP1_ind","HAP2_ind","UNPH_ind")
 	FULL[[name]] <- array( ,c(n.samps,length(colnames.FL)) )
 	rownames(FULL[[name]]) <- shared.samps
 	colnames(FULL[[name]]) <- colnames.FL
   	n.snp <- ncol(snp.raw)-6
 	n.ind <- ncol(ind.raw)-6
	 # How many Hets?
 	if ( n.snp>0 ) {
 		if ( n.snp==1 ) {
 			FULL[[name]][,"HET_snp"] <- as.numeric( snp.raw[order(snp.raw[,"IID"]),7]==1 )
 		}else{ FULL[[name]][,"HET_snp"] <- apply( snp.raw[order(snp.raw[,"IID"]),7:ncol(snp.raw)], 1, function(x) length(which(x==1)) ) }
 	}else{ FULL[[name]][,"HET_snp"] <- 0 }
 	if ( n.ind>0 ) {
 		if ( n.ind==1 ) {
 			FULL[[name]][,"HET_ind"] <- as.numeric( ind.raw[order(ind.raw[,"IID"]),7]==1 )
 		}else{ FULL[[name]][,"HET_ind"] <- apply( ind.raw[order(ind.raw[,"IID"]),7:ncol(ind.raw)], 1, function(x) length(which(x==1)) ) }
 	}else{ FULL[[name]][,"HET_ind"] <- 0 }
 	 # How many Hom_Vars?
 	if ( n.snp>0 ) {
 		if ( n.snp==1 ) {
 			FULL[[name]][,"HOM_VAR_snp"] <- as.numeric( snp.raw[order(snp.raw[,"IID"]),7]==2 )
 		}else{ FULL[[name]][,"HOM_VAR_snp"] <- apply( snp.raw[order(snp.raw[,"IID"]),7:ncol(snp.raw)], 1, function(x) length(which(x==2)) ) }
 	}else{ FULL[[name]][,"HOM_VAR_snp"] <- 0 }
 	if ( n.ind>0 ) {
 		if ( n.ind==1 ) {
 			FULL[[name]][,"HOM_VAR_ind"] <- as.numeric( ind.raw[order(ind.raw[,"IID"]),7]==2 )
 		}else{ FULL[[name]][,"HOM_VAR_ind"] <- apply( ind.raw[order(ind.raw[,"IID"]),7:ncol(ind.raw)], 1, function(x) length(which(x==2)) ) }
 	}else{ FULL[[name]][,"HOM_VAR_ind"] <- 0 }
 	 # What percent of variants were phased?
 	hap.1.cols <- seq( which(colnames(hap)=="ALT")+1,ncol(hap),2)
 	hap.2.cols <- seq( which(colnames(hap)=="ALT")+2,ncol(hap),2)
 	 # SNP
 	hap.snp <- hap[ which(hap$TYPE=="snp") , ]
 	hap.snp.diffs <- hap.snp[,hap.1.cols] - hap.snp[,hap.2.cols]
 	colnames(hap.snp.diffs) <- gsub( "_1","", colnames(hap.snp.diffs) )
 	hap.snp.n.hets <- apply( hap.snp.diffs, 2, function(x) length(which(x!=0)) )
 	FULL[[name]][,"PRC_PHAS_snp"] <- hap.snp.n.hets / FULL[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.n.hets / FULL[[name]][,"HET_snp"]
 	 # Indel
 	hap.ind <- hap[ which(hap$TYPE=="ind") , ]
 	hap.ind.diffs <- hap.ind[,hap.1.cols] - hap.ind[,hap.2.cols]
 	colnames(hap.ind.diffs) <- gsub( "_1","", colnames(hap.ind.diffs) )
 	hap.ind.n.hets <- apply( hap.ind.diffs, 2, function(x) length(which(x!=0)) )
 	FULL[[name]][,"PRC_PHAS_ind"] <- hap.ind.n.hets / FULL[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.n.hets / FULL[[name]][,"HET_ind"]
 	 # Any
 	hap.diffs <- hap[,hap.1.cols] - hap[,hap.2.cols]
 	colnames(hap.diffs) <- gsub( "_1","", colnames(hap.diffs) )
 	hap.n.hets <- apply( hap.diffs, 2, function(x) length(which(x!=0)) )
 	 # How many Variants are on each Strand and Unphased?
 	hap.snp.counts <- colSums(hap.snp[,8:ncol(hap.snp)])
 	FULL[[name]][,"HAP1_snp"] <- hap.snp.counts[seq(1,2*n.samps,2)]
 	FULL[[name]][,"HAP2_snp"] <- hap.snp.counts[seq(2,2*n.samps,2)]
 	FULL[[name]][,"UNPH_snp"] <- FULL[[name]][,"HET_snp"] - hap.snp.n.hets
 	hap.ind.counts <- colSums(hap.ind[,8:ncol(hap.ind)])
 	FULL[[name]][,"HAP1_ind"] <- hap.ind.counts[seq(1,2*n.samps,2)]
 	FULL[[name]][,"HAP2_ind"] <- hap.ind.counts[seq(2,2*n.samps,2)]
 	FULL[[name]][,"UNPH_ind"] <- FULL[[name]][,"HET_ind"] - hap.ind.n.hets
 	# SAMP <- "B012326"
	# SAMP.hets.hap <- as.character( hap.snp[ which( hap.snp.diffs[,SAMP]!=0 ), "TAG" ])
	# SAMP.hets.raw <- colnames(snp.raw)[ which( snp.raw[ which(snp.raw[,"IID"]==SAMP), ]==1 ) ]
	# SAMP.hets.bim <- as.character( GTX[ which(GTX[,"RAW_TAG"] %in% SAMP.hets.raw), "TAG" ] )
 	 # Which Samples have a Compound Het?
 	hap.snp.comp.hets <- as.numeric( apply( hap.snp.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	hap.ind.comp.hets <- as.numeric( apply( hap.ind.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
  	hap.comp.hets <- as.numeric( apply( hap.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
	FULL[[name]][,"COMP_HET_snp"] <- hap.snp.comp.hets
 	FULL[[name]][,"COMP_HET_ind"] <- hap.ind.comp.hets
	FULL[[name]][,"COMP_HET_any"] <- hap.comp.hets

 	## Cohort - Exon ###################################
 	print( "Compiling Cohort Stats (Exon)" )
 	colnames.EX <- c("HET_snp","HOM_VAR_snp","PRC_PHAS_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","PRC_PHAS_ind","COMP_HET_ind","COMP_HET_any","HAP1_snp","HAP2_snp","UNPH_snp","HAP1_ind","HAP2_ind","UNPH_ind")
 	EXON[[name]] <- array( ,c(n.samps,length(colnames.EX)) )
 	rownames(EXON[[name]]) <- shared.samps
 	colnames(EXON[[name]]) <- colnames.EX
 	exon.tag <- as.character( GTX[ which(GTX$LOC==4),"TAG"] )
 	exon.tag.raw <- as.character( GTX[ which(GTX$LOC==4),"RAW_TAG"] )
 	snp.raw.exon <- snp.raw[ ,c( 1:6,which(colnames(snp.raw)%in%exon.tag.raw) ) ]
 	ind.raw.exon <- ind.raw[ ,c( 1:6,which(colnames(ind.raw)%in%exon.tag.raw) ) ]
 	n.snp.exon <- ncol(snp.raw.exon)-6
 	n.ind.exon <- ncol(ind.raw.exon)-6
 	 # How many Hets?
 	if ( n.snp.exon>0 ) {
 		if ( n.snp.exon==1 ) {
 			EXON[[name]][,"HET_snp"] <- as.numeric( snp.raw.exon[order(snp.raw.exon[,"IID"]),7]==1 )
 		}else{ EXON[[name]][,"HET_snp"] <- apply( snp.raw.exon[order(snp.raw.exon[,"IID"]),7:ncol(snp.raw.exon)], 1, function(x) length(which(x==1)) ) }
 	}else{ EXON[[name]][,"HET_snp"] <- 0 }
 	if ( n.ind.exon>0 ) {
 		if ( n.ind.exon==1 ) {
 			EXON[[name]][,"HET_ind"] <- as.numeric( ind.raw.exon[order(ind.raw.exon[,"IID"]),7]==1 )
 		}else{ EXON[[name]][,"HET_ind"] <- apply( ind.raw.exon[order(ind.raw.exon[,"IID"]),7:ncol(ind.raw.exon)], 1, function(x) length(which(x==1)) ) }
 	}else{ EXON[[name]][,"HET_ind"] <- 0 }
 	 # How many Hom_Vars?
 	if ( n.snp.exon>0 ) {
 		if ( n.snp.exon==1 ) {
 			EXON[[name]][,"HOM_VAR_snp"] <- as.numeric( snp.raw.exon[order(snp.raw.exon[,"IID"]),7]==2 )
 		}else{ EXON[[name]][,"HOM_VAR_snp"] <- apply( snp.raw.exon[order(snp.raw.exon[,"IID"]),7:ncol(snp.raw.exon)], 1, function(x) length(which(x==2)) ) }
 	}else{ EXON[[name]][,"HOM_VAR_snp"] <- 0 }
 	if ( n.ind.exon>0 ) {
 		if ( n.ind.exon==1 ) {
 			EXON[[name]][,"HOM_VAR_ind"] <- as.numeric( ind.raw.exon[order(ind.raw.exon[,"IID"]),7]==2 )
 		}else{ EXON[[name]][,"HOM_VAR_ind"] <- apply( ind.raw.exon[order(ind.raw.exon[,"IID"]),7:ncol(ind.raw.exon)], 1, function(x) length(which(x==2)) ) }
 	}else{ EXON[[name]][,"HOM_VAR_ind"] <- 0 }
  	 # What percent of variants were phased?
 	hap.1.cols <- seq( which(colnames(hap)=="ALT")+1,ncol(hap),2)
 	hap.2.cols <- seq( which(colnames(hap)=="ALT")+2,ncol(hap),2)
 	 # SNP
 	hap.snp <- hap[ which(hap$TYPE=="snp") , ]
 	hap.snp.exon <- hap.snp[ which(hap.snp[,"TAG"]%in%exon.tag), ]
 	hap.snp.exon.diffs <- hap.snp.exon[,hap.1.cols] - hap.snp.exon[,hap.2.cols]
 	colnames(hap.snp.exon.diffs) <- gsub( "_1","", colnames(hap.snp.exon.diffs) )
 	hap.snp.exon.n.hets <- apply( hap.snp.exon.diffs, 2, function(x) length(which(x!=0)) )
 	EXON[[name]][,"PRC_PHAS_snp"] <- hap.snp.exon.n.hets / EXON[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.exon.n.hets / EXON[[name]][,"HET_snp"]
 	 # Indel
  	hap.ind <- hap[ which(hap$TYPE=="ind") , ]
 	hap.ind.exon <- hap.ind[ which(hap.ind[,"TAG"]%in%exon.tag), ]
 	hap.ind.exon.diffs <- hap.ind.exon[,hap.1.cols] - hap.ind.exon[,hap.2.cols]
 	colnames(hap.ind.exon.diffs) <- gsub( "_1","", colnames(hap.ind.exon.diffs) )
 	hap.ind.exon.n.hets <- apply( hap.ind.exon.diffs, 2, function(x) length(which(x!=0)) )
 	EXON[[name]][,"PRC_PHAS_ind"] <- hap.ind.exon.n.hets / EXON[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.exon.n.hets / EXON[[name]][,"HET_ind"]
 	 # Any
 	hap.exon.diffs <- hap.exon[,hap.1.cols] - hap.exon[,hap.2.cols]
 	colnames(hap.exon.diffs) <- gsub( "_1","", colnames(hap.exon.diffs) )
 	hap.exon.n.hets <- apply( hap.exon.diffs, 2, function(x) length(which(x!=0)) )
 	 # How many Variants are on each Strand and Unphased?
 	hap.snp.exon.counts <- colSums(hap.snp.exon[,8:ncol(hap.snp.exon)])
 	EXON[[name]][,"HAP1_snp"] <- hap.snp.exon.counts[seq(1,2*n.samps,2)]
 	EXON[[name]][,"HAP2_snp"] <- hap.snp.exon.counts[seq(2,2*n.samps,2)]
 	EXON[[name]][,"UNPH_snp"] <- EXON[[name]][,"HET_snp"] - hap.snp.exon.n.hets
 	hap.ind.exon.counts <- colSums(hap.ind.exon[,8:ncol(hap.ind.exon)])
 	EXON[[name]][,"HAP1_ind"] <- hap.ind.exon.counts[seq(1,2*n.samps,2)]
 	EXON[[name]][,"HAP2_ind"] <- hap.ind.exon.counts[seq(2,2*n.samps,2)]
 	EXON[[name]][,"UNPH_ind"] <- EXON[[name]][,"HET_ind"] - hap.ind.exon.n.hets
 	 # Which Samples have a Compound Het?
 	hap.snp.exon.comp.hets <- as.numeric( apply( hap.snp.exon.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	hap.ind.exon.comp.hets <- as.numeric( apply( hap.ind.exon.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	hap.exon.comp.hets <- as.numeric( apply( hap.exon.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	EXON[[name]][,"COMP_HET_snp"] <- hap.snp.exon.comp.hets
 	EXON[[name]][,"COMP_HET_ind"] <- hap.ind.exon.comp.hets
 	EXON[[name]][,"COMP_HET_any"] <- hap.exon.comp.hets

	####################################################
	## PLOT GENE & PHASING STATS #######################
	 # FULL: Boxplot HET_snp & HOM_VAR_snp & HET_ind & HOM_VAR_ind all in one
	 # EXON: Boxplot HET_snp & HOM_VAR_snp & HET_ind & HOM_VAR_ind all in one
	 # FULL: PRC_PHAS_snp & PRC_PHAS_ind & PRC_COMP_HET_snp & PRC_COMP_HET_ind
	 # EXON: PRC_PHAS_snp & PRC_PHAS_ind & PRC_COMP_HET_snp & PRC_COMP_HET_ind
	print( "Compiling Some Stats" )
	## Calculate % Compound Hets
	 # Full
	PRC.dat <- FULL[[name]][,c("PRC_PHAS_snp","PRC_PHAS_ind")]
	PRC_COMP_HET <- array(,c(2,2)) ; rownames(PRC_COMP_HET) <- c(0,1)
	PRC_COMP_HET_snp <- table( FULL[[name]][,"COMP_HET_snp"] ) / nrow(FULL[[name]])
	PRC_COMP_HET_ind <- table( FULL[[name]][,"COMP_HET_ind"] ) / nrow(FULL[[name]])
	PRC_COMP_HET[names(PRC_COMP_HET_snp),1] <- PRC_COMP_HET_snp[names(PRC_COMP_HET_snp)]
	PRC_COMP_HET[names(PRC_COMP_HET_ind),2] <- PRC_COMP_HET_ind[names(PRC_COMP_HET_ind)]
	PRC.dat.full <- PRC.dat
	PRC_COMP_HET.full <- PRC_COMP_HET
	 # Exon
	PRC.dat <- EXON[[name]][,c("PRC_PHAS_snp","PRC_PHAS_ind")]
	PRC_COMP_HET <- array(,c(2,2)) ; rownames(PRC_COMP_HET) <- c(0,1)
	PRC_COMP_HET_snp <- table( EXON[[name]][,"COMP_HET_snp"] ) / nrow(EXON[[name]])
	PRC_COMP_HET_ind <- table( EXON[[name]][,"COMP_HET_ind"] ) / nrow(EXON[[name]])
	PRC_COMP_HET[names(PRC_COMP_HET_snp),1] <- PRC_COMP_HET_snp[names(PRC_COMP_HET_snp)]
	PRC_COMP_HET[names(PRC_COMP_HET_ind),2] <- PRC_COMP_HET_ind[names(PRC_COMP_HET_ind)]
	PRC.dat.exon <- PRC.dat
	PRC_COMP_HET.exon <- PRC_COMP_HET
	## Compile Percent w/ Compound Het
	GENE[gtx,"p.COMP_HET_snp"] <- PRC_COMP_HET.full["1",1]
	GENE[gtx,"p.COMP_HET_ind"] <- PRC_COMP_HET.full["1",2]
	GENE[gtx,"p.COMP_HET_snp.ex"] <- PRC_COMP_HET.exon["1",1]
	GENE[gtx,"p.COMP_HET_ind.ex"] <- PRC_COMP_HET.exon["1",2]

	print( "Plotting Stats" )	
	if ( LIM[2]>4 | PLOT_RUNIF<PLOT_FRACTION ) {
		jpeg( paste(PathToOut,"/Plots/Gene_Stats_",name,".jpeg",sep=""), height=1600,width=2400, pointsize=32)
		par(mfrow=c(2,2))
		 # FULL: Boxplot
		COLS.cnt <- c("springgreen3","springgreen1","gold3","gold1")
		CNT.dat <- FULL[[name]][,c("HET_snp","HOM_VAR_snp","HET_ind","HOM_VAR_ind")]
		boxplot( CNT.dat, main=paste("Number Variants:",name), xlab="Category",ylab="# Vars in Individual", col=COLS.cnt, pch="" )
		for ( i in 1:4 ) { points( jitter(rep(i,nrow(CNT.dat)),amount=.1), CNT.dat[,i], pch="+" ) }
		 # EXON: Boxplot
		CNT.dat <- EXON[[name]][,c("HET_snp","HOM_VAR_snp","HET_ind","HOM_VAR_ind")]
		boxplot( CNT.dat, main=paste("Number Exonic Variants:",name), xlab="Category",ylab="# Exonic Vars in Individual", col=COLS.cnt, pch="" )
		for ( i in 1:4 ) { points( jitter(rep(i,nrow(CNT.dat)),amount=.1), CNT.dat[,i], pch="+" ) }
		 # FULL: Perc
		COLS.prc <- c("slateblue3","chocolate2","steelblue2")
		plot(0,0,type="n", xlim=c(0,4.5),ylim=c(0,1), main=paste("Percent Phased & Compound Hets:",name), xlab="",ylab="", xaxt="n",yaxt="n")
		abline( h=seq(0,1,.1), lty=2,lwd=1,col="grey50" )
		axis( 2, at=seq(0,1,.2) )
		axis( 4, at=seq(0,1,.2) )
		axis( 1, at=c(.5,1.5,3,4), labels=c("%Ph-SNP","%Ph-IND","%CH-SNP","%CH-IND") )
		boxplot( PRC.dat.full, at=c(.5,1.5), xaxt="n",yaxt="n", add=T, col=COLS.cnt[c(1,3)], pch="" )
		for ( i in 1:2 ) { points( jitter(rep(i-.5,nrow(PRC.dat.full)),amount=.1), PRC.dat.full[,i], pch="+" ) }
		barplot( PRC_COMP_HET.full, width=.8, space=c(3.25,.25), add=T, xaxt="n",yaxt="n", col=COLS.prc[2:3] )
		abline( v=2.25 )
		 # EXON: Perc
		plot(0,0,type="n", xlim=c(0,4.5),ylim=c(0,1), main=paste("Percent Phased & Compound Hets:",name), xlab="",ylab="", xaxt="n",yaxt="n")
		abline( h=seq(0,1,.1), lty=2,lwd=1,col="grey50" )
		axis( 2, at=seq(0,1,.2) )
		axis( 4, at=seq(0,1,.2) )
		axis( 1, at=c(.5,1.5,3,4), labels=c("%Ph-SNP","%Ph-IND","%CH-SNP","%CH-IND") )
		boxplot( PRC.dat.exon, at=c(.5,1.5), xaxt="n",yaxt="n", add=T, col=COLS.cnt[c(1,3)], pch="" )
		for ( i in 1:2 ) { points( jitter(rep(i-.5,nrow(PRC.dat.exon)),amount=.1), PRC.dat.exon[,i], pch="+" ) }
		barplot( PRC_COMP_HET.exon, width=.8, space=c(3.25,.25), add=T, xaxt="n",yaxt="n", col=COLS.prc[2:3] )
		abline( v=2.25 )
		dev.off()
	} # Close Compile Plot "IF"

	####################################################
	## UNIQUE HAPLOTYPES ###############################
	print( "Compiling Unique Haplotypes" )
	## Full
	N_hap <- ncol(hap)
	MAF <- rowMeans( data.matrix(hap[8:N_hap]) )
	UNIQ.all <- apply( hap[8:N_hap], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.all <- table( UNIQ.all )
	# sort( unname( TAB.uniq.all) ,decreasing=T )
	# hist( TAB.uniq.all, plot=F )$counts
	MAF.which.01 <- which( MAF > .01 )
	 # MAF .01
	hap.maf.01 <- hap[ MAF.which.01, ]
	UNIQ.maf.01 <- apply( hap.maf.01[8:N_hap], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.maf.01 <- table( UNIQ.maf.01 )
	 # MAF .05
	MAF.which.05 <- which( MAF > .05 )
	hap.maf.05 <- hap[ MAF.which.05, ]
	UNIQ.maf.05 <- apply( hap.maf.05[8:N_hap], 2, function(x) paste( x, collapse="" ) )
	TAB.uniq.maf.05 <- table( UNIQ.maf.05 )
	## Exon
	N_hap.exon <- ncol(hap.exon)
	MAF.exon <- rowMeans( data.matrix(hap.exon[8:N_hap.exon]) )
	VARS.exon.all <- hap.exon$TAG
	UNIQ.exon.all <- apply( hap.exon[8:N_hap.exon], 2, function(x) paste( x, collapse="" ) )
	TAB.exon.uniq.all <- table( UNIQ.exon.all )
	 # MAF.exon .01
	MAF.exon.which.01 <- which( MAF.exon > .01 )
	hap.exon.maf.01 <- hap.exon[ MAF.exon.which.01, ]
	UNIQ.exon.maf.01 <- apply( hap.exon.maf.01[8:N_hap.exon], 2, function(x) paste( x, collapse="" ) )
	TAB.exon.uniq.maf.01 <- table( UNIQ.exon.maf.01 )
	 # MAF.exon .05
	MAF.exon.which.05 <- which( MAF.exon > .05 )
	hap.exon.maf.05 <- hap.exon[ MAF.exon.which.05, ]
	UNIQ.exon.maf.05 <- apply( hap.exon.maf.05[8:N_hap.exon], 2, function(x) paste( x, collapse="" ) )
	TAB.exon.uniq.maf.05 <- table( UNIQ.exon.maf.05 )

	## Plot it
	COLS.full <- paste("steelblue",1:3,sep="")
	COLS.exon <- paste("slateblue",1:3,sep="") # c("tomato2","navyblue","gold1")
	if ( PLOT_RUNIF<PLOT_FRACTION ) {
		print( "Plotting Unique Haplotypes" )
		## Barplot Haplotype Frequencies
		jpeg( paste(PathToOut,"/Plots/Gene_HaploDistrib_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
		par(mfrow=c(2,3))
		 # Full
		barplot( sort(TAB.uniq.all,decreasing=T), main="Unique Haplotype Freq (Full/All)",xlab="Haplotype",las=2,col=COLS.full[1],border=NA,xaxt="n" )
		barplot( sort(TAB.uniq.maf.01,decreasing=T), main="Unique Haplotype Freq (Full/MAF>1%)",xlab="Haplotype",las=2,col=COLS.full[2],border=NA,xaxt="n" )
		barplot( sort(TAB.uniq.maf.05,decreasing=T), main="Unique Haplotype Freq (Full/MAF>5%)",xlab="Haplotype",las=2,col=COLS.full[3],border=NA,xaxt="n" )
		 # Exon
		barplot( sort(TAB.exon.uniq.all,decreasing=T), main="Unique Haplotype Freq (Exon/All)",xlab="Haplotype",las=2,col=COLS.exon[1],border=NA )
		barplot( sort(TAB.exon.uniq.maf.01,decreasing=T), main="Unique Haplotype Freq (Exon/1%)",xlab="Haplotype",las=2,col=COLS.exon[2],border=NA )
		barplot( sort(TAB.exon.uniq.maf.05,decreasing=T), main="Unique Haplotype Freq (Exon/5%)",xlab="Haplotype",las=2,col=COLS.exon[3],border=NA )
		dev.off()
		## Heatmap Showing Unique Haplotypes
		HAPLOS.arr <- matrix(as.numeric(unlist(sapply( names(TAB.exon.uniq.all), function(x) strsplit( x, ""), simplify="array" ))),nrow=length(TAB.exon.uniq.all),byrow=T )
		if ( all(dim(HAPLOS.arr)>=2) ) {
			colnames(HAPLOS.arr) <- VARS.exon.all
			MAFS.cols <- colorRampPalette(c("white","deepskyblue3"))(100)[ceiling(100*MAF.exon)]
			HAPLOS.cols <- colorRampPalette(c("white","chartreuse3"))(max(TAB.exon.uniq.all))[TAB.exon.uniq.all]
			jpeg( paste(PathToOut,"/Plots/Gene_HaploUniq_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
			heatmap.2( HAPLOS.arr, main="Unique Exonic Haplotypes",xlab="Variant Position",ylab="Haplotype",RowSideColors=HAPLOS.cols,ColSideColors=MAFS.cols,scale="none",trace="none",Rowv=T,Colv=F,dendrogram="row",col=c("black",COLS.exon[1]), lhei=c(1,6),lwid=c(1,6),margins=c(7,5) )
			dev.off()	
		}
	}
	
	####################################################
	## ASSOCIATION w/ PHENOTYPE ########################
	print( "Analyzing Associations w/ Phenotype" )
	if ( file.exists(PathToPheno) ) {
		MOD.covs <- lm( Pheno ~ . , data=PC.2[,-1] )
		RES.covs <- resid( MOD.covs )
		RES.covs.2 <- data.frame( IID=PC.2[,1],RES=RES.covs )

		## Unique Haplotypes vs Phenotype ##
		if ( nrow(hap.exon)>1 ) {
			Samps <- unique( PC.2$IID )
			N.Samps <- length(Samps)

			## Get Haplotypes for each Person
			HAPLOS.samp <- array( ,c(N.Samps,2) )
			colnames(HAPLOS.samp) <- paste("Hap",1:2,sep="_")
			rownames(HAPLOS.samp) <- Samps
			for ( s in 1:N.Samps ) {
				samp <- Samps[s]
				temp_columns <- grep( samp, colnames(hap.exon) )
				HAPLOS.samp[s,1] <- paste( hap.exon[,temp_columns[1]],collapse="" )
				HAPLOS.samp[s,2] <- paste( hap.exon[,temp_columns[2]],collapse="" )
				HAPLOS.samp[s,] <- sort(HAPLOS.samp[s,],decreasing=F)
			}
			# Put Clinical & Haplotype Data Together
			HAPLOS.samp.2 <- data.frame( IID=rep(rownames(HAPLOS.samp),2), HAP=c(HAPLOS.samp[,1],HAPLOS.samp[,2]), stringsAsFactors=F )
			PC.haps <- merge( PC.2, HAPLOS.samp.2, by="IID" )
			RES.covs.3 <- merge( RES.covs.2, HAPLOS.samp.2, by="IID" )
			# Specify "Rare" Haplotypes
			HAPLOS.rare <- names( which(TAB.exon.uniq.all<=5) )
			HAPLOS.samp.r <- HAPLOS.samp
			HAPLOS.samp.r[which(HAPLOS.samp.r[,"Hap_1"] %in% HAPLOS.rare),"Hap_1"] <- "Rare"
			HAPLOS.samp.r[which(HAPLOS.samp.r[,"Hap_2"] %in% HAPLOS.rare),"Hap_2"] <- "Rare"
			RES.covs.3[which(RES.covs.3[,"HAP"] %in% HAPLOS.rare),"HAP"] <- "Rare"
			# Run Analyses
			MOD.haps <- lm( Pheno ~ . , data=PC.haps[,-1] )
			P.haps <- anova(MOD.haps)["HAP","Pr(>F)"]
			GENE[gtx,"P_HAP"] <- P.haps
			MOD.haps.res <- lm( RES ~ HAP , data=RES.covs.3[,-1] )
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
		PC.ch <- merge( PC.2, EXON[[name]], by.x="IID",by.y="row.names" )
		PC.ch.2 <- PC.ch[,c(colnames(PC.2),"COMP_HET_any")]
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
		PC.burd.full <- merge( PC.2, FULL[[name]], by.x="IID",by.y="row.names" )
		SCORE.var <- PC.burd.full[,"HET_snp"] + PC.burd.full[,"HOM_VAR_snp"] + PC.burd.full[,"HET_ind"] + PC.burd.full[,"HOM_VAR_ind"]
		PC.burd.full.2 <- cbind( PC.burd.full[,colnames(PC.2)], SCORE.var )
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
				text( quantile(RES.burd.full.3$SCORE.var,.2),quantile(RES.burd.full.3$RES,.1), label=paste("P=",formatC(P.burd.full,format="e",digits=2)) )
				dev.off()
			}
		}
		 # Exon
		PC.burd.exon <- merge( PC.2, EXON[[name]], by.x="IID",by.y="row.names" )
		SCORE.var <- PC.burd.exon[,"HET_snp"] + PC.burd.exon[,"HOM_VAR_snp"] + PC.burd.exon[,"HET_ind"] + PC.burd.exon[,"HOM_VAR_ind"]
		PC.burd.exon.2 <- cbind( PC.burd.exon[,colnames(PC.2)], SCORE.var )
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
				text( quantile(RES.burd.exon.3$SCORE.var,.2),quantile(RES.burd.exon.3$RES,.1), label=paste("P=",formatC(P.burd.exon,format="e",digits=2)) )
				dev.off()
			}
		}
	}

	####################################################
	## EVERY FEW ITERATIONS, SAVE TABLES/DATA ##########
	if ( gtx%%50==0 ) {
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

print("Writing Tables")

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


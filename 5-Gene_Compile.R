## Script to Loop Through Genes and Compile/Plots Things ##
## Called by Gene_Score Pipeline ##
## February 16, 2015 ##
## Kristopher Standish ##


## Usage ##
# Rscript 5-Gene_Compile.R <Path/To/Out_Dir> <Path/To/Pheno> <Covs_Command>


###############################################################
## PARSE COMMAND LINE #########################################
###############################################################


LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( "/projects/janssen/Phased/20150216_Testing","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
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

###############################################################
## LOAD DATA ##################################################
###############################################################

## Load Phenotype/Covariate Files (if exists)
if ( file.exists( PathToPheno ) ) {
	# Load Phenotype/Covariate/Assoc Files
	PHENO <- read.table( PathToPheno, sep="\t",header=T )
	COVS.l <- read.table( PathToCovFile, sep="\t",header=T )
	# Load SNP Association Results
	SNP.bim <- read.table( paste(PathToOut,"/SNP_Vars.bim",sep=""), sep="\t",header=F)
	colnames(SNP.bim) <- c("CHR","SNP","XXX","BP","REF","ALT")
	SNP.P <- read.table( paste(PathToAssoc,"SNP/SNP_Assoc.P",sep=""), sep="\t",header=T)
	SNP.HWE <- read.table( paste(PathToAssoc,"SNP/SNP_Assoc.hwe",sep=""), sep="",header=T)
	# Load IND Association Results
	IND.bim <- read.table( paste(PathToOut,"/IND_Vars.bim",sep=""), sep="\t",header=F)
	colnames(IND.bim) <- c("CHR","SNP","XXX","BP","REF","ALT")
	IND.P <- read.table( paste(PathToAssoc,"IND/IND_Assoc.P",sep=""), sep="\t",header=T)
	IND.HWE <- read.table( paste(PathToAssoc,"IND/IND_Assoc.hwe",sep=""), sep="",header=T)
	# Reformat Covariate Table
	Cov_List.sp <- strsplit( Cov_List, "," )[[1]]
	COVS <- COVS.l[, c("IID",Cov_List.sp) ]
	# Merge Phenotype/Covariate Files
	PC <- merge( x=PHENO[,c("IID","Pheno")], y=COVS, by="IID" )
}else{ print("No Phenotype Provided") }

## Load Gene Coords
GENE_COORDS <- read.table( PathToGeneCoords, sep="\t", header=T, comment.char="" )
colnames(GENE_COORDS) <- gsub("X.","", colnames(GENE_COORDS), fixed=T )
colnames(GENE_COORDS) <- gsub("hg19.knownGene.","", colnames(GENE_COORDS), fixed=T )
colnames(GENE_COORDS) <- gsub("hg19.kgXref.","", colnames(GENE_COORDS), fixed=T )

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
hap.samps <- as.character( read.table( paste(PathToOut,"/Phased.sample",sep=""), sep="",header=T )[,1] )
n.samps <- length(hap.samps)
hap.colnames <- c("CHR","SNP","BP","REF","ALT", paste( rep(hap.samps, rep(2,n.samps)), 1:2, sep="_" ) )
## Loop Through Genes
for ( gtx in 1:n.gtx ) {
# for ( gtx in 10:14 ) {
	## Compile Info on Gene_Transcript
	name <- GTX_LIST[gtx]
	chr <- gsub( "chr","", as.character( CHR[gtx] ) )
	rng <- as.numeric( RNG_TAB[gtx,] )
	tx_rng <- as.numeric( TX_TAB[gtx,] )
	cd_rng <- as.numeric( CD_TAB[gtx,] )
	ex_b <- as.numeric( EX_B_LIST[[gtx]] )
	ex_e <- as.numeric( EX_E_LIST[[gtx]] )
	ex_cnt <- as.numeric( EX_COUNT[gtx] )

	####################################################
	## Load & Organize Transcript Data/Files
	
	## Load Genotype/HW Files
	snp.raw <- read.table( paste(PathToGenes,name,"/SNP_Vars.raw",sep=""), sep="",header=T )
	snp.hwe <- read.table( paste(PathToGenes,name,"/SNP_Vars.hwe",sep=""), sep="",header=T )
	ind.raw <- read.table( paste(PathToGenes,name,"/IND_Vars.raw",sep=""), sep="",header=T )
	ind.hwe <- read.table( paste(PathToGenes,name,"/IND_Vars.hwe",sep=""), sep="",header=T )
	gtx.hwe <- rbind( snp.hwe, ind.hwe )
	## Load Phased Haplotype File (& Rename Columns)
	hap <- read.table( paste(PathToGenes,name,"/Phased.haps",sep=""), sep="",header=F )
	colnames(hap) <- hap.colnames
	## Pull out Variant Positions from BIM files (& Compile to 1 file)
	gtx.snp.bim <- SNP.bim[ which(SNP.bim$CHR==chr & SNP.bim$BP>rng[1] & SNP.bim$BP<rng[2] ), ]
	gtx.snp.bim <- data.frame( gtx.snp.bim, TYPE="snp" )
	gtx.ind.bim <- IND.bim[ which(IND.bim$CHR==chr & IND.bim$BP>rng[1] & IND.bim$BP<rng[2] ), ]
	gtx.ind.bim <- data.frame( gtx.ind.bim, TYPE="ind" )
	gtx.bim <- rbind( gtx.snp.bim, gtx.ind.bim )
	gtx.bim <- gtx.bim[ order(gtx.bim[,"BP"]), ]
	## Merge BIM table w/ HWE table
	gtx.mg <- merge( gtx.bim[,c("SNP","BP","TYPE")], gtx.hwe, by="SNP", all=T )
	gtx.mg <- data.frame( TAG=paste(gtx.mg[,"CHR"],gtx.mg[,"BP"],sep="_"), gtx.mg )
	gtx.mg <- gtx.mg[ ,c("TAG","CHR","BP","SNP","TYPE","A1","A2","GENO","P") ]
	colnames(gtx.mg)[ncol(gtx.mg)] <- "P_HW"
	GTX <- gtx.mg
	## Compile which variants are in Exons
	TAG.exon <- c()
	for ( e in 1:ex_cnt ) {
		WHICH <- which( gtx.mg[,"BP"]>ex_b[e] & gtx.mg[,"BP"]<ex_e[e] )
		if ( length(WHICH)>0) { TAG.exon <- c( TAG.exon, as.character(gtx.mg[WHICH,"TAG"]) ) }
	}
	## Pull out Single-Locus Results
	if ( file.exists(PathToPheno) ) {
		snp.p <- SNP.P[ which(SNP.P$CHR==chr & SNP.P$BP>tx_rng[1] & SNP.P$BP<tx_rng[2] ), ]
		ind.p <- IND.P[ which(IND.P$CHR==chr & IND.P$BP>tx_rng[1] & IND.P$BP<tx_rng[2] ), ]
		gtx.p <- rbind( snp.p, ind.p )
		gtx.p <- gtx.p[order(gtx.p[,"BP"]),]
		gtx.p <- data.frame( TAG=paste(gtx.p[,"CHR"],gtx.p[,"BP"],sep="_"), gtx.p )
		# Merge with HWE Data (compile all PLink data)
		gtx.mg.p <- merge( gtx.mg, gtx.p[,c("SNP","P")], by="SNP", all=T )
		colnames(gtx.mg.p)[which(colnames(gtx.mg.p)=="P")] <- "P_Assoc"
		gtx.mg.p <- gtx.mg.p[ order(gtx.mg.p[,"BP"]), ]
		gtx.mg.p <- gtx.mg.p[ ,c("TAG","CHR","BP","SNP","A1","A2","P_Assoc","GENO","P_HW") ]
		GTX <- gtx.mg.p
	}
	## Specify Location of Variant (relative to Exon/Intron/Etc)
	GTX <- data.frame( GTX, LOC=rep(1,nrow(GTX)) )
	GTX[ which( GTX[,"BP"]>tx_rng[1] & GTX[,"BP"]<tx_rng[2] ), "LOC" ] <- 2
	GTX[ which( GTX[,"BP"]>cd_rng[1] & GTX[,"BP"]<cd_rng[2] ), "LOC" ] <- 3
	GTX[ which( GTX[,"TAG"] %in% TAG.exon ), "LOC" ] <- 4
	## Calculate Allele Frequencies
	gtx.geno.arr <- t(sapply( strsplit( as.character(GTX[,"GENO"]), "/" ), "[", 1:3 ))
	gtx.af <- matrix( as.numeric(gtx.geno.arr), ncol=3 ) %*% matrix(2:0,c(3,1)) / (2*n.samps)
	GTX <- data.frame( GTX, MAF=gtx.af )

	####################################################
	## Plot this Shiz

	## Set Plotting Parameters
	NUM_VARS <- nrow(GTX)
	COLS.gene <- c("slateblue3","cadetblue1","steelblue2","firebrick2","goldenrod2","springgreen3")
	COLS.P.list <- c("firebrick4","slateblue3","cadetblue3","steelblue2")
	XLIM <- c(rng[1],rng[2])
	if ( file.exists(PathToPheno) ) {
 		YLIM <- c( -1, max(4,max(GTX[,"P_Assoc"],na.rm=T)) )
 	}else{ YLIM <- c(-1,1) }
	Y <- 0
	jpeg( paste(PathToOut,"/Plots/Gene_Map_",name,".jpeg",sep=""), height=1400,width=2000, pointsize=30)
	plot(0,0,type="n", xlim=XLIM,ylim=YLIM, xlab=paste("Chromosome",chr,"Position"),ylab="-log10(p)", main=paste("Map & Single-Locus Results of:",name), yaxt="n" )
	## Plot Gene Borders
	arrows( tx_rng[1],Y,tx_rng[2],Y, code=3,angle=90, lwd=5,col=COLS.gene[1] )
	abline( h=Y, col="grey50",lwd=1 )
	polygon( cd_rng[c(1,2,2,1)],c(-.1,-.1,.1,.1)+Y, col=COLS.gene[2],lwd=1 )
	for ( e in 1:ex_cnt ) {
		X_COORDS <- c( ex_b[e],ex_e[e],ex_e[e],ex_b[e] )
		polygon( X_COORDS,c(-.2,-.2,.2,.2)+Y, col=COLS.gene[3],lwd=1 )	
	}
	arrows( tx_rng[1],Y,tx_rng[2],Y, code=3,angle=90,length=0, lwd=5,col=COLS.gene[1] )
	## Plot P-Values
	if ( file.exists(PathToPheno) ) {
		WHICH <- which(!is.na( GTX$P_Assoc))
		COLS.P <- COLS.P.list[ GTX$LOC ]
		abline( h=seq(0,YLIM[2],1), lty=2, col="grey50" )
		abline( h=seq(-1,0,.1), lty=2, col="grey50" )
		axis(2, at=seq(-1,YLIM[2],1) )
		points( GTX$BP[WHICH], -log10(GTX$P_Assoc[WHICH]), col=COLS.P[WHICH], pch="+", cex=2*(GTX$MAF[WHICH])^(1/6) )	
	}
	## Plot Allele Frequency
	arrows( GTX$BP, 0, GTX$BP, -GTX$MAF, col=COLS.P, angle=90 )
	## Plot Individual Variants on Haplotypes
	# if ( length(HAP_1)>0 ) { arrows( HAP_1,Y-.3,HAP_1,Y, lwd=2,col=COLS[5],length=.1 ) }
	# if ( length(HAP_2)>0 ) { arrows( HAP_2,Y+.3,HAP_2,Y, lwd=2,col=COLS[4],length=.1 ) }
	# if ( length(intersect(HAP_1,HAP_2))>0 ) {
	# 	arrows( intersect(HAP_1,HAP_2),Y-.3,intersect(HAP_1,HAP_2),Y, lwd=2,col=COLS[6],length=.1 )
	# 	arrows( intersect(HAP_1,HAP_2),Y+.3,intersect(HAP_1,HAP_2),Y, lwd=2,col=COLS[6],length=.1 )
	# }
	# text( X_LIM[1], Y-.4, labels=name, pos=4, cex=.7)
	dev.off()

	####################################################
	## Compile Stats about Variants/Cohort
	


	print(paste("Done with",gtx,"of",n.gtx,"-",name))
}



























###############################################################
## END OF DOC #################################################
###############################################################

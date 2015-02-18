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
# LINE <- c( "/projects/janssen/Phased/20150217_Testing","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt","DAS_BL_MN,PC1,PC2" )
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
GENE_COORDS <- GENE_COORDS[ which(GENE_COORDS[,"chrom"]!="chrX"), ]

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

## Set up Objects to Compile
 # Basic Gene Info
GENE.cats <- c("GTX","TX_S","TX_E","CD_S","CD_E","EX_S","EX_E","n.EX","n.VAR","n.SNP","n.IND","n.VAR.ph","n.SNP.ph","n.IND.ph","n.VAR.ex","n.SNP.ex","n.IND.ex","n.VAR.ph.ex","n.SNP.ph.ex","n.IND.ph.ex","AREA","AREA.ex","BEST_P","BEST_P.ex","p.COMP_HET_snp","p.COMP_HET_ind","p.COMP_HET_snp.ex","p.COMP_HET_ind.ex")
GENE <- array( , c(n.gtx,length(GENE.cats)) )
colnames(GENE) <- GENE.cats
GENE[,"GTX"] <- GTX_LIST
GENE[,"TX_S"] <- TX_TAB[,"TX_S"]
GENE[,"TX_E"] <- TX_TAB[,"TX_E"]
GENE[,"CD_S"] <- CD_TAB[,"CD_S"]
GENE[,"CD_E"] <- CD_TAB[,"CD_E"]
GENE[,"EX_S"] <- as.character( GENE_COORDS$exonStarts )
GENE[,"EX_E"] <- as.character( GENE_COORDS$exonEnds )
GENE[,"n.EX"] <- EX_COUNT
 # Cohort Data
FULL <- EXON <- DAMG <- list()

## Loop Through Genes
start_time <- proc.time()
for ( gtx in 1:n.gtx ) {
# for ( gtx in 1:50 ) {
	## Compile Info on Gene_Transcript
	name <- GTX_LIST[gtx]
	chr <- gsub( "chr","", as.character( CHR[gtx] ) )
	rng <- as.numeric( RNG_TAB[gtx,] )
	tx_rng <- as.numeric( TX_TAB[gtx,] )
	cd_rng <- as.numeric( CD_TAB[gtx,] )
	ex_b <- as.numeric( EX_B_LIST[[gtx]] )
	ex_e <- as.numeric( EX_E_LIST[[gtx]] )
	ex_cnt <- as.numeric( EX_COUNT[gtx] )

	print(paste( "### Starting Loop",gtx,"-",name,"#####",round(proc.time()-start_time,2)[3] ))
	####################################################
	## Load & Organize Transcript Data/Files
	
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
	## Pull out Variant Positions from BIM files (& Compile to 1 file)
	print( "Sorting Bim Files" )
	gtx.snp.bim <- SNP.bim[ which(SNP.bim$CHR==chr & SNP.bim$BP>rng[1] & SNP.bim$BP<rng[2] ), ]
	gtx.snp.bim <- data.frame( gtx.snp.bim, TYPE=rep("snp",nrow(gtx.snp.bim)), RAW_TAG=snp.raw.tag )
	gtx.ind.bim <- IND.bim[ which(IND.bim$CHR==chr & IND.bim$BP>rng[1] & IND.bim$BP<rng[2] ), ]
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
	## Compile which variants are in Exons
	print( "Pulling out Exonic Variants" )
	TAG.exon <- c()
	for ( e in 1:ex_cnt ) {
		WHICH <- which( gtx.mg[,"BP"]>ex_b[e] & gtx.mg[,"BP"]<ex_e[e] )
		if ( length(WHICH)>0) { TAG.exon <- c( TAG.exon, as.character(gtx.mg[WHICH,"TAG"]) ) }
	}
	## Pull out Single-Locus Results
	print( "Pulling out Single-Locus Results" )
	if ( file.exists(PathToPheno) ) {
		snp.p <- SNP.P[ which(SNP.P$CHR==chr & SNP.P$BP>tx_rng[1]-5000 & SNP.P$BP<tx_rng[2]+5000 ), ]
		ind.p <- IND.P[ which(IND.P$CHR==chr & IND.P$BP>tx_rng[1]-5000 & IND.P$BP<tx_rng[2]+5000 ), ]
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
	GTX[ which( GTX[,"BP"]>tx_rng[1] & GTX[,"BP"]<tx_rng[2] ), "LOC" ] <- 2
	GTX[ which( GTX[,"BP"]>cd_rng[1] & GTX[,"BP"]<cd_rng[2] ), "LOC" ] <- 3
	GTX[ which( GTX[,"TAG"] %in% TAG.exon ), "LOC" ] <- 4
	## Calculate Allele Frequencies
	print( "Calculating MAF" )
	gtx.geno.arr <- t(sapply( strsplit( as.character(GTX[,"GENO"]), "/" ), "[", 1:3 ))
	gtx.af <- matrix( as.numeric(gtx.geno.arr), ncol=3 ) %*% matrix(2:0,c(3,1)) / (2*n.samps)
	GTX <- data.frame( GTX, MAF=gtx.af )
	## Merge with Haplotype File
	# GTX.phased <- merge( GTX, hap, by="TAG", all=F )

	####################################################
	## Plot This Shiz

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
	## Plot Individual Variants on Haplotypes
	# if ( length(HAP_1)>0 ) { arrows( HAP_1,Y-.3,HAP_1,Y, lwd=2,col=COLS[5],length=.1 ) }
	# if ( length(HAP_2)>0 ) { arrows( HAP_2,Y+.3,HAP_2,Y, lwd=2,col=COLS[4],length=.1 ) }
	# if ( length(intersect(HAP_1,HAP_2))>0 ) {
	# 	arrows( intersect(HAP_1,HAP_2),Y-.3,intersect(HAP_1,HAP_2),Y, lwd=2,col=COLS[6],length=.1 )
	# 	arrows( intersect(HAP_1,HAP_2),Y+.3,intersect(HAP_1,HAP_2),Y, lwd=2,col=COLS[6],length=.1 )
	# }
	# text( X_LIM[1], Y-.4, labels=name, pos=4, cex=.7)
	dev.off()

	## QQ PLOTS ##
	print( "Plotting QQ Plot" )
	if ( file.exists(PathToPheno) ) {
		NUM_VARS <- length(which( !is.na(GTX$P_Assoc) ))
		NUM_EXONIC <- length(which( !is.na(GTX$P_Assoc) & GTX$LOC==4 ))
		COLS <- c("steelblue1","slateblue1")
		COLS.4 <- gsub("1","3",COLS)
		LIM <- c( 0, max( 4, max(-log10(GTX[,"P_Assoc"]),na.rm=T) ) )
		## Make Plot
		jpeg( paste(PathToOut,"/Plots/Gene_QQ_",name,".jpeg",sep=""), height=1500,width=1500, pointsize=32)
		plot(0,0,type="n", xlim=LIM,ylim=LIM, xlab="Expected -log10(p)",ylab="Observed -log10(p)", main=paste("QQ-Plot for:",name) )
		## Add Lines
		abline( h=seq(0,LIM[2],1), lty=2,lwd=1,col="grey50")
		abline( v=seq(0,LIM[2],1), lty=2,lwd=1,col="grey50")
		abline( 0,1, lty=1,lwd=2,col="black" )
		## Plot Full Set
		EXP <- -log10( 1:NUM_VARS / NUM_VARS )
		OBS <- -log10( sort( GTX$P_Assoc ) )
		IND <- c(1:length(EXP),length(EXP),1)
		polygon( EXP[IND], c(OBS,EXP[c(NUM_VARS,1)]), col=COLS[1], border=COLS.4[1], density=20,angle=45 )
		points( EXP, OBS, pch="+", col=COLS.4[1], type="o" )
		 # Calculate Area - Trapezoid Rule
		H_VALS <- EXP[2:NUM_VARS-1] - EXP[2:NUM_VARS]
		B_SUM <- ( OBS[2:NUM_VARS-1]-EXP[2:NUM_VARS-1] ) + ( OBS[2:NUM_VARS]-EXP[2:NUM_VARS] )
		AREA.traps <- .5*H_VALS*B_SUM
		AREA <- sum( AREA.traps ) # - .5*max(EXP)^2
		text( quantile(LIM,.8),quantile(LIM,.1), label=paste("Area:",round(AREA,3)), col=COLS.4[1], cex=1.2 )
		GENE[gtx,"AREA"] <- AREA
		GENE[gtx,"BEST_P"] <- 10^(min(-OBS,na.rm=T))
		## Plot Exon Set
		if ( NUM_EXONIC > 0 ) {
			EXP <- -log10( 1:NUM_EXONIC / NUM_EXONIC )
			OBS <- -log10( sort( GTX$P_Assoc[which(GTX$LOC==4)] ) )
			IND <- c(1:length(EXP),length(EXP),1)
			polygon( EXP[IND], c(OBS,EXP[c(NUM_EXONIC,1)]), col=COLS[2], border=COLS.4[2], density=20,angle=-45 )
			points( EXP, OBS, pch="+", col=COLS.4[2], type="o" )	
			 # Calculate Area - Trapezoid Rule
			H_VALS <- EXP[2:NUM_EXONIC-1] - EXP[2:NUM_EXONIC]
			B_SUM <- ( OBS[2:NUM_EXONIC-1]-EXP[2:NUM_EXONIC-1] ) + ( OBS[2:NUM_EXONIC]-EXP[2:NUM_EXONIC] )
			AREA.traps <- .5*H_VALS*B_SUM
			AREA <- sum( AREA.traps ) # - .5*max(EXP)^2
			text( quantile(LIM,.8),quantile(LIM,.07), label=paste("Area:",round(AREA,3)), col=COLS.4[2], cex=1.2 )
			GENE[gtx,"AREA.ex"] <- AREA
			GENE[gtx,"BEST_P.ex"] <- 10^(min(-OBS,na.rm=T))
		}
		dev.off()
	}

	####################################################
	## Compile Stats about Variants/Cohort
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
 	GENE[gtx,"n.SNP.ph"] <- length(which( hap[,"TYPE"]=="ind" ))
 	GENE[gtx,"n.IND.ph"] <- length(which( hap[,"TYPE"]=="ind" ))
 	 # Exon
 	GTX.exon <- GTX[ which(GTX$LOC==4), ]
 	hap.exon <- hap[ which(hap$TAG %in% TAG.exon), ]
 	GENE[gtx,"n.VAR.ex"] <- nrow(GTX.exon)
 	GENE[gtx,"n.SNP.ex"] <- length(which( GTX.exon$TYPE=="snp" ))
 	GENE[gtx,"n.IND.ex"] <- length(which( GTX.exon$TYPE=="ind" ))
 	GENE[gtx,"n.VAR.ph.ex"] <- nrow(hap.exon)
 	GENE[gtx,"n.SNP.ph.ex"] <- length(which( hap.exon[,"TYPE"]=="snp" ))
 	GENE[gtx,"n.IND.ph.ex"] <- length(which( hap.exon[,"TYPE"]=="ind" ))

 	## Cohort - Full ##
 	print( "Compiling Cohort Stats (Full)" )
 	FULL[[name]] <- array( ,c(n.samps,8) )
 	rownames(FULL[[name]]) <- shared.samps
 	colnames(FULL[[name]]) <- c("HET_snp","HOM_VAR_snp","PRC_PHAS_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","PRC_PHAS_ind","COMP_HET_ind")
 	 # How many Hets?
 	FULL[[name]][,"HET_snp"] <- apply( snp.raw[order(snp.raw[,"IID"]),7:ncol(snp.raw)], 1, function(x) length(which(x==1)) )
 	FULL[[name]][,"HET_ind"] <- apply( ind.raw[order(ind.raw[,"IID"]),7:ncol(ind.raw)], 1, function(x) length(which(x==1)) )
 	 # How many Hom_Vars?
 	FULL[[name]][,"HOM_VAR_snp"] <- apply( snp.raw[order(snp.raw[,"IID"]),7:ncol(snp.raw)], 1, function(x) length(which(x==2)) )
 	FULL[[name]][,"HOM_VAR_ind"] <- apply( ind.raw[order(ind.raw[,"IID"]),7:ncol(ind.raw)], 1, function(x) length(which(x==2)) )
 	 # What percent of variants were phased?
 	hap.1.cols <- seq( which(colnames(hap)=="ALT")+1,ncol(hap),2)
 	hap.2.cols <- seq( which(colnames(hap)=="ALT")+2,ncol(hap),2)
 	hap.snp <- hap[ which(hap$TYPE=="snp") , ]
 	hap.snp.diffs <- hap.snp[,hap.1.cols] - hap.snp[,hap.2.cols]
 	colnames(hap.snp.diffs) <- gsub( "_1","", colnames(hap.snp.diffs) )
 	hap.snp.n.hets <- apply( hap.snp.diffs, 2, function(x) length(which(x!=0)) )
 	FULL[[name]][,"PRC_PHAS_snp"] <- hap.snp.n.hets / FULL[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.n.hets / FULL[[name]][,"HET_snp"]
 	hap.ind <- hap[ which(hap$TYPE=="ind") , ]
 	hap.ind.diffs <- hap.ind[,hap.1.cols] - hap.ind[,hap.2.cols]
 	colnames(hap.ind.diffs) <- gsub( "_1","", colnames(hap.ind.diffs) )
 	hap.ind.n.hets <- apply( hap.ind.diffs, 2, function(x) length(which(x!=0)) )
 	FULL[[name]][,"PRC_PHAS_ind"] <- hap.ind.n.hets / FULL[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.n.hets / FULL[[name]][,"HET_ind"]
 	# SAMP <- "B012326"
	# SAMP.hets.hap <- as.character( hap.snp[ which( hap.snp.diffs[,SAMP]!=0 ), "TAG" ])
	# SAMP.hets.raw <- colnames(snp.raw)[ which( snp.raw[ which(snp.raw[,"IID"]==SAMP), ]==1 ) ]
	# SAMP.hets.bim <- as.character( GTX[ which(GTX[,"RAW_TAG"] %in% SAMP.hets.raw), "TAG" ] )
 	 # Which Samples have a Compound Het?
 	hap.snp.comp.hets <- as.numeric( apply( hap.snp.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	hap.ind.comp.hets <- as.numeric( apply( hap.ind.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	FULL[[name]][,"COMP_HET_snp"] <- hap.snp.comp.hets
 	FULL[[name]][,"COMP_HET_ind"] <- hap.ind.comp.hets

 	## Cohort - Exon ##
 	print( "Compiling Cohort Stats (Exon)" )
 	EXON[[name]] <- array( ,c(n.samps,8) )
 	rownames(EXON[[name]]) <- shared.samps
 	colnames(EXON[[name]]) <- c("HET_snp","HOM_VAR_snp","PRC_PHAS_snp","COMP_HET_snp","HET_ind","HOM_VAR_ind","PRC_PHAS_ind","COMP_HET_ind")
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
 	hap.snp <- hap[ which(hap$TYPE=="snp") , ]
 	hap.snp.exon <- hap.snp[ which(hap.snp[,"TAG"]%in%exon.tag), ]
 	hap.snp.exon.diffs <- hap.snp.exon[,hap.1.cols] - hap.snp.exon[,hap.2.cols]
 	colnames(hap.snp.exon.diffs) <- gsub( "_1","", colnames(hap.snp.exon.diffs) )
 	hap.snp.exon.n.hets <- apply( hap.snp.exon.diffs, 2, function(x) length(which(x!=0)) )
 	EXON[[name]][,"PRC_PHAS_snp"] <- hap.snp.exon.n.hets / EXON[[name]][,"HET_snp"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.snp.exon.n.hets / EXON[[name]][,"HET_snp"]
  	hap.ind <- hap[ which(hap$TYPE=="ind") , ]
 	hap.ind.exon <- hap.ind[ which(hap.ind[,"TAG"]%in%exon.tag), ]
 	hap.ind.exon.diffs <- hap.ind.exon[,hap.1.cols] - hap.ind.exon[,hap.2.cols]
 	colnames(hap.ind.exon.diffs) <- gsub( "_1","", colnames(hap.ind.exon.diffs) )
 	hap.ind.exon.n.hets <- apply( hap.ind.exon.diffs, 2, function(x) length(which(x!=0)) )
 	EXON[[name]][,"PRC_PHAS_ind"] <- hap.ind.exon.n.hets / EXON[[name]][,"HET_ind"] # length(intersect(SAMP.hets.hap,SAMP.hets.bim)) / length(SAMP.hets.bim) # hap.ind.exon.n.hets / EXON[[name]][,"HET_ind"]
 	 # Which Samples have a Compound Het?
 	hap.snp.exon.comp.hets <- as.numeric( apply( hap.snp.exon.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	hap.ind.exon.comp.hets <- as.numeric( apply( hap.ind.exon.diffs, 2, function(x) length(which( -1%in%x & 1%in%x ))>0 ) )
 	EXON[[name]][,"COMP_HET_snp"] <- hap.snp.exon.comp.hets
 	EXON[[name]][,"COMP_HET_ind"] <- hap.ind.exon.comp.hets

	####################################################
	## Plot This Shiz
	 # FULL: Boxplot HET_snp & HOM_VAR_snp & HET_ind & HOM_VAR_ind all in one
	 # EXON: Boxplot HET_snp & HOM_VAR_snp & HET_ind & HOM_VAR_ind all in one
	 # FULL: PRC_PHAS_snp & PRC_PHAS_ind & PRC_COMP_HET_snp & PRC_COMP_HET_ind
	 # EXON: PRC_PHAS_snp & PRC_PHAS_ind & PRC_COMP_HET_snp & PRC_COMP_HET_ind
	print( "Plotting Stats" )
	## Try all in one
	jpeg( paste(PathToOut,"/Plots/Gene_Stats_",name,".jpeg",sep=""), height=1600,width=2400, pointsize=32)
	par(mfrow=c(2,2))
	 # FULL: Boxplot
	COLS.cnt <- c("springgreen3","springgreen1","gold3","gold1")
	CNT.dat <- FULL[[name]][,c("HET_snp","HOM_VAR_snp","HET_ind","HOM_VAR_ind")]
	boxplot( CNT.dat, main=paste("Number Variants:",name), xlab="Category",ylab="# Vars in Individual", col=COLS.cnt )
	for ( i in 1:4 ) { points( jitter(rep(i,nrow(CNT.dat)),amount=.1), CNT.dat[,i], pch="+" ) }
	 # EXON: Boxplot
	CNT.dat <- EXON[[name]][,c("HET_snp","HOM_VAR_snp","HET_ind","HOM_VAR_ind")]
	boxplot( CNT.dat, main=paste("Number Exonic Variants:",name), xlab="Category",ylab="# Exonic Vars in Individual", col=COLS.cnt )
	for ( i in 1:4 ) { points( jitter(rep(i,nrow(CNT.dat)),amount=.1), CNT.dat[,i], pch="+" ) }
	 # FULL: Perc
	COLS.prc <- c("slateblue3","chocolate2","steelblue2")
	plot(0,0,type="n", xlim=c(0,4.5),ylim=c(0,1), main=paste("Percent Phased & Compound Hets:",name), xlab="",ylab="", xaxt="n",yaxt="n")
	abline( h=seq(0,1,.1), lty=2,lwd=1,col="grey50" )
	axis( 2, at=seq(0,1,.2) )
	axis( 4, at=seq(0,1,.2) )
	axis( 1, at=c(.5,1.5,3,4), labels=c("%Ph-SNP","%Ph-IND","%CH-SNP","%CH-IND") )
	PRC.dat <- FULL[[name]][,c("PRC_PHAS_snp","PRC_PHAS_ind")]
	boxplot( PRC.dat, at=c(.5,1.5), xaxt="n",yaxt="n", add=T, col=COLS.prc[1] )
	for ( i in 1:2 ) { points( jitter(rep(i-.5,nrow(PRC.dat)),amount=.1), PRC.dat[,i], pch="+" ) }
	PRC_COMP_HET <- array(,c(2,2)) ; rownames(PRC_COMP_HET) <- c(0,1)
	PRC_COMP_HET_snp <- table( FULL[[name]][,"COMP_HET_snp"] ) / nrow(FULL[[name]])
	PRC_COMP_HET_ind <- table( FULL[[name]][,"COMP_HET_ind"] ) / nrow(FULL[[name]])
	PRC_COMP_HET[names(PRC_COMP_HET_snp),1] <- PRC_COMP_HET_snp[names(PRC_COMP_HET_snp)]
	PRC_COMP_HET[names(PRC_COMP_HET_ind),2] <- PRC_COMP_HET_ind[names(PRC_COMP_HET_ind)]
	barplot( PRC_COMP_HET, width=.8, space=c(3.25,.25), add=T, xaxt="n",yaxt="n", col=COLS.prc[2:3] )
	abline( v=2.25 )
	## Compile Percept w/ Compound Het
	GENE[gtx,"p.COMP_HET_snp"] <- PRC_COMP_HET["1",1]
	GENE[gtx,"p.COMP_HET_ind"] <- PRC_COMP_HET["1",2]
	 # EXON: Perc
	plot(0,0,type="n", xlim=c(0,4.5),ylim=c(0,1), main=paste("Percent Phased & Compound Hets:",name), xlab="",ylab="", xaxt="n",yaxt="n")
	abline( h=seq(0,1,.1), lty=2,lwd=1,col="grey50" )
	axis( 2, at=seq(0,1,.2) )
	axis( 4, at=seq(0,1,.2) )
	axis( 1, at=c(.5,1.5,3,4), labels=c("%Ph-SNP","%Ph-IND","%CH-SNP","%CH-IND") )
	PRC.dat <- EXON[[name]][,c("PRC_PHAS_snp","PRC_PHAS_ind")]
	boxplot( PRC.dat, at=c(.5,1.5), xaxt="n",yaxt="n", add=T, col=COLS.prc[1] )
	for ( i in 1:2 ) { points( jitter(rep(i-.5,nrow(PRC.dat)),amount=.1), PRC.dat[,i], pch="+" ) }
	PRC_COMP_HET <- array(,c(2,2)) ; rownames(PRC_COMP_HET) <- c(0,1)
	PRC_COMP_HET_snp <- table( EXON[[name]][,"COMP_HET_snp"] ) / nrow(EXON[[name]])
	PRC_COMP_HET_ind <- table( EXON[[name]][,"COMP_HET_ind"] ) / nrow(EXON[[name]])
	PRC_COMP_HET[names(PRC_COMP_HET_snp),1] <- PRC_COMP_HET_snp[names(PRC_COMP_HET_snp)]
	PRC_COMP_HET[names(PRC_COMP_HET_ind),2] <- PRC_COMP_HET_ind[names(PRC_COMP_HET_ind)]
	barplot( PRC_COMP_HET, width=.8, space=c(3.25,.25), add=T, xaxt="n",yaxt="n", col=COLS.prc[2:3] )
	abline( v=2.25 )
	dev.off()
	## Compile Percept w/ Compound Het
	GENE[gtx,"p.COMP_HET_snp.ex"] <- PRC_COMP_HET["1",1]
	GENE[gtx,"p.COMP_HET_ind.ex"] <- PRC_COMP_HET["1",2]
	## Next...
	print(paste("Done with",gtx,"of",n.gtx,"-",name))
}

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
   # 
   # Area
   # % Samples that are Compound Hets

library(gplots)

## Heatmap Amongst Gene Variables
COR.dat <- cor( matrix(as.numeric(GENE[,8:28]),ncol=21), use="pairwise.complete.obs",method="spearman" )
colnames(COR.dat) <- rownames(COR.dat) <- colnames(GENE)[8:28]
COLS.list <- c("gold1","chocolate2","firebrick2","black","slateblue3","steelblue2","springgreen1")
COLS <- colorRampPalette(COLS.list)(100)
jpeg( paste(PathToOut,"/Plots/Compile_Gene_Heat.jpeg",sep=""), height=2000,width=2000, pointsize=36)
heatmap.2( COR.dat, col=COLS, trace="none" )
dev.off()

## Plot Size of Gene vs Size of Exons
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
COLS <- c("steelblue1","slateblue1","springgreen3","gold3")
BIN <- 2
jpeg( paste(PathToOut,"/Plots/Compile_Num_Vars.jpeg",sep=""), height=1200,width=2400, pointsize=32)
par(mfrow=c(1,2))
 # Full
VAR.dat <- matrix( as.numeric(GENE[,c("n.VAR","n.SNP","n.IND")]), ncol=3 )
BRKS <- seq( 0,as.numeric(max(VAR.dat))+250, 250)
hist( VAR.dat[,3], breaks=BRKS, col=COLS[4], density=20,angle=-30, main="Distribution of Variant Count (Full)", xlab="Variant Count", ylab="# Genes")
hist( VAR.dat[,1], breaks=BRKS, col=COLS[1], density=20,angle=30, main="Distribution of Variant Count (Full)", xlab="Variant Count", ylab="# Genes", add=T)
hist( VAR.dat[,2], breaks=BRKS, col=COLS[3], density=20,angle=60, main="Distribution of Variant Count (Full)", xlab="Variant Count", ylab="# Genes", add=T)
 # Exon
VAR.dat <- matrix( as.numeric(GENE[,c("n.VAR.ex","n.SNP.ex","n.IND.ex")]), ncol=3 )
BRKS <- seq( 0,as.numeric(max(VAR.dat))+5, 5)
hist( VAR.dat[,3], breaks=BRKS, col=COLS[4], density=20,angle=-30, main="Distribution of Variant Count (Exon)", xlab="Variant Count", ylab="# Genes")
hist( VAR.dat[,1], breaks=BRKS, col=COLS[2], density=20,angle=30, main="Distribution of Variant Count (Exon)", xlab="Variant Count", ylab="# Genes", add=T)
hist( VAR.dat[,2], breaks=BRKS, col=COLS[3], density=20,angle=60, main="Distribution of Variant Count (Exon)", xlab="Variant Count", ylab="# Genes", add=T)
dev.off()

## Distribution of AREA
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

## Distribution of % Samples that are Compound Het
COLS <- c("springgreen3","gold3")
BIN <- .05
XLIM <- c( 0,1 )
BRKS <- seq(0,1,BIN)
jpeg( paste(PathToOut,"/Plots/Compile_CH_Dist.jpeg",sep=""), height=1200,width=2400, pointsize=32)
par(mfrow=c(1,2))
 # Full
CH.dat.snp <- as.numeric( GENE[,"p.COMP_HET_snp"] )
CH.dat.ind <- as.numeric( GENE[,"p.COMP_HET_ind"] )
hist( CH.dat.snp, breaks=BRKS, col=COLS[1], main="Distribution of % Compound Hets (Full)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=45 )
hist( CH.dat.ind, breaks=BRKS, col=COLS[2], main="Distribution of % Compound Hets (Full)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=-45, add=T )
 # Exon
CH.dat <- as.numeric( GENE[,"p.COMP_HET_snp.ex"] )
hist( CH.dat, breaks=BRKS, col=COLS[1], main="Distribution of % Compound Hets (Exon)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=45 )
CH.dat <- as.numeric( GENE[,"p.COMP_HET_ind.ex"] )
hist( CH.dat, breaks=BRKS, col=COLS[2], main="Distribution of % Compound Hets (Exon)", xlab="% Samples w/ Compound Het in Gene", ylab="# Genes", density=20,angle=-45, add=T )
legend( "topright", fill=COLS, density=20, legend=c("SNP","IND") )
dev.off()

###############################################################
## SAVE COMPILED DATA #########################################
###############################################################

## Save Table that Gene Data are Compiled In
write.table( GENE, paste(PathToOut,"/Gene_Stats.txt",sep=""), sep="\t",row.names=F,col.names=T,quote=F )

## Save List of FULL/EXON Cohort data
COMPILE <- list( FULL, EXON )
names(COMPILE) <- c("Full","Exon")
save( COMPILE, file=paste(PathToOut,"/Gene_Stats.Rdata",sep="") )

###############################################################
## END OF DOC #################################################
###############################################################

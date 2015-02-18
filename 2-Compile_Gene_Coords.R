## Compile/Filter Gene Coordinates ##
## Called by Gene_Score Pipeline ##
## February 18, 2015 ##
## Kristopher Standish ##


## Usage ##
# Rscript 2-Compile_Gene_Coords.R <Path/To/Out_Dir>

###############################################################
## PARSE COMMAND LINE #########################################
###############################################################


LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- "/projects/janssen/Phased/20150218_Test_Genes"
PathToOut <- LINE[1]

# Specify other Paths
PathToGeneCoords <- paste( PathToOut, "Gene_Info.raw.txt", sep="/" )

# Check for proper parsing
print(paste( "Output:", PathToOut ))
print(paste( "Gene Coords:", PathToGeneCoords ))

###############################################################
## LOAD DATA ##################################################
###############################################################

## Load Gene Coordinate File
TAB <- read.table( PathToGeneCoords, sep="\t", header=T, comment.char="",quote="" )

###############################################################
## EDIT FOR NEW TABLES ########################################
###############################################################

## Combine Gene_Transcript Columns
GTX <- paste( TAB[,15],TAB[,1],sep="_")

## Rename Chromosome Column
CHR.range <- TAB[,2]
 # Remove "chr" tag
CHR.plink <- gsub( "chr","", CHR.range )

## Buffer Gene Coordinates
 # TX_S
TX_S.range <- TAB[,3]
TX_S.plink <- TAB[,3] - 5000
 # TX_E
TX_E.range <- TAB[,4]
TX_E.plink <- TAB[,4] + 5000

###############################################################
## FILTER OUT TRANSCRIPTS #####################################
###############################################################

## Non-Autosomal
RM.chr <- which( !(CHR.plink %in% 1:22) )

## Transcripts with same coordinates
 # String everything together: GENE_CHR_TXs_TXe_EXs_EXe_nEX
STR <- paste( TAB[,15],TAB[,2],TAB[,3],TAB[,4],TAB[,7],TAB[,8],TAB[,9], sep="_" )
RM.rpt <- which(duplicated( STR ))

## Combine Removal Criteria
RM <- Reduce( union, list(RM.chr,RM.rpt) )

## Create new Tables w/ Transcripts Removed
 # New "Gene_Info" table
TAB.rm <- TAB[ -RM, ]
 # "Gene_Range" table
RNG <- data.frame( GTX=GTX, CHR=CHR.range, TX_S=TX_S.range, TX_E=TX_E.range )
RNG.rm <- RNG[ -RM, ]
PLNK <- data.frame( CHR=CHR.plink, TX_S=TX_S.plink, TX_E=TX_E.plink, GTX=GTX )
PLNK.rm <- PLNK[ -RM, ]

## Write Tables
write.table( TAB.rm, paste(PathToOut,"Gene_Info.txt",sep="/" ), sep="\t",row.names=F,col.names=T,quote=F )
write.table( RNG.rm, paste(PathToOut,"Gene_Range.txt",sep="/" ), sep="\t",row.names=F,col.names=T,quote=F )
write.table( PLNK.rm, paste(PathToOut,"Gene_Range.plink.txt",sep="/" ), sep="\t",row.names=F,col.names=F,quote=F )






###############################################################
## END OF DOC #################################################
###############################################################
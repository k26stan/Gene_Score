## Pull out only JnJ Samples from Phased Genotypes Files ##
## February 13, 2015 ##
## Kristopher Standish ##

####################################################
## R CODE ##########################################
####################################################

## Load Haplotype File and Sample Name File
DAT.1 <- read.table( gzfile("/projects/janssen/ychoi/phasing_results/chr1.phased.haps.gz"), header=F, nrow=3 )
SMP.1 <- read.table( "/projects/janssen/ychoi/phasing_results/chr1.phased.sample", header=T )
dim(DAT.1) # 5 columns prior to first sample
dim(SMP.1) # 1 null row. Get rid of it

## Specify colnames for DAT.1
SMP.1.names <- as.character( SMP.1[2:nrow(SMP.1),1] )
SMP.1.names.x2.a <- rep( SMP.1.names, rep(2,length(SMP.1.names)) )
SMP.1.names.x2 <- paste( SMP.1.names.x2.a, c(1,2), sep="_")
DAT.1.colnames <- c( "CHROM","ID","POS","REF","ALT", SMP.1.names.x2 )
colnames(DAT.1) <- DAT.1.colnames

## Get out JnJ Samples only
RM.HG.names <- grep("HG", SMP.1.names.x2 )
RM.NA.names <- grep("NA", SMP.1.names.x2 )
RM.pub <- union( RM.HG.names, RM.NA.names )
JNJ.1.names <- SMP.1.names.x2[ -RM.pub ]

WHICH_COLS <- which( colnames(DAT.1) %in% JNJ.1.names )

### BOTTOM LINE ###
 ## First 5 columns are SNP information
 ## Next 874 columns are JnJ Samples
 ## Keep first 879 columns

####################################################
## BASH ############################################
####################################################

## Pull out JnJ Samples for each File
 # Unzip | pull out columns > write to separate file
for chrom in `seq 2 6`
# for chrom in `seq 22 22`
do
echo \### Running Chromosome ${chrom} - `date`
## Sample File
HAP_PATH=/projects/janssen/ychoi/phasing_results/chr${chrom}.phased.sample
OUT_PATH=/projects/janssen/Phased/Data/chr${chrom}.phased.sample
head -1 ${HAP_PATH} > ${OUT_PATH}
head -439 ${HAP_PATH} | tail -437 >> ${OUT_PATH}
## Hap File
HAP_PATH=/projects/janssen/ychoi/phasing_results/chr${chrom}.phased.haps.gz
OUT_PATH=/projects/janssen/Phased/Data/chr${chrom}.phased.haps
echo Writing Hap File - `date`
zcat ${HAP_PATH} | cut -d ' ' -f1-879 | sed 's/ /\t/g' > ${OUT_PATH}
echo BGZipping - `date`
bgzip ${OUT_PATH}
echo Tabixing - `date`
tabix ${OUT_PATH}.gz -s 1 -b 3 -e 3
done
echo DONE - `date`




# ## Load Haplotype File and Sample Name File
# DAT.16 <- read.table( gzfile("/projects/janssen/ychoi/phasing_results/chr16.phased.haps.gz"), header=F, nrow=3 )
# SMP.16 <- read.table( "/projects/janssen/ychoi/phasing_results/chr16.phased.sample", header=T )
# dim(DAT.16) # 5 columns prior to first sample
# dim(SMP.16) # 1 null row. Get rid of it

# ## Specify colnames for DAT.16
# SMP.16.names <- as.character( SMP.16[2:nrow(SMP.16),1] )
# SMP.16.names.x2.a <- rep( SMP.16.names, rep(2,length(SMP.16.names)) )
# SMP.16.names.x2 <- paste( SMP.16.names.x2.a, c(1,2), sep="_")
# DAT.16.colnames <- c( "CHROM","ID","POS","REF","ALT", SMP.16.names.x2 )
# colnames(DAT.16) <- DAT.16.colnames


































####################################################
## END OF DOC ######################################
####################################################
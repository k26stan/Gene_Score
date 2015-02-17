## Start compiling stats for Phased Janssen Genomes ##
## October 13, 2014 ##
## Kristopher Standish ##

######################################################
## OVERVIEW/GOALS ####################################
######################################################

## Notes ##
# accommodate non-coding RNA
#   - e.g., TOP1P2 (chr22)
# figure out what's going wrong w/ LOC/MRG in "Compile_Stats.R"
#   - e.g., TRAF3IP2-AS1 (chr6)

## Steps ##
#1) Pull out single chromosome sample from single sample of phased data
#2) Filter out non-variant positions
#3) Pull annotations for remaining positions
#4) Pull out list of genes for each chromosome
#5) Go thru variants and get counts of variant by functional position and gene for each homologue

##########################################################################
## 1 ## Set up Paths #####################################################
##########################################################################
 # Use Bash
 # Take in arguments, set up directories/paths for files/tools
echo \### 1 - `date` \###
echo \### Define Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

JANS_DIR=/projects/janssen

## Names/Paths for Output
DATE=20150213
HOME_DIR=${JANS_DIR}/Phased/${DATE}_Testing

## Parameters

## Files
ANNOTS=${JANS_DIR}/ANNOTATE/JnJ_121613_all_annotations.txt.bgz
PHASED_DIR=${JANS_DIR}/ychoi/phasing_results
PHASED_FILENAME=chrWHICH_CHROM.phased.haps.gz
SAMPLE_LIST=Path/To/Sample_List.txt

###########################################################
## Constant Paths ##

## Public Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt

## Custom Scripts

###########################################################
## Pull some Info out of Parameters ##






















###########################################################
## 1 - Pull out sample/chrom ##############################
###########################################################

## Loop Thru Chromsomes
for WHICH_CHROM in `seq 1 22`
do

## Set Path To Phased Chromosome
PHASED_FILENAME_CHR=`echo ${PHASED_FILENAME} | sed 's/WHICH_CHROM/9/g'`
PHASED_PATH=${PHASED_DIR}/${PHASED_FILENAME_CHR}
 # And Sample List that goes with it
PHASED_FILENAME_CHR_SAMP=`echo ${PHASED_FILENAME_CHR} | sed 's/haps\.gz/sample/g'`
PHASED_PATH_SAMP=${PHASED_DIR}/${PHASED_FILENAME_CHR_SAMP}

## Loop Thru Samples
for SAMPLE in `cat ${SAMPLE_LSIT}`
do

## Set Sample
WHICH_SAMP=${SAMPLE}
SAMPLE_NUM=$( expr `cat ${PHASED_PATH_SAMP} | grep ${WHICH_SAMP} -n | cut -d : -f 1` - 2 )
 # Specify Which Columns of Phased Data
COL_1=$( expr 4 + ${SAMPLE_NUM} \* 2)
COL_2=$( expr ${COL_1} + 1)

## Pull out haplotypes for desired sample and chromosome
OUT_FILE=${HOME_DIR}/OUT.${WHICH_SAMP}.${WHICH_CHROM}.phased
zcat ${PHASED_PATH} | cut -d ' ' -f1-5,${COL_1},${COL_2} > ${OUT_FILE}
# zcat ${JANS_DIR}/ychoi/phasing_results/chr${WHICH_CHROM}.phased.haps.gz | cut -d ' ' -f1-5,${COL_1},${COL_2} | head -100 > ${OUT_FILE}

###########################################################
## 2 - Filter out Non-Variant sites #######################
###########################################################

awk ' !( $6==0 && $7==0 ) ' ${OUT_FILE} > ${OUT_FILE}.2

###########################################################
## 3 - Pull Annotations for Variant Positions #############
###########################################################

## Make string of genomic positions
POSITS=`cat ${OUT_FILE}.2 | awk '{print "chr"$1":"$3"-"$3}'`

## Determine which columns to print for annotations
WHICH_COLS=`zcat ${ANNOTS} | head -1 | sed "s/\t/\n/g" | grep -nrx 'Gene\|Gene_Type\|Location\|Coding_Impact' | awk -F: '{print $1}'`
CUT_COLS=`echo ${WHICH_COLS} | sed 's/ /,/g'`

## Pull out annotations for desired locations
zcat ${ANNOTS} | head -1 | cut -f2-7,${CUT_COLS} > ${OUT_FILE}.annot
tabix ${ANNOTS} ${POSITS} | cut -f2-7,${CUT_COLS} >> ${OUT_FILE}.annot

###########################################################
## 4 - Pull Gene List for Chromosome ######################
###########################################################

## Copy Gene List into "HOME_DIR"
cp ${JANS_DIR}/RV_Rare_Variant/GENES_PER_CHROM/CHR_${WHICH_CHROM}_UNIQ.list ${HOME_DIR}/GENE_LIST_CHR${WHICH_CHROM}.list

###########################################################
## 5 - Get Counts by Gene/Genome/Functional Element #######
###########################################################

## Use Rscript to compile stats
Rscript ${JANS_DIR}/Phased/20141013_Compile_Stats.R ${HOME_DIR} ${WHICH_SAMP} ${WHICH_CHROM}

done # Close WHICH_CHROM loop






###########################################################
## END OF DOC #############################################
###########################################################

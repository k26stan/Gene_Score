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

## New Idea 2/13 ##
# ) Specify which genes I want to look at
# ) Collect Summary Statistics across Cohort
  # Number of Variants per person
  # Number of exonic/intronic/etc...
  # Number of hets vs homozygous
# ) Plot Summary Statistics
# ) Plot Gene Map (maybe w/ some variants?)
  # Plot common haplotypes??

## Specifics
#1) Compile Paths & whatnot
  # Specify list of genes to look at
#2) Pull out those genes into individual variant files & annotations
  # Include 5KB buffer on each end
  # Phased file and indel raw file and snp raw file
#3) Compile Stats about Gene across cohort
  # % Phased Variants (per person [pp])
  # % Hom/Het (pp?)
  # Num Exonic & Intronic & UTR vars
  # Num Compound Hets (intronic & exonic?)
  # 

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
DATE=$1 # 20150218
OUT_DIR=$2 # ${JANS_DIR}/Phased/${DATE}_Test_Genes
mkdir ${OUT_DIR}
cd ${OUT_DIR}

## Files
GENE_LIST=$3 # ${JANS_DIR}/Phased/20150217_Test_Genes.txt
ANNOTS=$4 # ${JANS_DIR}/ANNOTATE/JnJ_121613_all_annotations.txt.bgz
PHASED_PATH=$5 # ${JANS_DIR}/Phased/Data/chrWHICH_CHROM.phased.haps.gz
SNP_PATH=$6 # ${JANS_DIR}/VCFs/PLINK/FULL_BED_SNP/BED_FULL.SNP
IND_PATH=$7 # ${JANS_DIR}/VCFs/PLINK/FULL_BED_IND/BED_FULL.IND

## Association Testing Files/Parameters
PHENO_PATH=$8 # ${JANS_DIR}/ASSOCIATION/PH-PHENOTYPES/LT8_DEL_MNe_MN.txt # Which Phenotype File are you using?
PHENO_TYPE=$9 # C # Is phenotype (B)inary or (C)ontinuous?
COV_FILE=${10} # ${JANS_DIR}/ASSOCIATION/PH-PHENOTYPES/COV.txt # Path to Covariate File or "F"
COVS=${11} # `echo DAS_BL_MN` # Which Covariates to Include?
EIG_VEC=${12} # ${JANS_DIR}/ASSOCIATION/EIGEN/HC_FULL.eigenvec # Output from Plink's --pca command (MAF>1%) of "F"
PC_COUNT=${13} # 2 # How many PCs to Include as Covariates?
START_STEP=${14} # 5

###########################################################
## Constant Paths ##

## Public Tools
GATK_JAR=/projects/janssen/Tools/gatk2.7-2/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
VCF_TOOLS=/projects/janssen/Tools/vcftools_0.1.11/bin/vcftools
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink 
GENE_TABLE=/home/kstandis/HandyStuff/GG-Gene_Names_DB.txt

## Custom Scripts
s2_COMPILE_GENE_COORDS_R=/projects/janssen/Phased/SCRIPTS/2-Compile_Gene_Coords.R
s4_MAKE_COV_TAB_R=/projects/janssen/Phased/SCRIPTS/4-Make_Cov_Tab.R
s4_MANH_PLOT_R=/projects/janssen/Phased/SCRIPTS/4-Manhat_Plot.R
s5_GENE_COMPILE_R=/projects/janssen/Phased/SCRIPTS/5-Gene_Compile.R

###########################################################
## Pull some Info out of Parameters ##

# Get Names of Specific Files
DIRS=(${PHENO_PATH//\// })
PHENO=${DIRS[${#DIRS[@]} - 1]} # Get Name of Phenotype File

# Specify list of Covariates to include (for command ad for filename)
if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

# Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi

# Specify commands and extensions for Cont vs Bin Phenotype
if [ $PHENO_TYPE = "C" ]
then
SUFFIX=linear
else
SUFFIX=logistic
fi

## Specify a File to which to Write Updates
UPDATE_FILE=${OUT_DIR}/Update.txt

## Done
if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Define Set Variables and Paths - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 2 - Pull Gene Coordinates #############################################
##########################################################################
if [ "$START_STEP" -le 2 ]; then
echo \### 2 - `date` \###
echo \### Pull Gene Coordinates \###
echo `date` "2 - Pull Gene Coordinates" >> ${UPDATE_FILE}

## Get Gene Coordinates from Ref Table
head -1 ${GENE_TABLE} > ${OUT_DIR}/Gene_Info.raw.txt
cat ${GENE_TABLE} | grep -f ${GENE_LIST} >> ${OUT_DIR}/Gene_Info.raw.txt
# cat ${OUT_DIR}/Gene_Info.txt | awk -F\\t '{print $15"_"$1"\t"$2"\t"$3"\t"$4}' > ${OUT_DIR}/Gene_Range.txt

## Compile Coordinates into one file (w/ Buffer)
Rscript ${s2_COMPILE_GENE_COORDS_R} ${OUT_DIR}

# IFSo=$IFS
# IFS=$'\n' # Makes it so each line is read whole (not separated by tabs)
# for line in `tail -n+2 ${OUT_DIR}/Gene_Range.txt`
# do
#  # Pull out Name/Coordinates for each Transcripts
# tag=`echo ${line} | awk '{print $1}'`
# chr=`echo ${line} | awk '{print $2}' | sed 's/chr//g'`
# start=$( expr `echo ${line} | awk '{print $3}'` - 5000 )
# stop=$( expr `echo ${line} | awk '{print $4}'` + 5000 )
# echo -e ${chr}'\t'${start}"\t"${stop}"\t"${tag} >> ${OUT_DIR}/Gene_Range.plink.txt
# done
# IFS=$IFSo # Reset


## Done
echo `date` "2 - Pull Gene Coordinates - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 - Pull Variants & Annotations #######################################
##########################################################################
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Pull Variants and Annotations \###
echo `date` "3 - Pull Variants and Annotations" >> ${UPDATE_FILE}

echo \### Pulling All Variants to one BED \###
mkdir ${OUT_DIR}/Genes
## Pull out Variants for All Genes first
 # (quicker when looping through later)
 # Indels
${PLINK} \
--bfile ${IND_PATH} \
--make-bed \
--extract range ${OUT_DIR}/Gene_Range.plink.txt \
--out ${OUT_DIR}/IND_Vars
 # SNPs
${PLINK} \
--bfile ${SNP_PATH} \
--make-bed \
--extract range ${OUT_DIR}/Gene_Range.plink.txt \
--out ${OUT_DIR}/SNP_Vars

## Cycle Through Gene_Transcripts
echo \### Looping Through Gene_Transcripts \###
IFSo=$IFS
IFS=$'\n'
for line in `tail -n+2 ${OUT_DIR}/Gene_Range.txt`
do
 # Pull out Name/Coordinates for each Transcripts
tag=`echo ${line} | awk '{print $1}'`
chr=`echo ${line} | awk '{print $2}' | sed 's/chr//g'`
start=$( expr `echo ${line} | awk '{print $3}'` - 5000 )
stop=$( expr `echo ${line} | awk '{print $4}'` + 5000 )
 # Make Directory for each Gene_Transcript
OUT_PATH_GENE=${OUT_DIR}/Genes/${tag}
mkdir ${OUT_PATH_GENE}
echo ${line} > ${OUT_PATH_GENE}/Gene_Range.txt
echo -e ${tag}"\t"${chr}':'${start}"-"${stop} > ${OUT_PATH_GENE}/Coords_w_Buffer.txt
 # Pull out Phased Variants
PHASED_PATH_CHR=`echo ${PHASED_PATH} | sed "s/WHICH_CHROM/${chr}/g"`
tabix ${PHASED_PATH_CHR} ${chr}:${start}-${stop} > ${OUT_PATH_GENE}/Phased.haps
cp ${PHASED_PATH_CHR%%.haps.gz}.sample ${OUT_DIR}/Phased.sample
 # Pull out INDs to Raw File
${PLINK} \
--bfile ${OUT_DIR}/IND_Vars \
--recode A \
--silent \
--hardy midp \
--chr ${chr} \
--from-bp ${start} \
--to-bp ${stop} \
--out ${OUT_PATH_GENE}/IND_Vars
 # Pull out SNPs to Raw File
${PLINK} \
--bfile ${OUT_DIR}/SNP_Vars \
--recode A \
--silent \
--hardy midp \
--chr ${chr} \
--from-bp ${start} \
--to-bp ${stop} \
--out ${OUT_PATH_GENE}/SNP_Vars
 # Pull Gene Annotations from CYPHER file
if [ -e ${ANNOTS} ]
then
 # Determine which columns to print for annotations
WHICH_COLS=`zcat ${ANNOTS} | head -1 | sed "s/\t/\n/g" | grep -nrx 'Gene\|Gene_Type\|Location\|Coding_Impact' | awk -F: '{print $1}'`
CUT_COLS=`echo ${WHICH_COLS} | sed 's/ /,/g'`
 # Pull out annotations for desired locations
zcat ${ANNOTS} | head -1 > ${OUT_PATH_GENE}/Annots.txt
tabix ${ANNOTS} chr${chr}:${start}-${stop} >> ${OUT_PATH_GENE}/Annots.txt
cat ${OUT_PATH_GENE}/Annots.txt | cut -f2-7,${CUT_COLS} >> ${OUT_PATH_GENE}/Annots_Short.txt
fi

done # Close Gene_Transcript Loop
IFS=$IFSo # Reset

## Done
echo `date` "3 - Pull Variants and Annotations - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 4 - Single-Locus Association ##########################################
##########################################################################
if [ "$START_STEP" -le 4 ]; then
echo \### 4 - `date` \###
echo \### Single-Locus Association \###
echo `date` "4 - Single-Locus Association" >> ${UPDATE_FILE}

## Make directory for Plots
mkdir ${OUT_DIR}/Plots

## If a Phenotype is provided, Run Single Locus Testing
if [ -e ${PHENO_PATH} ]
then

## Make Directory for Association Files
ASSOC_DIR=${OUT_DIR}/Assoc
mkdir ${ASSOC_DIR}

echo \### Dealing with PC and Covariate Files \###
## If No Principal Components Exist, Make Them
if [ $EIG_VEC = "F" ]
then
 # Use BED to run PCA
${PLINK} --bfile ${SNP_PATH} \
--pca header \
--silent \
--allow-no-sex \
--out ${ASSOC_DIR}/PCs
EIG_VEC=${ASSOC_DIR}/PCs.eigenvec
fi

## Make new Covariate File
if [ $COV_FILE = "F" ] ; then
cp ${EIG_VEC} ${NEW_COV_FILE}
else
Rscript ${s4_MAKE_COV_TAB_R} ${EIG_VEC} ${COV_FILE} ${ASSOC_DIR}/Cov_w_PCs.txt
fi

echo \### Running Single-Locus Association \###
## Run Single-Locus Association on All Variants in Region
 # SNPs
${PLINK} \
--bfile ${OUT_DIR}/SNP_Vars \
--pheno ${PHENO_PATH} \
--covar ${ASSOC_DIR}/Cov_w_PCs.txt \
--covar-name ${COVS_COMMAND} \
--${SUFFIX} hide-covar \
--adjust qq-plot \
--allow-no-sex \
--maf 0.01 \
--keep-allele-order \
--freqx \
--hardy midp \
--silent \
--out ${ASSOC_DIR}/SNP_Assoc
 # Indels
${PLINK} \
--bfile ${OUT_DIR}/IND_Vars \
--pheno ${PHENO_PATH} \
--covar ${ASSOC_DIR}/Cov_w_PCs.txt \
--covar-name ${COVS_COMMAND} \
--${SUFFIX} hide-covar \
--adjust qq-plot \
--allow-no-sex \
--maf 0.01 \
--keep-allele-order \
--freqx \
--hardy midp \
--silent \
--out ${ASSOC_DIR}/IND_Assoc

echo \### Making Manhattan Plot \###
## Manhattan Plot of Single-Locus for these Genes
 # SNPs
mkdir ${ASSOC_DIR}/SNP/
cat ${ASSOC_DIR}/SNP_Assoc.assoc.${SUFFIX} | awk '{print $1"\t"$2"\t"$3"\t"$9}' > ${ASSOC_DIR}/SNP_Assoc.P
mv ${ASSOC_DIR}/SNP_Assoc* ${ASSOC_DIR}/SNP/
Rscript ${s4_MANH_PLOT_R} ${ASSOC_DIR}/SNP/SNP_Assoc.P ${PHENO} ${COVS_FILENAME}
 # Indels
mkdir ${ASSOC_DIR}/IND/
cat ${ASSOC_DIR}/IND_Assoc.assoc.${SUFFIX} | awk '{print $1"\t"$2"\t"$3"\t"$9}' > ${ASSOC_DIR}/IND_Assoc.P
mv ${ASSOC_DIR}/IND_Assoc* ${ASSOC_DIR}/IND/
Rscript ${s4_MANH_PLOT_R} ${ASSOC_DIR}/IND/IND_Assoc.P ${PHENO} ${COVS_FILENAME}

## Move Plots to deisgnated directory
for si in `echo SNP IND`
do
cd ${ASSOC_DIR}/${si}/
for file in `ls *jpeg`
do
mv ${file} ${OUT_DIR}/Plots/${si}_${file}
done
done
cd ${OUT_DIR}

fi # Closes if PHENO exists

## Done
echo `date` "4 - Single-Locus Association - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 5 - Compile Gene Stats & Plot #########################################
##########################################################################
if [ "$START_STEP" -le 5 ]; then
echo \### 5 - `date` \###
echo \### Compile Gene Stats and Plot \###
echo `date` "5 - Compile Gene Stats and Plot" >> ${UPDATE_FILE}

echo \### Doing Final Gene Stat Compilation and Plotting \###
## Rscript to Cycle Through Genes and Compile Stats & Make some Plots
Rscript ${s5_GENE_COMPILE_R} ${OUT_DIR} ${PHENO_PATH} ${COVS_COMMAND}


### Look at unique Haplotypes for each Gene




## Done
echo `date` "5 - Compile Gene Stats and Plot - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 6 - Clean Up Directory ################################################
##########################################################################
if [ "$START_STEP" -le 6 ]; then
echo \### 6 - `date` \###
echo \### Clean Up Directory \###
echo `date` "6 - Clean Up Directory" >> ${UPDATE_FILE}

### Look into making list of candidate Genes & keeping files for those genes
  # but deleting a bunch from the rest
    # To delete: Annots.txt, Phased.haps (?), *_Vars.raw (?)
    # Check how much space would be saved just by deleting the full Annots.txt









## Done
echo `date` "6 - Clean Up Directory - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
echo \### DONE \###






###########################################################
## END OF DOC #############################################
###########################################################
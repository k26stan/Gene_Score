# Make Covariate File w/ EigenVectors #

# Make_Cov_Tab.R <Path/To/EigenVectors> <Path/To/Covariate_Table> <Set>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c(""/)
PathToVec <- LINE[1]
PathToCov <- LINE[2]
PathToNewCov <- LINE[3]

## Load Tables
VEC <- read.table(PathToVec,header=T)
COV <- read.table(PathToCov,sep="\t",header=T)

## Merge Tables
VEC <- VEC[,c(1,3:ncol(VEC))]
MRG <- merge(x=COV,y=VEC,by="FID")

## Write Table
write.table(MRG,PathToNewCov,sep="\t",row.names=F,col.names=T,quote=F)

######### END OF DOC ###########

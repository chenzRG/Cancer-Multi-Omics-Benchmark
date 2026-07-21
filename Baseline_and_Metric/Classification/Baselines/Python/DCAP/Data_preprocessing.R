# MLOmics: accept command-line arguments so this script is not hardcoded to Windows paths
# Usage: Rscript Data_preprocessing.R <work_dir> <input_csv> <output_csv>
# Example: Rscript Data_preprocessing.R /data/brca brca_multitest.csv brcatest_go.csv
args <- commandArgs(trailingOnly = TRUE)
work_dir   <- if (length(args) >= 1) args[1] else 'D:/brca'           # MLOmics: default kept for backward compat
input_csv  <- if (length(args) >= 2) args[2] else 'brca_multitest.csv'  # MLOmics
output_csv <- if (length(args) >= 3) args[3] else 'brcatest_go.csv'     # MLOmics

setwd(work_dir)  # MLOmics: was hardcoded 'D:/brca'
library('imputeMissings')
data=read.csv(input_csv, row.names = 1)  # MLOmics: was hardcoded "brca_multitest.csv"

#################delete missing#######################
sum(is.na(data))
miss=c()
for (i in 1:nrow(data)){
  miss=c(miss,sum(is.na(data[i,])))
}
missrate=miss/ncol(data)
data1=data[missrate<0.2,]
data1[is.na(data1)] <- 0
###################delete zero#################
nz=c()
for (i in 1:nrow(data1)){
  nz=c(nz,sum(data1[i,]==0))
}
zerorate=nz/ncol(data1)
data2=data1[zerorate<0.2,]
###############impute######################

data3=t(data2)
data3=data.frame(data3)
data3[data3==0]=NA
data4<-impute(data3)
##########normalize##############
data5=t(data4)
data6=data5
data7=t(data6)
data8=scale(data7, center = T, scale = T)
data9=t(data8)
write.csv(data9, output_csv)  # MLOmics: was hardcoded 'brcatest_go.csv'


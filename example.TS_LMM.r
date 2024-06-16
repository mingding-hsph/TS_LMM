setwd("/Users/mingding/upload to github")

##download and read in TS_LMM macro
source('./TS_LMM.macro.r')

##read in dataset
new_data<-read.csv("./data.csv", row.names = NULL)

##read in corr_snps, which is the correlation between snps
##ld can be obtained using LDmatrix Tool in LDlink. Please refer to our paper how to obtain correlation between snps based on LD

corr_snps_readin<-read.csv("./corr_snps.csv", row.names = NULL)
corr_snps<- corr_snps_readin[,-1]  ##exclude the first column which is a character variable

##check and make sure that data are numeric and has 20 rows and 20 columns
##where 20 is the number of genetic variants
summary(corr_snps)  
nrow(corr_snps) 
ncol(corr_snps) 

##calculate correlation between summary statistics of risk factors 

corr_x<-cor(new_data[, c('bx1', 'bx2', 'bx3')])
  
##apply TS_LMM
  
  model_output<-TS_LMM(
    betaX=cbind(new_data$bx1, new_data$bx2, new_data$bx3),
    betaY=new_data$by, 
    betaX_se=cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
    betaY_se=new_data$by_se,
    corr_snps_tslmm =corr_snps,
    corr_X_tslmm =corr_x,
    loop_rem=500,
    cutoff_rem=0.00001 
  )
  
  model_output
  

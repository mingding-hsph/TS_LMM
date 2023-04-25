setwd("/Users/mingding/upload to github")

##download and read in TS_LMM macro
source('./TS_LMM.macro.r')

##read in dataset
new_data<-read.csv("./data.csv", row.names = NULL)

##read in corr_snps, which is the correlation between snps
##can be obtained using LDmatrix Tool in LDlink

ld_snps_readin<-read.csv("./corr_snps.csv", row.names = NULL)
ld_snps<- ld_snps_readin[,-1]  ##exclude the first column which is a character variable

##check and make sure that data are numeric and has 20 rows and 20 columns
##where 20 is the number of genetic variants
summary(ld_snps)  
nrow(ld_snps) 
ncol(ld_snps) 

##calculate correlation between summary statistics of risk factors 

ld_x<-cor(new_data[, c('bx1', 'bx2', 'bx3')])
  
##apply TS_LMM
  
  model_output<-TS_LMM(
    betaX=cbind(new_data$bx1, new_data$bx2, new_data$bx3),
    betaY=new_data$by, 
    betaX_se=cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
    betaY_se=new_data$by_se,
    ld_snps_tslmm =ld_snps,
    ld_X_tslmm =ld_x,
    loop_rem=50,
    cutoff_rem=0.00001 
  )
  
  model_output
  

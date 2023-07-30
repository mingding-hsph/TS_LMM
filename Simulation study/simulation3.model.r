setwd("/Users")

##read in our macro
source('/Users/TS_LMM.macro.r')

##read in ld_snps

data_snps<-read.csv("./simulated_snps_LD_low.csv", row.names = NULL)

ld_snps<-round(cor(data_snps, method = "pearson"),digits=3)

##start simulation

loop<-1000

m<-3

MRMV_model<-matrix(NA,nrow=loop, ncol=5*m) ##m is the number of exposures
MRMV_IVW_model<-matrix(NA,nrow=loop, ncol=5*m) ##m is the number of exposures
MRMV_MEDIAN_model<-matrix(NA,nrow=loop, ncol=5*m) ##m is the number of exposures
MMR_2stage_model<-matrix(NA,nrow=loop, ncol=5*m) ##m is the number of exposures
presso_model<-matrix(NA,nrow=loop, ncol=5*m) ##m is the number of exposures

IV_str<-matrix(NA,nrow=loop, ncol=m) ##m is the number of exposures

E_cor<-matrix(NA,nrow=loop, ncol=m*(m-1)/2)



colnames(MRMV_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                            "pvalue_x1", "pvalue_x2","pvalue_x3",
                            "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(MRMV_IVW_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                "pvalue_x1", "pvalue_x2","pvalue_x3",
                                "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(MRMV_MEDIAN_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                   "pvalue_x1", "pvalue_x2","pvalue_x3",
                                   "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(MMR_2stage_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                                  "pvalue_x1", "pvalue_x2","pvalue_x3",
                                  "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")
colnames(presso_model)<-rbind("mean_x1","mean_x2","mean_x3","se_x1","se_x2","se_x3",
                              "pvalue_x1", "pvalue_x2","pvalue_x3",
                              "ll_x1","ll_x2","ll_x3", "ul_x1",  "ul_x2", "ul_x3")

colnames(IV_str)<-rbind("IV_x1","IV_x2","IV_x3")

colnames(E_cor)<-rbind("x1_x2","x1_x3","x2_x3")

for (simu_loop in 1:loop)
  
{
  
source('/Users/simulation2.data.r')

  IV_str[simu_loop,1]<-table(new_data_f_di[,1])[2]
  IV_str[simu_loop,2]<-table(new_data_f_di[,2])[2]
  IV_str[simu_loop,3]<-table(new_data_f_di[,3])[2]
  
##install.packages("MendelianRandomization")

library(MendelianRandomization)

MRMV<-mr_mvegger(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                           bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                           by = new_data$by, byse = new_data$by_se, correlation=ld_snps))

MRMV_IVW<-mr_mvivw(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                    bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                    by = new_data$by, byse = new_data$by_se, correlation=ld_snps))

MRMV_MEDIAN<-mr_mvmedian(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                       bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                       by = new_data$by, byse = new_data$by_se, correlation=ld_snps), 
                       iterations = 1000)

##output

MRMV_model[simu_loop,1:3]<-MRMV$Estimate
MRMV_model[simu_loop,(1:3)+3]<-MRMV$StdError.Est
MRMV_model[simu_loop,(1:3)+2*3]<-MRMV$Pvalue.Est
MRMV_model[simu_loop,(1:3)+3*3]<-MRMV$CILower.Est
MRMV_model[simu_loop,(1:3)+4*3]<-MRMV$CIUpper.Est

MRMV_IVW_model[simu_loop,1:3]<-MRMV_IVW$Estimate
MRMV_IVW_model[simu_loop,(1:3)+3]<-MRMV_IVW$StdError
MRMV_IVW_model[simu_loop,(1:3)+2*3]<-MRMV_IVW$Pvalue
MRMV_IVW_model[simu_loop,(1:3)+3*3]<-MRMV_IVW$CILower
MRMV_IVW_model[simu_loop,(1:3)+4*3]<-MRMV_IVW$CIUpper

MRMV_MEDIAN_model[simu_loop,1:3]<-MRMV_MEDIAN$Estimate
MRMV_MEDIAN_model[simu_loop,(1:3)+3]<-MRMV_MEDIAN$StdError
MRMV_MEDIAN_model[simu_loop,(1:3)+2*3]<-MRMV_MEDIAN$Pvalue
MRMV_MEDIAN_model[simu_loop,(1:3)+3*3]<-MRMV_MEDIAN$CILower
MRMV_MEDIAN_model[simu_loop,(1:3)+4*3]<-MRMV_MEDIAN$CIUpper

###MMR-2stage

##read in correlation of exposures
ld_x<-cor(new_data[, c('bx1', 'bx2', 'bx3')])

E_cor[simu_loop,1]<-ld_x[1,2]
E_cor[simu_loop,2]<-ld_x[1,3]
E_cor[simu_loop,3]<-ld_x[2,3]

MMR<-TS_LMM(
  betaX=cbind(new_data$bx1, new_data$bx2, new_data$bx3),
  betaY=new_data$by, 
  betaX_se=cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
  betaY_se=new_data$by_se,
  ld_snps_tslmm =ld_snps,
  ld_X_tslmm =ld_x,
  loop_rem=50,
  cutoff_rem=0.00001 
)



MMR_2stage_model[simu_loop,1:3]<-t(MMR[,1])
MMR_2stage_model[simu_loop,(1:3)+3]<-t(MMR[,2])
MMR_2stage_model[simu_loop,(1:3)+2*3]<-t(MMR[,3])
MMR_2stage_model[simu_loop,(1:3)+3*3]<-t(MMR[,4])
MMR_2stage_model[simu_loop,(1:3)+4*3]<-t(MMR[,5])



##install mrpresso
##if (!require("devtools")) { install.packages("devtools") } else {}
##devtools::install_github("rondolab/MR-PRESSO")

library(MRPRESSO)

# Run MR-PRESSO on a multi-variable MR (MMR) model specifying several exposures

new_data_press<-data.frame(new_data)

MRPRESSO=mr_presso(BetaOutcome = "by", BetaExposure = c("bx1", "bx2", "bx3"), 
                   SdOutcome = "by_se", SdExposure = c("bx_se1", "bx_se2", "bx_se3"), 
                   OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = new_data_press, 
                   NbDistribution = 1000,  SignifThreshold = 0.05)

##checked, the warnings() are fine

presso_model[simu_loop,1:3]<-MRPRESSO$`Main MR results`[1:3,3]  ##estimate
presso_model[simu_loop,(1:3)+3]<-MRPRESSO$`Main MR results`[1:3,4]  ##although write sd, this shoudl be se
presso_model[simu_loop,(1:3)+2*3]<-MRPRESSO$`Main MR results`[1:3,6]  ##p value
presso_model[simu_loop,(1:3)+3*3]<-MRPRESSO$`Main MR results`[1:3,3]-1.96*MRPRESSO$`Main MR results`[1:3,4]
presso_model[simu_loop,(1:3)+4*3]<-MRPRESSO$`Main MR results`[1:3,3]+1.96*MRPRESSO$`Main MR results`[1:3,4]

write.csv(MRMV_model,"./table1_MRMV.csv")
write.csv(MRMV_IVW_model,"./table1_MRMV_IVW.csv")
write.csv(MRMV_MEDIAN_model,"./table1_MRMV_MEDIAN.csv")
write.csv(MMR_2stage_model,"./table1_TS_LMM.csv")
write.csv(presso_model,"./table1_presso.csv")
write.csv(IV_str,"./table1_IV.csv")
write.csv(E_cor,"./table1_E_cor.csv")

}






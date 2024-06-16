
##read in our macro
source('./TS_LMM.macro.r')


loop<-1000

qhet_model<-matrix(NA,nrow=loop, ncol=5*3) 
MRMV_model<-matrix(NA,nrow=loop, ncol=5*3) 
MRMV_IVW_model<-matrix(NA,nrow=loop, ncol=5*3) 
MRMV_MEDIAN_model<-matrix(NA,nrow=loop, ncol=5*3) 
MMR_2stage_model<-matrix(NA,nrow=loop, ncol=5*3) 
presso_model<-matrix(NA,nrow=loop, ncol=5*3) 
IV_strength<-matrix(NA,nrow=loop, ncol=3) 
IV_con_strength<-matrix(NA,nrow=loop, ncol=3) ##this include conditional IV
incept_pvalue<-matrix(NA,nrow=loop, ncol=1)

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

colnames(incept_pvalue)<-rbind("Pvalue")
                        
for (simu_loop in 1:loop)
  
{
 
source('./simulation.v2/simulation1.data.UBP.r')

  ##estimate IV strength
  
  ##estimate mean strength of IV
  
  library(tidyverse)
  
  F1=mean((new_data$bx1)^2/(new_data$bx_se1)^2)
  
  F2=mean((new_data$bx2)^2/(new_data$bx_se2)^2)
  
  F3=mean((new_data$bx3)^2/(new_data$bx_se3)^2)
  
  IV_strength[simu_loop,1]<-F1
  IV_strength[simu_loop,2]<-F2
  IV_strength[simu_loop,3]<-F3
  
  ##check pleiotropy
  
  incept_pvalue[simu_loop,1]=summary(lm(by~bx1+bx2+bx3,data=new_data))$coeff[1,4]  ##p value of intercept
  
  
##install.packages("MendelianRandomization")

library(MendelianRandomization)

MRMV<-mr_mvegger(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                           bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                           by = new_data$by, byse = new_data$by_se, correlation=corr_snps))

MRMV_IVW<-mr_mvivw(mr_mvinput(bx = cbind(new_data$bx1, new_data$bx2, new_data$bx3), 
                    bxse = cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
                    by = new_data$by, byse = new_data$by_se, correlation=corr_snps))


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


###MMR-2stage

##estimate conditional IV strength

library(MVMR)

con_F <- format_mvmr(BXGs =t(rbind(new_data$bx1,new_data$bx2,new_data$bx3)),
                     BYG = new_data$by,
                     seBXGs =t(rbind(new_data$bx_se1,new_data$bx_se2,new_data$bx_se3)),
                     seBYG = new_data$by_se,
                     RSID = seq(1:nrow(new_data)))

mvmrcovmatrix<-corr_x

Xcovmat<-phenocov_mvmr(mvmrcovmatrix, con_F[,7:9])

con_F_s <- strength_mvmr(r_input = con_F, gencov = Xcovmat)

IV_con_strength[simu_loop,1]=as.numeric(con_F_s[1])
IV_con_strength[simu_loop,2]=as.numeric(con_F_s[2])
IV_con_strength[simu_loop,3]=as.numeric(con_F_s[3])


##TS_LMM

MMR<-TS_LMM(
  betaX=cbind(new_data$bx1, new_data$bx2, new_data$bx3),
  betaY=new_data$by, 
  betaX_se=cbind(new_data$bx_se1, new_data$bx_se2, new_data$bx_se3),
  betaY_se=new_data$by_se,
  corr_snps_tslmm =corr_snps,
  corr_X_tslmm =corr_x,
  loop_rem=500,
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

write.csv(MRMV_model,"./simulation.v2/results/UBP_table1_MRMV.csv")
write.csv(MRMV_IVW_model,"./simulation.v2/results/UBP_table1_MRMV_IVW.csv")
write.csv(MMR_2stage_model,"./simulation.v2/results/UBP_table1_TS_LMM.csv")

write.csv(IV_con_strength,"./simulation.v2/results/UBP_table1_IV_con.csv")
write.csv(IV_strength,"./simulation.v2/results/UBP_table1_IV.csv")
write.csv(incept_pvalue,"./simulation.v2/results/UBP_incept_pvalue.csv")

}


###qhet model

##qhet<-qhet_mvmr(con_F, mvmrcovmatrix, CI = T, iterations = 100)
##qhet_model[simu_loop,1:3]<-t(matrix(qhet[,1]))   ##estimate
##qhet_model[simu_loop,(1:3)+3]<-t(matrix(qhet[,2]))  ##se
#write.csv(qhet_model,"./table1_qhet.csv")



setwd("/Users")

MRMV_model=read.csv("./uni_table1_MRMV.csv", row.names = NULL)
MRMV_IVW_model=read.csv("./uni_table1_MRMV_IVW.csv", row.names = NULL)
##MRMV_MEDIAN_model=read.csv("./uni_table1_MRMV_MEDIAN.csv", row.names = NULL)
MMR_2stage_model=read.csv("./uni_table1_TS_LMM.csv", row.names = NULL)
presso_model=read.csv("./uni_table1_presso.csv", row.names = NULL)

##output the main table

scenario=cbind(lp[1])  ##change risk factor here

##change the numbers here

MRMV_model$lp1=scenario[1,1]

MRMV_model$coverage1=NA
MRMV_model$power1=NA
MRMV_model$FPR1=NA
MRMV_model$mean1=NA

for (i in 1:nrow(MRMV_model))
{  
  
  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
  if (MRMV_model$lp1[i]>=MRMV_model$ll_x1[i] & MRMV_model$lp1[i]<=MRMV_model$ul_x1[i])
  { MRMV_model$coverage1[i]=1
  } else {
    MRMV_model$coverage1[i]=0
  }
  
  
  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
  if (MRMV_model$pvalue_x1[i]<0.05 & MRMV_model$lp1[i]==0) 
  { MRMV_model$FPR1[i]=1
  } else if (MRMV_model$pvalue_x1[i]>=0.05 & MRMV_model$lp1[i]==0) {
    MRMV_model$FPR1[i]=0  
  }
  
  
  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
  if (MRMV_model$pvalue_x1[i]<0.05 & MRMV_model$lp1[i]!=0) 
  { MRMV_model$power1[i]=1
  } else if (MRMV_model$pvalue_x1[i]>=0.05 & MRMV_model$lp1[i]!=0)
  { MRMV_model$power1[i]=0
  }

  
}

mean(MRMV_model$coverage1)
mean(MRMV_model$power1)
mean(MRMV_model$FPR1)


##ivw

##change the numbers here

MRMV_IVW_model$lp1=scenario[1,1]


MRMV_IVW_model$coverage1=NA
MRMV_IVW_model$power1=NA
MRMV_IVW_model$FPR1=NA
MRMV_IVW_model$mean1=NA


for (i in 1:nrow(MRMV_IVW_model))
{  
  
  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
  if (MRMV_IVW_model$lp1[i]>=MRMV_IVW_model$ll_x1[i] & MRMV_IVW_model$lp1[i]<=MRMV_IVW_model$ul_x1[i])
  { MRMV_IVW_model$coverage1[i]=1
  } else {
    MRMV_IVW_model$coverage1[i]=0
  }
  
  
  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
  if (MRMV_IVW_model$pvalue_x1[i]<0.05 & MRMV_IVW_model$lp1[i]==0) 
  { MRMV_IVW_model$FPR1[i]=1
  } else if (MRMV_IVW_model$pvalue_x1[i]>=0.05 & MRMV_IVW_model$lp1[i]==0) {
    MRMV_IVW_model$FPR1[i]=0  
  }

  
  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
  if (MRMV_IVW_model$pvalue_x1[i]<0.05 & MRMV_IVW_model$lp1[i]!=0) 
  { MRMV_IVW_model$power1[i]=1
  } else if (MRMV_IVW_model$pvalue_x1[i]>=0.05 & MRMV_IVW_model$lp1[i]!=0)
  { MRMV_IVW_model$power1[i]=0
  }
  
  
}

mean(MRMV_IVW_model$coverage1)
mean(MRMV_IVW_model$power1)
mean(MRMV_IVW_model$FPR1)

##median

#MRMV_MEDIAN_model$lp1=scenario[1,1]

#MRMV_MEDIAN_model$coverage1=NA
#MRMV_MEDIAN_model$power1=NA
#MRMV_MEDIAN_model$FPR1=NA
#MRMV_MEDIAN_model$mean1=NA


#for (i in 1:nrow(MRMV_MEDIAN_model))
#{  
  
  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
#  if (MRMV_MEDIAN_model$lp1[i]>=MRMV_MEDIAN_model$ll_x1[i] & MRMV_MEDIAN_model$lp1[i]<=MRMV_MEDIAN_model$ul_x1[i])
#  { MRMV_MEDIAN_model$coverage1[i]=1
#  } else {
#    MRMV_MEDIAN_model$coverage1[i]=0
    #  }


  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
    #  if (MRMV_MEDIAN_model$pvalue_x1[i]<0.05 & MRMV_MEDIAN_model$lp1[i]==0) 
    #  { MRMV_MEDIAN_model$FPR1[i]=1
    #  } else if (MRMV_MEDIAN_model$pvalue_x1[i]>=0.05 & MRMV_MEDIAN_model$lp1[i]==0) {
    #    MRMV_MEDIAN_model$FPR1[i]=0  
    #  }
  
  
  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
    #  if (MRMV_MEDIAN_model$pvalue_x1[i]<0.05 & MRMV_MEDIAN_model$lp1[i]!=0) 
    #  { MRMV_MEDIAN_model$power1[i]=1
    # } else if (MRMV_MEDIAN_model$pvalue_x1[i]>=0.05 & MRMV_MEDIAN_model$lp1[i]!=0)
    #  { MRMV_MEDIAN_model$power1[i]=0
    #  }
  
  
    #}

    #mean(MRMV_MEDIAN_model$coverage1)
    #mean(MRMV_MEDIAN_model$power1)
    #mean(MRMV_MEDIAN_model$FPR1)

###our model

MMR_2stage_model$lp1=scenario[1,1]

MMR_2stage_model$coverage1=NA
MMR_2stage_model$power1=NA
MMR_2stage_model$FPR1=NA


for (i in 1:nrow(MMR_2stage_model))
{  
  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
  if (MMR_2stage_model$lp1[i]>=MMR_2stage_model$ll_x1[i] & MMR_2stage_model$lp1[i]<=MMR_2stage_model$ul_x1[i])
  { MMR_2stage_model$coverage1[i]=1
  } else {
    MMR_2stage_model$coverage1[i]=0
  }
  
  
  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
  if (MMR_2stage_model$pvalue_x1[i]<0.05 & MMR_2stage_model$lp1[i]==0) 
  { MMR_2stage_model$FPR1[i]=1
  } else if (MMR_2stage_model$pvalue_x1[i]>=0.05 & MMR_2stage_model$lp1[i]==0) {
    MMR_2stage_model$FPR1[i]=0  
  }
  
  
  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
  if (MMR_2stage_model$pvalue_x1[i]<0.05 & MMR_2stage_model$lp1[i]!=0) 
  { MMR_2stage_model$power1[i]=1
  } else if (MMR_2stage_model$pvalue_x1[i]>=0.05 & MMR_2stage_model$lp1[i]!=0)
  { MMR_2stage_model$power1[i]=0
  }
  
  
}

mean(MMR_2stage_model$coverage1)
mean(MMR_2stage_model$power1)
mean(MMR_2stage_model$FPR1)

##presso model


presso_model$lp1=scenario[1,1]

presso_model$coverage1=NA
presso_model$power1=NA
presso_model$FPR1=NA
presso_model$mean1=NA

for (i in 1:nrow(presso_model))
{  
  
  ##Coverage: proportion of 95% CIs including the true causal effect over all simulated datasets; 
  
  if (presso_model$lp1[i]>=presso_model$ll_x1[i] & presso_model$lp1[i]<=presso_model$ul_x1[i])
  { presso_model$coverage1[i]=1
  } else {
    presso_model$coverage1[i]=0
  }
  

  ##false-positive rate: the proportion of significant causal estimates in settings with no causal effect. 
  
  if (presso_model$pvalue_x1[i]<0.05 & presso_model$lp1[i]==0) 
  { presso_model$FPR1[i]=1
  } else if (presso_model$pvalue_x1[i]>=0.05 & presso_model$lp1[i]==0) {
    presso_model$FPR1[i]=0  
  }
  

  
  ##power is estimated by the proportion of 95% CIs excluding zero in settings with significant causal effect; 
  
  if (presso_model$pvalue_x1[i]<0.05 & presso_model$lp1[i]!=0) 
  { presso_model$power1[i]=1
  } else if (presso_model$pvalue_x1[i]>=0.05 & presso_model$lp1[i]!=0)
  { presso_model$power1[i]=0
  }
  

  
}


##output into tables

##install.packages("plotrix")                           # Install plotrix R package
library("plotrix")                                    # to calculate standard error


MRMV_v1_mean=mean(MRMV_model$mean_x1)
MRMV_v1_std=mean(MRMV_model$se_x1)
MRMV_v1_coverage=mean(MRMV_model$coverage1)
MRMV_v1_power=mean(MRMV_model$power1)
MRMV_v1_fpr=mean(MRMV_model$FPR1)



MRMV_IVW_v1_mean=mean(MRMV_IVW_model$mean_x1)
MRMV_IVW_v1_std=mean(MRMV_IVW_model$se_x1)
MRMV_IVW_v1_coverage=mean(MRMV_IVW_model$coverage1)
MRMV_IVW_v1_power=mean(MRMV_IVW_model$power1)
MRMV_IVW_v1_fpr=mean(MRMV_IVW_model$FPR1)


#MRMV_MEDIAN_v1_mean=mean(MRMV_MEDIAN_model$mean_x1)
#MRMV_MEDIAN_v1_std=mean(MRMV_MEDIAN_model$se_x1)
#MRMV_MEDIAN_v1_coverage=mean(MRMV_MEDIAN_model$coverage1)
#MRMV_MEDIAN_v1_power=mean(MRMV_MEDIAN_model$power1)
#MRMV_MEDIAN_v1_fpr=mean(MRMV_MEDIAN_model$FPR1)


MMR_v1_mean=mean(MMR_2stage_model$mean_x1)
MMR_v1_std=mean(MMR_2stage_model$se_x1)
MMR_v1_coverage=mean(MMR_2stage_model$coverage1)
MMR_v1_power=mean(MMR_2stage_model$power1)
MMR_v1_fpr=mean(MMR_2stage_model$FPR1)

presso_v1_mean=mean(presso_model$mean_x1)
presso_v1_std=mean(presso_model$se_x1)
presso_v1_coverage=mean(presso_model$coverage1)
presso_v1_power=mean(presso_model$power1)
presso_v1_fpr=mean(presso_model$FPR1)


v1_mean=round(rbind(MMR_v1_mean, MRMV_v1_mean, MRMV_IVW_v1_mean, presso_v1_mean),digits=2)
v1_std=round(rbind(MMR_v1_std,MRMV_v1_std, MRMV_IVW_v1_std,  presso_v1_std),digits=2)
v1_coverage=100*rbind(MMR_v1_coverage,MRMV_v1_coverage, MRMV_IVW_v1_coverage, presso_v1_coverage)
v1_power=100*rbind(MMR_v1_power,MRMV_v1_power, MRMV_IVW_v1_power,presso_v1_power)
v1_FPR=100*rbind(MMR_v1_fpr,MRMV_v1_fpr, MRMV_IVW_v1_fpr, presso_v1_fpr)



combined_table=round(cbind(rbind(v1_mean,v1_std,v1_coverage,v1_power,v1_FPR)),digits=2)

write.csv(combined_table,"./uni_table1_output.table.csv")

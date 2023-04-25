setwd("/Users")

##change the parameters in the following lines for a variety of scenarios

##strength of IV: strong: strength_IV=0; weak: strength_IV=1

##correlation between summary statistics of risk factors:
##strong correlation:   lmk[i,j]<-lmk_half1*lmk_half2+rnorm(1, mean = 0, sd = 0.1)  
##weak correlation:   lmk[i,j]<-lmk_half1*lmk_half2+rnorm(1, mean = 0, sd = 1)  


##read in simulated genetic variants

data_snps<-read.csv("./simulated_snps_LD_low.csv", row.names = NULL)
colnames(data_snps)
summary(data_snps)
data<-data_snps

##general information of simulated data

n=1000      ##number of participants
k=20         ##number of genetic variants
m=3          ##number of exposures
strength_IV=1 ##1  

##simulate coefficients

##lmk, a k*(m+1) matrix, represents coefficients of k genetic variants with m exposures and Y
##u, a n*(m+1) matrix, also create correlations among m exposures and Y
##lp, a m*1 matrix, represents coefficients of m exposures with outcome Y

##simulate lmk

lmk<-matrix(nrow=m+1, ncol=k)  

for (j in 1:k)
  {
  
lmk_half1<- sample(cbind(-1,1),1)*runif(1, min = -1, max = 1)  ##sample(matrix(cbind(-1,1)),1)*
lmk_half2<- runif(1, min = -1, max = 1)  ##sample(matrix(cbind(-1,1)),1)*

    for (i in 1:(m+1))
    {
      lmk[i,j]<-lmk_half1*lmk_half2+rnorm(1, mean = 0, sd = 1)  
      ##strong correlation 0.1, moderate, 0.5, weak 1
  }
}

##simulate u for exposures and outcome

library(mvtnorm)

u_mean=c(1,1,1,1) 

cov_u <- matrix(nrow = m+1, ncol = m+1)

for (i in 1:(m+1))
{
  for (j in 1:(m+1))
  {
    if (i==j) {
      cov_u[i,j]=1  
    } 
    else if (i<j){
      cov_u[i,j]=runif(1,min=0.5, max=0.7)   ##runif(1,min=0.5, max=0.7) 
    }
    else if (i>j){
      cov_u[i,j]=cov_u[j,i]
    }  
  }
}

u<-rmvnorm(n, mean=u_mean, sigma=cov_u)  ##sigma is the covariance for each variable


##generate exposure a

for (i in 1:m)
{
  for (j in 1:k)
  {
    data[,paste0('a_sub_snp', i,j)]<-(data[,paste0('snp', j)]+rnorm(n, mean = 0, sd = strength_IV))*lmk[i,j]
  }
}

for (i in 1:m)
{
data[,paste0('a',i)]<-
  rowSums(data[,c(paste0("a_sub_snp", i, 1:k))])  ##sum up all snps*lmk for exposure i
+ u[,i]   
+ rnorm(n, mean = 0, sd = 1) 
}

##simulate lp, p between exposure m and outcome

lp<-matrix(nrow=m)  

for (i in 1:m)
{
  lp[1]=1
  lp[2]=0
  lp[3]=-1
  
  data[,paste0('y_sub_a', i)]<-data[,paste0('a', i)]*lp[i]
}

##generate outcome y

  for (j in 1:k)
  {
    data[,paste0('y_sub_snp', j)]<-lmk[m+1,j]*(data[,paste0('snp', j)]+rnorm(n, mean = 0, sd = strength_IV))
  }


data$y<-rowSums(data[,c(paste0("y_sub_a", 1:m))])+rowSums(data[,c(paste0("y_sub_snp", 1:k))])
+ u[,m]  
+ rnorm(n, mean = 0, sd = 1)

## obtain summary statistics

new_data_b<-data.frame(matrix(nrow = k,ncol = m+1))
new_data_se<-data.frame(matrix(nrow = k,ncol = m+1))
new_data_pvalue<-data.frame(matrix(nrow = k,ncol = m+1))
new_data_r2<-data.frame(matrix(nrow = k,ncol = m+1))
new_data_f<-data.frame(matrix(nrow = k,ncol = m+1))
new_data_f_di<-data.frame(matrix(nrow = k,ncol = m+1))

for (i in 1:m)
{
colnames(new_data_b)[i]<-paste0('bx', i)
colnames(new_data_se)[i]<-paste0('bx_se', i)
colnames(new_data_pvalue)[i]<-paste0('bx_p', i)
colnames(new_data_r2)[i]<-paste0('bx_r', i)
colnames(new_data_f)[i]<-paste0('bx_f', i)
colnames(new_data_f_di)[i]<-paste0('bx_f_di', i)

}

colnames(new_data_b)[m+1]<-paste0('by')
colnames(new_data_se)[m+1]<-paste0('by_se')
colnames(new_data_pvalue)[m+1]<-paste0('by_p')
colnames(new_data_r2)[m+1]<-paste0('by_r')
colnames(new_data_f)[m+1]<-paste0('by_f')
colnames(new_data_f_di)[m+1]<-paste0('by_f_di')

for (j in 1:k)
{
for (i in 1:m)
{
new_data_b[j,i]<-summary(lm(data[,c(paste0('a', i))]~data[,c(paste0('snp', j))]))$coeff[2,1]

new_data_se[j,i]<-summary(lm(data[,c(paste0('a', i))]~data[,c(paste0('snp', j))]))$coeff[2,2]

new_data_pvalue[j,i]<-summary(lm(data[,c(paste0('a', i))]~data[,c(paste0('snp', j))]))$coeff[2,4]

new_data_r2[j,i]<-summary(lm(data[,c(paste0('a', i))]~data[,c(paste0('snp', j))]))$adj.r.squared

new_data_f[j,i]<-summary(lm(data[,c(paste0('a', i))]~data[,c(paste0('snp', j))]))$fstatistic[1]

if (new_data_f[j,i]<10) {new_data_f_di[j,i]=0} else {
new_data_f_di[j,i]=1
}

}
}



for (j in 1:k)
{
  ##obtain beta coefficient
  new_data_b[j,m+1]<-summary(lm(data$y~data[,c(paste0('snp', j))]))$coeff[2,1]
  
  ##obtain se of beta coefficient
  new_data_se[j,m+1]<-summary(lm(data$y~data[,c(paste0('snp', j))]))$coeff[2,2]

  new_data_pvalue[j,m+1]<-summary(lm(data$y~data[,c(paste0('snp', j))]))$coeff[2,4]
  
  new_data_r2[j,m+1]<-summary(lm(data$y~data[,c(paste0('snp', j))]))$adj.r.squared
  
  new_data_f[j,m+1]<-summary(lm(data$y~data[,c(paste0('snp', j))]))$fstatistic[1]
  
  if (new_data_f[j,m+1]<10) {new_data_f_di[j,m+1]=0} else {
    new_data_f_di[j,m+1]=1
  }
  
  }

new_data<-cbind(new_data_b,new_data_se,new_data_pvalue,new_data_r2,new_data_f)

##output correlation between summary statistics of genetic variants and x

cor(new_data_b[, c('bx1', 'bx2', 'bx3')])



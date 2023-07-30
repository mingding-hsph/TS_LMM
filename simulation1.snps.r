setwd("/Users")

#In case the simulation of genetic variants is too slow, 
##another way is to simulate 20 genetic variants in 2 linkage disequilibrium (LD) blocks 
##with each block including 10 genetic variants. 

##change parameters in line LD[i,j]=runif(1,min=0.1, max=0.5) to simulate low and moderate lD
##simulate low LD: LD[i,j]=runif(1,min=0, max=0.1)
##simulate  moderate LD: LD[i,j]=runif(1,min=0.1, max=0.5)

n=1000    ##number of participants

data_snps<-matrix(NA, nrow=n, ncol=20)

k=20     ##number of genetic variants within each LD block

##simulate correlations between snps 

LD <- matrix(nrow = k, ncol = k)

for (i in 1:k)
{
  for (j in 1:k)
  {
    if (i==j) {
      LD[i,j]=1  
    } 
    else if (i<j){
      LD[i,j]=runif(1,min=0.1, max=0.5)  ##low to moderate ld 
    }
    else if (i>j){
      LD[i,j]=LD[j,i]
    }  
  }
}


# marginal probabilities to represent minor allele frequency of genetic variants

p_snps <- matrix(nrow = k, ncol = 1)

for (i in 1:k)
{
  p_snps[i]=runif(1,min=0.1, max=0.5)  ##p is between 0.1-0.5 to exclude rare alleles
}


# estimate the joint-distribution

library(mipfp)
p.joint <- ObtainMultBinaryDist(corr = LD, marg.probs = p_snps)

# simulate n draws from the obtained joint-distribution


data_snps[1:n,1:k] <- 
  RMultBinary(n = n, mult.bin.dist = p.joint)$binary.sequences


##generate data

colnames(data_snps)=paste0('snp', 1:20)

write.csv(data_snps,"./simulated_snps_LD_new_aa.csv", row.names = F)

##check data

summary(data_snps)  ##check if mean values are in the range of 0.1-0.5

check_cor<-matrix(NA, ncol=20, nrow=20)

for (m in 1:20){
  for (n in 1:20) {
check_cor[m,n]<-cor(data_snps[,m], data_snps[,n])
##check correlation between snps
}
}

summary(check_cor)

######################################################################
               ##                MR-CAUSE           ##
######################################################################
# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("/projects/MRC-IEU/research/projects/icep1/wp2/040/working/data/CAUSE") 

library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(tidyverse)
library(dplyr)
library(readr)
library(ieugwasr)
library(cause)
library(genetics.binaRies)

#---------------------------------------------------------------------#
#                            Read exposure                             #----
#---------------------------------------------------------------------#
#read exposure
AFS <- read.table("afs.tsv", header = T, sep="\t")
NSP <- read.table("NUMBER_SEXUAL_PARTNERS_GWAS.txt?dl=0", header = T, sep="\t")
#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#

AllSites <- read.table("data.AllSites.txt", header = T)
OPC <- read.table("data.Oropharynx.txt", header = T)  
  
#---------------------------------------------------------------------#
#                            AFS -> AllSites                        #----
#---------------------------------------------------------------------#
  
#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(AFS, AllSites, 
                snp_name_cols = c("variant_id", "RSID"), 
                beta_hat_cols = c("beta", "ESTIMATE"), 
                se_cols = c("standard_error", "SE"), 
                A1_cols = c("effect_allele", "EFFECT_ALL"), 
                A2_cols = c("other_allele", "OTHER_ALL"))

#Formatting X1
#There are  9851866  variants.
#Removing  14738  duplicated variants leaving  9822416 variants.
#No variants have illegal alleles.
#Removed  1481612  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  8340804  variants.
#Formatting X2
#There are  7528473  variants.
#Removing  0  duplicated variants leaving  7528473 variants.
#No variants have illegal alleles.
#Removed  1138139  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  6390334  variants.
#After merging and removing variants with inconsistent alleles,  there are  6309519  variants that are present in both studies and can be used with CAUSE.

#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#
X$p_value <- 2*pnorm(abs(X$beta_hat_1/X$seb1), lower.tail=FALSE)
X_clump <- X %>% rename(rsid = snp,
                        pval = p_value) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.001,
                     clump_p = 5e-08,
                     # plink_bin = genetics.binaRies::get_plink_binary(),
                     #bfile = "~/EUR"
  )
keep_snps <- as.character(X_clump$rsid) 
#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#
res <- cause(X=X, variants = keep_snps, param_ests = params)
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('CAUSE_AFS_AllSites1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('CAUSE_AFS_AllSites2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()


#---------------------------------------------------------------------#
#                            NSP -> AllSites                        #----
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(NSP, AllSites, 
                snp_name_cols = c("MarkerName", "RSID"), 
                beta_hat_cols = c("Beta", "ESTIMATE"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("A1", "EFFECT_ALL"), 
                A2_cols = c("A2", "OTHER_ALL"))

#Formatting X1
#There are  11515109  variants.
#Removing  97101  duplicated variants leaving  11418007 variants.
#No variants have illegal alleles.
#Removed  1707624  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  9710383  variants.
#Formatting X2
#There are  7528473  variants.
#Removing  0  duplicated variants leaving  7528473 variants.
#No variants have illegal alleles.
#Removed  1138139  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  6390334  variants.
#After merging and removing variants with inconsistent alleles,  there are  6307835  variants that are present in both studies and can be used with CAUSE
#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#
X$p_value <- 2*pnorm(abs(X$beta_hat_1/X$seb1), lower.tail=FALSE)
X_clump <- X %>% rename(rsid = snp,
                        pval = p_value) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.001,
                     clump_p = 5e-08,
                     # plink_bin = genetics.binaRies::get_plink_binary(),
                     #bfile = "~/EUR"
  )
keep_snps <- as.character(X_clump$rsid) 
#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#
res <- cause(X=X, variants = keep_snps, param_ests = params)
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('CAUSE_NSP_AllSites1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('CAUSE_NSP_AllSites2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()


#---------------------------------------------------------------------#
#                            AFS -> OPC                       #----
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(AFS, OPC, 
                snp_name_cols = c("variant_id", "CHR_POS_EFF_OTH"), 
                beta_hat_cols = c("beta", "ESTIMATE"), 
                se_cols = c("standard_error", "SE"), 
                A1_cols = c("effect_allele", "EFFECT_ALL"), 
                A2_cols = c("other_allele", "OTHER_ALL"))
#There are  9851866  variants.
#Removing  14738  duplicated variants leaving  9822416 variants.
#No variants have illegal alleles.
#Removed  1481612  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  8340804  variants.
#Formatting X2
#There are  7562369  variants.
#Removing  38129  duplicated variants leaving  7486115 variants.
#Removing  1  variants with illegal allelse leaving  7486115 variants.
#Removed  1119573  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  6366542  variants.
#After merging and removing variants with inconsistent alleles,  there are  6291355  variants that are present in both studies and can be used with CAUSE.

#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#
X$p_value <- 2*pnorm(abs(X$beta_hat_1/X$seb1), lower.tail=FALSE)
X_clump <- X %>% rename(rsid = snp,
                        pval = p_value) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.001,
                     clump_p = 5e-08,
                     # plink_bin = genetics.binaRies::get_plink_binary(),
                     #bfile = "~/EUR"
  )
keep_snps <- as.character(X_clump$rsid) 
#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#
res <- cause(X=X, variants = keep_snps, param_ests = params)
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('CAUSE_AFS_OPC1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('CAUSE_AFS_OPC2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()

#---------------------------------------------------------------------#
#                            NSP -> OPC                       #----
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#                            Merge GWAS data                          #----
#---------------------------------------------------------------------#

X <- gwas_merge(NSP, OPC, 
                snp_name_cols = c("MarkerName", "CHR_POS_EFF_OTH"), 
                beta_hat_cols = c("Beta", "ESTIMATE"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("A1", "EFFECT_ALL"), 
                A2_cols = c("A2", "OTHER_ALL"))

#There are  11515109  variants.
#Removing  97101  duplicated variants leaving  11418007 variants.
#No variants have illegal alleles.
#Removed  1707624  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  9710383  variants.
#Formatting X2
#There are  7562369  variants.
#Removing  38129  duplicated variants leaving  7486115 variants.
#Removing  1  variants with illegal allelse leaving  7486115 variants.
#Removed  1119573  variants with ambiguous strand.
#Flipping strand and effect allele so A1 is always A
#Returning  6366542  variants.
#After merging and removing variants with inconsistent alleles,  there are  6289967  variants that are present in both studies and can be used with CAUSE.

#---------------------------------------------------------------------#
#                    Calculate nuisance parameters                    #----
#---------------------------------------------------------------------#
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
head(params$mix_grid)

#---------------------------------------------------------------------#
#                                Clump data                           #----
#---------------------------------------------------------------------#
X$p_value <- 2*pnorm(abs(X$beta_hat_1/X$seb1), lower.tail=FALSE)
X_clump <- X %>% rename(rsid = snp,
                        pval = p_value) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.001,
                     clump_p = 5e-08,
                     # plink_bin = genetics.binaRies::get_plink_binary(),
                     #bfile = "~/EUR"
  )
keep_snps <- as.character(X_clump$rsid) 
#---------------------------------------------------------------------#
#                    MR-CAUSE analysis                                #----
#---------------------------------------------------------------------#
res <- cause(X=X, variants = keep_snps, param_ests = params)
plot(res$sharing)
plot(res$causal)
summary(res, ci_size=0.95)
plot(res)
plot(res, type="data")

png('CAUSE_NSP_OPC1.png', res=300, height=2000, width=3500)
plot(res)
dev.off()

png('CAUSE_NSP_OPC2.png', res=300, height=2000, width=3500)
plot(res, type="data")
dev.off()

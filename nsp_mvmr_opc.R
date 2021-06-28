
#library(digest, lib.loc = "C:/Program Files/R/R-3.6.1/library")
#install.packages("devtools")
library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#devtools::install_github("MRCIEU/MRInstruments")
library(MRInstruments)

library(MendelianRandomization)

factors=c("csi", 
          "dpw", 
          "rt", 
          "si")
df_bmis=as.matrix(cbind(factors))
for (row in 1:nrow(df_bmis)) {
  print(row)
  factor=df_bmis[row, "factors"]
  print(factor)
  setwd("/filepath")  

  #Read in top SNPs in exposure 1 
  exposure1 <- read.table(paste(factor, "_snp_list.txt", sep=""), header=T)
  exposure1$exposure=c(factor)
  #Read in / extract summary data for top SNPs in exposure 2 
  exposure2 <- read.table("nsp_snp_list.txt", header=T)
  exposure2$exposure=c("nsp")
  #Exposures
  exposure1SNPs=as.character(exposure1$SNP)
  exposure2SNPs=as.character(exposure2$SNP)
  
  exposures <- rbind(exposure1, exposure2)
  exposure_SNPs <- data.frame(c(exposure1SNPs, exposure2SNPs))
  colnames(exposure_SNPs)=c("SNP")
  check <- data.frame(table(exposure_SNPs$SNP))
  check <- data.frame(check[check$Freq > 1,])

  #Extract exposure 1 SNPs from exposure 1 
  exposure1 <- read_exposure_data(
    filename=paste(factor, "_allsnps.txt", sep=""),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    eaf_col = "eaf"
  )
  exposure1 <- exposure1[exposure1$SNP %in% exposure1SNPs,]
  exposure1$exposure=c(factor)
  n_occur <- data.frame(table(exposure1$SNP))
  n_occur <- data.frame(n_occur[n_occur$Freq > 1,])
  exposure1 = exposure1[!(exposure1$SNP %in% n_occur$Var1),]
  
  #Extract exposure 1 SNPs from exposure 2 
  outcome2 <- read_outcome_data(
    filename="nsp_allsnps.txt",
    snps = exposure1$SNP, 
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    eaf_col = "eaf"
  )
  outcome2$outcome=c("nsp")
  n_occur <- data.frame(table(outcome2$SNP))
  n_occur <- data.frame(n_occur[n_occur$Freq > 1,])
  outcome2 = outcome2[!(outcome2$SNP %in% n_occur$Var1),]
  

  #Harmonize datasets 
  dat1 <- harmonise_data(exposure1, outcome2)
  dat1<-dat1[!(dat1$remove=="TRUE"),]
  
  
  #Extract exposure 2 SNPs from exposure 2
  exposure2 <- read_exposure_data(
    filename=paste(factor, "_allsnps.txt", sep=""),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    eaf_col = "eaf"
  )
  exposure2 <- exposure2[exposure2$SNP %in% exposure2SNPs,]
  exposure2$exposure=c(factor)
  n_occur <- data.frame(table(exposure2$SNP))
  n_occur <- data.frame(n_occur[n_occur$Freq > 1,])
  exposure2 = exposure2[!(exposure2$SNP %in% n_occur$Var1),]
  
  #Extract exposure 2 SNPs from exposure 1
  outcome1 <- read_outcome_data(
    filename=paste(factor, "_allsnps.txt", sep=""),
    snps = exposure2$SNP, 
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    eaf_col = "eaf"
  )
  outcome1$outcome=c(factor)
  n_occur <- data.frame(table(outcome1$SNP))
  n_occur <- data.frame(n_occur[n_occur$Freq > 1,])
  outcome1 = outcome1[!(outcome1$SNP %in% n_occur$Var1),]
  #Harmonize datasets 
  dat2 <- harmonise_data(exposure2, outcome1)
  dat2<-dat2[!(dat2$remove=="TRUE"),]
  
  #Extract exposure 1 SNPs from exposure 1 
  outcome3 <- read_outcome_data(
    filename= paste(factor, "_allsnps.txt", sep=""),
    snps = exposure1$SNP,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    eaf_col = "eaf"
  )
  outcome3$outcome=c(factor)
  n_occur <- data.frame(table(outcome3$SNP))
  n_occur <- data.frame(n_occur[n_occur$Freq > 1,])
  outcome3 = outcome3[!(outcome3$SNP %in% n_occur$Var1),]
  #Harmonize datasets 
  dat3 <- harmonise_data(exposure1, outcome3)
  dat3<-dat3[!(dat3$remove=="TRUE"),]
  
  #Extract exposure 2 SNPs from exposure 2
  outcome4 <- read_outcome_data(
    filename="nsp_allsnps.txt",
    snps = exposure2$SNP, 
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    eaf_col = "eaf"
  )
  outcome4$outcome=c("nsp")
  n_occur <- data.frame(table(outcome4$SNP))
  n_occur <- data.frame(n_occur[n_occur$Freq > 1,])
  outcome4 = outcome4[!(outcome4$SNP %in% n_occur$Var1),]
  #Harmonize datasets 
  dat4 <- harmonise_data(exposure2, outcome4)
  dat4<-dat4[!(dat4$remove=="TRUE"),]
  
  #Append datasets 
  dat <- rbind(dat1, dat2, dat3, dat4)
  dat = dat[!(dat$SNP %in% check$Var1),]
  dat$rsid = dat$SNP
  dat <- ieugwasr::ld_clump(dat)
  dat$rsid = NULL
  dat=dat[,c(1,4,5,7,9, 13:19,28:30)]
  #double check ordering of column names 
  names(dat)=c("SNP", "effect_allele.exposure", "other_allele.exposure",
               "beta.exposure",  "eaf.exposure","id.exposure","se.exposure",
               "pval.exposure" ,"exposure","mr_keep.exposure","pval_origin.exposure",
               "data_source.exposure", "samplesize.exposure","pval", "id")
  write.csv(dat, paste(factor, "_nsp_exposures.csv", sep=""))
  exposure_SNPs =unique(dat$SNP)
  exposures=exposures[exposures$SNP %in% exposure_SNPs, ]
  #Extract SNPs for both exposures from outcome dataset 
  outcome_dat <- read_outcome_data("opc_allsnps.txt", 
                                   sep = "\t",
                                   snps = exposure_SNPs, 
                                   snp_col="SNP", 
                                   beta_col="beta", 
                                   se_col="se", 
                                   eaf_col="eaf", 
                                   effect_allele_col="effect_allele", 
                                   other_allele_col="other_allele")
  #if (is.null(outcome_dat)) {next}
  #Harmonize datasets 
  #dat <- dat[dat$SNP %in% outcome_dat$SNP,]
  dat <- harmonise_data(dat, outcome_dat, action=2)
  dat = dat[dat$ambiguous=="FALSE",]
  write.csv(dat, paste(factor, "_nsp_exposures_outcome_harmonised_opc.csv", sep=""))

  #use MVMR package
  library(devtools)
  #install.packages("remotes")
  library(remotes)
  #install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = T)
  library(MVMR)
  
  XGs <- read.csv(paste(factor, "_nsp_exposures.csv", sep=""))
  YG <- read.csv(paste(factor, "_nsp_exposures_outcome_harmonised_opc.csv", sep=""))
  
  
  XGs <- XGs[XGs$SNP %in% exposure_SNPs,]
  YG <- YG[YG$SNP %in% exposure_SNPs,]
  
  XGs = XGs[,c(2,5,8,10)]
  colnames(XGs) = c("SNP", "xg", "xgse", "Exposure")
  
  
  YG = YG[,c(2,8,15)]
  colnames(YG) = c("SNP", "yg", "ygse")
  YG=unique(YG)
  
  library(tidyr)
  XGs_betas <- XGs[,c(1:2,4)]
  XGs_betas <- spread(XGs_betas, Exposure, xg)
  XGs_se <- XGs[,c(1,3:4)]
  XGs_se <- spread(XGs_se, Exposure, xgse)
  
  #Remove NAs 
  XGs_betas <- na.omit(XGs_betas)
  XGs_se <- na.omit(XGs_se)
  
  YG <- YG[YG$SNP %in% XGs_betas$SNP,]
  XGs_betas <- XGs_betas[XGs_betas$SNP %in% YG$SNP,]
  XGs_se <- XGs_se[XGs_se$SNP %in% YG$SNP,]
  
  XGs_betas <- XGs_betas[order(XGs_betas$SNP),]
  XGs_se <- XGs_se[order(XGs_se$SNP),]
  YG <- YG[order(YG$SNP),]
  
  mvmr <- format_mvmr(XGs_betas[,c(2:3)], YG$yg, XGs_se[,c(2:3)], YG$ygse, XGs_betas$SNP)
 
  print(names(XGs_betas))  
  
  #MVMR-IVW
  mr_mvivw <- mr_mvivw(mr_mvinput(bx = as.matrix(cbind(XGs_betas[2], XGs_betas[3])), bxse = as.matrix(cbind(XGs_se[2], XGs_se[3])), by = YG$yg, byse =YG$ygse))
  mr_mvegger <- mr_mvegger(mr_mvinput(bx = as.matrix(cbind(XGs_betas[2], XGs_betas[3])), bxse = as.matrix(cbind(XGs_se[2], XGs_se[3])), by = YG$yg, byse =YG$ygse), orientate = 1)
  strength_mvmr <- strength_mvmr(mvmr)
  pleiotropy_mvmr <- pleiotropy_mvmr(mvmr)
  
  #exposure1
  #exposure1
  df <- data.frame(matrix(ncol = 14, nrow = 0))
  colnames(df)=c("b", "se", "lci", "uci", "pval", 'or', 'or_lci95', 'or_uci95', "nsnps", "exposure", "model", "F-stat", "Qstat", "Qpval")
  dfEstimate <-as.matrix(mr_mvivw@Estimate)
  df[1,1] = dfEstimate[1,1]
  dfSE = as.matrix(mr_mvivw@StdError)
  df[1,2] = dfSE[1,1]
  dfCILower <-as.matrix(mr_mvivw@CILower)
  df[1,3] = dfCILower[1,1]
  dfCIUpper <-as.matrix(mr_mvivw@CIUpper)
  df[1,4] = dfCIUpper[1,1]
  dfpval <-as.matrix(mr_mvivw@Pvalue)
  df[1,5] = dfpval[1,1]
  df[1,6] = exp(df[1,1])
  df[1,7] = exp(df[1,3])
  df[1,8] = exp(df[1,4])
  df[1,9] = as.matrix(mr_mvivw@SNPs)
  df[1,10] = names(XGs_betas[2])
  df[1,11] = "IVW" 
  df[1,12] = strength_mvmr$exposure1
  df[1,13] = pleiotropy_mvmr$Qstat
  df[1,14] = pleiotropy_mvmr$Qpval 
  df[2,1] = dfEstimate[2,1]
  df[2,2] = dfSE[2,1]
  df[2,3] = dfCILower[2,1]
  df[2,4] = dfCIUpper[2,1]
  df[2,5] = dfpval[2,1]
  df[2,6] = exp(df[2,1])
  df[2,7] = exp(df[2,3])
  df[2,8] = exp(df[2,4])
  df[2,9] = as.matrix(mr_mvivw@SNPs)
  df[2,10] = names(XGs_betas[3])
  df[2,11] = "IVW"
  df[2,12] = strength_mvmr$exposure2
  df[2,13] = pleiotropy_mvmr$Qstat
  df[2,14] = pleiotropy_mvmr$Qpval 
  dfEstimate <-as.matrix(mr_mvegger@Estimate)
  df[3,1] = dfEstimate[1,1]
  dfSE = as.matrix(mr_mvegger@StdError.Est)
  df[3,2] = dfSE[1,1]
  dfCILower <-as.matrix(mr_mvegger@CILower.Est)
  df[3,3] = dfCILower[1,1]
  dfCIUpper <-as.matrix(mr_mvegger@CIUpper.Est)
  df[3,4] = dfCIUpper[1,1]
  dfpval <-as.matrix(mr_mvegger@Pvalue.Est)
  df[3,5] = dfpval[1,1]
  df[3,6] = exp(df[3,1])
  df[3,7] = exp(df[3,3])
  df[3,8] = exp(df[3,4])
  df[3,9] = as.matrix(mr_mvivw@SNPs)
  df[3,10] = names(XGs_betas[2])
  df[3,11] = "Egger" 
  df[4,1] = dfEstimate[2,1]
  df[4,2] = dfSE[2,1]
  df[4,3] = dfCILower[2,1]
  df[4,4] = dfCIUpper[2,1]
  df[4,5] = dfpval[2,1]
  df[4,6] = exp(df[4,1])
  df[4,7] = exp(df[4,3])
  df[4,8] = exp(df[4,4])
  df[4,9] = as.matrix(mr_mvivw@SNPs)
  df[4,10] = names(XGs_betas[3])
  df[4,11] = "Egger" 
  dfIntercept <-as.matrix(mr_mvegger@Intercept)
  df[5,1] = dfEstimate[1,1]
  dfSE = as.matrix(mr_mvegger@StdError.Int)
  df[5,2] = dfSE[1,1]
  dfCILower <-as.matrix(mr_mvegger@CILower.Int)
  df[5,3] = dfCILower[1,1]
  dfCIUpper <-as.matrix(mr_mvegger@CIUpper.Int)
  df[5,4] = dfCIUpper[1,1]
  dfpval <-as.matrix(mr_mvegger@Pvalue.Int)
  df[5,5] = dfpval[1,1]
  df[5,6] = exp(df[5,1])
  df[5,7] = exp(df[5,3])
  df[5,8] = exp(df[5,4])
  df[5,9] = as.matrix(mr_mvivw@SNPs)
  df[5,10] = "intercept" 
  df[5,11] = "Egger"
  write.table(df, file=paste("nsp_", factor, "_opc.txt", sep=""), quote=F, row.names=F, sep="\t")
  
}  
setwd("")

install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")

install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")
install.packages("devtools")
devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
devtools::install_github('MRCIEU/TwoSampleMR')
install.packages("MendelianRandomization", force = TRUE)
install.packages("LDlinkR")
install.packages("plyr")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("simex")
devtools::install_github("rondolab/MR-PRESSO")

library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(MendelianRandomization)
library(gridExtra)
library(grid)
library(lattice)
library(LDlinkR)
library(ggpubr)
library(simex)
library(MRPRESSO)
ao <- available_outcomes()

#Oropharyngeal cancer analysis
#Outcome data extracted from Lesseur et al. GWAS available in IEU OpenGWAS project and dbGaP.

exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- read_outcome_data("nsp_opc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oropharyngeal cancer"
dat$exposure <- "Number of Sexual Partners"
write.csv(dat, "nsp_opc_harmonised.csv")    
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Number of sexual partners","Oropharyngeal cancer",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nsp_opc_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "nsp_opc_heterogeneity.txt")
write.table(Plt, "nsp_opc_pleiotropy.txt")
write.table(Sin, "nsp_opc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
#0 outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#2 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "nsp_opc_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nsp_opc_nooutliers.txt")

#MR presso 
outcome_dat <- read_outcome_data("hnc_snps.txt", sep="\t", snp_col="rs", beta_col = "est", se_col="SE", effect_allele_col="eff_all",other_allele_col = "other_all", pval_col = "P")
dat <- harmonise_data(exposure_dat, outcome_dat)

devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "nsp_opc_mr_presso.txt")

#Calculate I-squared, r-squared and F-statistics
#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(exposure_dat, outcome_dat, action = 1)
str(dat)

#F-statistic
dat$samplesize.exposure <- 370711
dat$samplesize.outcome <- 2641 
dat <- steiger_filtering(dat) 

N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure) 
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         
# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F = BXG^2/seBetaXG^2
mF = mean(F)
mF
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="regression_dilution_isq_weighted_nsp.csv", row.names = FALSE)

#8. SIMEX correction
#Create empty dataframe to store output
simexegger<-c()

#Run SIMEX
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG

BYG <- BetaYG*sign(BetaXG)#Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

#MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

#MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

#Simulation extrapolation 
mod.sim1 <- simex(Fit1,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
#Print results
mod1
mod2
#Use BXG result and exponentiate the estimates
OR = exp(1.288376) 
CIL = exp(1.288376  - 1.96 * 1.112615) 
CIU = exp(1.288376  + 1.96 * 1.112615) 


#Regional analysis of Oropharyngeal cancer
exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-b-96'))
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oropharyngeal cancer (European)"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
write.table(or_results, "opc_nsp_europe_results.txt")

exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-b-97'))
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oropharyngeal cancer (North American)"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
write.table(or_results, "opc_nsp_northamerican_results.txt")

exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-b-98'))
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oropharyngeal cancer (South American)"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
write.table(or_results, "opc_nsp_southamerican_results.txt")


#Positive control analysis

#Cervical cancer
#Outcome data from a GWAS generated using UKBB data 
exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- read_outcome_data("nsp_cervical_ukbb_lk.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Cervical cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Number of sexual partners","Cervical cancer",sep=""))
write.table(or_results, "nsp_cervical_cancer_results_lk.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "nsp_cervical_heterogeneity.txt")
write.table(Plt, "nsp_cervical_pleiotropy.txt")
write.table(Sin, "nsp_cervical_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 

#0 outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#4 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#5 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "nsp_cervical_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nsp_cervical_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "nsp_cervical_mr_presso.txt")


#Chlamydia infection
#Outcome data extracted from a GWAS using UKBB data
exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- read_outcome_data("nsp_chlamydia.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Chlamydia infection"
dat$exposure <- "Number of Sexual Partners"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Number of sexual partners","Chlamydia infection",sep=""))
write.table(or_results, "nsp_chlamydia_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "nsp_chlamydia_heterogeneity.txt")
write.table(Plt, "nsp_chlamydia_pleiotropy.txt")
write.table(Sin, "nsp_chlamydia_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
#0 outliers at bonf 

ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#11 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#11 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "nsp_chlamydia_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nsp_chlamydia_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "nsp_chlamydia_mr_presso.txt") 


#Negative control analysis

#Lung cancer 
#Outcome data extracted from Wang et al. GWAS (2014) available in IEU OpenGWAS project
exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-a-966'))
mu <- 11348/27209
outcome_dat$beta.outcome <- outcome_dat$beta.outcome/(mu*(1-mu))
outcome_dat$se.outcome <- outcome_dat$se.outcome/(mu*(1-mu))
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Lung cancer"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Number of sexual partners","Lung cancer",sep=""))
write.table(or_results, "nsp_lung_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "nsp_lung_heterogeneity.txt")
write.table(Plt, "nsp_lung_pleiotropy.txt")
write.table(Sin, "nsp_lung_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
#0 outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#2 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "nsp_lung_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nsp_lung_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "nsp_lung_mr_presso.txt")

#Oral cancer
#Outcome data extracted from Lesseur et al. GWAS available in IEU OpenGWAS project

exposure_dat <- read_exposure_data("nsp_exposure.txt", sep="\t")
outcome_dat <- read_outcome_data("nsp_oc.txt", sep="\t")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "Oral cancer"
dat$exposure <- "Number of Sexual Partners"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Number of sexual partners","Oral cancer",sep=""))
write.table(or_results, "nsp_oc_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
write.table(Het, "nsp_oc_heterogeneity.txt")
write.table(Plt, "nsp_oc_pleiotropy.txt")
write.table(Sin, "nsp_oc_single_snp.txt")
p1 <- mr_scatter_plot(mr_results, dat)
p1[[1]]

res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

p4 <- mr_funnel_plot(res_single)
p4[[1]]

#Assess outliers
#Radial plots 
devtools::install_github("WSpiller/RadialMR")
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 
#0 outliers at bonf 

ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#11 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
#11 outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "nsp_oc_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "nsp_oc_nooutliers.txt")

#MR presso 
devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
mr_presso
write.table(mr_presso, "nsp_oc_mr_presso.txt") 

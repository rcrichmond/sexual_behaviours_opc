#setwd("")

library(devtools)
library(TwoSampleMR)
library(MRInstruments)

#Outcome data from HPV GWAS conducted in UK Biobank

#HPV16 L1 seroprevalence (Definition 1)
exposure_dat <- read_exposure_data("./exposure_AFS.csv", sep=",")
outcome_dat <- read.csv("./AFS_HPVL1_seroprev_out.csv", header=TRUE)
outcome_dat <- read_outcome_data("AFS_HPVL1_seroprev_out.csv", exposure_dat$SNP, sep=",")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "HPV16 L1 Seroprevalence"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Age at first sex","HPV16 L1 Seroprevalence",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPVL1_seroprev_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
all_res <- combine_all_mrresults(mr_results,Het,Plt,Sin,ao_slc=F,Exp=T)
write.table(all_res, "AFS_HPVL1_seroprev_allres.txt")
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

# outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
# outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "AFS_HPVL1_seroprev_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPVL1_seroprev_nooutliers.txt")

#MR presso 
outcome_dat <- read_outcome_data("AFS_HPVL1_seroprev_out.csv", sep="\t", snp_col="rs", beta_col = "est", se_col="SE", effect_allele_col="eff_all",other_allele_col = "other_all", pval_col = "P")
dat <- harmonise_data(exposure_dat, outcome_dat)

devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)


#HPV E7 Seroprevalence (Definition 2)
exposure_dat <- read_exposure_data("./exposure_AFS.csv", sep=",")
outcome_dat <- read.csv("./AFS_HPVE7_seroprev_out.csv", header=TRUE)
outcome_dat <- read_outcome_data("AFS_HPVE7_seroprev_out.csv", exposure_dat$SNP, sep=",")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "HPV16 E7 Seroprevalence"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Age at first sex","HPV16 E7 Seroprevalence",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPVE7_seroprev_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
all_res <- combine_all_mrresults(mr_results,Het,Plt,Sin,ao_slc=F,Exp=T)
write.table(all_res, "AFS_HPVE7_seroprev_allres.txt")
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
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 

# outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
# outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "AFS_HPVE7_seroprev_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPVE7_seroprev_nooutliers.txt")

#MR presso 
outcome_dat <- read_outcome_data("AFS_HPVE7_seroprev_out.csv", sep="\t", snp_col="rs", beta_col = "est", se_col="SE", effect_allele_col="eff_all",other_allele_col = "other_all", pval_col = "P")
dat <- harmonise_data(exposure_dat, outcome_dat)

devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)


#HPV E6 Seroprevalence (Definition 2)
exposure_dat <- read_exposure_data("./exposure_AFS.csv", sep=",")
outcome_dat <- read.csv("./AFS_HPVE6_seroprev_out.csv", header=TRUE)
outcome_dat <- read_outcome_data("AFS_HPVE6_seroprev_out.csv", exposure_dat$SNP, sep=",")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "HPV16 E6 Seroprevalence"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Age at first sex","HPV16 E6 Seroprevalence",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPVE6_seroprev_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
all_res <- combine_all_mrresults(mr_results,Het,Plt,Sin,ao_slc=F,Exp=T)
write.table(all_res, "AFS_HPVE6_seroprev_allres.txt")
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
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 

# outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
# outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "AFS_HPVE6_seroprev_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPVE6_seroprev_nooutliers.txt")

#MR presso 
outcome_dat <- read_outcome_data("AFS_HPVE6_seroprev_out.csv", sep="\t", snp_col="rs", beta_col = "est", se_col="SE", effect_allele_col="eff_all",other_allele_col = "other_all", pval_col = "P")
dat <- harmonise_data(exposure_dat, outcome_dat)

devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)


#HPV18 Seroprevalence
exposure_dat <- read_exposure_data("./exposure_AFS.csv", sep=",")
outcome_dat <- read.csv("./AFS_HPV18_seroprev_out.csv", header=TRUE)
outcome_dat <- read_outcome_data("AFS_HPV18_seroprev_out.csv", exposure_dat$SNP, sep=",")
dat <- harmonise_data(exposure_dat, outcome_dat)
dat$outcome <- "HPV18 Seroprevalence"
mr_results <- mr(dat)
mr_results
or_results <- generate_odds_ratios(mr_results)
or_results
mr_report(dat, output_path = "~/filepath", author="Mark Gormley", study = paste("Age at first sex","HPV18 Seroprevalence",sep=""))
results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPV18_seroprev_results.txt")

Het<-mr_heterogeneity(dat)
Plt<-mr_pleiotropy_test(dat)
Sin<-mr_singlesnp(dat)
all_res <- combine_all_mrresults(mr_results,Het,Plt,Sin,ao_slc=F,Exp=T)
write.table(all_res, "AFS_HPV18_seroprev_allres.txt")
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
library(RadialMR)
library(MendelianRandomization)

dat <- dat[dat$SNP%in%res_single$SNP,]

raddat <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
ivwrad <- ivw_radial(raddat, alpha=0.05/63, weights=3)
dim(ivwrad$outliers)[1] 

# outliers at bonf 
ivwrad <- ivw_radial(raddat, alpha=0.05, weights=3)
dim(ivwrad$outliers)[1] 
#2 outliers at 0.05

eggrad <- egger_radial(raddat, alpha=0.05, weights=3)
eggrad$coef 
dim(eggrad$outliers)[1] 
# outliers at 0.05 

#plot_radial(ivwrad, TRUE, FALSE, TRUE)
plot_radial(c(ivwrad,eggrad), TRUE, FALSE, TRUE)

ivwrad$qstatistic 
ivwrad$sortoutliers <- ivwrad$outliers[order(ivwrad$outliers$p.value),]
ivwrad$sortoutliers$Qsum <- cumsum(ivwrad$sortoutliers$Q_statistic)
ivwrad$sortoutliers$Qdif <- ivwrad$sortoutliers$Qsum - ivwrad$qstatistic
write.csv(ivwrad$sortoutliers, "AFS_HPV18_seroprev_outliers.csv", row.names=F, quote=F)

#Remove top outliers
dat2 <- dat[!dat$SNP %in% ivwrad$outliers$SNP,]
mr_results2 <- mr(dat2)
or_results <- generate_odds_ratios(mr_results2)

results<-cbind.data.frame(or_results$outcome,or_results$nsnp,or_results$method,or_results$b,or_results$se,or_results$pval,or_results$or,or_results$or_lci95,or_results$or_uci95)
write.table(results, "AFS_HPV18_seroprev_nooutliers.txt")

#MR presso 
outcome_dat <- read_outcome_data("AFS_HPV18_seroprev_out.csv", sep="\t", snp_col="rs", beta_col = "est", se_col="SE", effect_allele_col="eff_all",other_allele_col = "other_all", pval_col = "P")
dat <- harmonise_data(exposure_dat, outcome_dat)

devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)


# title: "BirdScan analysis - model fitting"
# author: "Mikko Jimenez"
# date: "2023-10-12"

# load packages
library(tidyverse)
library(mgcv)
library(lme4)
library(AICcmodavg)
library(flextable)
library(officer)
library(broom)
library(MuMIn)
library(merTools)
library(lmerTest)
library(dplyr)
library(viridis)

# read in tot_df_dens
tot_df_dens <- read.csv("/Users/mikkojimenez/Downloads/26828305/tot_df_dens.csv")
tot_df_dens_meannight <- read.csv("/Users/mikkojimenez/Downloads/26828305/tot_df_dens_meannight.csv")

# convert site to factor
tot_df_dens_meannight$site<-as.factor(tot_df_dens_meannight$site)
tot_df_dens_meannight$date<-as.factor(tot_df_dens_meannight$date)
tot_df_dens$site<-as.factor(tot_df_dens$site)
tot_df_dens$date<-as.factor(tot_df_dens$date)

# get mean and SD of nightly dens as sum stat
nightlymean_sum <- tot_df_dens_meannight %>%
  group_by(site) %>%
  summarise(
    mean_dens_mean = mean(mean_dens, na.rm = TRUE),
    mean_dens_sd = sd(mean_dens, na.rm = TRUE)
  )

# Perform ANOVA and print summary
anova_result <- aov(mean_dens ~ site, data = tot_df_dens_meannight)
summary(anova_result) 

# calculate correlation coefficient for temporal scale x dataset
cor.test(tot_df_dens_meannight$mean_dens, tot_df_dens_meannight$mean_uncorrected)
cor.test(tot_df_dens_meannight$mean_dens, tot_df_dens_meannight$mean_corrected)
cor.test(tot_df_dens$dens, tot_df_dens$mean_eta_uncorrected)
cor.test(tot_df_dens$dens, tot_df_dens$mean_eta_corrected)

# ACROSS NIGHT MODELS
# fit models
## uncorrected
meannight_uncorr_null_lm<-lm(mean_dens~mean_uncorrected, data=tot_df_dens_meannight)
meannight_uncorr_lm<-lm(mean_dens~mean_uncorrected + dist*mean_uncorrected + block*mean_uncorrected, data=tot_df_dens_meannight)
meannight_uncorr_int_lm<-lm(mean_dens~mean_uncorrected + dist*block*mean_uncorrected, data=tot_df_dens_meannight)
## corrected
meannight_corr_null_lm<-lm(mean_dens~mean_corrected, data=tot_df_dens_meannight)
meannight_corr_lm<-lm(mean_dens~mean_corrected + dist*mean_corrected + block*mean_corrected, data=tot_df_dens_meannight)
meannight_corr_int_lm<-lm(mean_dens~mean_corrected + dist*block*mean_corrected, data=tot_df_dens_meannight)

## adjusted r2s
meannight_uncorr_null_sum<-summary(meannight_uncorr_null_lm)
meannight_uncorr_null_sum$adj.r.squared
meannight_uncorr_lm_sum<-summary(meannight_uncorr_lm)
meannight_uncorr_lm_sum$adj.r.squared
meannight_uncorr_int_lm_sum<-summary(meannight_uncorr_int_lm)
meannight_uncorr_int_lm_sum$adj.r.squared

meannight_corr_null_sum<-summary(meannight_corr_null_lm)
meannight_corr_null_sum$adj.r.squared
meannight_corr_sum<-summary(meannight_corr_lm)
meannight_corr_sum$adj.r.squared
meannight_corr_int_sum<-summary(meannight_corr_int_lm)
meannight_corr_int_sum$adj.r.squared

# likelihood ratio test
anova(meannight_corr_null_lm, meannight_corr_lm)
# summary of top model
summary(meannight_corr_lm)

# WITHIN NIGHT MODELS
## fit models
night_uncorr_null_lmm=lmer(dens~mean_eta_uncorrected + (1|date), data=tot_df_dens)
night_uncorr_lmm=lmer(dens~mean_eta_uncorrected + block_scaled*mean_eta_uncorrected + dist_scaled*mean_eta_uncorrected + (1|date), data=tot_df_dens)
night_uncorr_int_lmm=lmer(dens~mean_eta_uncorrected + block_scaled*dist_scaled*mean_eta_uncorrected + (1|date), data=tot_df_dens)

night_corr_null_lmm=lmer(dens~mean_eta_corrected + (1|date), data=tot_df_dens)
night_corr_lmm=lmer(dens~mean_eta_corrected + block_scaled*mean_eta_corrected + dist_scaled*mean_eta_corrected + (1|date), data=tot_df_dens)
night_corr_int_lmm=lmer(dens~mean_eta_corrected + block_scaled*dist_scaled*mean_eta_corrected + (1|date), data=tot_df_dens)

## summarize models
# get r2 values
MuMIn::r.squaredGLMM(night_uncorr_null_lmm)
MuMIn::r.squaredGLMM(night_uncorr_lmm)
MuMIn::r.squaredGLMM(night_uncorr_int_lmm)
MuMIn::r.squaredGLMM(night_corr_null_lmm)
MuMIn::r.squaredGLMM(night_corr_lmm)
MuMIn::r.squaredGLMM(night_corr_int_lmm)

# likelihood ratio test
anova(night_uncorr_null_lmm, night_uncorr_int_lmm)
# summary of top model
summary(night_uncorr_int_lmm)

# Create a vector of thresholds
thresholds <- c(30, 45, 65, 85, 100)

# Initialize an empty list to store subsets
subset_list <- list()

# Create subsets using for loop
for (thresh in thresholds) {
  # Subset the dataframe
  subset_name <- paste("tot_df_dens_", thresh, sep = "")
  subset_df <- tot_df_dens[tot_df_dens$dist <= thresh, ]
  
  # Assign the subset to a list element
  subset_list[[subset_name]] <- subset_df
}

# create distance threshold datasets from tot_df_dens
tot_df_dens_30 <- subset_list[["tot_df_dens_30"]]
tot_df_dens_45 <- subset_list[["tot_df_dens_45"]]
tot_df_dens_65 <- subset_list[["tot_df_dens_65"]]
tot_df_dens_85 <- subset_list[["tot_df_dens_85"]]

# Loop over each suffix and create the models
for (threshold in thresholds) {
  # Dynamically get the dataset from the list
  dataset <- subset_list[[paste0("tot_df_dens_", threshold)]]
  
  # Dynamically create model names
  uncorrected_model_name <- paste0("night_uncorr_int_", threshold)
  corrected_model_name <- paste0("night_corr_int_", threshold)
  
  # Create models and assign them to the dynamically created names
  assign(uncorrected_model_name, 
         lmer(dens ~ mean_eta_uncorrected + block_scaled * dist_scaled * mean_eta_uncorrected + (1 | date), 
              data = dataset))
  
  assign(corrected_model_name, 
         lmer(dens ~ mean_eta_corrected + block_scaled * dist_scaled * mean_eta_corrected + (1 | date), 
              data = dataset))
}

# get r2 values
MuMIn::r.squaredGLMM(night_uncorr_int_30)
MuMIn::r.squaredGLMM(night_corr_int_30)
MuMIn::r.squaredGLMM(night_uncorr_int_45)
MuMIn::r.squaredGLMM(night_corr_int_45)
MuMIn::r.squaredGLMM(night_uncorr_int_65)
MuMIn::r.squaredGLMM(night_corr_int_65)
MuMIn::r.squaredGLMM(night_uncorr_int_85)
MuMIn::r.squaredGLMM(night_corr_int_85)
MuMIn::r.squaredGLMM(night_uncorr_int_lmm)
MuMIn::r.squaredGLMM(night_corr_int_lmm)
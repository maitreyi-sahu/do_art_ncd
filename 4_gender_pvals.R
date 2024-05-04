# MSahu
# Feb 24, 2024

# Chi-squared and t-tests for gender comparisons

# ------------------------------------------------------------------------------


rm(list=ls())

library(tidyverse, ggplot2)

# Load data

dir <- "C:/Users/msahu/OneDrive - UW/Documents/Research/DGH+MGH/DO_ART_NCD/"
in_dir <- paste0(dir, "0_data/3_cleaned_data/")
out_dir <-  paste0(dir, "/3_plots/")

# Subset to South Africa only!!
# NOTE: we combine the community and hybrid arm (community follow-up), compared to clinic

data <- readRDS(paste0(in_dir, "ncd_merged_subset.rds")) %>% 
  filter(country == "South Africa") %>% 
  mutate(site = factor (site,
                        levels = c("HSRC",
                                   "AHRI"))) %>% 
  mutate(armR = case_when(arm == "Clinic arm" ~ "Clinic follow-up",
                          arm == "Community arm" ~ "Community follow-up",
                          arm == "Hybrid arm" ~ "Community follow-up" )) %>% 
  filter(!is.na(exit_viral_load_suppressed)) 


# ------------------------------------------------------------------------------

bmi_test <- t.test(table(data$genderR, data$bmi))
bmi_test$p.value

bmi_cat_test <- chisq.test(table(data$genderR, data$bmi_cat_exit))
bmi_cat_test$p.value

smoking_test <- chisq.test(table(data$genderR, data$smokingR))
smoking_test$p.value

bp_test <- chisq.test(table(data$genderR, data$bp_cat_exit))
bp_test$p.value

lipid_test <- t.test(table(data$genderR, data$lipid_result_exitR))
lipid_test$p.value


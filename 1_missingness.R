# From class project (SAP code -- RMD file)
# MSahu
# Updated March 5, 2024

# ------------------------------------------------------------------------------

rm(list=ls())

# Packages
pacman::p_load(naniar, ggplot2, dplyr)

# Directories
dir <- "C:/Users/msahu/OneDrive - UW/Documents/Research/DGH+MGH/DO_ART_NCD/"
in_dir <- paste0(dir, "0_data/3_cleaned_data/")
out_dir <-  paste0(dir, "/3_plots/")

# Load data
ncd_merged <- readRDS(paste0(in_dir,"ncd_merged.rds"))
ncd_merged_sa <- readRDS(paste0(in_dir,"ncd_merged_subset.rds")) %>% 
  filter(country == "South Africa")

# ------------------------------------------------------------------------------

# Visualize Missingness - by site, for outcome vars

ncd_outcomes <- ncd_merged_sa %>% 
  select(siteR,
         exit_viral_load_suppressed,
         #baseline_bmi, 
         bmi, 
         hemoglobin_a1c_result, 
         lipid_result,
         #sbp_mean_baseline, 
         sbp_mean_exit,
         #dbp_mean_baseline, 
         dbp_mean_exit)


vis_miss(ncd_outcomes) 
gg_miss_fct(x = ncd_outcomes, fct = exit_viral_load_suppressed)  + labs(title = "Missing data by viral suppression status")
ggplot(ncd_outcomes,
       aes(x = siteR,
           y = lipid_result)) +
  geom_miss_point() + 
  facet_wrap(~exit_viral_load_suppressed)

# Replace missing values with NA
ncd_outcomes[ncd_outcomes==""]<-NA
ncd_outcomes[ncd_outcomes=="Unknown"]<-NA
ncd_outcomes[ncd_outcomes=="Don't Know"]<-NA

# Save pdf - by site
pdf(file = paste0(out_dir, "supplement_missing_data_by_site.pdf"),
    width = 8, 
    height = 4) 

gg_miss_fct(x = ncd_outcomes, fct = siteR) + labs(title = "Missing data by trial site")

dev.off()

# By viral suppression

pdf(file = paste0(out_dir, "supplement_missing_data_by_vs.pdf"),
    width = 8, 
    height = 4) 

gg_miss_fct(ncd_outcomes, fct = exit_viral_load_suppressed) + labs(title = "Missing data by viral suppression status at exit")

dev.off()



# ------------------------------------------------------------------------------

# Visualize missingness - ALL VARS

# NCD
ncd_outcomes <- ncd_merged_sa %>% 
  
  select(
         # site
        pid, 
        siteR,
         
         # demographics  
         age,
         gender,
         education,
         occupation,
         
         # risk 
         smoking_status,
         days_exercise,
         eat_vegetables,
         ever_had_stroke,
         
         # baseline
         baseline_bmi,
         sbp_mean_baseline,
         dbp_mean_baseline,
         
         #outcomes
         exit_viral_load_suppressed,
         bmi, 
         hemoglobin_a1c_result, 
         lipid_result,
         sbp_mean_exit,
         dbp_mean_exit)


# Replace missing values with NA
ncd_outcomes[ncd_outcomes==""]<-NA
ncd_outcomes[ncd_outcomes=="Unknown"]<-NA
ncd_outcomes[ncd_outcomes=="Don't Know"]<-NA


# Save pdf - overall missingness
pdf(file = paste0(out_dir, "supplement_missing_data.pdf"),
    width = 12, 
    height = 8) 

missmap(ncd_outcomes,
        main = "Missingness of CVD risk data in the DO ART Study, South Africa (N = 1010)",
        rank.order = T,
        margins = c(10, 5)
        ) 

dev.off()

# ------------------------------------------------------------------------------

# IMPUTATION -- ABANDONED BECAUSE TOO MUCH COLLINEARITY BETWEEN CVD VARIABLES

str(ncd_outcomes)
imputation_model_types <- c("polyreg", 
                            "norm.nob", "logreg", "polyreg", "polyreg", 
                            "polyreg", "polyreg", "polyreg", "logreg",
                            "norm.nob", "norm.nob", "norm.nob", 
                            "logreg", "norm.nob", "norm.nob", "norm.nob", "norm.nob", "norm.nob") 

# MCMC imputation using Amelia
set.seed(12345)
mice_ncd <- mice(ncd_outcomes, m=10, maxit=50, seed=12345, method=imputation_model_types)
summary(mice_ncd)

# 10 imputed datasets

# First imputation and second imputation have different values for BMI for person 1 because each was imputed
mice_ncd$imp$ever_had_stroke[1:2,]
mice_ncd$imp$ever_had_stroke[1:2,]
mheart5[1:2,]


mice::complete(mice_ncd, 1)[1:50,]
mice::complete(mice_imputed, 2)[1:2,]

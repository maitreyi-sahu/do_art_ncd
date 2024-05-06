# DO ART NCD DATA CLEANING
# MSahu
# April 14, 2022

# ==============================================================================================

rm(list=ls())

# Load packages

library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)

# WHO CVD Risk Calculator: 

# install.packages("remotes")
# remotes::install_github("DylanRJCollins/whoishRisk")

library(whoishRisk)

# Notes:
# https://pubmed-ncbi-nlm-nih-gov.offcampus.lib.washington.edu/28357040/
# https://rdrr.io/github/DylanRJCollins/whoishRisk/

# ================================================================================================

# Directories

dir <- "C:/Users/msahu/OneDrive - UW/Documents/Research/DGH+MGH/DO_ART_NCD/"
in_dir <- paste0(dir, "0_data/2_raw_data/")
out_dir <- paste0(dir, "0_data/3_cleaned_data/")

# Load data

ncd <- read.csv(paste0(in_dir,"doart_ncd_MS_v3.csv"), stringsAsFactors = T) %>% setDT() 
bp <- read.csv(paste0(in_dir, "doart_blood_pressure_MS.csv")) %>%  setDT() 
health_cond <- read.csv(paste0(in_dir, "doart_health_conditions_MS.csv")) %>%  setDT() 
meds <- read.csv(paste0(in_dir, "doart_conmeds_MS.csv")) %>%  setDT() 

# ================================================================================================

# REFORMAT DATA

# Keep BP results from enrollment and exit. Take average of SBP and DBP results
bp$sbp_mean <- apply(bp[ , c("sbp_1", "sbp_2", "sbp_3")], MARGIN = 1, FUN = mean)
bp$dbp_mean <- apply(bp[ , c("dbp_1", "dbp_2", "dbp_3")], MARGIN = 1, FUN = mean)    

bpR <- bp %>% mutate(visit2 = case_when(visit=="Enrollment" ~ "baseline",
                                        visit=="Exit" ~ "exit")) %>% 
  filter(visit2 %in% c("baseline", "exit")) %>% 
  dcast(pid ~ visit2, value.var = c("sbp_mean", "dbp_mean"))

ncd_merged <- ncd %>% left_join(bpR, by = "pid") 

# ================================================================================================

# Recode outcomes

ncd_merged <- ncd_merged %>% 
  
  mutate(
    
    # -------------------------------------------------------------------------------------------
    
    # BLOOD PRESSURE
    
    # BP - categorical
    # Mutually exclusive and collectively exhaustive categories:
    # https://www.health.harvard.edu/heart-health/reading-the-new-blood-pressure-guidelines
    # Code from: https://www.r4epi.com/conditional-operations.html
    
    bp_cat_baseline = factor( 
      
      case_when( 
        
        sbp_mean_baseline < 120 & dbp_mean_baseline < 80 ~ "Normal",
        
        sbp_mean_baseline >= 120 & sbp_mean_baseline< 130 & dbp_mean_baseline < 80 ~ "Elevated",
        
        sbp_mean_baseline  >= 130 & sbp_mean_baseline  < 140 | dbp_mean_baseline >= 80 & dbp_mean_baseline < 90 ~ "Hypertension Stage 1",
        
        sbp_mean_baseline  >= 140 & sbp_mean_baseline <= 180 | dbp_mean_baseline >= 90 & dbp_mean_baseline <= 120 ~ "Hypertension Stage 2",
        
        sbp_mean_baseline > 180 | dbp_mean_baseline > 120 ~ "Crisis", 
        
        is.na(sbp_mean_baseline) | is.na(dbp_mean_baseline) ~ "Missing"),
      
      levels = c("Normal", "Elevated", "Hypertension Stage 1", "Hypertension Stage 2", "Missing")), # No clients with crisis
    
    bp_cat_exit = factor( 
      
      case_when( 
      
        sbp_mean_exit < 120 & dbp_mean_exit < 80 ~ "Normal",
        
        sbp_mean_exit >= 120 & sbp_mean_exit< 130 & dbp_mean_exit < 80 ~ "Elevated",
        
        sbp_mean_exit  >= 130 & sbp_mean_exit  < 140 | dbp_mean_exit >= 80 & dbp_mean_exit < 90 ~ "Hypertension Stage 1",
        
        sbp_mean_exit >= 140 & sbp_mean_exit <= 180 | dbp_mean_exit >= 90 & dbp_mean_exit <= 120 ~ "Hypertension Stage 2",
        
        sbp_mean_exit > 180 | dbp_mean_exit > 120 ~ "Crisis",
        
        is.na(sbp_mean_exit) | is.na(dbp_mean_exit) ~ "Missing"),
      
      levels = c("Normal", "Elevated", "Hypertension Stage 1", "Hypertension Stage 2", "Missing")),  # No clients with crisis
    
    
    # BP: differences from baseline to endline
    
    sbp_diff = sbp_mean_baseline - sbp_mean_exit,
    dbp_diff = dbp_mean_baseline - dbp_mean_exit,
    
    # BP -- elevated (binary)
    
    bp_hi_baseline = ifelse(sbp_mean_baseline >= 120, 1, 0),
    bp_hi_exit = ifelse(sbp_mean_exit >= 120, 1, 0),
    
    # ----------------------------------------------------------------------------------------------------------
    
    # BMI
    
    bmi_diff = bmi - baseline_bmi,
    bmi_any = ifelse(is.na( baseline_bmi), bmi, baseline_bmi),
    
    # BMI categories 
    
    bmi_cat_baseline = case_when(
                                  baseline_bmi < 18.5 ~ "<18.5",
                                  baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "18.5-24.9",
                                  baseline_bmi >= 25 & baseline_bmi < 30 ~ "25-29.9",
                                  baseline_bmi >= 30 ~ "30+",
                                  is.na(baseline_bmi) ~ "Missing"),
    
    bmi_cat_exit = case_when( bmi < 18.5 ~ "<18.5",
                                  bmi >= 18.5 & bmi < 25 ~ "18.5-24.9",
                                  bmi >= 25 & bmi < 30 ~ "25-29.9",
                                  bmi >= 30 ~ "30+",
                                  is.na(bmi) ~ "Missing"),
    
    
    # BMI - elevated (binary)
    
    overwt_baseline = ifelse( baseline_bmi >= 25, 1, 0),
    overwt_exit = ifelse( bmi >= 25, 1, 0),
    
    # -----------------------------------------------------------------------------------------------------------
    
    # BLOOD SUGAR 
    
    # A1C - categorical
    
    a1c_cat_exit = factor(case_when( hemoglobin_a1c_result < 5.7 ~ "Normal",
                                     hemoglobin_a1c_result >= 5.7 & hemoglobin_a1c_result < 6.5 ~ "Prediabetes",
                                     hemoglobin_a1c_result >= 6.5 ~ "Diabetes",
                                     is.na(hemoglobin_a1c_result) ~ "Missing"),
                          levels = c("Normal", "Prediabetes", "Diabetes", "Missing")),
    
    # A1C (binary)
    
    a1c_hi_exit = ifelse(hemoglobin_a1c_result >= 5.7, 1, 0),
    a1c_diabetes_exit = ifelse(hemoglobin_a1c_result >= 6.5, 1, 0),
    
    # -------------------------------------------------------------------------------------------------------
    
    # LIPIDS 
    
    # Total cholesterol - categorical
    # https://medlineplus.gov/lab-tests/cholesterol-levels/
    
    total_cholesterol_exit = factor(case_when( lipid_result < 200 ~ "Normal",
                                               lipid_result >= 200 & lipid_result < 240 ~ "Borderline High",
                                               lipid_result >= 240 ~ "High",
                                               is.na(lipid_result) ~ "Missing"),
                    levels = c("Normal", "Borderline High", "High", "Missing")),
    
    # Total cholesterol (binary)
    
    total_cholesterol_hi_exit = ifelse(lipid_result >= 200, 1, 0),
    
    # Recode lipid result in mmol/ L for WHO package
    
    # https://www.omnicalculator.com/health/cholesterol-units#:~:text=The%20rules%20for%20converting%20cholesterol,mmol%2FL%20multiply%20by%200.02586%20.
    # To get from mg/dL to mmol/L multiply by 0.02586
    
    lipid_result_exitR = ifelse(is.na(lipid_result), 0, lipid_result*0.02586) 
)


# ================================================================================================

# WHO RISK SCORE -- Exit only

# Calculated using "whoishRisk" package
# Arguments are described here:
# https://www-ncbi-nlm-nih-gov.offcampus.lib.washington.edu/pmc/articles/PMC5345772/

ncd_merged <- ncd_merged %>% 
  
  # Recode gender
  
  mutate(genderB = ifelse(gender=="Female", 0, 1)) %>% 
  
  # Recode smoking
  
  mutate(smokingR = case_when ( smoking_status %in% c("Not at all", 
                                                      "I used to smoke, but I quit") ~ 0,
                                smoking_status %in% c("Yes, occasionally",
                                                      "Yes, daily or most days") ~ 1)) 

ncd_cvd_subset <- ncd_merged %>%   
  
  # Package will give error for missing values - need to drop if N/A
  
  filter(!is.na(sbp_mean_exit)) %>% 
  filter(!is.na(a1c_diabetes_exit))  

ncd_cvd_subset <- ncd_cvd_subset %>%

# RUN WHO CVD Calculator  
    
  mutate(
    who_cvd_risk_exit = 
           WHO_ISH_Risk(
             age = ncd_cvd_subset$age, # numeric age
             gdr = ncd_cvd_subset$genderB, # binary (0 - female, 1 = male)
             smk = ncd_cvd_subset$smokingR, # binary (0 = non-smoker, 1 = smoker)
             sbp = ncd_cvd_subset$sbp_mean_exit, # continuous (mmHg)
             dm = ncd_cvd_subset$a1c_diabetes_exit, # binary (0=not diabetic; 1 = diabetic)
             chl = ncd_cvd_subset$lipid_result_exitR, # Continuous (mmol/ L); 0=unknown cholesterol)
             
             # Code WHO subregion: 
             # https://plos.figshare.com/articles/dataset/World_Health_Organization_WHO_Member_States_by_subregion_/2579569
             subregion = "AFR_E"))  # South Africa and Uganda both in Africa E subregion 

# Merge back with main dataset, recode NA as "unknown"

ncd_merged <- ncd_merged %>% 
  left_join(ncd_cvd_subset[ , c("pid", "who_cvd_risk_exit")], by = "pid") %>%  
  mutate(who_cvd_risk_exit = ifelse(is.na(who_cvd_risk_exit), "Missing", who_cvd_risk_exit))
  
# Create binary "elevated risk" or not  

ncd_merged <- ncd_merged %>%   
  mutate(who_cvd_hi_exit = 
           case_when(
             who_cvd_risk_exit == "<10%" ~ 0,
             who_cvd_risk_exit == "10% to <20%" ~ 1,
             who_cvd_risk_exit == "20% to <30%" ~ 1,
             who_cvd_risk_exit == "30% to <40%" ~ 1
             ))

rm(ncd_cvd_subset)
  
# =================================================================================================

# Recode categorical vars

ncd_merged <- ncd_merged %>% 
 
   mutate(
        
        # version of trial arms with only 2 levels
        armR = factor(case_when(arm == "Clinic arm" ~ "Clinic follow-up",
                                    arm == "Community arm" ~ "Community follow-up",
                                    arm == "Hybrid arm" ~ "Community follow-up" ),
                          levels = c("Clinic follow-up", "Community follow-up")),
     
        site = factor(site,
                      levels = c("ICOBI", "HSRC", "AHRI")),
     
        # Site
         siteR = factor(case_when(site == "ICOBI" ~ "SW Uganda",
                          site == "AHRI" ~ "Northern KZN SA",
                          site == "HSRC" ~ "Midlands KZN SA"),
                        levels = c("SW Uganda", "Midlands KZN SA", "Northern KZN SA")),
         
         # Country
         country = case_when(site == "ICOBI" ~ "Uganda",
                           site == "AHRI" ~ "South Africa",
                           site == "HSRC" ~ "South Africa"),
         
         # Age
         age_catR = factor(case_when(age >= 18 & age <40 ~ "18-39",
                                   age >= 40 & age < 60 ~ "40-59", 
                                   age >= 60 ~ "60+"),
                         levels = c("18-39", "40-59", "60+")),
        
         # age_catR = factor(case_when(age >= 18 & age <30 ~ "18-29",
         #                     age >= 30 & age < 45 ~ "30-44",        
         #                     age >= 45 & age < 60 ~ "45-59",
         #                     age >= 60 ~ "60+"),
         #                  levels = c("18-29", "30-44", "45-59", "60+")),
         
         # Gender
         genderR = factor(ifelse(gender=="Female", "Women", "Men"),
                          levels = c("Women", "Men")), 
        
        # Education
        educationR = factor(case_when(education == "Primary" ~ "Primary",
                                      education == "Secondary" ~ "Secondary",
                                      education == "Tertiary and above" ~ "Tertiary+",
                                      education == "Unknown" ~ "Missing",
                                      education == "Declined to answer" ~ "Missing"),
        levels = c("Primary", "Secondary", "Tertiary+", "Missing")),
         
         # Viral suppression at exit 
         
         exit_viral_load_suppressedR = factor(case_when(
           exit_viral_load_suppressed == T ~ "Yes",
           exit_viral_load_suppressed == F ~ "No",
           is.na(exit_viral_load_suppressed) ~ "Missing"), 
           levels = c("Yes", "No", "Missing")),
         
         # Education
         edu_catR = case_when(education == "Primary" ~ 0,
                             education == "Secondary" ~ 1,
                             education == "Tertiary and above" ~ 2,
                             education == "Declined to answer" | education == "Unknown" ~ 3),
         
         # Occupation
         occ_catR = case_when(occupation == "Unemployed" ~ 0,
                             occupation == "Farming/animal raising" ~ 1,
                             occupation == "Labourer/semi skilled" ~ 2,
                             occupation == "Trade/sales" ~ 3,
                             occupation == "Housewife" ~ 4, 
                             occupation == "Student" ~ 5,
                             occupation == "Professional" ~ 6,
                             occupation == "Other" | occupation == "Declined to answer" ~ 7),
         
         # Smoking 
         smoking_catR = factor(case_when( smoking_status == "Not at all" ~ "Not at all",
                                  smoking_status == "I used to smoke, but I quit" ~ "Former",
                                  smoking_status == "Yes, occasionally" ~ "Current - occasional",
                                  smoking_status == "Yes, daily or most days" ~ "Current - frequent"),
                       levels = c("Not at all", "Former", "Current - occasional", "Current - frequent")),
                         
         
         smokes_day_catR = case_when(smokes_per_day >=1 & smokes_per_day <=4 ~ 0,
                                    smokes_per_day>=5 & smokes_per_day <=10 ~ 1,
                                    smokes_per_day >=11 & smokes_per_day <= 20 ~ 2),
         
         # Stroke or heart attack in the past
         stroke_catR = factor( case_when(ever_had_stroke =="No" ~ "No",
                                ever_had_stroke == "Yes" ~ "Yes", 
                                T ~ "Missing"),
                              levels = c("No", "Yes", "Missing")),
         
         # Vegetables
         veg_catR = factor(case_when(eat_vegetables == "Rarely" | eat_vegetables =="Never" ~ "Never or rarely",
                             eat_vegetables == "Sometimes" ~ "Sometimes",
                             eat_vegetables =="Always" | eat_vegetables =="Usually" ~ "Always or usually",
                             eat_vegetables == "Refused" ~ "Missing",
                             T ~ "Missing"),
                          levels = c("Always or usually", "Sometimes", "Never or rarely", "Missing")),
         
         # Exercise
         exer_catR = factor(case_when(days_exercise == "0" ~ "0",
                                     days_exercise == "1-2" ~ "1-2",
                                     days_exercise == "3-4" ~ "3-4",
                                     days_exercise == "5-7" ~ "5-7",
                                     days_exercise == "Don't Know" ~ "Missing",
                                     T ~ "Missing"),
                           levels = c("5-7","3-4","1-2", "0", "Missing"))
)

# ==========================================================================================================

length(ncd_merged$baseline_viral_load_suppressed[ncd_merged$baseline_viral_load_suppressed==T]) #216 

ncd_merged_subset <- ncd_merged %>% 
  filter(baseline_viral_load_suppressed == F | is.na(baseline_viral_load_suppressed)) 

ncd_merged_subset %>%  nrow() # 1315

# Check missingness
ncd_merged_subset %>% filter(is.na(who_cvd_risk_exit)) %>%  nrow() # 170 missing; 13% of data
ncd_merged_subset %>% filter(sbp_mean_exit == "Missing") %>%  nrow() # 18
ncd_merged_subset %>% filter(a1c_diabetes_exit == "Missing") %>%  nrow() # 170
ncd_merged_subset %>% filter(sbp_mean_exit == "Missing" & !(a1c_diabetes_exit == "Missing")) %>%  nrow() # 0

# SAVE ------------------------------------------------------------------------------------------------------

saveRDS(ncd_merged, file = paste0(out_dir,"ncd_merged.rds"))
saveRDS(ncd_merged_subset, file = paste0(out_dir,"ncd_merged_subset.rds"))

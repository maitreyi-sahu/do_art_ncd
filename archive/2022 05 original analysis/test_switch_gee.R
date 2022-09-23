
## NEED TO ASK ADAM

library(geepack)
library(sandwich)
library(lme4) # note, NLME doesn't include GLM (Poisson)
library(clubSandwich)
library(merDeriv)

ncd_noHybrid <- ncd_merged_subset %>% 
  filter(arm!= "Hybrid arm") %>% 
  mutate(arm = factor (arm,
                       levels = c("Clinic arm",
                                  "Community arm")))

  
  # Adjusted for site only
  m <- geeglm(sbp_mean_exit ~ arm + site, 
              family = poisson, data = ncd_noHybrid,
              id = hhid, corstr = "independence") 
 
 m.rr <- exp(coef(m)[["armCommunity arm"]]) # Relative Risk
 m.se <- summary(m)$coef["armCommunity arm",2] # SE
 m.rr_ci <- exp(c(coef(m)[["armCommunity arm"]] - qnorm(0.975)*m.se, 
                  coef(m)[["armCommunity arm"]] + qnorm(0.975)*m.se)) # RR: 95% CI with robust SE
 m.p <- summary(m)$coef["armCommunity arm",4]
 
 m <- geeglm(bp_hi_exit ~ site + exit_viral_load_suppressed, 
             family = poisson, data = ncd_merged_subset,
             id = hhid, corstr = "independence") 
 
 
 
   m <- geeglm(sbp_mean_exit ~ arm + site, 
              family = gaussian(), data = ncd_noHybrid,
              id = hhid, corstr = "independence") 
  
# ----------------------------------------------------------------------------

rm(list=ls())

# Load packages

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

# Table formatting packages
library(janitor)
library(tables)
library(flextable)

# Regression packages
library(geepack)
library(sandwich)

# Load data

dir <- "C:/Users/msahu/Documents/Other_Research/DO_ART_NCD/"
in_dir <- paste0(dir, "0_data/3_cleaned_data/")

ncd_merged_subset <- readRDS(paste0(in_dir, "ncd_merged_subset.rds"))
ncd_noHybrid <- ncd_merged_subset %>% filter(arm!= "Hybrid arm")

# Set up --------------------------------------------------------------------------

headings <- data.frame(
  section = c("A", "B", "C"),
  title = c( "", 
             "ii) Lifestyle and other risk",
             "iii) Cardiovascular Risk")
)

# A) Demographics
A_vars = c(
  "age_cat",
  "genderR"
)
A_var_labels = c(
  "Age",
  "Gender"
)

# B) Lifestyle and Other Risk

B_vars = c(  
  "smoking_cat",
  "exer_cat",
  "veg_cat",
  "stroke_cat")
B_var_labels = c(
  "Smoking Status  \n [Baseline]",
  "Days of exercise (per week) \n [Endline]",
  "Vegetable intake \n [Endline]",
  "Prior stroke or heart attack  \n [Endline]")


# C) Cardiovascular risk
C_vars = c(
  "bp_cat_exit", 
  "bmi_cat_exit", 
  "a1c_cat_exit", 
  "total_cholesterol_exit",
  "who_cvd_risk_exit")
C_var_labels = c(
  "Blood pressure \n [Endline]",
  "BMI (kg/m^2) \n [Endline]",
  "Hemoglobin A1C (%) \n [Endline]",
  "Total cholesterol (mg/dL) \n [Endline]",
  "10-year CVD risk score \n [Endline]")

table_vars <- data.frame(
  
  var_name = c(A_vars, B_vars, C_vars),
  
  var_label = c(A_var_labels, B_var_labels, C_var_labels)) %>% 
  
  mutate( 
    
    section = case_when(
      var_name %in% c(A_vars) ~ "A",
      var_name %in% c(B_vars) ~ "B",
      var_name %in% c(C_vars) ~ "C")
  )

# Descriptive Proportion Table ----------------------------------------------------

col_headings_freq <- c(
  "Category",
  "Clinic Arm \n n (%)", 
  "Community Arm \n n (%)",
  "Hybrid Arm \n n (%)",
  "Women \n n (%)",
  "Men \n n (%)",
  "Total \n n (%)"
)

for (i in unique(headings$section)) {
  
  table_subset <- table_vars[table_vars$section==i, ]
  
  tableDF <- data.frame(
    variable = character(),
    category = character(),
    clinic_arm = character(),
    community_arm = character(),
    hybrid_arm = character(),
    Women = character(),
    Men = character(),
    Total = character()
  )
  
  for (v in unique(table_subset$var_name)) {
    
    data = ncd_merged_subset
    
    temp1 <- data %>% 
      
      # Format table using Janitor package 
      tabyl(!!sym(v), arm) %>% 
      adorn_totals(where = "col") %>% 
      adorn_percentages(denominator = "col") %>% 
      adorn_pct_formatting(digits = 0, affix_sign = TRUE) %>% 
      adorn_ns(position = "front") %>%
      as.data.frame() %>% 
      
      # Clean up headers
      rename(category = 1,
             clinic_arm = `Clinic arm`,
             community_arm = `Community arm`,
             hybrid_arm = `Hybrid arm`) %>% 
      mutate(variable = table_vars[table_vars$var_name == v, "var_label"] ) %>% 
      select(variable, category, clinic_arm, community_arm, hybrid_arm, Total)
    
    if (v != "genderR") {
      
      temp2 <- data %>% 
        
        # Format table using Janitor package 
        tabyl(!!sym(v), genderR) %>% 
        adorn_totals(where = "col") %>% 
        adorn_percentages(denominator = "col") %>% 
        adorn_pct_formatting(digits = 0, affix_sign = TRUE) %>% 
        adorn_ns(position = "front") %>%
        as.data.frame() %>% 
        
        # Clean up headers
        select(Women, Men)
      
    }
    
    if (v == "genderR") {
      
      temp2 <- data.frame(
        Women = "--",
        Men = "--"
      )
      
    }
    
    temp <- cbind(temp1, temp2) %>% 
      select(variable, category, clinic_arm, community_arm, hybrid_arm, Women, Men, Total)
    
    tableDF <- rbind(tableDF, temp)
    
  }

}


results_list <- vector(mode = "list")  # BIND RESULTS FOR PLOT

# Poisson REGRESSIONS --------------------------------------------------------------

reg_vars <- data.frame(
  
  var_name = c(
    "bp_hi_exit", 
    "overwt_exit", 
    "a1c_hi_exit",
    "total_cholesterol_hi_exit"),
  
  var_label = c(
    "Elevated BP",
    "Overweight (BMI >= 25)",
    "Elevated blood sugar (a1c >= 5.7)",
    "Elevated lipids (total cholesterol >= 200)"))


regressionDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_ci = character(),
  rr_p = numeric(),
  adj.rr = character(),
  adj.rr_ci = character(), 
  adj.rr_p = character(),
  ext.rr = character(),
  ext.rr_ci = character(), 
  ext.rr_p = character() 
)

plotDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_ci = character(),
  rr_p = numeric(),
  adj.rr = character(),
  adj.rr_ci = character(), 
  adj.rr_p = character(),
  ext.rr = character(),
  ext.rr_ci = character(), 
  ext.rr_p = character() 
)

library(lme4) # note, NLME doesn't include GLM (Poisson)
library(clubSandwich)
library(merDeriv)
for (v in unique(reg_vars$var_name)) {
  
  # Adjusted for site only
  m <- geeglm(sbp_mean_exit ~ arm + site, 
              family = poisson, data = ncd_noHybrid,
              id = hhid, corstr = "exchangeable") 
  m <- glmer(sbp_mean_exit ~ arm + site + (1 | hhid), #lme4
              family = poisson, data = ncd_noHybrid) 
  m.rr <- exp(fixef(m)[["armCommunity arm"]]) # Relative Risk
  m.se <- sqrt(diag( vcovHC(m, type = "HC0")))["armCommunity arm"] # Robust SE (unexponentiated)
  m.rr_ci <- exp(c(coef(m)[["armCommunity arm"]] - qnorm(0.975)*m.se, 
                   coef(m)[["armCommunity arm"]] + qnorm(0.975)*m.se)) # RR: 95% CI with robust SE
  m.p <- 2*pnorm(abs(coef(m)[["armCommunity arm"]]/m.se), lower.tail = F)
  
  # Adjusted for site, age, gender, smoking 
  
  if (v=="who_cvd_hi_exit") { # 3-way interaction term for cvd risk only
    
    adj.m <- glm(get(v) ~ arm + site + age_cat * genderR * smoking_status, 
                 data = ncd_noHybrid,
                 family = quasipoisson(link = "log")) }
  
  if (v!="who_cvd_hi_exit") {
    
    adj.m <- glm(get(v) ~ arm + site + age_cat + genderR + smoking_status, 
                 data = ncd_noHybrid,
                 family = quasipoisson(link = "log")) }
  
  adj.rr <- exp(coef(adj.m)[["armCommunity arm"]]) # Relative Risk
  adj.se <- sqrt(diag(vcovHC(adj.m, type = "HC0")))["armCommunity arm"] # Robust SE (unexponentiated)
  adj.rr_ci <- exp(c(coef(adj.m)[["armCommunity arm"]] - qnorm(0.975)*adj.se, 
                     coef(adj.m)[["armCommunity arm"]] + qnorm(0.975)*adj.se)) # RR: 95% CI with robust SE
  adj.p <- 2*pnorm(abs(coef(adj.m)[["armCommunity arm"]]/adj.se), lower.tail = F)
  
  # Adjusted for VS and site, age, gender, smoking 
  
  if (v=="who_cvd_hi_exit") { # 3-way interaction term for cvd risk only
    
    ext.m <- glm( get(v) ~ arm + site + age_cat * genderR * smoking_status + exit_viral_load_suppressedR, 
                  data = ncd_noHybrid,
                  family = quasipoisson(link = "log")) }
  
  if (v!="who_cvd_hi_exit") { # 3-way interaction term for cvd risk only
    
    ext.m <- glm( get(v) ~ arm + site + age_cat + genderR + smoking_status + exit_viral_load_suppressedR, 
                  data = ncd_noHybrid,
                  family = quasipoisson(link = "log")) }
  
  ext.rr <- exp(coef(ext.m)[["armCommunity arm"]]) # Relative Risk
  ext.se <- sqrt(diag(vcovHC(ext.m, type = "HC0")))["armCommunity arm"] # Robust SE (unexponentiated)
  ext.rr_ci <- exp(c(coef(ext.m)[["armCommunity arm"]] - qnorm(0.975)*ext.se, 
                     coef(ext.m)[["armCommunity arm"]] + qnorm(0.975)*ext.se)) # RR: 95% CI with robust SE
  ext.p <- 2*pnorm(abs(coef(ext.m)[["armCommunity arm"]]/ext.se), lower.tail = F)
  
  temp_row <- c(reg_vars[reg_vars$var_name == v, "var_label"], 
                as.numeric(round(m.rr, round_digits)), 
                paste0("(",
                       as.numeric(round(m.rr_ci, round_digits))[1], 
                       ",",
                       as.numeric(round(m.rr_ci, round_digits))[2],
                       ")"),
                as.numeric(round(m.p, p_digits)),
                as.numeric(round(adj.rr, round_digits)), 
                paste0("(",
                       as.numeric(round(adj.rr_ci, round_digits))[1], 
                       ",",
                       as.numeric(round(adj.rr_ci, round_digits))[2],
                       ")"),
                as.numeric(round(adj.p, p_digits)),
                as.numeric(round(ext.rr, round_digits)), 
                paste0("(",
                       as.numeric(round(ext.rr_ci, round_digits))[1], 
                       ",",
                       as.numeric(round(ext.rr_ci, round_digits))[2],
                       ")"),
                as.numeric(round(ext.p, p_digits))
  )
  
  regressionDF <- rbind(regressionDF, temp_row)
  
  temp_row_plot <- c(reg_vars[reg_vars$var_name == v, "var_label"], 
                     as.numeric(m.rr), 
                     as.numeric(m.rr_ci[1]),
                     as.numeric(m.rr_ci[2]),
                     as.numeric(round(m.p, p_digits)),
                     as.numeric(adj.rr), 
                     as.numeric(adj.rr_ci[1]),
                     as.numeric(adj.rr_ci[2]),
                     as.numeric(round(adj.p, p_digits)),
                     as.numeric(ext.rr), 
                     as.numeric(ext.rr_ci[1]),
                     as.numeric(ext.rr_ci[2]),
                     as.numeric(round(ext.p, p_digits))
  )
  
  plotDF <- rbind(plotDF, temp_row_plot)
}

names(regressionDF) =
  c("Variable",
    "Relative Risk",
    "95% CI",
    "p-value",
    "Adj. Relative Risk",
    "Adj. 95% CI",
    "Adj. p-value",
    "Ext. Relative Risk",
    "Ext. 95% CI",
    "Ext. p-value")

names(plotDF) =
  c("Variable",
    "main.estimate",
    "main.lower",
    "main.upper",
    "main.p",
    "adj.estimate",
    "adj.lower",
    "adj.upper",
    "adj.p",
    "ext.estimate",
    "ext.lower",
    "ext.upper",
    "ext.p")

plotDF <-  plotDF %>% 
  gather(label, value, 2:13) %>% 
  separate(label, into = c("model", "estimate")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = estimate, value = value) %>% 
  mutate(country = "Overall",
         outcome = "relative_risk")

results_list[[length(results_list) + 1]] <- plotDF
names(results_list)[length(results_list)] <- "rr.Overall"

# Word Doc output
ft <- flextable(regressionDF) %>% 
  colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
  add_header_row(values = c("", "Baseline Analysis", "Adjusted Analysis", "Extended Analysis \n (including Viral Suppression)"), 
                 top = T, colwidths = c(1, 3, 3, 3)) %>% 
  width(width = .94) %>%
  width(j = 1, width = 1.6) %>%
  theme_vanilla() %>% 
  vline(j = c(1, 4, 7), border = fp_border_default(), part = "all") %>% 
  align(i =1, align = "center", part = "header") 

ft

# LINEAR REGRESSIONS --------------------------------------------------------------

lin_reg_vars <- data.frame(
  
  var_name = c(
    "sbp_mean_exit",
    "dbp_mean_exit", 
    "bmi", 
    "hemoglobin_a1c_result",
    "lipid_result"),
  
  var_label = c(
    "Systolic BP (mmHg)",
    "Diastolic BP (mmHg)",
    "BMI (kg/m^2)",
    "Hemoglobin A1c (%)",
    "Total cholesterol (mg/dL)") 
)


regressionDF <- data.frame(
  variable = character(),
  rd = numeric(),
  rd_ci = character(),
  rd_p = numeric(),
  adj.rd = character(),
  adj.rd_ci = character(), 
  adj.rd_p = character(), 
  ext.rd = character(),
  ext.rd_ci = character(), 
  ext.rd_p = character()
)

plotDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_lower = numeric(),
  rr_upper = numeric(),
  rr_p = numeric(),
  adj.rr = numeric(),
  adj.rr_lower = numeric(),
  adj.rr_upper = numeric(),
  adj.rr_p = numeric(),
  ext.rr = numeric(),
  ext.rr_lower = numeric(),
  ext.rr_upper = numeric(),
  ext.rr_p = numeric()
)

for (v in unique(lin_reg_vars$var_name)) {
  
  # Adjusted for site only
  lm <- lm(get(v) ~ arm + site, data = ncd_noHybrid)
  rd <- coef(lm)[["armCommunity arm"]] # Risk Difference
  lm.se <- sqrt(diag( vcovHC(lm, type = "HC0")))["armCommunity arm"] # Robust SE (unexponentiated)
  rd_ci <- c(coef(lm)[["armCommunity arm"]] - qnorm(0.975)*lm.se, 
             coef(lm)[["armCommunity arm"]] + qnorm(0.975)*lm.se) # RD: 95% CI with robust SE
  rd_p <- 2*pnorm(abs(coef(lm)[["armCommunity arm"]]/lm.se), lower.tail = F)
  
  # Adjusted for site, age, gender, smoking (no 3-way interaction effects)
  adj.lm <- lm(get(v) ~ arm + site + age_cat + genderR + smoking_status, 
               data = ncd_noHybrid)
  adj.rd <- coef(adj.lm)[["armCommunity arm"]] # Risk Difference
  adj.lm.se <- sqrt(diag(vcovHC(adj.lm, type = "HC0")))["armCommunity arm"] # Robust SE (unexponentiated)
  adj.rd_ci <- c(coef(adj.lm)[["armCommunity arm"]] - qnorm(0.975)*adj.lm.se, 
                 coef(adj.lm)[["armCommunity arm"]] + qnorm(0.975)*adj.lm.se) # RD: 95% CI with robust SE
  adj.rd_p <- 2*pnorm(abs(coef(adj.lm)[["armCommunity arm"]]/adj.lm.se), lower.tail = F)
  
  # Adjusted for VS and site, age, gender, smoking (with 3-way interaction effects for all vars)
  ext.lm <- lm(get(v) ~ arm + site + age_cat + genderR + smoking_status + exit_viral_load_suppressedR, 
               data = ncd_noHybrid)
  ext.rd <- coef(ext.lm)[["armCommunity arm"]] # Risk Difference
  ext.lm.se <- sqrt(diag(vcovHC(ext.lm, type = "HC0")))["armCommunity arm"] # Robust SE (unexponentiated)
  ext.rd_ci <- c(coef(ext.lm)[["armCommunity arm"]] - qnorm(0.975)*ext.lm.se, 
                 coef(ext.lm)[["armCommunity arm"]] + qnorm(0.975)*ext.lm.se) # RD: 95% CI with robust SE
  ext.rd_p <- 2*pnorm(abs(coef(ext.lm)[["armCommunity arm"]]/ext.lm.se), lower.tail = F)
  
  
  temp_row <- c(lin_reg_vars[lin_reg_vars$var_name == v, "var_label"], 
                as.numeric(round(rd, round_digits)), 
                paste0("(",
                       as.numeric(round(rd_ci, round_digits))[1], 
                       ",",
                       as.numeric(round(rd_ci, round_digits))[2],
                       ")"),
                as.numeric(round(rd_p, p_digits)),
                as.numeric(round(adj.rd, round_digits)), 
                paste0("(",
                       as.numeric(round(adj.rd_ci, round_digits))[1], 
                       ",",
                       as.numeric(round(adj.rd_ci, round_digits))[2],
                       ")"),
                as.numeric(round(adj.rd_p, p_digits)),
                as.numeric(round(ext.rd, round_digits)), 
                paste0("(",
                       as.numeric(round(ext.rd_ci, round_digits))[1], 
                       ",",
                       as.numeric(round(ext.rd_ci, round_digits))[2],
                       ")"),
                as.numeric(round(ext.rd_p, p_digits))
  )
  
  regressionDF <- rbind(regressionDF, temp_row)
  
  temp_row_plot <- c(lin_reg_vars[lin_reg_vars$var_name == v, "var_label"], 
                     as.numeric(rd), 
                     as.numeric(rd_ci[1]),
                     as.numeric(rd_ci[2]),
                     as.numeric(round(rd_p, p_digits)),
                     as.numeric(adj.rd), 
                     as.numeric(adj.rd_ci[1]),
                     as.numeric(adj.rd_ci[2]),
                     as.numeric(round(adj.rd_p, p_digits)),
                     as.numeric(ext.rd), 
                     as.numeric(ext.rd_ci[1]),
                     as.numeric(ext.rd_ci[2]),
                     as.numeric(round(ext.rd_p, p_digits))
  )
  
  plotDF <- rbind(plotDF, temp_row_plot)
}

names(regressionDF) =
  c("Variable",
    "Risk Difference",
    "95% CI",
    "p-value",
    "Adj. Risk Difference",
    "Adj. 95% CI",
    "Adj. p-value",
    "Ext. Risk Difference",
    "Ext. 95% CI",
    "Ext. p-value")

names(plotDF) =
  c("Variable",
    "main.estimate",
    "main.lower",
    "main.upper",
    "main.p",
    "adj.estimate",
    "adj.lower",
    "adj.upper",
    "adj.p",
    "ext.estimate",
    "ext.lower",
    "ext.upper",
    "ext.p")

plotDF <-  plotDF %>% 
  gather(label, value, 2:13) %>% 
  separate(label, into = c("model", "estimate")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = estimate, value = value) %>% 
  mutate(country = "Overall",
         outcome = "risk_difference")

results_list[[length(results_list) + 1]] <- plotDF
names(results_list)[length(results_list)] <- "rd.Overall"




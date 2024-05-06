# Functions for GEE
# 5/4/2024

# ==============================================================================

# GEE FORMULA
run_gee <- function(formula, family, df) {geeglm(formula = formula, 
                                                 family = family, 
                                                 data = df, 
                                                 id = hhid, 
                                                 corstr = "independence") }

# Set formula and adjustment variables
unadj_vars <- c("site")
adj_vars <- c(unadj_vars, "age_catR", "genderR", "smokingR", "educationR") # UPDATE HERE FOR ANY CHANGES TO CONTROL VARS
ext_var <- "baseline_cd4"

# Settings
round_digits = 2
p_digits = 3

# ==============================================================================

# CVD VARS + VAR LABELS

# Poisson (Binary)

poisson_reg_vars <- data.frame(
  
  var_name = c(
    "bp_hi_exit", 
    "overwt_exit", 
    "a1c_hi_exit",
    #  "total_cholesterol_hi_exit",
    "smokingR"),
  
  var_label = c(
    "Elevated BP",
    "Overweight (BMI >= 25)",
    "Elevated blood sugar (a1c >= 5.7)",
    #   "Elevated lipids (total cholesterol >= 200)",
    "Current smoker [baseline]"))

# Linear (Continuous)

lin_reg_vars <- data.frame(
  
  var_name = c(
    "sbp_mean_exit",
    "dbp_mean_exit", 
    "bmi", 
    "hemoglobin_a1c_result",
    "lipid_result"
   ),
  
  var_label = c(
    "Systolic BP (mmHg)",
    "Diastolic BP (mmHg)",
    "BMI (kg/m^2)",
    "Hemoglobin A1c (%)",
    "Total cholesterol (mg/dL)"
    ) 
)

# ==============================================================================

# Col headings for descriptive tables


headings <- data.frame(
  section = c("A", "B", "C"),
  title = c( "", 
             "ii) Lifestyle and other risk",
             "iii) Cardiovascular Risk")
)

# A) Demographics
A_vars = c(
  "age_catR",
  "genderR",
  "educationR"
)
A_var_labels = c(
  "Age",
  "Gender",
  "Education"
)

# B) Lifestyle and Other Risk

B_vars = c(  
  "smoking_catR",
  "exer_catR",
  "veg_catR",
  "stroke_catR")
B_var_labels = c(
  "Smoking Status  \n [Baseline]",
  "Days of exercise (per week) \n [Exit]",
  "Vegetable intake \n [Exit]",
  "Prior stroke or heart attack  \n [Exit]")


# C) Cardiovascular risk
C_vars = c(
  "bp_cat_exit", 
  "bmi_cat_exit", 
  "a1c_cat_exit", 
  "total_cholesterol_exit",
  "who_cvd_risk_exit")
C_var_labels = c(
  "Blood pressure \n [Exit]",
  "BMI (kg/m^2) \n [Exit]",
  "Hemoglobin A1C (%) \n [Exit]",
  "Total cholesterol (mg/dL) \n [Exit]",
  "10-year CVD risk score \n [Exit]")

# D) Change in CVD risk
D_vars = c(  
  "bmi_diff",
  "sbp_diff")
D_var_labels = c(
  "Change in BMI, baseline to endline  \n [AHRI only]",
  "Change in Systolic Blood pressure, baseline to endline  \n [AHRI only]" )


# ==============================================================================

# Get results for single CVD var!

# ext = T
# df <- ncd_merged_subset
# var_to_pull = "armRCommunity follow-up"
# cvd_var = poisson_reg_vars$var_name[1]
# family = "poisson"

get_gee_results <- function(df = df, 
                            family, # poisson or gaussian
                            cvd_var,
                            var_of_interest, # e.g. "armR"
                            var_to_pull, # coefficient to pull, e.g. "armRCommunity follow-up"
                            ext = ext) { # is this the extended analysis? defaults to no
  
  make_formula <- function(vars) {as.formula(paste0(cvd_var, " ~ ", paste(vars, collapse = " + ")))}
  
  # Get var labels
  if (family == "poisson") { reg_vars = poisson_reg_vars} 
  if (family == "gaussian") { reg_vars = lin_reg_vars}
  
  # Run GEE and get coeffs - unadjusted
  vars <-  c(var_of_interest, unadj_vars)
  formula <- make_formula(vars)
  n <- df %>% select(c(cvd_var, vars, hhid)) %>% filter(complete.cases(.)) %>% nrow() # N obs
  
  m <- run_gee(formula, family, df)
  if (family == "poisson") {
    beta <- exp(coef(m)[[var_to_pull]]) # Relative Risk 
    se <- summary(m)$coef[var_to_pull,2] # SE
    ci <- exp(c(coef(m)[[var_to_pull]] - qnorm(0.975)*se, 
              coef(m)[[var_to_pull]] + qnorm(0.975)*se)) # RR: 95% CI with robust SE
    } 
  if (family == "gaussian") {
    beta <- coef(m)[[var_to_pull]] # Relative Difference
    se <- summary(m)$coef[var_to_pull,2] # SE
    ci <- c(coef(m)[[var_to_pull]] - qnorm(0.975)*se, 
                coef(m)[[var_to_pull]] + qnorm(0.975)*se) # RD: 95% CI with robust SE
  } 
  p <- summary(m)$coef[var_to_pull,4]
  
  # Run adjusted GEE and get coefs
  if (cvd_var != "smokingR") { 
    vars <- c(var_of_interest, adj_vars)
    formula <- make_formula(vars)
  }
  if (cvd_var == "smokingR") {
    vars <- c(var_of_interest, adj_vars[adj_vars != "smokingR"])
    formula <- make_formula(vars)
  }
  adj.m <- run_gee(formula, family, df)
  if (family == "poisson") {
    adj.beta <- exp(coef(adj.m)[[var_to_pull]]) # Relative Risk 
    adj.se <- summary(adj.m)$coef[var_to_pull,2] # SE
    adj.ci <- exp(c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*se, 
                coef(adj.m)[[var_to_pull]] + qnorm(0.975)*se)) # RR: 95% CI with robust SE
  } 
  if (family == "gaussian") {
    adj.beta <- coef(adj.m)[[var_to_pull]] # Relative Difference
    adj.se <- summary(adj.m)$coef[var_to_pull,2] # SE
    adj.ci <- c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*se, 
            coef(adj.m)[[var_to_pull]] + qnorm(0.975)*se) # RD: 95% CI with robust SE
  } 
  adj.p <- summary(adj.m)$coef[var_to_pull,4]

  # Format row
  temp_row <- c(reg_vars[reg_vars$var_name == cvd_var, "var_label"], 
                paste0(n),
                format(round(beta, round_digits), nsmall = round_digits), 
                paste0("(",
                       format(round(ci, round_digits), nsmall = round_digits)[1], 
                       ",",
                       format(round(ci, round_digits), nsmall = round_digits)[2],
                       ")"),
                format(round(p, p_digits), nsmall = p_digits),
                format(round(adj.beta, round_digits), nsmall = round_digits), 
                paste0("(",
                       format(round(adj.ci, round_digits), nsmall = round_digits)[1], 
                       ",",
                       format(round(adj.ci, round_digits), nsmall = round_digits)[2],
                       ")"),
                format(round(adj.p, p_digits), nsmall = p_digits)
                )
  
  # Run extended GEE and get coefs [CONTROLLING FOR BASELINE CD4+]
  if (ext == T) {
  
    if (cvd_var != "smokingR") { 
      vars <- c(var_of_interest, adj_vars, ext_var)
      formula <- make_formula(vars)
    }
    if (cvd_var == "smokingR") {
      vars <- c(var_of_interest, adj_vars[adj_vars != "smokingR"], ext_var)
      formula <- make_formula(vars)
    }
    
    ext.m <- run_gee(formula, family, df)
    if (family == "poisson") {
      ext.beta <- exp(coef(ext.m)[[var_to_pull]]) # Relative Risk 
      ext.se <- summary(ext.m)$coef[var_to_pull,2] # SE
      ext.ci <- exp(c(coef(ext.m)[[var_to_pull]] - qnorm(0.975)*se, 
                  coef(ext.m)[[var_to_pull]] + qnorm(0.975)*se)) # RR: 95% CI with robust SE
    } 
    if (family == "gaussian") {
      ext.beta <- coef(ext.m)[[var_to_pull]] # Relative Difference
      ext.se <- summary(ext.m)$coef[var_to_pull,2] # SE
      ext.ci <- c(coef(ext.m)[[var_to_pull]] - qnorm(0.975)*se, 
              coef(ext.m)[[var_to_pull]] + qnorm(0.975)*se) # RD: 95% CI with robust SE
    } 
    ext.p <- summary(ext.m)$coef[var_to_pull,4]
  
    temp_row <- c(temp_row,
                  format(round(ext.beta, round_digits), nsmall = round_digits), 
                  paste0("(",
                         format(round(ext.ci, round_digits), nsmall = round_digits)[1], 
                         ",",
                         format(round(ext.ci, round_digits), nsmall = round_digits)[2],
                          ")"),
                  format(round(ext.p, p_digits), nsmall = p_digits)
    )
    
    }
  
    return(temp_row)
}

# ==============================================================================

# Regression results for ALL CVD vars- unadjusted and adjusted

get_regression_results <- function(df, var_of_interest, var_to_pull, family, ext = F) {
  
  if (family == "poisson") {
    
    reg_vars = poisson_reg_vars
    var_names <- unique(reg_vars$var_name)
    beta_label = "Relative Risk"
  }
  
  if (family == "gaussian") {
    
    reg_vars = lin_reg_vars
    var_names <- unique(reg_vars$var_name)
    beta_label = "Mean Difference"
  }
  
  # Loop through CVD vars and run GEE
  
  results_list <- list()
  
  for (v in unique(var_names)) {
    
    cvd_var <- v  
    
    temp_row <- get_gee_results (df = df, family = family, cvd_var = cvd_var, var_of_interest = var_of_interest, var_to_pull = var_to_pull, ext = ext)
   
    results_list[[v]] <- temp_row
    regressionDF <- as.data.frame(do.call(rbind, results_list))
  }
  
  # ADD COL NAMES --> FORMATTED WORD TABLES  

  if (ext == F) {
    
    names(regressionDF) =
      c("Variable", "n",
        beta_label, "95% CI", "p-value",
        paste0("Adj. ", beta_label), "Adj. 95% CI", "Adj. p-value"
      )
    
    ft <- flextable(regressionDF) %>% 
      colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
      add_header_row(values = c("", "", "Unadjusted Analysis", "Adjusted Analysis"), 
                     top = T, colwidths = c(1, 1, 3, 3)) %>% 
      width(width = 1.1) %>%
      width(j = 1, width = 2.5) %>%
      theme_vanilla() %>% 
      vline(j = c(1, 2, 5), border = fp_border_default(), part = "all") %>% 
      align(i = 1, align = "center", part = "header") %>% 
      fontsize(size = 10)
  }
  
  if (ext == T) {
    
    names(regressionDF) = c(
      "Variable", "n",
      beta_label, "95% CI", "p-value",
      paste0("Adj. ", beta_label), "Adj. 95% CI", "Adj. p-value",
      paste0("Ext. ", beta_label), "Ext. 95% CI", "Ext. p-value")
    
    ft <- flextable(regressionDF) %>% 
      colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
      add_header_row(values = c("", "", "Unadjusted Analysis", "Adjusted Analysis", "Extended Analysis"), 
                     top = T, colwidths = c(1, 1, 3, 3, 3)) %>% 
      width(width = 0.9) %>%
      width(j = 1, width = 1.3) %>%
      theme_vanilla() %>% 
      vline(j = c(1, 2, 5, 8), border = fp_border_default(), part = "all") %>% 
      align(i =1, align = "center", part = "header") %>% 
      fontsize(size = 10)
  }
      
  return(ft)    
}

#regressionDF <- get_regression_results(df = ncd_merged_subset, var_of_interest = "armR", var_to_pull = "armRCommunity follow-up", ext = F, family = "gaussian") 

# ==============================================================================

get_gee_strat_results <- function(df = df, 
                                  family, # poisson or gaussian
                                  cvd_var,
                                  var_of_interest,  # e.g. "exit_viral_load_suppressed"
                                  var_to_pull, # coefficient to pull, e.g. "exit_viral_load_suppressedTRUE"
                                  strat_var, # needs to be a binary var = e.g. "genderR", "armR"
                                  strat_var_to_pull, # e.g. "armRCommunity follow-up"
                                  ext = ext) { # is this the extended analysis? defaults to no
  
  # Get var labels
  if (family == "poisson") { reg_vars = poisson_reg_vars} 
  if (family == "gaussian") { reg_vars = lin_reg_vars}
  
  # Get formula
  make_formula <- function(vars) {as.formula(paste0(cvd_var, " ~ ", paste(vars, collapse = " + ")))}
  
  if (cvd_var != "smokingR") { 
    vars <- c(var_of_interest, adj_vars)
    formula <- make_formula(vars)
  }
  if (cvd_var == "smokingR") {
    vars <- c(var_of_interest, adj_vars[adj_vars != "smokingR"])
    formula <- make_formula(vars)
  }
  
  # Run GEE and get coefs - first level
  df1 <- df[df[[strat_var]] == levels(df[[strat_var]])[1], ]
  n.l1 <- df1 %>% select(c(cvd_var, vars, hhid)) %>% filter(complete.cases(.)) %>% nrow() # N obs
  adj.m <- run_gee(formula, family, df1)
  if (family == "poisson") {
    beta.l1 <- exp(coef(adj.m)[[var_to_pull]]) # Relative Risk 
    se <- summary(adj.m)$coef[var_to_pull,2] # SE
    ci.l1 <- exp(c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*se, 
                coef(adj.m)[[var_to_pull]] + qnorm(0.975)*se)) # RR: 95% CI with robust SE
  } 
  if (family == "gaussian") {
    beta.l1 <- coef(adj.m)[[var_to_pull]] # Relative Difference
    se <- summary(adj.m)$coef[var_to_pull,2] # SE
    ci.l1 <- c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*se, 
            coef(adj.m)[[var_to_pull]] + qnorm(0.975)*se) # RD: 95% CI with robust SE
  } 
  p.l1 <- summary(adj.m)$coef[var_to_pull,4]
  
  # Run GEE and get coefs - 2nd level
  df2 <- df[df[[strat_var]] == levels(df[[strat_var]])[2], ]
  n.l2 <- df2 %>% select(c(cvd_var, vars, hhid)) %>% filter(complete.cases(.)) %>% nrow() # N obs
  adj.m <- run_gee(formula, family, df2)
  if (family == "poisson") {
    beta.l2 <- exp(coef(adj.m)[[var_to_pull]]) # Relative Risk 
    se <- summary(adj.m)$coef[var_to_pull,2] # SE
    ci.l2 <- exp(c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*se, 
                   coef(adj.m)[[var_to_pull]] + qnorm(0.975)*se)) # RR: 95% CI with robust SE
  } 
  if (family == "gaussian") {
    beta.l2 <- coef(adj.m)[[var_to_pull]] # Relative Difference
    se <- summary(adj.m)$coef[var_to_pull,2] # SE
    ci.l2 <- c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*se, 
               coef(adj.m)[[var_to_pull]] + qnorm(0.975)*se) # RD: 95% CI with robust SE
  } 
  p.l2 <- summary(adj.m)$coef[var_to_pull,4]
  
  # Run GEE with interaction term and get p-value only
  
  if (cvd_var != "smokingR") { 
    vars <- c(var_of_interest, adj_vars, paste0(var_of_interest, "*", strat_var))
    formula <- make_formula(vars)
  }
  if (cvd_var == "smokingR") {
    vars <- c(var_of_interest, adj_vars[adj_vars != "smokingR"], paste0(var_of_interest, "*", strat_var))
    formula <- make_formula(vars)
  }
  
  int.m <- run_gee(formula, family, df)
  int.p <- summary(int.m)$coef[paste0(var_to_pull, ":", strat_var_to_pull),4]
  
  # Format row
  temp_row <- c(reg_vars[reg_vars$var_name == cvd_var, "var_label"], 
                paste0(n.l1),
                format(round(beta.l1, round_digits), nsmall = round_digits), 
                paste0("(",
                       format(round(ci.l1, round_digits), nsmall = round_digits)[1], 
                       ",",
                       format(round(ci.l1, round_digits), nsmall = round_digits)[2],
                       ")"),
                format(round(p.l1, p_digits), nsmall = p_digits),
                paste0(n.l2),
                format(round(beta.l2, round_digits), nsmall = round_digits), 
                paste0("(",
                       format(round(ci.l2, round_digits), nsmall = round_digits)[1], 
                       ",",
                       format(round(ci.l2, round_digits), nsmall = round_digits)[2],
                       ")"),
                format(round(p.l2, p_digits), nsmall = p_digits),
                as.character(round(int.p, p_digits))
  )
  
  return(temp_row)
}


# ==============================================================================

# Regression results for ALL CVD vars- stratified analysis

get_reg_strat_results <- function(df, 
                                   family, 
                                   cvd_var,
                                   var_of_interest, var_to_pull, 
                                   strat_var, strat_var_to_pull,
                                   ext = F) {
  
  if (family == "poisson") {
    
    reg_vars = poisson_reg_vars
    var_names <- unique(reg_vars$var_name)
    beta_label = "Relative Risk"
  }
  
  if (family == "gaussian") {
    
    reg_vars = lin_reg_vars
    var_names <- unique(reg_vars$var_name)
    beta_label = "Mean Difference"
  }
  
  # Loop through CVD vars and run GEE
  
  results_list <- list()
  
  for (v in unique(var_names)) {
    
    cvd_var <- v  
    
    temp_row <- get_gee_strat_results (df = df, family = family, cvd_var = cvd_var, 
                                       var_of_interest = var_of_interest, var_to_pull = var_to_pull, 
                                       strat_var = strat_var, strat_var_to_pull = strat_var_to_pull,
                                       ext = ext)
    
    results_list[[v]] <- temp_row
    regressionDF <- as.data.frame(do.call(rbind, results_list))
  }
  
  # Col names
  names(regressionDF) = c(
      "Variable", 
      "n", paste0("Adj. ", beta_label), "Adj. 95% CI", "Adj. p-value",
      "n ", paste0("Adj. ", beta_label, " "), "Adj. 95% CI ", "Adj. p-value ",
      "Adj. p-value  "
      )
  
  # Formatted Word Table
  ft <- flextable(regressionDF) %>%
    colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
    add_header_row(values = c("", levels(df[[strat_var]])[1], levels(df[[strat_var]])[2], "Interaction"),
                   top = T, colwidths = c(1, 4, 4, 1)) %>%
    width(width = 1.1) %>%
    width(j = 1, width = 1.75) %>%
    theme_vanilla()  %>%
    vline(j = c(1, 5, 9), border = fp_border_default(), part = "all") %>%
    align(i =1, align = "center", part = "header")  %>%
    fontsize(size = 10)

    return(ft) 
}

# regressionDF <- get_reg_strat_results(df <- ncd_merged_subset %>% filter(!is.na(exit_viral_load_suppressed)),
#                                       var_of_interest = "exit_viral_load_suppressed",
#                                       var_to_pull = "exit_viral_load_suppressedTRUE",
#                                       cvd_var = poisson_reg_vars$var_name[1],
#                                       family = "poisson",
#                                       strat_var = "armR",
#                                       strat_var_to_pull = "armRCommunity follow-up",
#                                       ext = T)

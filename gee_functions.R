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
unadj_vars <- c("armR", "site")
adj_vars <- c(unadj_vars, "age_catR", "genderR", "smokingR", "educationR") # UPDATE HERE FOR ANY CHANGES TO CONTROL VARS

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

# Get results for single CVD var!

#  ext = T
#  df <- ncd_merged_subset
#  var_to_pull = "armRCommunity follow-up"
# cvd_var = poisson_reg_vars$var_name[1]
# family = "poisson"
#  
get_gee_results <- function(df = df, 
                            family, # poisson or gaussian
                            cvd_var,
                            var_to_pull, # coefficient to pull, e.g. "armRCommunity follow-up"
                            ext = ext) { # is this the extended analysis? defaults to no
  
  # Get var labels
  if (family == "poisson") { reg_vars = poisson_reg_vars} 
  if (family == "gaussian") { reg_vars = lin_reg_vars}
  
  # Run GEE and get coeffs - unadjusted
  vars <- unadj_vars
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
    vars <- adj_vars
    formula <- make_formula(vars)
  }
  if (cvd_var == "smokingR") {
    vars <- adj_vars[adj_vars != "smokingR"]
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
      vars <- c(adj_vars, "baseline_cd4")
      formula <- make_formula(vars)
    }
    if (cvd_var == "smokingR") {
      vars <- c(adj_vars[adj_vars != "smokingR"], "baseline_cd4")
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

# Make Word table with Poisson and Linear regressions

format_regression_table <- function(df, ext = F) {
  
  if (ext == F) {
    
    ft <- flextable(df) %>% 
      colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
      add_header_row(values = c("", "", "Unadjusted Analysis", "Adjusted Analysis"), 
                     top = T, colwidths = c(1, 1, 3, 3)) %>% 
      width(width = 1.1) %>%
      width(j = 1, width = 2.5) %>%
      theme_vanilla() %>% 
      vline(j = c(1, 2, 5), border = fp_border_default(), part = "all") %>% 
      align(i =1, align = "center", part = "header") 
  
  }
  
  if (ext == T) {
    
    ft <- flextable(df) %>% 
      colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
      add_header_row(values = c("", "", "Unadjusted Analysis", "Adjusted Analysis", "Extended Analysis"), 
                     top = T, colwidths = c(1, 1, 3, 3, 3)) %>% 
      width(width = 1.1) %>%
      width(j = 1, width = 2.5) %>%
      theme_vanilla() %>% 
      vline(j = c(1, 2, 5, 8), border = fp_border_default(), part = "all") %>% 
      align(i =1, align = "center", part = "header") 
    
  }
  
  return(ft)
}

# ==============================================================================

get_regression_results <- function(df, var_to_pull, family, ext = F) {
  
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
  
  # LIST OF RESULTS FOR SET OF CVD VARS
  
  results_list <- list()
  
  for (v in unique(var_names)) {
    
      cvd_var <- v  
      make_formula <- function(vars) {as.formula(paste0(cvd_var, " ~ ", paste(vars, collapse = " + ")))}
      
      # Run GEE and get coeffs - unadjusted
      vars <- unadj_vars
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
        vars <- adj_vars
        formula <- make_formula(vars)
      }
      if (cvd_var == "smokingR") {
        vars <- adj_vars[adj_vars != "smokingR"]
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
          vars <- c(adj_vars, "baseline_cd4")
          formula <- make_formula(vars)
        }
        if (cvd_var == "smokingR") {
          vars <- c(adj_vars[adj_vars != "smokingR"], "baseline_cd4")
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
                      format(round(ext.p, p_digits), nsmall = p_digits))
      }
        
    # Store temp_row in results_list
    results_list[[v]] <- temp_row
    # Combine the results from the list into a single data frame
    regressionDF <- as.data.frame(do.call(rbind, results_list))
        
  }
  
  # Names
  if (ext == F) {
    
    names(regressionDF) =
      c("Variable", "n",
        beta_label, "95% CI", "p-value",
        paste0("Adj. ", beta_label), "Adj. 95% CI", "Adj. p-value"
      )
  }
  
  if (ext == T) {
    
    names(regressionDF) = c(
      "Variable", "n",
      beta_label, "95% CI", "p-value",
      paste0("Adj. ", beta_label), "Adj. 95% CI", "Adj. p-value",
      paste0("Ext. ", beta_label), "Ext. 95% CI", "Ext. p-value")
  }
      
  return(regressionDF)    
}

regressionDF <- get_regression_results(df = ncd_merged_subset, var_to_pull = "armRCommunity follow-up", ext = F, family = "gaussian") 

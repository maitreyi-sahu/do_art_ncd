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
make_formula <- function(vars) {as.formula(paste0(v, " ~ ", paste(vars, collapse = " + ")))}
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
    "lipid_result"),
  
  var_label = c(
    "Systolic BP (mmHg)",
    "Diastolic BP (mmHg)",
    "BMI (kg/m^2)",
    "Hemoglobin A1c (%)",
    "Total cholesterol (mg/dL)") 
)

# ==============================================================================

# Get results for all vars!

get_empty_reg_df <- function(ext = F) {
  
  regressionDF <- data.frame(
    
    variable = character(),
    beta = numeric(),
    ci = character(),
    p = numeric(),
    n = integer(),
    adj.beta = character(),
    adj.ci = character(), 
    adj.p = character(),
    adj.n = integer()
  )
  
  if (ext == T) {
    
    regressionDF$ext.beta = character()
    regressionDF$ext.ci = character() 
    regressionDF$ext.p = character()
    regressionDF$ext.n = integer()
  }
  
  return(regressionDF)
}

get_gee_results <- function(df = ncd_merged_subset, 
                            family, # poisson or gaussian
                            cvd_var,
                            var_to_pull, # coefficient to pull, e.g. "armRCommunity follow-up"
                            ext = F) { # is this the extended analysis? defaults to no

  # Get var labels
  if (family == "poisson") { reg_vars = poisson_reg_vars} 
  else if (family == "gaussian") { reg_vars = lin_reg_vars}
  
  # Run GEE and get coeffs - unadjusted
  vars <- unadj_vars
  formula <- make_formula(vars)
  m <- run_gee(formula, family, df)
  beta <- exp(coef(m)[[var_to_pull]]) # Relative Risk / Relative Difference
  se <- summary(m)$coef[var_to_pull,2] # SE
  ci <- exp(c(coef(m)[[var_to_pull]] - qnorm(0.975)*se, 
              coef(m)[[var_to_pull]] + qnorm(0.975)*se)) # RR: 95% CI with robust SE
  p <- summary(m)$coef[var_to_pull,4]
  n <- df %>% select(c(v, vars)) %>% filter(complete.cases(.)) %>% nrow() # N obs
  
  # Run adjusted GEE and get coefs
  if (v != "smokingR") { 
    vars <- adj_vars
    formula <- make_formula(vars)
  }
  if (v == "smokingR") {
    vars <- adj_vars[adj_vars != "smokingR"]
    formula <- make_formula(vars)
  }
  adj.m <- run_gee(formula, family, df)
  adj.beta <- exp(coef(adj.m)[[var_to_pull]]) # Relative Risk / Relative Difference
  adj.se <- summary(adj.m)$coef[var_to_pull,2] # SE
  adj.ci <- exp(c(coef(adj.m)[[var_to_pull]] - qnorm(0.975)*adj.se, 
              coef(adj.m)[[var_to_pull]] + qnorm(0.975)*adj.se)) # RR: 95% CI with robust SE
  adj.p <- summary(adj.m)$coef[var_to_pull,4]
  adj.n <- df %>% select(c(v, vars)) %>% filter(complete.cases(.)) %>% nrow() # N obs

  # Format row
  temp_row <- c(reg_vars[reg_vars$var_name == v, "var_label"], 
                format(round(beta, round_digits), nsmall = round_digits), 
                paste0("(",
                       format(round(ci, round_digits), nsmall = round_digits)[1], 
                       ",",
                       format(round(ci, round_digits), nsmall = round_digits)[2],
                       ")"),
                format(round(p, p_digits), nsmall = p_digits),
                paste0(n),
                format(round(adj.beta, round_digits), nsmall = round_digits), 
                paste0("(",
                       format(round(adj.ci, round_digits), nsmall = round_digits)[1], 
                       ",",
                       format(round(adj.ci, round_digits), nsmall = round_digits)[2],
                       ")"),
                format(round(adj.p, p_digits), nsmall = p_digits),
                paste0(adj.n)
                )
  
  # Column names
  if (family == "poisson") { beta_label = "Relative Risk"} 
  else if (family == "gaussian") { beta_label = "Mean Difference"}
  
  names(regressionDF) =
    c("Variable", 
      beta_label, "95% CI", "p-value", "n",
      paste0("Adj. ", beta_label), "Adj. 95% CI", "Adj. p-value", "Adj. n"
    )
  
  # Run extended GEE and get coefs [CONTROLLING FOR BASELINE CD4+]
  if (ext == T) {
  
    if (v != "smokingR") { 
      vars <- c(adj_vars, "baseline_cd4")
      formula <- make_formula(vars)
    }
    if (v == "smokingR") {
      vars <- c(adj_vars[adj_vars != "smokingR"], "baseline_cd4")
      formula <- make_formula(vars)
    }
    ext.m <- run_gee(formula, family, df)
    ext.beta <- exp(coef(ext.m)[[var_to_pull]]) # Relative Risk / Relative Difference
    ext.se <- summary(ext.m)$coef[var_to_pull,2] # SE
    ext.ci <- exp(c(coef(ext.m)[[var_to_pull]] - qnorm(0.975)*m.se, 
                    coef(ext.m)[[var_to_pull]] + qnorm(0.975)*m.se)) # RR: 95% CI with robust SE
    ext.p <- summary(ext.m)$coef[var_to_pull,4]
    ext.n <- df %>% filter(!is.na(get(v))) %>% filter_all(any_vars(!is.na(.)), vars = vars) %>% nrow() # number of observations
  
    temp_row <- c(temp_row,
                  format(round(ext.rr, round_digits), nsmall = round_digits), 
                  paste0("(",
                         format(round(ext.rr_ci, round_digits), nsmall = round_digits)[1], 
                         ",",
                         format(round(ext.rr_ci, round_digits), nsmall = round_digits)[2],
                          ")"),
                  format(round(ext.p, p_digits), nsmall = p_digits),
                  paste0(ext.n)
    )
    
    }
  
    return(temp_row)
}

# 
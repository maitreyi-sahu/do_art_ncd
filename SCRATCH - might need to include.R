
### Page Break

# Table 4: Cardiovascular risk for virally suppressed versus not, by treatment arm (N = 983)

** Note, we drop 27 participants without endline viral load available.

&nbsp;

```{r, echo = F, warning = F, message = F, ft.split = T, ft.keepnext = F, results = 'asis'}

vs_subset <- ncd_merged_subset %>% 
  filter(!is.na(exit_viral_load_result))

for (s in unique(vs_subset$armR)) {
  
  vs_arm_subset <- vs_subset %>% filter(armR == s)  
  
  N = nrow(vs_arm_subset)
  
  # POISSON REGRESSIONS -------------------------------------------------------------
  
  reg_vars <- data.frame(
    
    var_name = c(
      "bp_hi_exit", 
      "overwt_exit", 
      "a1c_hi_exit",
      "total_cholesterol_hi_exit",
      "smokingR"),
    
    var_label = c(
      "Elevated BP",
      "Overweight (BMI >= 25)",
      "Elevated blood sugar (a1c >= 5.7)",
      "Elevated lipids (total cholesterol >= 200)",
      "Current smoker [baseline]"))
  
  
  regressionDF <- data.frame(
    variable = character(),
    rr = numeric(),
    rr_ci = character(),
    rr_p = numeric(),
    adj.rr =  numeric(),
    adj.rr_ci = character(), 
    adj.rr_p = character()
  )
  
  for (v in unique(reg_vars$var_name)) {
    
    # Unadjusted
    
    m <- geeglm(get(v) ~ exit_viral_load_suppressed, 
                family = poisson, data = vs_arm_subset,
                id = hhid, corstr = "independence") 
    
    m.rr <- exp(coef(m)[["exit_viral_load_suppressedTRUE"]]) # Relative Risk
    m.se <- summary(m)$coef["exit_viral_load_suppressedTRUE",2] # SE
    m.rr_ci <- exp(c(coef(m)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*m.se, 
                     coef(m)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*m.se)) # RR: 95% CI 
    m.p <- summary(m)$coef["exit_viral_load_suppressedTRUE",4]
    
    # Adjusted for site, age, gender, smoking 
    
    if (v != "smokingR") {
      
      adj.m <- geeglm(get(v) ~ exit_viral_load_suppressed + site + age_cat + genderR + smokingR, 
                      family = poisson, data = vs_arm_subset,
                      id = hhid, corstr = "independence") 
    }
    
    if (v == "smokingR") {
      
      adj.m <- geeglm(get(v) ~ exit_viral_load_suppressed + site + age_cat + genderR, 
                      family = poisson, data = vs_arm_subset,
                      id = hhid, corstr = "independence") 
    }
    
    adj.rr <- exp(coef(adj.m)[["exit_viral_load_suppressedTRUE"]]) # Relative Risk
    adj.se <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE",2] # SE
    adj.rr_ci <- exp(c(coef(adj.m)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*adj.se, 
                       coef(adj.m)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*adj.se)) # RR: 95% CI 
    adj.p <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE",4]
    
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
                  as.numeric(round(adj.p, p_digits)))
    
    regressionDF <- rbind(regressionDF, temp_row)
    
  }
  
  names(regressionDF) =
    c("Variable",
      "Relative Risk",
      "95% CI",
      "p-value",
      "Adj. Relative Risk",
      "Adj. 95% CI",
      "Adj. p-value")
  
  # Word Doc output
  ft <- flextable(regressionDF) %>% 
    colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
    add_header_row(values = c(paste0(s, " (n = " , N, ")"), "Baseline Analysis", "Adjusted Analysis"), 
                   top = T, colwidths = c(1, 3, 3)) %>% 
    width(width = 1.1) %>%
    width(j = 1, width = 2.5) %>%
    theme_vanilla() %>% 
    vline(j = c(1, 4), border = fp_border_default(), part = "all") %>% 
    align(i =1, align = "center", part = "header") 
  
  flextable_to_rmd(ft) 
  
  # LINEAR REGRESSIONS -----------------------------------------------------------------
  
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
    adj.rd_p = character() 
  )
  
  
  for (v in unique(lin_reg_vars$var_name)) {
    
    # Unadjusted 
    
    lm <- geeglm(get(v) ~ exit_viral_load_suppressed, 
                 family = gaussian, data = vs_arm_subset,
                 id = hhid, corstr = "independence") 
    rd <- coef(lm)[["exit_viral_load_suppressedTRUE"]] # Mean Difference
    lm.se <- summary(lm)$coef["exit_viral_load_suppressedTRUE",2] # SE
    rd_ci <- c(coef(lm)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*lm.se, 
               coef(lm)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*lm.se) # 95% CI
    rd_p <- summary(lm)$coef["exit_viral_load_suppressedTRUE",4]
    
    # Adjusted for site, age, gender, smoking 
    
    adj.lm <- geeglm(get(v) ~ exit_viral_load_suppressed + site + age_cat + genderR + smokingR, # NOTE THIS SMOKING VARIABLE IS DIFFERENT
                     family = gaussian, data = vs_arm_subset,
                     id = hhid, corstr = "independence") 
    adj.rd <- coef(adj.lm)[["exit_viral_load_suppressedTRUE"]] # Mean Difference
    adj.lm.se <- summary(adj.lm)$coef["exit_viral_load_suppressedTRUE",2] # SE
    adj.rd_ci <- c(coef(adj.lm)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*adj.lm.se, 
                   coef(adj.lm)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*adj.lm.se) # 95% CI
    adj.rd_p <- summary(adj.lm)$coef["exit_viral_load_suppressedTRUE",4]
    
    
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
                  as.numeric(round(adj.rd_p, p_digits)))
    
    regressionDF <- rbind(regressionDF, temp_row)
    
  }
  
  names(regressionDF) =
    c("Variable",
      "Mean Difference",
      "95% CI",
      "p-value",
      "Adj. Mean Difference",
      "Adj. 95% CI",
      "Adj. p-value")
  
  # Word Doc output
  ft <- flextable(regressionDF) %>% 
    colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
    add_header_row(values = c(paste0(s, " (n = " , N, ")"), "Baseline Analysis", "Adjusted Analysis"), 
                   top = T, colwidths = c(1, 3, 3)) %>% 
    width(width = 1.1) %>%
    width(j = 1, width = 2.5) %>%
    theme_vanilla()  %>% 
    vline(j = c(1, 4), border = fp_border_default(), part = "all") %>% 
    align(i =1, align = "center", part = "header") 
  
  flextable_to_rmd(ft)  
  
}

```
+ NEED TO DROP CHOLESTEROL FROM THIS TABLE?
  + FIX GLITCH WITH HEADINGS on pg 2
+ This analysis includes smoking as a binary variable because otherwise getting error

### Page Break

Table 5: Cardiovascular risk among VS versus not VS at endline, adjusting for treatment arm (N = 983)

```{r, echo = F, warning = F, message = F, ft.split = T, ft.keepnext = F, results = 'asis'}

# DROP NA for GEE to run

vs_subset <- ncd_merged_subset %>% 
  filter(!is.na(exit_viral_load_result))

# POISSON REGRESSIONS --------------------------------------------------------------------------

reg_vars <- data.frame(
  
  var_name = c(
    "bp_hi_exit", 
    "overwt_exit", 
    "a1c_hi_exit",
    "total_cholesterol_hi_exit",
    "smokingR"
  ),
  
  var_label = c(
    "Elevated BP",
    "Overweight (BMI >= 25)",
    "Elevated blood sugar (a1c >= 5.7)",
    "Elevated lipids (total cholesterol >= 200)",
    "Current smoker [baseline]"
  ))

regressionDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_ci = character(),
  rr_p = numeric(),
  adj.rr = character(),
  adj.rr_ci = character(), 
  adj.rr_p = character() 
)


plotDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_ci = character(),
  rr_p = numeric(),
  adj.rr = character(),
  adj.rr_ci = character(), 
  adj.rr_p = character()
)


for (v in unique(reg_vars$var_name)) {
  
  # Adjusted for site only
  m <- geeglm(get(v) ~ site + exit_viral_load_suppressed + armR, 
              family = poisson, data = vs_subset,
              id = hhid, corstr = "independence") 
  m.rr <- exp(coef(m)[["exit_viral_load_suppressedTRUE"]]) # Relative Risk
  m.se <- summary(m)$coef["exit_viral_load_suppressedTRUE",2] # SE
  m.rr_ci <- exp(c(coef(m)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*m.se, 
                   coef(m)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*m.se)) # RR: 95% CI 
  m.p <- summary(m)$coef["exit_viral_load_suppressedTRUE",4]
  
  # Adjusted for site, age, gender, smoking 
  
  if (v != "smokingR") {
    
    adj.m <- geeglm(get(v) ~ site + age_cat + genderR + smoking_status + exit_viral_load_suppressed + armR, 
                    family = poisson, data = vs_subset,
                    id = hhid, corstr = "independence") 
    
    adj.rr <- exp(coef(adj.m)[["exit_viral_load_suppressedTRUE"]]) # Relative Risk
    adj.se <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE",2] # SE
    adj.rr_ci <- exp(c(coef(adj.m)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*adj.se, 
                       coef(adj.m)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*adj.se)) # RR: 95% CI with robust SE
    adj.p <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE",4]
    
  }
  
  if (v == "smokingR") {
    
    adj.m <- geeglm(get(v) ~ site + age_cat + genderR + exit_viral_load_suppressed + armR, 
                    family = poisson, data = vs_subset,
                    id = hhid, corstr = "independence") 
    
    adj.rr <- exp(coef(adj.m)[["exit_viral_load_suppressedTRUE"]]) # Relative Risk
    adj.se <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE",2] # SE
    adj.rr_ci <- exp(c(coef(adj.m)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*adj.se, 
                       coef(adj.m)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*adj.se)) # RR: 95% CI with robust SE
    adj.p <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE",4]
    
  }
  
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
                as.numeric(round(adj.p, p_digits)))
  
  regressionDF <- rbind(regressionDF, temp_row)
  
  temp_row_plot <- c(reg_vars[reg_vars$var_name == v, "var_label"], 
                     as.numeric(m.rr), 
                     as.numeric(m.rr_ci[1]),
                     as.numeric(m.rr_ci[2]),
                     as.numeric(round(m.p, p_digits)),
                     as.numeric(adj.rr), 
                     as.numeric(adj.rr_ci[1]),
                     as.numeric(adj.rr_ci[2]),
                     as.numeric(round(adj.p, p_digits))
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
    "Adj. p-value")

names(plotDF) =
  c("Variable",
    "main.estimate",
    "main.lower",
    "main.upper",
    "main.p",
    "adj.estimate",
    "adj.lower",
    "adj.upper",
    "adj.p")

plotDF <-  plotDF %>% 
  gather(label, value, 2:9) %>% 
  separate(label, into = c("model", "estimate")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = estimate, value = value) %>% 
  mutate(country = "South Africa",
         outcome = "relative_risk",
         version = "vs")

results_list[[length(results_list) + 1]] <- plotDF
names(results_list)[length(results_list)] <- "vs.rr.south_africa"

# Word Doc output
ft <- flextable(regressionDF) %>% 
  colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
  add_header_row(values = c("", "Baseline Analysis", "Adjusted Analysis"), 
                 top = T, colwidths = c(1, 3, 3)) %>% 
  width(width = 1.1) %>%
  width(j = 1, width = 2.5) %>%
  theme_vanilla() %>% 
  vline(j = c(1, 4), border = fp_border_default(), part = "all") %>% 
  align(i =1, align = "center", part = "header") 

ft

# LINEAR REGRESSIONS -------------------------------------------------------------------------------

lin_reg_vars <- data.frame(
  
  var_name = c(
    "sbp_mean_exit",
    "dbp_mean_exit", 
    "bmi", 
    "hemoglobin_a1c_result",
    "lipid_result"),
  
  var_label = c(
    "Systolic blood pressure (mmHg)",
    "Diastolic blood pressure (mmHg)",
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
  adj.rd_p = character() 
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
  adj.rr_p = numeric()
)


for (v in unique(lin_reg_vars$var_name)) {
  
  # Adjusted for site only
  lm <- geeglm(get(v) ~ site + exit_viral_load_suppressed + armR, 
               family = gaussian, data = vs_subset,
               id = hhid, corstr = "independence") 
  rd <- coef(lm)[["exit_viral_load_suppressedTRUE"]] # Mean Difference
  lm.se <- summary(lm)$coef["exit_viral_load_suppressedTRUE",2] # SE
  rd_ci <- c(coef(lm)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*lm.se, 
             coef(lm)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*lm.se) # 95% CI
  rd_p <- summary(lm)$coef["exit_viral_load_suppressedTRUE",4]
  
  # Adjusted for site, age, gender, smoking 
  
  adj.lm <- geeglm(get(v) ~  site + age_cat + genderR + smoking_status + exit_viral_load_suppressed + armR, 
                   family = gaussian, data = vs_subset,
                   id = hhid, corstr = "independence") 
  adj.rd <- coef(adj.lm)[["exit_viral_load_suppressedTRUE"]] # Mean Difference
  adj.lm.se <- summary(adj.lm)$coef["exit_viral_load_suppressedTRUE",2] # SE
  adj.rd_ci <- c(coef(adj.lm)[["exit_viral_load_suppressedTRUE"]] - qnorm(0.975)*adj.lm.se, 
                 coef(adj.lm)[["exit_viral_load_suppressedTRUE"]] + qnorm(0.975)*adj.lm.se) # 95% CI
  adj.rd_p <- summary(adj.lm)$coef["exit_viral_load_suppressedTRUE",4]
  
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
                as.numeric(round(adj.rd_p, p_digits)))
  
  regressionDF <- rbind(regressionDF, temp_row)
  
  
  temp_row_plot <- c(lin_reg_vars[lin_reg_vars$var_name == v, "var_label"], 
                     as.numeric(rd), 
                     as.numeric(rd_ci[1]),
                     as.numeric(rd_ci[2]),
                     as.numeric(round(rd_p, p_digits)),
                     as.numeric(adj.rd), 
                     as.numeric(adj.rd_ci[1]),
                     as.numeric(adj.rd_ci[2]),
                     as.numeric(round(adj.rd_p, p_digits))
  )
  
  plotDF <- rbind(plotDF, temp_row_plot)
}

names(regressionDF) =
  c("Variable",
    "Mean Difference",
    "95% CI",
    "p-value",
    "Adj. Mean Difference",
    "Adj. 95% CI",
    "Adj. p-value")


names(plotDF) =
  c("Variable",
    "main.estimate",
    "main.lower",
    "main.upper",
    "main.p",
    "adj.estimate",
    "adj.lower",
    "adj.upper",
    "adj.p")

plotDF <-  plotDF %>% 
  gather(label, value, 2:9) %>% 
  separate(label, into = c("model", "estimate")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = estimate, value = value) %>% 
  mutate(country = "South Africa",
         outcome = "mean_difference",
         version = "vs")

results_list[[length(results_list) + 1]] <- plotDF
names(results_list)[length(results_list)] <- "vs.rd.south_africa"


# Word Doc output
ft <- flextable(regressionDF) %>% 
  colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
  add_header_row(values = c("", "Baseline Analysis", "Adjusted Analysis"), 
                 top = T, colwidths = c(1, 3, 3)) %>% 
  width(width = 1.1) %>%
  width(j = 1, width = 2.5) %>%
  theme_vanilla()  %>% 
  vline(j = c(1, 4), border = fp_border_default(), part = "all") %>% 
  align(i =1, align = "center", part = "header") 

ft
```

### Page Break

## Table 6: Coefficients for interaction term for poisson + linear regressions with VS as predictor, with interaction VS * arm

** Note, we drop 27 participants without endline viral load available.

```{r, echo = F, warning = F, message = F, ft.split = T, ft.keepnext = F, results = 'asis'}

# DROP NA for GEE to run

vs_subset <- ncd_merged_subset %>% 
  filter(!is.na(exit_viral_load_result))

# POISSON REGRESSIONS --------------------------------------------------------------------------

reg_vars <- data.frame(
  
  var_name = c(
    "bp_hi_exit", 
    "overwt_exit", 
    "a1c_hi_exit",
    "total_cholesterol_hi_exit",
    "smokingR"
  ),
  
  var_label = c(
    "Elevated BP",
    "Overweight (BMI >= 25)",
    "Elevated blood sugar (a1c >= 5.7)",
    "Elevated lipids (total cholesterol >= 200)",
    "Current smoker [baseline]"
  ))

regressionDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_ci = character(),
  rr_p = numeric(),
  adj.rr = character(),
  adj.rr_ci = character(), 
  adj.rr_p = character() 
)


plotDF <- data.frame(
  variable = character(),
  rr = numeric(),
  rr_ci = character(),
  rr_p = numeric(),
  adj.rr = character(),
  adj.rr_ci = character(), 
  adj.rr_p = character()
)


for (v in unique(reg_vars$var_name)) {
  
  # Adjusted for site only
  m <- geeglm(get(v) ~ site + exit_viral_load_suppressed*armR, 
              family = poisson, data = vs_subset,
              id = hhid, corstr = "independence") 
  m.rr <- exp(coef(m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]]) # Relative Risk
  m.se <- summary(m)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",2] # SE
  m.rr_ci <- exp(c(coef(m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] - qnorm(0.975)*m.se, 
                   coef(m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] + qnorm(0.975)*m.se)) # RR: 95% CI 
  m.p <- summary(m)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",4]
  
  # Adjusted for site, age, gender, smoking 
  
  if (v != "smokingR") {
    
    adj.m <- geeglm(get(v) ~ site + age_cat + genderR + smoking_status + exit_viral_load_suppressed*armR, 
                    family = poisson, data = vs_subset,
                    id = hhid, corstr = "independence") 
    
    adj.rr <- exp(coef(adj.m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]]) # Relative Risk
    adj.se <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",2] # SE
    adj.rr_ci <- exp(c(coef(adj.m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] - qnorm(0.975)*adj.se, 
                       coef(adj.m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] + qnorm(0.975)*adj.se)) # RR: 95% CI with robust SE
    adj.p <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",4]
    
  }
  
  if (v == "smokingR") {
    
    adj.m <- geeglm(get(v) ~ site + age_cat + genderR + exit_viral_load_suppressed*armR, 
                    family = poisson, data = vs_subset,
                    id = hhid, corstr = "independence") 
    
    adj.rr <- exp(coef(adj.m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]]) # Relative Risk
    adj.se <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",2] # SE
    adj.rr_ci <- exp(c(coef(adj.m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] - qnorm(0.975)*adj.se, 
                       coef(adj.m)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] + qnorm(0.975)*adj.se)) # RR: 95% CI with robust SE
    adj.p <- summary(adj.m)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",4]
    
  }
  
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
                as.numeric(round(adj.p, p_digits)))
  
  regressionDF <- rbind(regressionDF, temp_row)
  
  temp_row_plot <- c(reg_vars[reg_vars$var_name == v, "var_label"], 
                     as.numeric(m.rr), 
                     as.numeric(m.rr_ci[1]),
                     as.numeric(m.rr_ci[2]),
                     as.numeric(round(m.p, p_digits)),
                     as.numeric(adj.rr), 
                     as.numeric(adj.rr_ci[1]),
                     as.numeric(adj.rr_ci[2]),
                     as.numeric(round(adj.p, p_digits))
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
    "Adj. p-value")

names(plotDF) =
  c("Variable",
    "main.estimate",
    "main.lower",
    "main.upper",
    "main.p",
    "adj.estimate",
    "adj.lower",
    "adj.upper",
    "adj.p")

plotDF <-  plotDF %>% 
  gather(label, value, 2:9) %>% 
  separate(label, into = c("model", "estimate")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = estimate, value = value) %>% 
  mutate(country = "South Africa",
         outcome = "relative_risk",
         version = "vs")

results_list[[length(results_list) + 1]] <- plotDF
names(results_list)[length(results_list)] <- "vs.rr.south_africa"

# Word Doc output
ft <- flextable(regressionDF) %>% 
  colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
  add_header_row(values = c("", "Baseline Analysis", "Adjusted Analysis"), 
                 top = T, colwidths = c(1, 3, 3)) %>% 
  width(width = 1.1) %>%
  width(j = 1, width = 2.5) %>%
  theme_vanilla() %>% 
  vline(j = c(1, 4), border = fp_border_default(), part = "all") %>% 
  align(i =1, align = "center", part = "header") 

ft

# LINEAR REGRESSIONS -------------------------------------------------------------------------------

lin_reg_vars <- data.frame(
  
  var_name = c(
    "sbp_mean_exit",
    "dbp_mean_exit", 
    "bmi", 
    "hemoglobin_a1c_result",
    "lipid_result"),
  
  var_label = c(
    "Systolic blood pressure (mmHg)",
    "Diastolic blood pressure (mmHg)",
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
  adj.rd_p = character() 
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
  adj.rr_p = numeric()
)


for (v in unique(lin_reg_vars$var_name)) {
  
  # Adjusted for site only
  lm <- geeglm(get(v) ~ site + exit_viral_load_suppressed*armR, 
               family = gaussian, data = vs_subset,
               id = hhid, corstr = "independence") 
  rd <- coef(lm)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] # Mean Difference
  lm.se <- summary(lm)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",2] # SE
  rd_ci <- c(coef(lm)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] - qnorm(0.975)*lm.se, 
             coef(lm)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] + qnorm(0.975)*lm.se) # 95% CI
  rd_p <- summary(lm)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",4]
  
  # Adjusted for site, age, gender, smoking 
  
  adj.lm <- geeglm(get(v) ~  site + age_cat + genderR + smoking_status + exit_viral_load_suppressed*armR, 
                   family = gaussian, data = vs_subset,
                   id = hhid, corstr = "independence") 
  adj.rd <- coef(adj.lm)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] # Mean Difference
  adj.lm.se <- summary(adj.lm)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",2] # SE
  adj.rd_ci <- c(coef(adj.lm)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] - qnorm(0.975)*adj.lm.se, 
                 coef(adj.lm)[["exit_viral_load_suppressedTRUE:armRCommunity follow-up"]] + qnorm(0.975)*adj.lm.se) # 95% CI
  adj.rd_p <- summary(adj.lm)$coef["exit_viral_load_suppressedTRUE:armRCommunity follow-up",4]
  
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
                as.numeric(round(adj.rd_p, p_digits)))
  
  regressionDF <- rbind(regressionDF, temp_row)
  
  
  temp_row_plot <- c(lin_reg_vars[lin_reg_vars$var_name == v, "var_label"], 
                     as.numeric(rd), 
                     as.numeric(rd_ci[1]),
                     as.numeric(rd_ci[2]),
                     as.numeric(round(rd_p, p_digits)),
                     as.numeric(adj.rd), 
                     as.numeric(adj.rd_ci[1]),
                     as.numeric(adj.rd_ci[2]),
                     as.numeric(round(adj.rd_p, p_digits))
  )
  
  plotDF <- rbind(plotDF, temp_row_plot)
}

names(regressionDF) =
  c("Variable",
    "Mean Difference",
    "95% CI",
    "p-value",
    "Adj. Mean Difference",
    "Adj. 95% CI",
    "Adj. p-value")


names(plotDF) =
  c("Variable",
    "main.estimate",
    "main.lower",
    "main.upper",
    "main.p",
    "adj.estimate",
    "adj.lower",
    "adj.upper",
    "adj.p")

plotDF <-  plotDF %>% 
  gather(label, value, 2:9) %>% 
  separate(label, into = c("model", "estimate")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = estimate, value = value) %>% 
  mutate(country = "South Africa",
         outcome = "mean_difference",
         version = "vs")

results_list[[length(results_list) + 1]] <- plotDF
names(results_list)[length(results_list)] <- "vs.rd.south_africa"


# Word Doc output
ft <- flextable(regressionDF) %>% 
  colformat_double(big.mark=",", digits = 1, na_str = "N/A") %>%
  add_header_row(values = c("", "Baseline Analysis", "Adjusted Analysis"), 
                 top = T, colwidths = c(1, 3, 3)) %>% 
  width(width = 1.1) %>%
  width(j = 1, width = 2.5) %>%
  theme_vanilla()  %>% 
  vline(j = c(1, 4), border = fp_border_default(), part = "all") %>% 
  align(i =1, align = "center", part = "header") 

ft
```


### Page Break

## Figure 1: Distribution of clinical CVD measures, by arm [dashed lines are group medians]

```{r, echo = F, fig.height = 3.4, fig.width = 10, warning = F, message = F, results = 'hide'}
# Graph, by arm

df <- ncd_merged_subset %>% 
  mutate(`Trial arm` = armR)

alpha_level = 0.35

lp = "bottom" # legend position, for all distribution figures

# SBP        

cvd_cutoff1 = 120
cvd_cutoff2 = NA
cvd_label = "Systolic BP (mmHg)"

mu <- plyr::ddply(df, "`Trial arm`", summarise, grp.median=median(sbp_mean_exit, na.rm = T))

p1 <-  ggplot(data = df, aes(x = sbp_mean_exit, fill = `Trial arm`)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values=c("#1b9e77", "#d95f02")) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=`Trial arm`),
             linetype="dashed") +
  scale_color_manual(values=c("#1b9e77", "#d95f02")) +
  theme_bw() +
  theme(legend.position = lp)

# DBP

cvd_cutoff1 = 80
cvd_cutoff2 = NA
cvd_label = "Diastolic BP (mmHg)"

mu <- plyr::ddply(df, "`Trial arm`", summarise, grp.median=median(dbp_mean_exit, na.rm = T))

p2 <-  ggplot(data = df, aes(x = dbp_mean_exit, fill = `Trial arm`)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values=c("#1b9e77", "#d95f02")) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=`Trial arm`),
             linetype="dashed") +
  scale_color_manual(values=c("#1b9e77", "#d95f02")) +
  theme_bw() +
  theme(legend.position = lp)

# BMI

cvd_cutoff1 = 25
cvd_cutoff2 = 30
cvd_label = "BMI (kg/m^2)"

mu <- plyr::ddply(df, "`Trial arm`", summarise, grp.median=median(bmi, na.rm = T))

p3 <-  ggplot(data = df, aes(x = bmi, fill = `Trial arm`)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values=c("#1b9e77", "#d95f02")) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=`Trial arm`),
             linetype="dashed") +
  scale_color_manual(values=c("#1b9e77", "#d95f02")) +
  theme_bw() +
  theme(legend.position = lp)

# A1c

cvd_cutoff1 = 5.7
cvd_cutoff2 = 6.5
cvd_label = "Hemoglobin A1c (%)"

mu <- plyr::ddply(df, "`Trial arm`", summarise, grp.median=median(hemoglobin_a1c_result, na.rm = T))

p4 <-  ggplot(data = df, aes(x = hemoglobin_a1c_result, fill = `Trial arm`)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values=c("#1b9e77", "#d95f02")) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=`Trial arm`),
             linetype="dashed") +
  scale_color_manual(values=c("#1b9e77", "#d95f02")) +
  theme_bw() +
  theme(legend.position = lp)

# Cholesterol


cvd_cutoff1 = 200
cvd_cutoff2 = 240
cvd_label = "Total cholesterol (mg/dL"

mu <- plyr::ddply(df, "`Trial arm`", summarise, grp.median=median(lipid_result, na.rm = T))

p5 <-  ggplot(data = df, aes(x = lipid_result, fill = `Trial arm`)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values=c("#1b9e77", "#d95f02")) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color= `Trial arm`),
             linetype="dashed") +
  scale_color_manual(values=c("#1b9e77", "#d95f02")) +
  theme_bw() +
  theme(legend.position = lp)

# Arrange!

ggarrange(p1, p3, p4, p5, # skip DBP
          ncol = 2)

```

### Page Break

## Figure 2: Distribution, by viral suppression status at endline [dashed lines are group medians]

```{r, echo = F, fig.height = 3.4, fig.width = 10, warning = F, message = F, results = 'hide'}
df = ncd_merged_subset %>%  filter(exit_viral_load_suppressedR != "Missing")
alpha_level = 0.5

# Excluding missing

vs_colors = c("#00A08A", "#FF0000")
vs_title = "Virally suppressed at exit"

# SBP        

cvd_cutoff1 = 120
cvd_cutoff2 = NA
cvd_label = "Systolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(sbp_mean_exit, na.rm = T))

p1 <-  ggplot(df, aes(x = sbp_mean_exit, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values= vs_colors) +
  guides(fill=guide_legend(title = vs_title)) +
  guides(color=guide_legend(title = vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# DBP

cvd_cutoff1 = 80
cvd_cutoff2 = NA
cvd_label = "Diastolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(dbp_mean_exit, na.rm = T))

p2 <-  ggplot(df, aes(x = dbp_mean_exit, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values= vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values= vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# BMI

cvd_cutoff1 = 25
cvd_cutoff2 = 30
cvd_label = "BMI (kg/m^2)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(bmi, na.rm = T))

p3 <-  ggplot(df, aes(x = bmi, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# A1c

cvd_cutoff1 = 5.7
cvd_cutoff2 = 6.5
cvd_label = "Hemoglobin A1c (%)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(hemoglobin_a1c_result, na.rm = T))

p4 <-  ggplot(df, aes(x = hemoglobin_a1c_result, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# Cholesterol


cvd_cutoff1 = 200
cvd_cutoff2 = 240
cvd_label = "Total cholesterol (mg/dL)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(lipid_result, na.rm = T))

p5 <-  ggplot(df, aes(x = lipid_result, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# Arrange!

ggarrange(p1, p3, p4, p5, # skip DBP
          ncol = 2)


```

### Page Break

## Figure 3: Distribution, by VS status at endline [dashed lines are group medians] - COMMUNITY ONLY

```{r, echo = F, fig.height = 3.4, fig.width = 10, warning = F, message = F, results = 'hide'}
df = ncd_merged_subset %>%  
  filter(exit_viral_load_suppressedR != "Missing") %>% 
  filter(armR == "Community follow-up")

alpha_level = 0.5

# Excluding missing

vs_colors = c("#00A08A", "#FF0000")
vs_title = "Virally suppressed at exit"

# SBP        

cvd_cutoff1 = 120
cvd_cutoff2 = NA
cvd_label = "Systolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(sbp_mean_exit, na.rm = T))

p1 <-  ggplot(df, aes(x = sbp_mean_exit, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values= vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# DBP

cvd_cutoff1 = 80
cvd_cutoff2 = NA
cvd_label = "Diastolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(dbp_mean_exit, na.rm = T))

p2 <-  ggplot(df, aes(x = dbp_mean_exit, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values= vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values= vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# BMI

cvd_cutoff1 = 25
cvd_cutoff2 = 30
cvd_label = "BMI (kg/m^2)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(bmi, na.rm = T))

p3 <-  ggplot(df, aes(x = bmi, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# A1c

cvd_cutoff1 = 5.7
cvd_cutoff2 = 6.5
cvd_label = "Hemoglobin A1c (%)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(hemoglobin_a1c_result, na.rm = T))

p4 <-  ggplot(df, aes(x = hemoglobin_a1c_result, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# Cholesterol


cvd_cutoff1 = 200
cvd_cutoff2 = 240
cvd_label = "Total cholesterol (mg/dL)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(lipid_result, na.rm = T))

p5 <-  ggplot(df, aes(x = lipid_result, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# Arrange!

ggarrange(p1, p3, p4, p5, # skip DBP
          ncol = 2)


```

### Page Break

## Figure 4: Distribution, by VS status at endline [dashed lines are group medians] - CLINIC ONLY

```{r, echo = F, fig.height = 3.4, fig.width = 10, warning = F, message = F, results = 'hide'}
df = ncd_merged_subset %>%  
  filter(exit_viral_load_suppressedR != "Missing") %>% 
  filter(armR == "Clinic follow-up")

alpha_level = 0.5

# Excluding missing

vs_colors = c("#00A08A", "#FF0000")
vs_title = "Virally suppressed at exit"

# SBP        

cvd_cutoff1 = 120
cvd_cutoff2 = NA
cvd_label = "Systolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(sbp_mean_exit, na.rm = T))

p1 <-  ggplot(df, aes(x = sbp_mean_exit, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values= vs_colors) +
  guides(fill=guide_legend(title= vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# DBP

cvd_cutoff1 = 80
cvd_cutoff2 = NA
cvd_label = "Diastolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(dbp_mean_exit, na.rm = T))

p2 <-  ggplot(df, aes(x = dbp_mean_exit, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values= vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values= vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# BMI

cvd_cutoff1 = 25
cvd_cutoff2 = 30
cvd_label = "BMI (kg/m^2)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(bmi, na.rm = T))

p3 <-  ggplot(df, aes(x = bmi, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# A1c

cvd_cutoff1 = 5.7
cvd_cutoff2 = 6.5
cvd_label = "Hemoglobin A1c (%)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(hemoglobin_a1c_result, na.rm = T))

p4 <-  ggplot(df, aes(x = hemoglobin_a1c_result, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# Cholesterol


cvd_cutoff1 = 200
cvd_cutoff2 = 240
cvd_label = "Total cholesterol (mg/dL)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(lipid_result, na.rm = T))

p5 <-  ggplot(df, aes(x = lipid_result, fill = exit_viral_load_suppressedR)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = vs_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# Arrange!

ggarrange(p1, p3, p4, p5, # skip DBP
          ncol = 2)

```

### Page Break

## Figure 5: Distribution of CVD variables by age

```{r, echo = F, fig.height = 3.4, fig.width = 10, warning = F, message = F, results = 'hide'}
df = ncd_merged_subset %>%
  filter(exit_viral_load_suppressedR != "Missing")

alpha_low = 0.7
alpha_hi = 0.8
point_size = 2
lp = "bottom" # legend position, for all distribution figures


# Excluding missing

vs_colors = c("#00A08A", "#FF0000")
vs_title = "Virally suppressed at exit"

# SBP        

cvd_cutoff1 = 120
cvd_cutoff2 = NA
cvd_label = "Systolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(sbp_mean_exit, na.rm = T))

p1 <- ggplot(df, aes(x = age, y = sbp_mean_exit, group = exit_viral_load_suppressedR)) + 
  geom_point(aes(color = exit_viral_load_suppressedR,
                 alpha = ifelse(sbp_mean_exit > 120, alpha_hi, alpha_low)), size = point_size) +   
  ylab(cvd_label) + xlab("Age") +
  geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_hline(data=mu, aes(yintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  scale_fill_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# DBP

cvd_cutoff1 = 80
cvd_cutoff2 = NA
cvd_label = "Diastolic BP (mmHg)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(dbp_mean_exit, na.rm = T))

p2 <-  ggplot(df, aes(x = age, y = dbp_mean_exit, group = exit_viral_load_suppressedR)) + 
  geom_point(aes(color = exit_viral_load_suppressedR,
                 alpha = ifelse(dbp_mean_exit > 80, alpha_hi, alpha_low)), size = point_size) +   
  ylab(cvd_label) + xlab("Age") +
  geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_hline(data=mu, aes(yintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  scale_fill_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# BMI

cvd_cutoff1 = 25
cvd_cutoff2 = 30
cvd_label = "BMI (kg/m^2)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(bmi, na.rm = T))

p3 <-  ggplot(df, aes(x = age, y = bmi, group = exit_viral_load_suppressedR)) + 
  geom_point(aes(color = exit_viral_load_suppressedR,
                 alpha = ifelse(bmi > 25, alpha_hi, alpha_low)), size = point_size) +     
  ylab(cvd_label) + xlab("Age") +
  geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_hline(data=mu, aes(yintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  scale_fill_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# A1c

cvd_cutoff1 = 5.7
cvd_cutoff2 = 6.5
cvd_label = "Hemoglobin A1c (%)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(hemoglobin_a1c_result, na.rm = T))

p4 <-  ggplot(df, aes(x = age, y = hemoglobin_a1c_result, group = exit_viral_load_suppressedR)) + 
  geom_point(aes(color = exit_viral_load_suppressedR,
                 alpha = ifelse(hemoglobin_a1c_result > 5.7, alpha_hi, alpha_low)), size = point_size) +    
  ylab(cvd_label) + xlab("Age") +
  geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_hline(data=mu, aes(yintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  scale_fill_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# Cholesterol


cvd_cutoff1 = 200
cvd_cutoff2 = 240
cvd_label = "Total cholesterol (mg/dL)"

mu <- plyr::ddply(df, "exit_viral_load_suppressedR", summarise, grp.median=median(lipid_result, na.rm = T))

p5 <-  ggplot(df, aes(x = age, y = lipid_result, group = exit_viral_load_suppressedR)) + 
  geom_point(aes(color = exit_viral_load_suppressedR,
                 alpha = ifelse(lipid_result > 200, alpha_hi, alpha_low)), size = point_size) +       
  ylab(cvd_label) + xlab("Age") +
  geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_hline(data=mu, aes(yintercept=grp.median, color=exit_viral_load_suppressedR),
             linetype="dashed") +
  scale_color_manual(values = vs_colors) +
  scale_fill_manual(values = vs_colors) +
  guides(fill=guide_legend(title=vs_title)) +
  guides(color=guide_legend(title=vs_title)) +
  theme_bw() +
  theme(legend.position = lp)

# Arrange!

ggarrange(p1, p3, p4, p5, # skip DBP
          ncol = 2)

```

### Page Break

## Figure 6: Distribution of CVD variables by age category

```{r, echo = F, fig.height = 3.4, fig.width = 10, warning = F, message = F, results = 'hide'}
alpha_level = 0.3

# Excluding missing

age_colors = c("#1b9e77", "#d95f02", "#00A08A", "#FF0000")
age_title = "Age category"

# SBP        

cvd_cutoff1 = 120
cvd_cutoff2 = NA
cvd_label = "Systolic BP (mmHg)"

mu <- plyr::ddply(df, "age_cat", summarise, grp.median=median(sbp_mean_exit, na.rm = T))

p1 <-  ggplot(df, aes(x = sbp_mean_exit, fill = age_cat)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = age_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=age_cat),
             linetype="dashed") +
  scale_color_manual(values= age_colors) +
  guides(fill=guide_legend(title= age_title)) +
  guides(color=guide_legend(title=age_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# DBP

cvd_cutoff1 = 80
cvd_cutoff2 = NA
cvd_label = "Diastolic BP (mmHg)"

mu <- plyr::ddply(df, "age_cat", summarise, grp.median=median(dbp_mean_exit, na.rm = T))

p2 <-  ggplot(df, aes(x = dbp_mean_exit, fill = age_cat)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values= age_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=age_cat),
             linetype="dashed") +
  scale_color_manual(values= age_colors) +
  guides(fill=guide_legend(title=age_title)) +
  guides(color=guide_legend(title=age_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# BMI

cvd_cutoff1 = 25
cvd_cutoff2 = 30
cvd_label = "BMI (kg/m^2)"

mu <- plyr::ddply(df, "age_cat", summarise, grp.median=median(bmi, na.rm = T))

p3 <-  ggplot(df, aes(x = bmi, fill = age_cat)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = age_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=age_cat),
             linetype="dashed") +
  scale_color_manual(values = age_colors) +
  guides(fill=guide_legend(title=age_title)) +
  guides(color=guide_legend(title=age_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# A1c

cvd_cutoff1 = 5.7
cvd_cutoff2 = 6.5
cvd_label = "Hemoglobin A1c (%)"

mu <- plyr::ddply(df, "age_cat", summarise, grp.median=median(hemoglobin_a1c_result, na.rm = T))

p4 <-  ggplot(df, aes(x = hemoglobin_a1c_result, fill = age_cat)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = age_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=age_cat),
             linetype="dashed") +
  scale_color_manual(values = age_colors) +
  guides(fill=guide_legend(title=age_title)) +
  guides(color=guide_legend(title=age_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# Cholesterol


cvd_cutoff1 = 200
cvd_cutoff2 = 240
cvd_label = "Total cholesterol (mg/dL)"

mu <- plyr::ddply(df, "age_cat", summarise, grp.median=median(lipid_result, na.rm = T))

p5 <-  ggplot(df, aes(x = lipid_result, fill = age_cat)) + 
  geom_density(alpha = alpha_level) +   
  scale_fill_manual(values = age_colors) +
  xlab(cvd_label) +
  geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
  geom_vline(data=mu, aes(xintercept=grp.median, color=age_cat),
             linetype="dashed") +
  scale_color_manual(values = age_colors) +
  guides(fill=guide_legend(title=age_title)) +
  guides(color=guide_legend(title=age_title)) +
  theme_bw()  +
  theme(legend.position = lp)

# Arrange!

ggarrange(p1, p3, p4, p5, # skip DBP
          ncol = 2)

```

### Page Break

## Figure S1: Relative risk of cardiovascular risk for community follow-up compared with clinic follow-up

```{r, echo = F, fig.width = 12,  fig.height = 8, warning = F, message = F}
plotDF <- rbind(results_list[["arm.rr.south_africa"]],
                results_list[["arm.rd.south_africa"]])

plotDF <- plotDF %>% 
  select(country, model, outcome, Variable, estimate, lower, upper, p) %>% 
  group_by(outcome) %>% 
  arrange(Variable) %>% 
  mutate(index = as.numeric(as.factor(Variable))) 

# Relative Risk

plot_subset <- plotDF %>% 
  filter(outcome == "relative_risk")

y_breaks <- -sort(unique(plot_subset$index))
y_labels <- unique(plot_subset$Variable)

graph_intercept = 1

ggplot() + 
  geom_vline(xintercept = graph_intercept, 
             linetype="dashed",size = .9) +
  geom_point(data = plot_subset[plot_subset$model == "main", ], 
             aes(x = estimate, y = -index, color = "1. Baseline model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "main", ], 
               aes(y= -index, yend=-index, x = lower, xend = upper,
                   color = "1. Baseline model"), 
               size=.95) +
  geom_point(data = plot_subset[plot_subset$model == "adj", ], 
             aes(x = estimate, y = -index + .15, color = "2. Primary (adjusted)\n model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "adj", ], 
               aes(y= -index + .15, yend=-index + .15, x = lower, xend = upper,
                   color = "2. Primary (adjusted)\n model"), 
               size=.95) +
  geom_point(data = plot_subset[plot_subset$model == "ext", ], 
             aes(x = estimate, y = -index + .3, color = "3. Extended model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "ext", ], 
               aes(y= -index + .3, yend=-index +.3, x = lower, xend = upper,
                   color = "3. Extended model"), 
               size=.95) +
  scale_color_manual(values = c("#374E55FF", #gray
                                "#DF8F44FF", #orange
                                "#79AF97FF" # green
  )) +
  facet_wrap(~country, ncol = 1) +
  
  # Labels
  # ggtitle("Estimated Relative Risk of Cardiovascular Risk") +
  xlab("Estimate and 95% Confidence Interval") +
  ylab("Variable") +
  scale_y_continuous(labels = y_labels, breaks=y_breaks) +
  
  # Theme
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 18, face="bold"), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text.y= element_text(face='plain',size=16),
        axis.text.x=element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y=element_text(size=16),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16)) 
```

## Figure S2: Mean difference of cardiovascular risk variables for community follow-up compared with clinic follow-up

```{r, echo = F, fig.width = 12,  fig.height = 8}

# RISK DIFF

plot_subset <- plotDF %>% 
  filter(outcome == "mean_difference")

y_breaks <- -sort(unique(plot_subset$index))
y_labels <- unique(plot_subset$Variable)

graph_intercept = 0

ggplot() + 
  geom_vline(xintercept = graph_intercept, 
             linetype="dashed", size = .9) +
  geom_point(data = plot_subset[plot_subset$model == "main", ], 
             aes(x = estimate, y = -index, color = "1. Baseline model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "main", ], 
               aes(y= -index, yend=-index, x = lower, xend = upper,
                   color = "1. Baseline model"), 
               size=.95) +
  geom_point(data = plot_subset[plot_subset$model == "adj", ], 
             aes(x = estimate, y = -index + .15, color = "2. Primary (adjusted)\n model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "adj", ], 
               aes(y= -index + .15, yend=-index + .15, x = lower, xend = upper,
                   color = "2. Primary (adjusted)\n model"), 
               size=.95) +
  geom_point(data = plot_subset[plot_subset$model == "ext", ], 
             aes(x = estimate, y = -index + .3, color = "3. Extended model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "ext", ], 
               aes(y= -index + .3, yend=-index +.3, x = lower, xend = upper,
                   color = "3. Extended model"), 
               size=.95) +
  scale_color_manual(values = c("#374E55FF", #gray
                                "#DF8F44FF", #orange
                                "#79AF97FF" # green
  )) +
  facet_wrap(~country, ncol = 1) +
  
  # Labels
  # ggtitle("Estimated Difference in Cardiovascular Risk") +
  xlab("Estimate and 95% Confidence Interval") +
  ylab("Variable") +
  scale_y_continuous(labels = y_labels, breaks=y_breaks) +
  
  # Theme
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 18, face="bold"), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text.y= element_text(face='plain',size=16),
        axis.text.x=element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y=element_text(size=16),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16)) 
```

## Figure S3: Relative risk of cardiovascular risk for people who are virally suppressed versus not

```{r, echo = F, fig.width = 12,  fig.height = 8, warning = F, message = F}
plotDF <- rbind(results_list[["vs.rr.south_africa"]],
                results_list[["vs.rd.south_africa"]])

plotDF <- plotDF %>% 
  select(country, model, outcome, Variable, estimate, lower, upper, p) %>% 
  group_by(outcome) %>% 
  arrange(Variable) %>% 
  mutate(index = as.numeric(as.factor(Variable))) 

# Relative Risk

plot_subset <- plotDF %>% 
  filter(outcome == "relative_risk")

y_breaks <- -sort(unique(plot_subset$index))
y_labels <- unique(plot_subset$Variable)

graph_intercept = 1

ggplot() + 
  geom_vline(xintercept = graph_intercept, 
             linetype="dashed",size = .9) +
  geom_point(data = plot_subset[plot_subset$model == "main", ], 
             aes(x = estimate, y = -index, color = "1. Baseline model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "main", ], 
               aes(y= -index, yend=-index, x = lower, xend = upper,
                   color = "1. Baseline model"), 
               size=.95) +
  geom_point(data = plot_subset[plot_subset$model == "adj", ], 
             aes(x = estimate, y = -index + .15, color = "2. Primary (adjusted)\n model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "adj", ], 
               aes(y= -index + .15, yend=-index + .15, x = lower, xend = upper,
                   color = "2. Primary (adjusted)\n model"), 
               size=.95) +
  scale_color_manual(values = c("#374E55FF", #gray
                                "#DF8F44FF", #orange
                                "#79AF97FF" # green
  )) +
  facet_wrap(~country, ncol = 1) +
  
  # Labels
  # ggtitle("Estimated Relative Risk of Cardiovascular Risk") +
  xlab("Estimate and 95% Confidence Interval") +
  ylab("Variable") +
  scale_y_continuous(labels = y_labels, breaks=y_breaks) +
  
  # Theme
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 18, face="bold"), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text.y= element_text(face='plain',size=16),
        axis.text.x=element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y=element_text(size=16),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16)) 
```

## Figure S4: Mean difference in cardiovascular risk variables for people who are virally suppressed versus not

```{r, echo = F, fig.width = 12,  fig.height = 8}

# RISK DIFF

plot_subset <- plotDF %>% 
  filter(outcome == "mean_difference")

y_breaks <- -sort(unique(plot_subset$index))
y_labels <- unique(plot_subset$Variable)

graph_intercept = 0

ggplot() + 
  geom_vline(xintercept = graph_intercept, 
             linetype="dashed", size = .9) +
  geom_point(data = plot_subset[plot_subset$model == "main", ], 
             aes(x = estimate, y = -index, color = "1. Baseline model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "main", ], 
               aes(y= -index, yend=-index, x = lower, xend = upper,
                   color = "1. Baseline model"), 
               size=.95) +
  geom_point(data = plot_subset[plot_subset$model == "adj", ], 
             aes(x = estimate, y = -index + .15, color = "2. Primary (adjusted)\n model"), 
             size=3.5) +
  geom_segment(data = plot_subset[plot_subset$model == "adj", ], 
               aes(y= -index + .15, yend=-index + .15, x = lower, xend = upper,
                   color = "2. Primary (adjusted)\n model"), 
               size=.95) +
  scale_color_manual(values = c("#374E55FF", #gray
                                "#DF8F44FF", #orange
                                "#79AF97FF" # green
  )) +
  facet_wrap(~country, ncol = 1) +
  
  # Labels
  # ggtitle("Estimated Difference in Cardiovascular Risk") +
  xlab("Estimate and 95% Confidence Interval") +
  ylab("Variable") +
  scale_y_continuous(labels = y_labels, breaks=y_breaks) +
  
  # Theme
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.title = element_text(size = 18, face="bold"), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.text.y= element_text(face='plain',size=16),
        axis.text.x=element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.title.y=element_text(size=16),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size = 16)) 
```
# NCD plot functions for appendix figures showing distribution of CVD vars
# 5/4/2019

# To do: should also print the N / missing afterwards

# ==============================================================================

# PLOT LABELS

# CVD 'at-risk' cutoffs
cvd_var_df <- data.frame(
  cvd_var = c("sbp_mean_exit", "dbp_mean_exit", "bmi", "hemoglobin_a1c_result", "lipid_result"),
  cvd_cutoff1 = c(120, 80, 25, 5.7, 200),
  cvd_cutoff2 = c(NA, NA, 30, 6.5, 240),
  cvd_label = c("Systolic BP (mmHg)", "Diastolic BP (mmHg)", "BMI (kg/m^2)", "Hemoglobin A1c (%)", "Total cholesterol (mg/dL)")
)

# ==============================================================================

#################
# DENSITY PLOTS #
#################

# Figure function for single plot
ggdensity <- function(df, fill_var, cvd_var) {
  
  # df <- ncd_merged_subset
  # fill_var = "exit_viral_load_suppressedR"
  # cvd_var = "sbp_mean_exit"
  
  plot_df = df %>%  filter( get(fill_var) != "Missing")
  lp = "bottom" # legend position
  
  if (fill_var == "exit_viral_load_suppressedR") {
    
    fill_colors = c("#00A08A", "#FF0000")
    legend_title = "Virally suppressed at exit"
    alpha_level = 0.5 
  
  } else if (fill_var == "armR") {
    
    fill_colors = c("#1b9e77", "#d95f02")
    legend_title = "Trial arm"
    alpha_level = 0.4 
    
  } else if (fill_var == "genderR") {
    
    fill_colors = c("darkgoldenrod1", "#00A08A")
    legend_title = "Gender"
    alpha_level = 0.5

  }
  
  # Cutoffs for 'at-risk'
  selected_var <- cvd_var_df[cvd_var_df$cvd_var == cvd_var, ]
  cvd_cutoff1 <- selected_var$cvd_cutoff1
  cvd_cutoff2 <- selected_var$cvd_cutoff2
  cvd_label <- selected_var$cvd_label
  
  # Dashed line for medians
  mu <- plyr::ddply(plot_df, fill_var, summarise, grp.median = median(get(cvd_var), na.rm = T))
  
  # Density plot!
  p <- ggplot(plot_df, aes_string(x = cvd_var, fill = fill_var)) + 
    geom_density(alpha = alpha_level) +   
    scale_fill_manual(values = fill_colors) + # for density plots
    xlab(cvd_label) +
    geom_vline(xintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
    geom_vline(data = mu, aes(xintercept = grp.median, color = get(fill_var)),
               linetype="dashed") +
    scale_color_manual(values = fill_colors) + # for lines
    guides(fill = guide_legend(title = legend_title)) +
    guides(color = guide_legend(title = legend_title)) +
    theme_bw() +
    theme(legend.position = lp)
  
  return(p)
}

# Figure function for 4 CVD plots together
make_density_figs <- function(df, fill_var) {
  
  p1 <- ggdensity(df, fill_var, "sbp_mean_exit")
  #p2 <-  ggdensity(fill_var, "dbp_mean_exit")
  p3 <- ggdensity(df, fill_var, "bmi")
  p4 <- ggdensity(df, fill_var, "hemoglobin_a1c_result")
  p5 <- ggdensity(df, fill_var, "lipid_result")
  
  # Arrange!
  ggarrange(p1, p3, p4, p5, # skip DBP
            ncol = 2)
}

# make_density_figs(ncd_merged_subset, "exit_viral_load_suppressedR")

# ==============================================================================

#######################################
# DISTRIBUTION BY ANOTHER VAR 
# If continuous --> scatterplot
# If categorical --> boxplot
#######################################

ggdistribution_byVar <- function(df, by_var, fill_var, cvd_var) {
  
  # df <- ncd_merged_subset
  # fill_var = "exit_viral_load_suppressedR"
  # cvd_var = "sbp_mean_exit"
  # by_var = "age"
  alpha_level = 0.8
  point_size = 2
  
  plot_df = df %>%  filter( get(fill_var) != "Missing")
  lp = "bottom" # legend position
  
  if (fill_var == "exit_viral_load_suppressedR") {
    
    fill_colors = c("#00A08A", "#FF0000")
    legend_title = "Virally suppressed at exit"
    
  } else if (fill_var == "armR") {
    
    fill_colors = c("#1b9e77", "#d95f02")
    legend_title = "Trial arm"
    
  } else if (fill_var == "genderR") {
    
    fill_colors = c("darkgoldenrod1", "#00A08A")
    legend_title = "Gender"
    
  }
  
  # Cutoffs for 'at-risk'
  selected_var <- cvd_var_df[cvd_var_df$cvd_var == cvd_var, ]
  cvd_cutoff1 <- selected_var$cvd_cutoff1
  cvd_cutoff2 <- selected_var$cvd_cutoff2
  cvd_label <- selected_var$cvd_label
  
  # SCATTERPLOT for numeric vars
  
  if (is.numeric(plot_df[[by_var]]) | is.integer(plot_df[[by_var]])) {
  
    # Dashed line for medians
    mu <- plyr::ddply(plot_df, fill_var, summarise, grp.median = median(get(cvd_var), na.rm = T))

    p <- ggplot(plot_df, aes_string(x = by_var, y = cvd_var, group = fill_var)) + 
      geom_point(aes_string(color = fill_var),
                 alpha = alpha_level, size = point_size) +   
      ylab(cvd_label) + 
      xlab(stringr::str_to_title(by_var)) + 
      geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), size = .6) +
      geom_hline(data=mu, aes(yintercept=grp.median, color=exit_viral_load_suppressedR),
                 linetype="dashed") +
      scale_color_manual(values = fill_colors) +
      scale_fill_manual(values = fill_colors) +
      guides(fill = guide_legend(title = legend_title)) +
      guides(color = guide_legend(title = legend_title)) +
      theme_bw() +
      theme(legend.position = lp)
  
  # BOX PLOT for categorical vars  
    
  } else if (is.character(plot_df[[by_var]]) | is.factor(plot_df[[by_var]])) {
    
    p <-  ggplot(plot_df, aes_string(x = by_var, y = cvd_var, fill = fill_var)) +
      geom_boxplot(varwidth = TRUE, alpha = 0.7)  +
      scale_fill_manual(values = fill_colors) +
      ylab(cvd_label) + 
      xlab(stringr::str_to_title(by_var)) + 
      geom_hline(yintercept = c(cvd_cutoff1, cvd_cutoff2), linetype = "dashed", size = .6) +
      scale_color_manual(values = fill_colors) +
      scale_fill_manual(values = fill_colors) +
      guides(fill = guide_legend(title = legend_title)) +
      guides(color = guide_legend(title = legend_title)) +
      theme_bw() +
      theme(legend.position = lp)
  }
    
  return(p)
}

# Figure function for 4 CVD plots together
make_distribution_figs <- function(df, by_var, fill_var) {
  
  p1 <- ggdistribution_byVar(df, by_var, fill_var, "sbp_mean_exit")
  #p2 <-  ggdistribution_byVar(df, by_var, fill_var, "dbp_mean_exit")
  p3 <- ggdistribution_byVar(df, by_var, fill_var, "bmi")
  p4 <- ggdistribution_byVar(df, by_var, fill_var, "hemoglobin_a1c_result")
  p5 <- ggdistribution_byVar(df, by_var, fill_var, "lipid_result")
  
  # Arrange!
  ggarrange(p1, p3, p4, p5, # skip DBP
            ncol = 2)
}

# make_distribution_figs(ncd_merged_subset, "age", "exit_viral_load_suppressedR")
# make_distribution_figs(ncd_merged_subset, "education", "exit_viral_load_suppressedR")
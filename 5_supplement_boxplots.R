# MSahu
# Feb 24, 2024

# Boxplots - change in BMI / SBP from baseline to exit

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

# Some values with very high BMI!

check <- data %>%  filter(baseline_bmi > 40 | bmi > 40)

prop_table(data$bmi_cat_exit, data$exit_viral_load_suppressedR, data$siteR)


# ------------------------------------------------------------------------------

ggplot(data = data, aes(x = genderR, y = baseline_bmi, fill = exit_viral_load_suppressedR)) + 
  geom_boxplot(varwidth = TRUE) +
  facet_wrap(~siteR + armR)

ggplot(data = data, aes(x = genderR, y = bmi, fill = exit_viral_load_suppressedR)) + 
  geom_bar(position="dodge", stat = "summary", fun.y = "mean") +
  facet_wrap(~siteR + armR)

ggplot(data = data, aes(x = genderR, y = bmi_diff, fill = exit_viral_load_suppressedR)) + 
  geom_bar(position="dodge", stat = "summary", fun.y = "mean") +
  facet_wrap(~siteR + armR)

# BOXPLOTS


pdf(file = paste0(out_dir, "supplement_boxplots_BMI_diff.pdf"),
    width = 10, 
    height = 10) 

ggplot(data = data, aes(x = genderR, y = bmi_diff, fill = exit_viral_load_suppressedR)) + 
  geom_boxplot(varwidth = F) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(aes(group = exit_viral_load_suppressedR), 
               fun.y=mean, shape = 18, geom="point", size=4, position = position_dodge(0.75), color="darkred", fill="darkred") +
  facet_wrap(~armR + siteR) +
  ylim(-20,20) +
  theme_bw(base_size = 15) + theme(strip.background = element_rect(fill="white")) +
  xlab("Gender") + ylab("Difference in BMI from baseline to exit (kg/m^2)") + labs(fill='Virally\nsuppressed\nat exit?') 

dev.off()

data %>%  filter (!is.na(bmi_diff)) %>%  nrow()




pdf(file = paste0(out_dir, "supplement_boxplots_sbp_diff.pdf"),
    width = 10, 
    height = 10) 

ggplot(data = data, aes(x = genderR, y = sbp_diff, fill = exit_viral_load_suppressedR)) + 
  geom_boxplot(varwidth = F) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(aes(group = exit_viral_load_suppressedR), 
               fun.y=mean, shape = 18, geom="point", size=4, position = position_dodge(0.75), color="darkred", fill="darkred") +
  facet_wrap(~armR + siteR) +
  ylim(-50,50) +
  theme_bw(base_size = 15) + theme(strip.background = element_rect(fill="white")) +
  xlab("Gender") + ylab("Difference in systolic blood pressure from baseline to exit (mmHg)") + labs(fill='Virally\nsuppressed\nat exit?') 

dev.off()

data %>%  filter (!is.na(sbp_diff)) %>%  nrow()



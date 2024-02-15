# From class project (SAP code -- RMD file)
# Need to update...

# ------------------------------------------------------------------------------

rm(list=ls())

# Directories
dir <- "C:/Users/msahu/OneDrive - UW/Documents/Research/DGH+MGH/DO_ART_NCD/"
in_dir <- paste0(dir, "0_data/3_cleaned_data/")
out_dir <-  paste0(dir, "/3_plots/")

# Load data
ncd_merged <- readRDS(paste0(in_dir,"ncd_merged.rds"))
ncd_merged_sa <- readRDS(paste0(in_dir,"ncd_merged_subset.rds")) %>% 
  filter(country == "South Africa")

# Packages
pacman::p_load(naniar, ggplot2)

# ------------------------------------------------------------------------------

# Visualize Missingness

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

# Save pdf - by site
pdf(file = paste0(out_dir, "supplement_missing_data_by_site.pdf"),
    width = 8, 
    height = 4) 

gg_miss_fct(x = ncd_outcomes, fct = siteR) + labs(title = "Missing data by trial site")

dev.off()

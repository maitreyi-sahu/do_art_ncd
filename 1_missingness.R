# From class project (SAP code -- RMD file)
# Need to update...


# Table One
ncd_subset <- ncd_merged[,c(2:4, 7:23, 26:29)]
Vars <- colnames(ncd_subset)
# you need to specify which variables are factors
factor_var <- c("arm", "country", "site", "gender","education","occupation","smoking_status", 
                "baseline_who_stage", "baseline_viral_load_result_cat","exit_viral_load_result_cat") 
tb1 <- CreateTableOne(vars = Vars, data = ncd_subset, factorVars = factor_var)
print(tb1, quote = F, noSpaces = TRUE)

ncd_subset <- ncd_merged[,c(2:4, 7:23, 26:29)]

ncd_outcomes <- ncd_merged %>% 
  select(arm, country, site,
         baseline_viral_load_result, exit_viral_load_result,
         baseline_bmi, bmi, 
         hemoglobin_a1c_result, lipid_result,
         sbp_mean_Enrollment, sbp_mean_Exit,
         dbp_mean_Enrollment, dbp_mean_Exit,
         sbp_diff, bmi_diff)

ncd_outcomes %>%  group_by(arm) %>% descr(stats = "common", transpose = F, headings = FALSE)
library(naniar)
vis_miss(ncd_outcomes) 
gg_miss_fct(x = ncd_outcomes, fct = arm)  + labs(title = "Missing Data by trial arm")
gg_miss_fct(x = ncd_outcomes, fct = site) + labs(title = "Missing Data by trial site")
gg_miss_fct(x = ncd_outcomes, fct = country) + labs(title = "Missing Data by Country")
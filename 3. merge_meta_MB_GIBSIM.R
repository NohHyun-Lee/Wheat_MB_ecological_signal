library(stringr)
### library --------------------------------------------------------------------
library(randomForest)
library(caret)
library(xgboost)
library(dplyr)
library(Metrics)
library(DescTools)
library(dismo)
library(ggplot2)
library(gridExtra)
library(forcats)
library(readxl)
library(writexl)
### data ------------------------------------------------------------------------------------------------------------
df_rslt_stat <- read_xlsx("E:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v3.xlsx") 
df_rslt_gz <- read_xlsx("E:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/1.Run_GIBSIM/MB_v3.1/MB_wheat_wth_data_adj_inc_MB_v3.1.xlsx")
names(df_rslt_stat)
df_rslt_gz_wthid_year_gz <- df_rslt_gz %>% dplyr::select("wth_ID","year", "GZ_mean_max", "GZ_mean_mean", "GZ_mean_median")

df_rslt_stat_join <- left_join(df_rslt_stat, df_rslt_gz_wthid_year_gz, by = c("wth_ID", "year"))

writexl::write_xlsx(df_rslt_stat_join, "E:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v3.1.xlsx")


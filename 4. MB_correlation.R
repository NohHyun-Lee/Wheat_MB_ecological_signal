### library --------------------------------------------------------------------
library(readxl)
library(writexl)
library(dplyr)

### setting --------------------------------------------------------------------
# wdir <- "C:/Users/B550/Desktop/work/9. 숭실대 MB"
# setwd("C:/Users/B550/Desktop/work/9. 숭실대 MB")
# getwd()

### import data -----------------------------------------------------------------------
df_rslt_stat <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v4.xlsx") 
MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_ITSfull_v4.1.xlsx")
MB_max_filtered_data$Taxa[1:21]
df_rslt_stat$wth_ID
names(df_rslt_stat)

tail(df_rslt_stat)

sum(is.na(df_rslt_stat)==T)
df_rslt_stat <- df_rslt_stat[complete.cases(df_rslt_stat),]
nrow(df_rslt_stat)
names(df_rslt_stat)


library(tidyverse)
max_filtered_taxa <- MB_max_filtered_data$Taxa[1:21]
df_rslt_stat2 <- df_rslt_stat %>% dplyr::select("incidence", all_of(max_filtered_taxa))
names(df_rslt_stat2)

head(as.data.frame(df_rslt_stat2))

res_cor <- map_dfr(names(df_rslt_stat2)[2:ncol(df_rslt_stat2)], function(taxa){
  test <- cor.test(df_rslt_stat2$incidence, df_rslt_stat2[[taxa]], method = "spearman")
  
  tibble(
    taxa = taxa,
    rho = test$estimate,
    p = test$p.value
  )
})
res_cor <- res_cor %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),   # BH correction
    
    sig = p_adj < 0.05,                  # 
    
    rho_plot = ifelse(sig, rho, NA)
  ) %>%
  arrange(rho) %>%
  mutate(taxa = factor(taxa, levels = taxa))

cor_m_inc_g <- ggplot(res_cor, aes(x = taxa, y = rho)) +
 
  geom_segment(aes(x = taxa, xend = taxa,
                   y = 0, yend = rho),
               color = "grey70", linewidth = 0.6) +
 
  geom_point(data = res_cor %>% filter(!sig),
             color = "grey80",
             size = 3) +
  
  geom_point(data = res_cor %>% filter(sig),
             aes(color = rho),
             size = 3) +
  scale_color_gradient2(
    # low = "#4575b4",
    # mid = "white",
    # high = "#d73027",
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", linewidth = 0.6)+
  geom_hline(yintercept = 0)+
  # coord_flip() +
  theme_bw() +
  labs(
    x = "",
    y = "Rho",
    color = "Rho"
  )+
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  ) + guides(
    color = guide_colorbar(
      direction = "horizontal",
      barwidth = 8,
      barheight = 0.5,
      title.position = "top"
    )
  ) + 
  theme(
    legend.position = c(0.99, 0.01),   
    legend.justification = c(1, 0),    
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  )

save_cor_version <- "v7.3"
if(!dir.exists(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", save_cor_version))){
  dir.create(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", save_cor_version))
}

ggsave(filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", 
                         save_cor_version, "/MB_incidence_cor_plot_low_pval_max_filtered_", save_cor_version, "_(Rho0.3).png"),
       plot = cor_m_inc_g, 
       width  = 6.4, 
       height  = 5)

res_cor_for_save <- res_cor %>% dplyr::select(taxa, rho, p, p_adj, sig)

writexl::write_xlsx(res_cor_for_save, paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", 
                                             save_cor_version, "/MB_incidence_cor_df_low_pval_max_filtered_", save_cor_version, "_(Rho0.3).xlsx"))

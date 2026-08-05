#v1 : group by loc
#v4 : group by loc, data = merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v2.xlsx
#v5 : group by loc, data = merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v2.1.xlsx (CLR version)
#v6 : group by loc, data = merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v3.xlsx (CLR version, genus level)
#v6.2 : group by loc, data = merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v3.xlsx (CLR version, genus level, select key_microbes which mean(CLR(rel_abu)) > 2)
#v7 : data = merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v3.xlsx (CLR version, genus level, select key_microbes which mean(CLR(rel_abu)) > 2)
#v7.2 : same as v7, but Benjamini???Hochberg performed
#v7.2 : same as v7.2, max filter is different

### library --------------------------------------------------------------------
library(readxl)
library(writexl)
library(dplyr)

### setting --------------------------------------------------------------------
# wdir <- "C:/Users/B550/Desktop/work/9. 숭실대 MB"
# setwd("C:/Users/B550/Desktop/work/9. 숭실대 MB")
# getwd()

### import data -----------------------------------------------------------------------
# df_rslt_stat <- read_xlsx("E:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v1.xlsx")
df_rslt_stat <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v4.xlsx") 
# MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_more_ITSfull_v3.2.xlsx")
MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_ITSfull_v4.1.xlsx")
MB_max_filtered_data$Taxa[1:21]
df_rslt_stat$wth_ID
names(df_rslt_stat)

tail(df_rslt_stat)

# df_rslt_stat <- read_xlsx("E:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2-1. Run_GIBSIM/Run_GIBSIM/output/for_MB_analysis_v1.3.3/MB_cor_wth_data/MB_wheat_cor_wth_data_v1.3.3.xlsx")
# df_rslt_stat <- read_xlsx("F:/서울대/microbiome_SSU/Rcode/FHB_analysis/2-1. Run_GIBSIM/Run_GIBSIM/output/for_MB_analysis_v1.3.3/MB_cor_wth_data/MB_wheat_cor_wth_data_v1.3.3.xlsx")
# df_rslt_stat <- read_xlsx("E:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2-1. Run_GIBSIM/Run_GIBSIM/output/for_MB_analysis_v1.3/MB_cor_wth_data/MB_wheat_cor_wth_data_v1.3.xlsx")

sum(is.na(df_rslt_stat)==T)
df_rslt_stat <- df_rslt_stat[complete.cases(df_rslt_stat),]
nrow(df_rslt_stat)
names(df_rslt_stat)


library(tidyverse)
max_filtered_taxa <- MB_max_filtered_data$Taxa[1:21]
df_rslt_stat2 <- df_rslt_stat %>% dplyr::select("incidence", all_of(max_filtered_taxa))
names(df_rslt_stat2)

head(as.data.frame(df_rslt_stat2))


# correlation 계산
res_cor <- map_dfr(names(df_rslt_stat2)[2:ncol(df_rslt_stat2)], function(taxa){
  test <- cor.test(df_rslt_stat2$incidence, df_rslt_stat2[[taxa]], method = "spearman")
  
  tibble(
    taxa = taxa,
    rho = test$estimate,
    p = test$p.value
  )
})

# res_cor <- res_cor %>%
#   mutate(
#     sig = p < 0.05,
#     
#     # plotting용 색상 변수
#     rho_plot = ifelse(sig, rho, NA)
#   ) %>%
#   arrange(rho) %>%  # 정렬 (중요)
#   mutate(taxa = factor(taxa, levels = taxa))
res_cor <- res_cor %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),   # BH correction
    
    sig = p_adj < 0.05,                  # 
    
    rho_plot = ifelse(sig, rho, NA)
  ) %>%
  arrange(rho) %>%
  mutate(taxa = factor(taxa, levels = taxa))

cor_m_inc_g <- ggplot(res_cor, aes(x = taxa, y = rho)) +
  # stem
  geom_segment(aes(x = taxa, xend = taxa,
                   y = 0, yend = rho),
               color = "grey70", linewidth = 0.6) +
  # non-significant (먼저 깔기)
  geom_point(data = res_cor %>% filter(!sig),
             color = "grey80",
             size = 3) +
  # significant (색상 적용)
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
    legend.position = c(0.99, 0.01),   # (x, y) in [0,1]
    legend.justification = c(1, 0),    # 오른쪽 아래 기준 anchor
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
### <old version> ===============================================================
# df_rslt_stat <- df_rslt_stat %>% filter(wth_ID != "ID_043") %>% filter(wth_ID != "ID_124") %>%
#   filter(wth_ID != "ID_022") %>% filter(wth_ID != "ID_021")# %>% filter(wth_ID != "ID_038")# %>% filter(wth_ID != "ID_027") %>%
# # names(df_rslt_stat)[which(names(df_rslt_stat) == "adj_inc_pred")] <- "incidence"
# nrow(df_rslt_stat)
# plot(df_rslt_stat$Alternaria, df_rslt_stat$incidence)
# plot(df_rslt_stat$Epicoccum, df_rslt_stat$incidence)
# plot(df_rslt_stat$Ustilago, df_rslt_stat$incidence)
# plot(df_rslt_stat$Cladosporium, df_rslt_stat$incidence)
# plot(df_rslt_stat$Papiliotrema, df_rslt_stat$incidence)
# plot(df_rslt_stat$Fusarium, df_rslt_stat$incidence)
# plot(df_rslt_stat$Botrytis, df_rslt_stat$incidence)
# plot(df_rslt_stat$Hannaella, df_rslt_stat$incidence)
# plot(df_rslt_stat$Cladosporium, df_rslt_stat$incidence)
# plot(df_rslt_stat$Periconia, df_rslt_stat$incidence)
# plot(df_rslt_stat$Alternaria, df_rslt_stat$incidence)
# plot(df_rslt_stat$Cercospora, df_rslt_stat$incidence)
# 
# ### choose data for analysis
# MB_intensity_wth_cor_df <- df_rslt_stat
# nrow(MB_intensity_wth_cor_df)
# 
# ### every data should be numeric------------------------------------------------
# # MB_intensity_wth_cor_df <- MB_intensity_wth_cor_df %>% select(-contains("wsd"))
# MB_intensity_wth_cor_df <- as.data.frame(MB_intensity_wth_cor_df)
# 
# # for(x in c(18:ncol(MB_intensity_wth_cor_df))){ #18:ncol(MB_intensity_wth_cor_df) ==> only numeric! (OTU column and wth column)
# #   #x=18
# #   MB_intensity_wth_cor_df[,x] <- as.numeric(MB_intensity_wth_cor_df[,x])
# # }
# 
# ### remove only 0 --------------------------------------------------------------
# ncol(MB_intensity_wth_cor_df)
# colnames(MB_intensity_wth_cor_df)
# MB_intensity_wth_cor_df_numeric <- MB_intensity_wth_cor_df[,-c(1,2)]
# max_vals <- sapply(MB_intensity_wth_cor_df_numeric, max, na.rm = TRUE)
# imsi_data <- MB_intensity_wth_cor_df_numeric[, max_vals != 0]
# MB_intensity_wth_cor_df <- cbind(MB_intensity_wth_cor_df[,c(1,2)], imsi_data)
# ncol(MB_intensity_wth_cor_df)
# # MB_intensity_wth_cor_df$OTU0067 <- NULL
# # MB_intensity_wth_cor_df$OTU0077 <- NULL
# 
# ### correlation analysis -------------------------------------------------------
# # x <- MB_intensity_wth_cor_df$Ustilago_nuda
# # y <- MB_intensity_wth_cor_df$incidence
# # 
# # cor.test(x,y,method = "spearman")[[4]][1] #rho
# # cor.test(x,y,method = "spearman")[[3]][1] #p-value
# # 
# # summary(MB_intensity_wth_cor_df)
# # 
# # colnames(MB_intensity_wth_cor_df)
# 
# # spearman correlation : MB - weather
# names(MB_intensity_wth_cor_df)
# MB_intensity_wth_cor_df <- as.data.frame(MB_intensity_wth_cor_df)
# 
# OTU_col_num <- c(8:390)#grep("OTU", names(MB_intensity_wth_cor_df))
# wth_col_num <- c(396:422)#names(MB_intensity_wth_cor_df)[1017:1043]
# 
# 
# MB_wth_cor_result_df <- data.frame()
# imsi_MB_wth_cor_result_df <- data.frame()
# imsi_3 <- data.frame()
# for(i in c(OTU_col_num[1]:OTU_col_num[length(OTU_col_num)])){ # 18~100 are MB column
#   for(j in c(wth_col_num[1]:wth_col_num[length(wth_col_num)])){ # 101 ~ 112 are wth column
#     #i=1
#     # j=32
#     x <- MB_intensity_wth_cor_df[,i]
#     y <- MB_intensity_wth_cor_df[,j]
#     
#     imsi_MB_wth_cor_result_df[1,1] <- names(MB_intensity_wth_cor_df)[i] 
#     imsi_MB_wth_cor_result_df[1,2] <- names(MB_intensity_wth_cor_df)[j] 
#     imsi_MB_wth_cor_result_df[1,3] <- cor.test(x,y,method = "spearman")[[4]][1] #rho
#     imsi_MB_wth_cor_result_df[1,4] <- cor.test(x,y,method = "spearman")[[3]][1]
#     imsi_3 <- rbind(imsi_3, imsi_MB_wth_cor_result_df)
#   }
#   MB_wth_cor_result_df <- rbind(MB_wth_cor_result_df, imsi_3)
#   imsi_3 <- data.frame()
#   imsi_MB_wth_cor_result_df <- data.frame()
# }
# names(MB_wth_cor_result_df) <- c("MB", "wth", "rho", "p_value")
# nrow(MB_wth_cor_result_df)
# MB_wth_cor_result_df
# 
# save_version <- "v6.2"
# writexl::write_xlsx(MB_wth_cor_result_df, path = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", save_version, "/MB_wth_cor_result_",save_version,".xlsx"))
# 
# ### heatmap : correlation between MB and wth -----------------------------------
# library(ggplot2)
# # 
# # if(sum(is.na(MB_wth_cor_result_df$rho)== T) == 0){
# #   MB_wth_cor_result_df_2 <- MB_wth_cor_result_df %>% filter(p_value <= 0.05) 
# # }else{
# #   print("find NA data")
# # }
# MB_wth_cor_result_df
# 
# MB_wth_cor_result_df$label_color <- ifelse(MB_wth_cor_result_df$p_val <= 0.05, "black", "lightgray")
# MB_wth_cor_result_df_no_after_flower <- MB_wth_cor_result_df[!grepl("after_flower", MB_wth_cor_result_df$wth), ]
# 
# heatmap_plot <- ggplot(MB_wth_cor_result_df_no_after_flower, aes(x = wth, y = MB, fill = rho)) +
#   geom_tile(color = "white") +
#   # ggtitle(Y_value)+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) +
#   theme_minimal() +
#   geom_text(aes(label = round(rho,2), color = label_color), size = 3) +
#   scale_color_identity() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 20),
#         axis.text.y = element_text(color = "black", size = 15)) +
#   labs(x = "", y = "", fill = "Rho") + 
#   coord_fixed(ratio = 0.5)
# heatmap_plot
# 
# high_cor_df <- MB_wth_cor_result_df_no_after_flower %>% filter(rho >= 0.3 | rho <= -0.3)
# high_cor_heatmap_plot <- ggplot(high_cor_df, aes(x = wth, y = MB, fill = rho)) +
#   geom_tile(color = "white") +
#   # ggtitle(Y_value)+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1, 1)) +
#   theme_minimal() +
#   geom_text(aes(label = round(rho,2), color = label_color), size = 3) +
#   scale_color_identity() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 20),
#         axis.text.y = element_text(color = "black", size = 15)) +
#   labs(x = "", y = "", fill = "Rho") + 
#   coord_fixed(ratio = 0.5)
# 
# high_cor_heatmap_plot
# 
# version_heatmap <- "v6.2"
# ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", 
#                                    save_version, "/","MB_wth_heatmap_",version_heatmap, ".png")),plot = heatmap_plot, width = 50, height = 150, bg = "white",limitsize = FALSE)
# ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", 
#                                    save_version, "/","MB_wth_heatmap_high_cor_",version_heatmap, ".png")),plot = high_cor_heatmap_plot, width = 50, height = 100, bg = "white",limitsize = FALSE)
# 
# ################################################################################
# ### MB - intensity
# # spearman correlation : MB - intensity
# colnames(MB_intensity_wth_cor_df)
# MB_incidence_cor_result_df <- data.frame()
# imsi_MB_incidence_cor_result_df <- data.frame()
# imsi_4 <- data.frame()
# for(i in c(8:390)){ # 18~100 are MB column c(OTU_col_num[1]:OTU_col_num[length(OTU_col_num)])
#   for(j in c(6:7)){ # 13 ~ 14 are intensity column
#     #i=1
#     # j=32
#     x <- MB_intensity_wth_cor_df[,i]
#     y <- MB_intensity_wth_cor_df[,j]
#     
#     imsi_MB_incidence_cor_result_df[1,1] <- names(MB_intensity_wth_cor_df)[i] 
#     imsi_MB_incidence_cor_result_df[1,2] <- names(MB_intensity_wth_cor_df)[j] 
#     imsi_MB_incidence_cor_result_df[1,3] <- cor.test(x,y,method = "spearman")[[4]][1] #rho
#     imsi_MB_incidence_cor_result_df[1,4] <- cor.test(x,y,method = "spearman")[[3]][1]
#     imsi_4 <- rbind(imsi_4, imsi_MB_incidence_cor_result_df)
#   }
#   MB_incidence_cor_result_df <- rbind(MB_incidence_cor_result_df, imsi_4)
#   imsi_4 <- data.frame()
#   imsi_MB_incidence_cor_result_df <- data.frame()
# }
# names(MB_incidence_cor_result_df) <- c("MB", "intensity", "rho", "p_value")
# nrow(MB_incidence_cor_result_df)
# MB_incidence_cor_result_df
# 
# save_version <- "v6.2"
# writexl::write_xlsx(MB_incidence_cor_result_df, path = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",
#                                                               save_version, "/MB_intensity_cor_result_",save_version,".xlsx"))
# 
# ### heatmap : correlation between MB and intensity------------------------------
# #MB_incidence_cor_result_df <- MB_incidence_cor_result_df[-(which(is.na(MB_incidence_cor_result_df$rho)== T)),]
# MB_intensity_heatmap_plot <- ggplot(MB_incidence_cor_result_df,aes(x=reorder(MB, rho), y=intensity, fill=rho)) + 
#   geom_tile() +
#   xlab("\nMB") + ylab("Intensity\n") +
#   theme_bw() + 
#   theme(plot.background = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         strip.background = element_rect(fill = "white", colour = "white"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-1,1))
# MB_intensity_heatmap_plot
# 
# MB_intensity_high_cor_df <- MB_incidence_cor_result_df %>% filter(rho >= 0.3 | rho <= -0.3)
# MB_intensity_heatmap_plot_high_cor <- ggplot(MB_intensity_high_cor_df,aes(x= reorder(MB, rho), y=intensity, fill=rho)) + 
#   geom_tile() + 
#   xlab("\nMB") + ylab("Intensity\n") +
#   theme_bw() + 
#   theme(plot.background = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         strip.background = element_rect(fill = "white", colour = "white"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-1,1))
# MB_intensity_heatmap_plot_high_cor
# 
# 
# ################################################################################
# ### wth - intensity
# # spearman correlation : wth - intensity
# colnames(MB_intensity_wth_cor_df)
# wth_intensity_cor_result_df <- data.frame()
# imsi_wth_intensity_cor_result_df <- data.frame()
# imsi_5 <- data.frame()
# for(i in c(396:422)){ # 28~58 are wth column #c(wth_col_num[1]:wth_col_num[length(wth_col_num)])
#   for(j in c(6:7)){ # 6 ~ 7 are intensity column
#     #i=1
#     # j=32
#     x <- MB_intensity_wth_cor_df[,i]
#     y <- MB_intensity_wth_cor_df[,j]
#     
#     imsi_wth_intensity_cor_result_df[1,1] <- names(MB_intensity_wth_cor_df)[i] 
#     imsi_wth_intensity_cor_result_df[1,2] <- names(MB_intensity_wth_cor_df)[j] 
#     imsi_wth_intensity_cor_result_df[1,3] <- cor.test(x,y,method = "spearman")[[4]][1] #rho
#     imsi_wth_intensity_cor_result_df[1,4] <- cor.test(x,y,method = "spearman")[[3]][1]
#     imsi_5 <- rbind(imsi_5, imsi_wth_intensity_cor_result_df)
#   }
#   wth_intensity_cor_result_df <- rbind(wth_intensity_cor_result_df, imsi_5)
#   imsi_5 <- data.frame()
#   imsi_wth_intensity_cor_result_df <- data.frame()
# }
# names(wth_intensity_cor_result_df) <- c("wth", "intensity", "rho", "p_value")
# nrow(wth_intensity_cor_result_df)
# wth_intensity_cor_result_df
# 
# save_version <- "v6.2"
# writexl::write_xlsx(wth_intensity_cor_result_df, path = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/", 
#                                                                save_version, "/wth_intensity_cor_result_",save_version,".xlsx"))
# 
# ### heatmap : correlation between wth and intensity------------------------------
# wth_intensity_cor_result_df <- wth_intensity_cor_result_df
# MB_intensity_heatmap_plot <- ggplot(wth_intensity_cor_result_df,aes(y=reorder(wth, rho), x=intensity, fill=rho)) + 
#   geom_tile() + 
#   xlab("\nIntensity") + ylab("wth\n") +
#   theme_bw() + 
#   theme(plot.background = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         strip.background = element_rect(fill = "white", colour = "white"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-1,1))
# MB_intensity_heatmap_plot
# 
# MB_intensity_high_cor_df <- wth_intensity_cor_result_df %>% filter(rho >= 0.5 | rho <= -0.5)
# MB_intensity_heatmap_plot_high_cor <- ggplot(MB_intensity_high_cor_df,aes(y= reorder(wth, rho), x=intensity, fill=rho)) + 
#   geom_tile() + 
#   xlab("\nIntensity") + ylab("wth\n") +
#   theme_bw() + 
#   theme(plot.background = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         strip.background = element_rect(fill = "white", colour = "white"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "white", limits=c(-1,1))
# MB_intensity_heatmap_plot_high_cor
# 
# 
# ### bar plot : correlation -----------------------------------------------------
# ## data 
# MB_wth_cor_result_df
# MB_incidence_cor_result_df
# wth_intensity_cor_result_df
# # version_cor_plot <- "v2"
# # what_do_you_want_to_compare_with_intensity <- "MB"
# make_cor_plot2 <- function(cor_result, line_threshold, what_do_you_want_to_compare_with_intensity, version_cor_plot){
#   # cor_result = MB_incidence_cor_result_df_max_filtered 
#   # what_do_you_want_to_compare_with_intensity = "MB" 
#   # line_threshold = 0.3
#   # version_cor_plot = "max_filtered_v6.2"
#   
#   cor_result <- cor_result
#   # important_result_order <- important_result_order
#   
#   # cor_result_ordered <- cor_result[match(important_result_order, cor_result$Wth),]
#   cor_result_ordered <- cor_result
#   
#   cor_result_ordered$color <- with(cor_result_ordered, 
#                                    ifelse(p_value > 0.05, "#dbd7d7",
#                                           ifelse(rho > 0, "#d11c19", "#23346e")))
#   
#   incidence_df <- cor_result_ordered %>% filter(intensity == "incidence")
#   severity_df <- cor_result_ordered %>% filter(intensity == "severity")
#   
#   incidence_df_ordered <- incidence_df[order(incidence_df$rho, decreasing = TRUE),]
#   severity_df_ordered <- severity_df[order(severity_df$rho, decreasing = TRUE),]
#   
#   
#   # cor_result_ordered$Wth <- factor(cor_result_ordered$MB, levels = rev(important_result_order))
#   if(names(cor_result)[1] == "MB"){
#     incidence_df_ordered$MB <- factor(incidence_df_ordered$MB, levels = rev(incidence_df_ordered$MB))
#     severity_df_ordered$MB <- factor(severity_df_ordered$MB, levels = rev(severity_df_ordered$MB))
#     
#     incidence_bar_cor_plot <- ggplot(incidence_df_ordered, aes(x = MB, y = rho, fill = color)) +
#       geom_bar(stat = "identity") +
#       scale_fill_identity() +
#       theme_minimal() +
#       coord_flip() + # x축 레이블을 세로로 보기 좋게 하기 위해 사용
#       labs(title = "FHB_incidence",x = "", y = "\nRho", fill = "p_val Color") +
#       theme(axis.text.x = element_text(angle = 360)) + 
#       theme_minimal() + 
#       theme(axis.text = element_text(size=17),
#             axis.title = element_text(size=17)) + 
#       theme(axis.title.x = element_text(size = 21,hjust = 0.5, face='bold')) + 
#       theme(axis.title.y = element_text(size = 21,hjust = 0.5, face='bold')) +
#       scale_y_continuous(limits = c(-1,1)) + 
#       geom_hline(yintercept = c(line_threshold, -line_threshold), linetype = "dashed", color = "black")
#     incidence_bar_cor_plot
#     severity_bar_cor_plot <- ggplot(severity_df_ordered, aes(x = MB, y = rho, fill = color)) +
#       geom_bar(stat = "identity") +
#       scale_fill_identity() +
#       theme_minimal() +
#       coord_flip() + # x축 레이블을 세로로 보기 좋게 하기 위해 사용
#       labs(title = "FHB_severity",x = "", y = "\nRho", fill = "p_val Color") +
#       theme(axis.text.x = element_text(angle = 360)) + 
#       theme_minimal() + 
#       theme(axis.text = element_text(size=17),
#             axis.title = element_text(size=23)) + 
#       theme(axis.title.x = element_text(size = 21,hjust = 0.5, face='bold')) + 
#       theme(axis.title.y = element_text(size = 21,hjust = 0.5, face='bold')) +
#       scale_y_continuous(limits = c(-1,1)) + 
#       geom_hline(yintercept = c(line_threshold, -line_threshold), linetype = "dashed", color = "black")
#     severity_bar_cor_plot
#     
#     ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",what_do_you_want_to_compare_with_intensity,"_incidence_bar_cor_plot_",version_cor_plot,"_(Rho ",line_threshold,  ").png")),plot = incidence_bar_cor_plot, width = 7, height = 25, bg = "white", dpi = 1500)
#     ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",what_do_you_want_to_compare_with_intensity,"_severity_bar_cor_plot_",version_cor_plot,"_(Rho ",line_threshold,  ").png")),plot = severity_bar_cor_plot, width = 7, height = 25, bg = "white", dpi = 1500)
#     
#     
#     #p-val
#     incidence_df_ordered_low_p_val <- incidence_df_ordered %>% filter(p_value <= 0.05)
#     severity_df_ordered_low_p_val <- severity_df_ordered %>% filter(p_value <= 0.05)
#     
#     incidence_bar_cor_plot_pval <- ggplot(incidence_df_ordered_low_p_val, aes(x = MB, y = rho, fill = color)) +
#       geom_bar(stat = "identity") +
#       scale_fill_identity() +
#       theme_minimal() +
#       coord_flip() + # x축 레이블을 세로로 보기 좋게 하기 위해 사용
#       labs(title = "FHB_incidence",x = "", y = "\nRho", fill = "p_val Color") +
#       theme(axis.text.x = element_text(angle = 360)) + 
#       theme_minimal() + 
#       theme(axis.text = element_text(size=17),
#             axis.title = element_text(size=17)) + 
#       theme(axis.title.x = element_text(size = 21,hjust = 0.5, face='bold')) + 
#       theme(axis.title.y = element_text(size = 21,hjust = 0.5, face='bold')) +
#       scale_y_continuous(limits = c(-1,1)) + 
#       geom_hline(yintercept = c(line_threshold, -line_threshold), linetype = "dashed", color = "black")
#     incidence_bar_cor_plot_pval
#     severity_bar_cor_plot_pval <- ggplot(severity_df_ordered_low_p_val, aes(x = MB, y = rho, fill = color)) +
#       geom_bar(stat = "identity") +
#       scale_fill_identity() +
#       theme_minimal() +
#       coord_flip() + # x축 레이블을 세로로 보기 좋게 하기 위해 사용
#       labs(title = "FHB_severity",x = "", y = "\nRho", fill = "p_val Color") +
#       theme(axis.text.x = element_text(angle = 360)) + 
#       theme_minimal() + 
#       theme(axis.text = element_text(size=17),
#             axis.title = element_text(size=23)) + 
#       theme(axis.title.x = element_text(size = 21,hjust = 0.5, face='bold')) + 
#       theme(axis.title.y = element_text(size = 21,hjust = 0.5, face='bold')) +
#       scale_y_continuous(limits = c(-1,1)) + 
#       geom_hline(yintercept = c(line_threshold, -line_threshold), linetype = "dashed", color = "black")
#     severity_bar_cor_plot_pval
#     
#     ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",what_do_you_want_to_compare_with_intensity,"_incidence_bar_cor_plot_low_pval_",version_cor_plot,"_(Rho ",line_threshold,  ").png")),plot = incidence_bar_cor_plot_pval, width = 7, height = 14, bg = "white", dpi = 1500)
#     ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",what_do_you_want_to_compare_with_intensity,"_severity_bar_cor_plot_low_pval_",version_cor_plot,"_(Rho ",line_threshold,  ").png")),plot = severity_bar_cor_plot_pval, width = 7, height = 14, bg = "white", dpi = 1500)
#     
#     
#     
#   }else if(names(cor_result)[1] == "wth"){
#     incidence_df_ordered$wth <- factor(incidence_df_ordered$wth, levels = rev(incidence_df_ordered$wth))
#     severity_df_ordered$wth <- factor(severity_df_ordered$wth, levels = rev(severity_df_ordered$wth))
#     
#     
#     incidence_bar_cor_plot <- ggplot(incidence_df_ordered, aes(x = wth, y = rho, fill = color)) +
#       geom_bar(stat = "identity") +
#       scale_fill_identity() +
#       theme_minimal() +
#       coord_flip() + # x축 레이블을 세로로 보기 좋게 하기 위해 사용
#       labs(title = "FHB_incidence",x = "", y = "\nRho", fill = "p_val Color") +
#       theme(axis.text.x = element_text(angle = 360)) + 
#       theme_minimal() + 
#       theme(axis.text = element_text(size=17),
#             axis.title = element_text(size=17)) + 
#       theme(axis.title.x = element_text(size = 21,hjust = 0.5, face='bold')) + 
#       theme(axis.title.y = element_text(size = 21,hjust = 0.5, face='bold')) +
#       scale_y_continuous(limits = c(-1,1)) + 
#       geom_hline(yintercept = c(line_threshold, -line_threshold), linetype = "dashed", color = "black")
#     incidence_bar_cor_plot
#     severity_bar_cor_plot <- ggplot(severity_df_ordered, aes(x = wth, y = rho, fill = color)) +
#       geom_bar(stat = "identity") +
#       scale_fill_identity() +
#       theme_minimal() +
#       coord_flip() + # x축 레이블을 세로로 보기 좋게 하기 위해 사용
#       labs(title = "FHB_severity",x = "", y = "\nRho", fill = "p_val Color") +
#       theme(axis.text.x = element_text(angle = 360)) + 
#       theme_minimal() + 
#       theme(axis.text = element_text(size=17),
#             axis.title = element_text(size=17)) + 
#       theme(axis.title.x = element_text(size = 21,hjust = 0.5, face='bold')) + 
#       theme(axis.title.y = element_text(size = 21,hjust = 0.5, face='bold')) +
#       scale_y_continuous(limits = c(-1,1)) + 
#       geom_hline(yintercept = c(line_threshold, -line_threshold), linetype = "dashed", color = "black")
#     severity_bar_cor_plot
#     
#     ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",what_do_you_want_to_compare_with_intensity,"_incidence_bar_cor_plot_",version_cor_plot,"_(Rho ",line_threshold,  ").png")),plot = incidence_bar_cor_plot, width = 7, height = 14, bg = "white", dpi = 1500)
#     ggsave(filename = file.path(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",what_do_you_want_to_compare_with_intensity,"_severity_bar_cor_plot_",version_cor_plot,"_(Rho ",line_threshold,  ").png")),plot = severity_bar_cor_plot, width = 7, height = 14, bg = "white", dpi = 1500)
#   }
#   
#   # version_cor_plot <- "v3"
#   
#   
# }
# 
# ## check save version
# make_cor_plot2(cor_result = MB_incidence_cor_result_df, 
#                what_do_you_want_to_compare_with_intensity = "MB", 
#                line_threshold = 0.4,
#                version_cor_plot = "v6.2")
# 
# make_cor_plot2(cor_result = wth_intensity_cor_result_df,
#                what_do_you_want_to_compare_with_intensity = "wth",
#                line_threshold = 0.4,
#                version_cor_plot = "v6.2")
# 
# ### filter based on relative abundance (max value >= 0.01) -------------------------------
# MB_for_filter_df <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_more_ITSfull_v3.2.xlsx")
# # MB_for_filter_df <- MB_for_filter_df %>% dplyr::select(!contains("AF"))
# # names(MB_for_filter_df)
# # 
# # MB_Taxa <- MB_for_filter_df$Taxa
# # MB_for_filter_df$Taxa <- NULL
# # 
# # # calculate prevalence
# # n_samples <- ncol(MB_for_filter_df)
# # prevalence <- rowSums(MB_for_filter_df > 0) / n_samples   # 각 OTU가 나타난 비율
# # 
# # # calculate relative abundance
# # rel_abund <- sweep(MB_for_filter_df, 2, colSums(MB_for_filter_df), "/")
# # 
# # # average of each OTU's abundance
# # max_abund <- apply(rel_abund, 1, max)
# # 
# # # filtering
# # otu_keep <- (prevalence >= 0.5) & (max_abund >= 0.01)
# # OTU_kept <- rownames(MB_for_filter_df)[otu_keep]
# # OTU_removed <- rownames(MB_for_filter_df)[!otu_keep]
# # MB_filtered <- MB_for_filter_df[otu_keep, ]
# # 
# # MB_Taxa[as.numeric(OTU_kept)]
# # 
# # # make data frame
# # MB_filter_df_all_results <- data.frame(
# #   OTU <- MB_Taxa,
# #   prevalence <- prevalence,
# #   max_val <- max_abund
# # )
# 
# 
# MB_filter_df_results <- MB_for_filter_df # MB_filter_df_all_results %>% filter(prevalence >= 0.5) %>% filter(max_val >= 0.01)
# # names(MB_filter_df_all_results) <- c("OTU", "prevalence", "max_val")
# MB_filter_df_results <- MB_filter_df_results[-which(MB_filter_df_results$Taxa == "unclassified"),]
# ### correlation plot -----------------------------------------------------
# MB_incidence_cor_result_df_max_filtered <- MB_incidence_cor_result_df %>% filter(MB %in% MB_filter_df_results$Taxa)
# 
# make_cor_plot2(cor_result = MB_incidence_cor_result_df_max_filtered, 
#                what_do_you_want_to_compare_with_intensity = "MB", 
#                line_threshold = 0.3,
#                version_cor_plot = "max_filtered_v6.2")
# 
# make_cor_plot2(cor_result = MB_incidence_cor_result_df_max_filtered, 
#                what_do_you_want_to_compare_with_intensity = "MB", 
#                line_threshold = 0.281,
#                version_cor_plot = "max_filtered_v6.2")
# 
# library(pwr)
# pwr.r.test(
#   n = 96,
#   # r = 0.30,
#   sig.level = 0.05,   # 네가 쓰려는 기준
#   power = 0.8,
#   alternative = "two.sided"
# )
# 
# pwr.r.test(
#   n = 96,
#   # r = 0.40,
#   sig.level = 0.05,   # 네가 쓰려는 기준
#   power = 0.8,
#   alternative = "two.sided"
# )
# 
# pwr.r.test(
#   n = 96,
#   r = 0.40,
#   sig.level = 0.005,
#   alternative = "two.sided"
# )
# 
# MB_incidence_cor_result_df_max_filtered %>% filter(intensity == "incidence")
# # # save
# # save_version
# # writexl::write_xlsx(MB_filter_df_all_results, path = paste0("E:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",
# #                                                             save_version, "/MB_filter_df_all_results_", save_version, ".xlsx"))
# # 
# # writexl::write_xlsx(MB_filter_df_results, path = paste0("E:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/",
# #                                                         save_version, "/MB_filter_df_results_", save_version, ".xlsx"))

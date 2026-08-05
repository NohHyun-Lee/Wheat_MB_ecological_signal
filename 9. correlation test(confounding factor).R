library(dplyr)
library(readxl)
library(ggplot2)

### data 
# GIBSIM_daily_res_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/1.Run_GIBSIM/GIBSIM_daily_results/MB_v3.3/"
# GIBSIM_daily_res_filelist <- list.files(GIBSIM_daily_res_dir)
df_for_more_analysis <- read_xlsx("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/df_for_more_analysis_good_20_bad_100.xlsx")
names(df_for_more_analysis)
df_for_CI <- df_for_more_analysis %>% dplyr::select(loc, wth_ID, incidence, 
                                                    Alternaria, Epicoccum, Hannaella, Periconia, 
                                                    Cladosporium, Papiliotrema, 
                                                    temp_before_hd_sampling_10days, rhum_before_hd_sampling_10days, prcp_before_hd_sampling_10days, 
                                                    temp_hd_sampling_10days, rhum_hd_sampling_10days, prcp_hd_sampling_10days, 
                                                    temp_flower_sampling_10days, rhum_flower_sampling_10days, prcp_flower_sampling_10days,
                                                    temp_after_flower_sampling_10days, rhum_after_flower_sampling_10days, prcp_after_flower_sampling_10days,
                                                    temp_after_sampling, rhum_after_sampling, prcp_after_sampling
)
writexl::write_xlsx(df_for_CI,
                    "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/data_for_causal_inference_2025_v12.xlsx")
### 
wth_var <- c("temp_before_hd_sampling_10days",
             "rhum_before_hd_sampling_10days",
             "prcp_before_hd_sampling_10days",
             "temp_hd_sampling_10days",
             "rhum_hd_sampling_10days",
             "prcp_hd_sampling_10days",
             "temp_flower_sampling_10days",
             "rhum_flower_sampling_10days",
             "prcp_flower_sampling_10days",
             "temp_after_flower_sampling_10days",
             "rhum_after_flower_sampling_10days",
             "prcp_after_flower_sampling_10days",
             "temp_after_sampling",
             "rhum_after_sampling",
             "prcp_after_sampling")
inc_and_m_var <- c("incidence","Alternaria", "Epicoccum", "Hannaella", "Periconia", "Cladosporium", "Papiliotrema")

wth_cor_res_df <- data.frame()
df_for_CI <- as.data.frame(df_for_CI)
for(im in c(1:length(inc_and_m_var))){
  # im = 1
  imsi_target_im <- inc_and_m_var[im]
  for(w in c(1:length(wth_var))){
    # w = 1
    imsi_target_w <- wth_var[w]
    imsi_cor <- stats::cor.test(df_for_CI[,imsi_target_im], df_for_CI[,imsi_target_w], method = "spearman")
    wth_cor_res_df_raw <- as.data.frame(matrix(nrow = 1, ncol = 4))  
    wth_cor_res_df_raw[1,1] <- inc_and_m_var[im]
    wth_cor_res_df_raw[1,2] <- wth_var[w]
    wth_cor_res_df_raw[1,3] <- as.vector(imsi_cor$estimate)
    wth_cor_res_df_raw[1,4] <- imsi_cor$p.value
    
    wth_cor_res_df <- rbind(wth_cor_res_df, wth_cor_res_df_raw)
  }
}
names(wth_cor_res_df) <- c("inc_or_M", "wth", "rho", "p_val")
wth_cor_res_df$p_val <- round(wth_cor_res_df$p_val, 3)

# cor graph
wth_cor_res_df2 <- wth_cor_res_df
wth_cor_res_df2 <- wth_cor_res_df %>%
  mutate(
    rho_plot = ifelse(p_val <= 0.05, rho, NA)
  )

wth_cor_res_df2 <- wth_cor_res_df %>%
  mutate(
    rho_plot = ifelse(p_val <= 0.05, rho, NA),
    sig = case_when(
      p_val <= 0.001 ~ "***",
      p_val <= 0.01 ~ "**",
      p_val <= 0.05 ~ "*",
      TRUE ~ ""
    )
  )
wth_cor_res_df2$wth <- factor(wth_cor_res_df2$wth,levels = wth_var)
wth_cor_res_df2$inc_or_M <- factor(wth_cor_res_df2$inc_or_M,levels = inc_and_m_var)

p_all <- ggplot(
  wth_cor_res_df2,
  aes(
    x = inc_or_M,
    y = wth
  )
) +
  # grid background
  geom_tile(
    fill = "white",
    color = "grey85",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      fill = rho_plot,
      size = abs(rho)
    ),
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-1,1),
    na.value = "#D9D9D9",
    name = "Rho"
  ) +
  
  scale_size(
    range = c(2,10),
    name = "|Rho|"
  ) +
  
  # geom_text(
  #   aes(label = sig),
  #   size = 4,
  #   vjust = 0.5,
  #   hjust = 0.5,
  #   color = "black"
  # ) +
  geom_text(
    aes(label = ifelse(p_val <= 0.05, round(rho, 2), "")),
    size = 2
  ) + 
  theme_minimal(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.title = element_blank(),
    
    legend.position = "right",
    
    panel.grid = element_blank()
  ) +
  
  scale_y_discrete(labels = c(
    "temp_before_hd_sampling_10days" = "temp_before_hd",
    "rhum_before_hd_sampling_10days" = "rhum_before_hd",
    "prcp_before_hd_sampling_10days" = "prcp_before_hd",
    "temp_hd_sampling_10days" = "temp_hd",
    "rhum_hd_sampling_10days" = "rhum_hd",
    "prcp_hd_sampling_10days" = "prcp_hd",
    "temp_flower_sampling_10days" = "temp_flower",
    "rhum_flower_sampling_10days" = "rhum_flower",
    "prcp_flower_sampling_10days" = "prcp_flower",
    "temp_after_flower_sampling_10days" = "temp_after_flower",
    "rhum_after_flower_sampling_10days" = "rhum_after_flower",
    "prcp_after_flower_sampling_10days" = "prcp_after_flower",
    "temp_after_sampling" = "temp_after_sampling",
    "rhum_after_sampling" = "rhum_after_sampling",
    "prcp_after_sampling" = "prcp_after_sampling"
  ))

p_all

wth_cor_res_df2$wth
library(dplyr)
library(stringr)

wth_cor_res_df_af_flower <- wth_cor_res_df2 %>%
  dplyr::filter(
    str_detect(wth, "flower_sampling") |
      str_detect(wth, "after_flower_sampling") |
      str_detect(wth, "_after_sampling")
  )

p_after_flower <- ggplot(
  wth_cor_res_df_af_flower,
  aes(
    x = inc_or_M,
    y = wth
  )
) +
  # grid background
  geom_tile(
    fill = "white",
    color = "grey85",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      fill = rho_plot,
      size = abs(rho)
    ),
    shape = 21,
    color = "black",
    stroke = 0.3
  ) +
  
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-1,1),
    na.value = "#D9D9D9",
    name = "Rho"
  ) +
  
  scale_size(
    range = c(2,10),
    name = "|Rho|"
  ) +
  
  # geom_text(
  #   aes(label = sig),
  #   size = 4,
  #   vjust = 0.5,
  #   hjust = 0.5,
  #   color = "black"
  # ) +
  geom_text(
    aes(label = ifelse(p_val <= 0.05, round(rho, 2), "")),
    size = 2
  ) + 
  theme_minimal(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.title = element_blank(),
    
    legend.position = "right",
    
    panel.grid = element_blank()
  ) +
  
  scale_y_discrete(labels = c(
    "temp_flower_sampling_10days" = "temp_flower",
    "rhum_flower_sampling_10days" = "rhum_flower",
    "prcp_flower_sampling_10days" = "prcp_flower",
    "temp_after_flower_sampling_10days" = "temp_after_flower",
    "rhum_after_flower_sampling_10days" = "rhum_after_flower",
    "prcp_after_flower_sampling_10days" = "prcp_after_flower",
    "temp_after_sampling" = "temp_after_sampling",
    "rhum_after_sampling" = "rhum_after_sampling",
    "prcp_after_sampling" = "prcp_after_sampling"
  ))

p_after_flower




cor_version <- "v1"
ggsave(filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/", 
                         version_pred_obs, 
                         "/causal_inference/correlation_plot_",
                         cor_version, ".png"),
       plot = p_all,
       width = 8,
       height = 7
)

ggsave(filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/", 
                         version_pred_obs, 
                         "/causal_inference/correlation_plot_only_after_sampling_",
                         cor_version, ".png"),
       plot = p_after_flower,
       width = 6,
       height = 5
)


#temp_before_hd_sampling_10days
library(tidyverse)
vars_microbe <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia", "incidence")
vars_weather <- c("temp_before_hd_sampling_10days",
                  "rhum_before_hd_sampling_10days")
df_for_CI2 <- df_for_CI
df_for_CI2$incidence <- df_for_CI2$incidence*100
df_long_sc <- df_for_CI2 %>%
  select(all_of(c(vars_microbe, vars_weather))) %>%
  pivot_longer(cols = all_of(vars_microbe),
               names_to = "microbe",
               values_to = "value") %>%
  pivot_longer(cols = all_of(vars_weather),
               names_to = "weather",
               values_to = "weather_value")

plot_scatter <- function(data, microbe_name, weather_name) {
  
  df_sub <- data %>%
    filter(microbe == microbe_name,
           weather == weather_name)
  
  # correlation °è»ê
  cor_test <- cor.test(df_sub$value, df_sub$weather_value, method = "spearman")
  
  rho <- round(cor_test$estimate, 2)
  pval <- signif(cor_test$p.value, 2)
  
  ggplot(df_sub, aes(x = value, y = weather_value)) +
    geom_point(size = 1, alpha = 0.7, color = "black") +
    geom_smooth(method = "lm",
                color = "#d73027",
                fill = "grey",
                alpha = 0.5) +
    
    theme_classic() +
    labs(
      x = paste0("\n",microbe_name),
      y = paste0(weather_name,"\n")
      # title = paste0(microbe_name, " vs ", weather_name)
    )    +
    theme(
      plot.title = element_text(size = 6, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 10, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                                 size = 12, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 12, face = 'bold', color = 'black')
    ) +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("Rho = ", rho, "\np = ", pval),
             hjust = 1.1, vjust = 1.5,
             size = 3)
}

if(!dir.exists(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference",
                      "/scatter_plots/"))){
  dir.create(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference",
                    "/scatter_plots/"))
}


for (m in vars_microbe) {
  for (w in vars_weather) {
    
    p <- plot_scatter(df_long_sc, m, w)
    
    ggsave(
      filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference",
                        "/scatter_plots/scatter_plots_", w, "_vs_", m, ".png"),
      plot = p,
      width = 2.5,
      height = 2.5,
      dpi = 300
    )
  }
}


### histogram
df_hist_temp <- df_for_CI2 %>%
  select(temp_before_hd_sampling_10days,
         temp_hd_sampling_10days) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value")

hist_temp_g <- ggplot(df_hist, aes(x = value, fill = variable)) +
  geom_histogram(aes(y = ..density..),
                 alpha = 0.3,
                 position = "identity",
                 bins = 30) +
  geom_density(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "\nTemp",
    y = "Density\n",
    fill = ""
  ) + scale_fill_manual(
    values = c(
      "temp_before_hd_sampling_10days" = "#FFA339",
      "temp_hd_sampling_10days" = "#FD3F3F"
    )
  )+
  theme(
    plot.title = element_text(size = 6, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 12, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 12, face = 'bold', color = 'black'),
    legend.position = "bottom"
  ) 


df_hist_rhum <- df_for_CI2 %>%
  select(rhum_before_hd_sampling_10days,
         rhum_hd_sampling_10days) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value")

hist_rhum_g <- ggplot(df_hist_rhum, aes(x = value, fill = variable)) +
  geom_histogram(aes(y = ..density..),
                 alpha = 0.3,
                 position = "identity",
                 bins = 30) +
  geom_density(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "\nRhum",
    y = "Density\n",
    fill = ""
  ) + scale_fill_manual(
    values = c(
      "rhum_before_hd_sampling_10days" = "#9A9400",
      "rhum_hd_sampling_10days" = "#006666"
    )
  )+
  theme(
    plot.title = element_text(size = 6, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 12, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 12, face = 'bold', color = 'black'),
    legend.position = "bottom"
  ) 
hist_rhum_g

df_hist_prcp <- df_for_CI2 %>%
  select(prcp_before_hd_sampling_10days,
         prcp_hd_sampling_10days) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value")

hist_prcp_g <- ggplot(df_hist_prcp, aes(x = value, fill = variable)) +
  geom_histogram(aes(y = ..density..),
                 alpha = 0.3,
                 position = "identity",
                 bins = 30) +
  geom_density(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "\nPrcp",
    y = "Density\n",
    fill = ""
  ) + scale_fill_manual(
    values = c(
      "prcp_before_hd_sampling_10days" = "#3596EB",
      "prcp_hd_sampling_10days" = "#202020"
    )
  )+
  theme(
    plot.title = element_text(size = 6, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 12, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 12, face = 'bold', color = 'black'),
    legend.position = "bottom"
  ) 
hist_prcp_g


ggsave(
  filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference",
                    "/hist_plots/hist_temp.png"),
  plot = hist_temp_g,
  width = 4,
  height = 4,
  dpi = 300
)

ggsave(
  filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference",
                    "/hist_plots/hist_rhum", ".png"),
  plot = hist_rhum_g,
  width = 4,
  height = 4,
  dpi = 300
)

ggsave(
  filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference",
                    "/hist_plots/hist_prcp", ".png"),
  plot = hist_prcp_g,
  width = 4,
  height = 4,
  dpi = 300
)
hist(df_for_CI2$temp_before_hd_sampling_10days)
hist(df_for_CI2$temp_hd_sampling_10days)

hist(df_for_CI2$rhum_before_hd_sampling_10days)
hist(df_for_CI2$rhum_hd_sampling_10days)

hist(df_for_CI2$prcp_before_hd_sampling_10days)
hist(df_for_CI2$prcp_hd_sampling_10days)


# p_r <- ggplot(wth_cor_res_df2,
#             aes(x = inc_or_M,
#                 y = wth,
#                 fill = rho_plot)) +
#   
#   geom_tile(color = "white", linewidth = 0.5) +
#   
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = c(-1, 1),
#     na.value = "#E1E1E1",
#     name = "Spearman rho"
#   ) +
#   geom_text(
#     aes(label = ifelse(p_val <= 0.05, round(rho, 2), "")),
#     size = 4
#   ) + 
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold"),
#     axis.text.y = element_text(size = 12, face = "bold"),
#     axis.title = element_blank()
#   ) + 
#   scale_y_discrete(labels = c(
#     "temp_before_hd_sampling_10days" = "temp_before_hd",
#     "rhum_before_hd_sampling_10days" = "rhum_before_hd",
#     "prcp_before_hd_sampling_10days" = "prcp_before_hd",
#     "temp_hd_sampling_10days" = "temp_hd",
#     "rhum_hd_sampling_10days" = "rhum_hd",
#     "prcp_hd_sampling_10days" = "prcp_hd",
#     "temp_flower_sampling_10days" = "temp_flower",
#     "rhum_flower_sampling_10days" = "rhum_flower",
#     "prcp_flower_sampling_10days" = "prcp_flower",
#     "temp_after_flower_sampling_10days" = "temp_after_flower",
#     "rhum_after_flower_sampling_10days" = "rhum_after_flower",
#     "prcp_after_flower_sampling_10days" = "prcp_after_flower",
#     "temp_after_sampling" = "temp_after_sampling",
#     "rhum_after_sampling" = "rhum_after_sampling",
#     "prcp_after_sampling" = "prcp_after_sampling"
#   ))
# 
# p_r
# ggsave(filename = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/", 
#                          version_pred_obs, 
#                          "/causal_inference/correlation_plot_Rectangle_",
#                          cor_version, ".png"),
#        plot = p_r,
#        width = 8,
#        height = 7
# )

# ### bipartite correlation network
# 
# cor_df <- wth_cor_res_df %>% dplyr::filter(p_val <= 0.05)
# cor_df <- cor_df %>%
#   mutate(
#     node_type = case_when(
#       inc_or_M == "incidence" ~ "incidence",
#       TRUE ~ "microbe"
#     )
#   )
# edges <- cor_df %>%
#   transmute(
#     from = wth,
#     to = inc_or_M,
#     rho = rho
#   )
# nodes <- data.frame(
#   name = unique(c(edges$from, edges$to))
# )
# 
# nodes$type <- case_when(
#   nodes$name %in% unique(wth_cor_res_df$wth) ~ "weather",
#   nodes$name == "incidence" ~ "incidence",
#   TRUE ~ "microbe"
# )
# library(igraph)
# library(ggraph)
# 
# graph <- graph_from_data_frame(edges, vertices = nodes)
# 
# ggraph(graph, layout = "fr") +
#   
#   geom_edge_link(aes(color = rho,
#                      width = abs(rho)),
#                  alpha = 0.8) +
#   
#   scale_edge_color_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     limits = c(-1,1)
#   ) +
#   
#   scale_edge_width(range = c(0.5,2)) +
#   
#   geom_node_point(aes(color = type),
#                   size = 5) +
#   
#   scale_color_manual(values = c(
#     weather = "#2e8b57",
#     microbe = "#ffa500",
#     incidence = "#404040"
#   )) +
#   
#   geom_node_text(aes(label = name),
#                  repel = TRUE,
#                  size = 4) +
#   
#   theme_void()
# 











cor.test(df_for_CI$temp_flower_sampling_10days, df_for_CI$temp_after_flower_sampling_10days)


stats::cor(df_for_CI$incidence, df_for_CI$temp_before_hd_sampling_10days) #
stats::cor(df_for_CI$incidence, df_for_CI$rhum_before_hd_sampling_10days) #
stats::cor(df_for_CI$incidence, df_for_CI$prcp_before_hd_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$temp_hd_sampling_10days) #
stats::cor(df_for_CI$incidence, df_for_CI$rhum_hd_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$prcp_hd_sampling_10days)
print("ss")
stats::cor(df_for_CI$incidence, df_for_CI$temp_flower_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$rhum_flower_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$prcp_flower_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$temp_after_flower_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$rhum_after_flower_sampling_10days)
stats::cor(df_for_CI$incidence, df_for_CI$prcp_after_flower_sampling_10days)

stats::cor(df_for_CI$incidence, df_for_CI$temp_after_sampling)
stats::cor(df_for_CI$incidence, df_for_CI$rhum_after_sampling)
stats::cor(df_for_CI$incidence, df_for_CI$prcp_after_sampling)


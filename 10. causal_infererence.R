### library
library(readxl)
library(brms)
library(ggplot2)
library(dplyr)
### data
df_for_more_analysis <- read_xlsx("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/df_for_more_analysis_good_20_bad_100.xlsx")
MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_more_ITSfull_v3.2.xlsx")
MB_max_filtered_data <- MB_max_filtered_data[-which(MB_max_filtered_data$Taxa == "unclassified"),]
MB_max_filtered_microbes <- MB_max_filtered_data$Taxa
names(df_for_more_analysis)
df_for_CI_selected <- df_for_more_analysis[,c("incidence", MB_max_filtered_microbes, 
                                              "GZ_mean_max",
                                              "temp_before_hd_sampling_10days", "rhum_before_hd_sampling_10days",
                                              "temp_flower_sampling_10days", "rhum_flower_sampling_10days", "prcp_after_sampling")]


head(as.data.frame(df_for_CI_selected))


### data preprocessing
# change 0 to 0.001 (because of beta distridution)
min(df_for_CI_selected$incidence[which(df_for_CI_selected$incidence!= 0)])
df_for_CI_selected$incidence <- ifelse(df_for_CI_selected$incidence == 0, 0.00001, df_for_CI_selected$incidence)
df_for_CI_selected$incidence <- ifelse(df_for_CI_selected$incidence == 1, 0.9999, df_for_CI_selected$incidence)
if(max(df_for_CI_selected$incidence) > 1.1){
  df_for_CI_selected$incidence <- df_for_CI_selected$incidence / 100
  print(df_for_CI_selected$incidence)
}else{
  print(df_for_CI_selected$incidence)
}

df_scaled <- df_for_CI_selected
vars_for_scale <- c(
  "temp_before_hd_sampling_10days",
  "rhum_before_hd_sampling_10days",
  "GZ_mean_max",
  "temp_flower_sampling_10days",
  "rhum_flower_sampling_10days" ,
  "prcp_after_sampling"
)

df_scaled[vars_for_scale] <- scale(df_scaled[vars_for_scale])

# for loop -----------------------------------------------------------------------------------------------------
library(brms)
print("1st")
# 분석할 미생물 변수들
microbes <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia")

#######====================================================================================================================================================
imsi_save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference/ITE_ATE_CATE"
ITE_save_version <- "v2.5.2"
ITE_save_version2 <- "GZ_interaction_v2.5.2"
high_weather_quantile <- 0.75
low_weather_quantile <- 0.25
microbes_q <- seq(0.1,0.9,0.1)
for(ITE_save_version in c(ITE_save_version)){
  print(ITE_save_version)
  print(ITE_save_version2)
  results <- list()
  ATE_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
  CATE_high_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
  CATE_low_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
  
  for (mi in c(1:length(microbes))) {
    # microbe <- microbes[1]
    
    microbe <- microbes[mi]
    print(microbe)

    if(microbe == "Alternaria"){
      formula_str <- as.formula(
        paste0("incidence ~ ", 
               "s(",microbe,", k=3) +",
               " + 
              t2(Alternaria,Papiliotrema, k=4) + 
              t2(Alternaria,Periconia, k=4) + 
              s(GZ_mean_max, k=3) +
              s(temp_flower_sampling_10days, k=3) +
              
              s(rhum_flower_sampling_10days, k=3)")
      )
    }else if(microbe == "Epicoccum"){
      formula_str <- as.formula(
        paste0("incidence ~ ",
               "s(",microbe,", k=3) +",
               " + 
              t2(Epicoccum, Cladosporium, k=4) + 
              s(GZ_mean_max, k=3) +
              s(temp_flower_sampling_10days, k=3) +
              
              s(rhum_flower_sampling_10days, k=3)")
      )
    }else if(microbe == "Periconia"){
      formula_str <- as.formula(
        paste0("incidence ~ ",
               "s(",microbe,", k=3) +",
               " + 
              t2(Periconia, Alternaria, k=4) + 
              s(GZ_mean_max, k=3) +
              s(temp_flower_sampling_10days, k=3) +
              
              s(rhum_flower_sampling_10days, k=3)")
      )
    }else{
      formula_str <- as.formula(
        paste0("incidence ~ ",
               
               "s(",microbe,", k=3) +",
               
               "s(GZ_mean_max, k=3) +",
               
               "s(temp_flower_sampling_10days, k=3) +
              
              s(rhum_flower_sampling_10days, k=3)")
      )
      
    }
    
    
    
    
    model <- brm(
      formula = formula_str,
      family = Beta(),
      data = df_scaled,
      prior = c(
        prior(normal(0,1), class="b"),
        prior(exponential(1), class="sds")
      ),
      chains = 4,
      iter = 6000,
      warmup = 3000,
      cores = 4,
      control = list(adapt_delta=0.995, max_treedepth=15),
      silent = 2
    )
    
    cat(paste0(microbe, " model finished"))
    
    ##### ATE #####
    
    new_high <- df_scaled
    new_low  <- df_scaled
    
    new_high[[microbe]] <- quantile(df_scaled[[microbe]],0.9)
    new_low[[microbe]]  <- quantile(df_scaled[[microbe]],0.1)
    
    pred_high <- posterior_epred(model,newdata=new_high)
    pred_low  <- posterior_epred(model,newdata=new_low)
    
    ATE_draws <- pred_high - pred_low
    ATE_draw_mean <- rowMeans(ATE_draws)
    ITE_draw_mean <- colMeans(ATE_draws)
    
    # save ITE
    pred_high_df <- as.data.frame(pred_high)
    pred_low_df <- as.data.frame(pred_low)
    pred_high_minus_low <- pred_high_df-pred_low_df
    
    names(pred_high_minus_low) <- names(pred_low_df) <- names(pred_high_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
    
    if(!dir.exists(paste0(imsi_save_dir, "/", ITE_save_version))){
      dir.create(paste0(imsi_save_dir, "/", ITE_save_version))
    }
    
    writexl::write_xlsx(pred_low_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_low_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_high_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_high_",microbe, "_",ITE_save_version2,".xlsx"))
    writexl::write_xlsx(pred_high_minus_low, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_high_minus_low_",microbe, "_",ITE_save_version2,".xlsx"))
    
    print(paste0(microbe, " ITE finished"))
    
    # # ATE 계산
    ATE <- colMeans(pred_high - pred_low)

    ATE_results <- cbind(ATE_results, ATE)
    names(ATE_results)[mi] <- paste0(microbe, "_ATE")

    print(paste0(microbe, " ATE finished"))
    
    
    results[[microbe]] <- list(
      
      model=model,
      
      ATE_draw=ATE_draw_mean,
      ITE_draw_mean = ITE_draw_mean
    )
    
    cat("// Done:", microbe, "\n")
    
  }
  writexl::write_xlsx(ATE_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ATE_results_",ITE_save_version2,".xlsx"))
  
}

ATE_plot_df <- data.frame()
for(microbe in microbes){
  ATE <- results[[microbe]]$ITE_draw_mean

  tmp <- data.frame(
    microbe = microbe,
    effect_type = c("ATE"), 
    mean = c(mean(ATE)),
    lower = c(quantile(ATE,0.025)),
    upper = c(quantile(ATE,0.975))
  )
  
  ATE_plot_df <- rbind(ATE_plot_df,tmp)
}
ATE_plot_df$mean <- ATE_plot_df$mean*100
ATE_plot_df$lower <- ATE_plot_df$lower*100
ATE_plot_df$upper <- ATE_plot_df$upper*100
ATE_plot_df$color_flag <- ifelse(ATE_plot_df$upper > 0, "pos", "default")
ATE_posterior_g <- ggplot(ATE_plot_df,
                          aes(x = mean,
                              y = microbe,
                              color = color_flag)) +
  
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey40") +
  
  geom_errorbarh(aes(xmin = lower,
                     xmax = upper),
                 height = 0.2,
                 size = 1.2) +
  
  geom_point(size = 4, shape = "x") +
  
  scale_color_manual(values = c(
    default = "#004C99",
    pos = "grey60"
  )) +
  
  labs(
    x = "\nEffect on FHB incidence (%)",
    y = "",
    color = ""
  ) +
  
  theme_bw(base_size = 15) + 
  theme(
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  )
version_pred_obs <- "2025_v12"
ggsave(
  plot = ATE_posterior_g,
  file = paste0(
    "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/causal_inference/ITE_ATE_CATE/",
    ITE_save_version,
    "/graph/",
    "ATE_posterior_",
    ITE_save_version2,
    ".png"
  ),
  width = 20,
  height = 10,
  units = c("cm")
)
getwd()

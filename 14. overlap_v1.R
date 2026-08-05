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

df_for_CI_selected_imsi <- df_for_CI_selected %>% dplyr::select("incidence", "Alternaria", "temp_flower_sampling_10days", "rhum_flower_sampling_10days", "GZ_mean_max")
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

# quantile(df_scaled[["Epicoccum"]],0.1)
# quantile(df_scaled[["Epicoccum"]],0.9)
# 
# quantile(df_scaled[["Alternaria"]],0.1)
# quantile(df_scaled[["Alternaria"]],0.9)
# 
# quantile(df_scaled[["Periconia"]],0.1)
# quantile(df_scaled[["Periconia"]],0.9)
# 
# quantile(df_scaled[["Hannaella"]],0.1)
# quantile(df_scaled[["Hannaella"]],0.9)
# 
# hist(df_scaled[["Alternaria"]])
# hist(df_scaled[["Epicoccum"]])
# hist(df_scaled[["Hannaella"]])
# hist(df_scaled[["Periconia"]])

### modelinglibrary(brms)
# detach("package:conflicted", unload = TRUE)
# prior <- prior(horseshoe(), class = "b")

### overlap
df_for_CI_selected2 <- df_for_CI_selected
quantile(df_for_CI_selected2$Alternaria, 0.2);quantile(df_for_CI_selected2$Alternaria, 0.8)
plot(df_for_CI_selected2$Alternaria, df_for_CI_selected2$temp_before_hd_sampling_10days)
plot(df_for_CI_selected2$Alternaria, df_for_CI_selected2$rhum_before_hd_sampling_10days)

quantile(df_for_CI_selected2$Epicoccum, 0.2);quantile(df_for_CI_selected2$Epicoccum, 0.8)
plot(df_for_CI_selected2$Epicoccum, df_for_CI_selected2$temp_before_hd_sampling_10days)
plot(df_for_CI_selected2$Epicoccum, df_for_CI_selected2$rhum_before_hd_sampling_10days)

quantile(df_for_CI_selected2$Periconia, 0.2);quantile(df_for_CI_selected2$Periconia, 0.8)
plot(df_for_CI_selected2$Periconia, df_for_CI_selected2$temp_before_hd_sampling_10days)
plot(df_for_CI_selected2$Periconia, df_for_CI_selected2$rhum_before_hd_sampling_10days)

quantile(df_for_CI_selected2$Hannaella, 0.2);quantile(df_for_CI_selected2$Hannaella, 0.8)
plot(df_for_CI_selected2$Hannaella, df_for_CI_selected2$temp_before_hd_sampling_10days)
plot(df_for_CI_selected2$Hannaella, df_for_CI_selected2$rhum_before_hd_sampling_10days)

df_for_CI_selected2$A_group <- ifelse(df_for_CI_selected2$Hannaella < mean(df_for_CI_selected2$Hannaella), "low","high")
ggplot(df_for_CI_selected2, aes(x=temp_before_hd_sampling_10days, fill=A_group)) +
  geom_density(alpha=0.4)
ggplot(df_for_CI_selected2, aes(x=rhum_before_hd_sampling_10days, fill=A_group)) +
  geom_density(alpha=0.4)

# library("tableone")
# 
# CreateTableOne(vars=c("temp_before_hd_sampling_10days",
#                       "rhum_before_hd_sampling_10days"),
#                strata="A_group",
#                data=df_for_CI_selected2,
#                test=FALSE)

# # for loop -----------------------------------------------------------------------------------------------------
# library(brms)
# print("1st")
# # 분석할 미생물 변수들
# microbes <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia")
# 
# #######====================================================================================================================================================
# imsi_save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/causal_inference/ITE_ATE_CATE"
# ITE_save_version <- "v2.5"
# ITE_save_version2 <- "GZ_interaction_v2.5"
# high_weather_quantile <- 0.75
# low_weather_quantile <- 0.25
# microbes_q <- seq(0.1,0.9,0.1)
# for(ITE_save_version in c(ITE_save_version)){
#   print(ITE_save_version)
#   print(ITE_save_version2)
#   results <- list()
#   ATE_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
#   CATE_high_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
#   CATE_low_results <- data.frame(matrix(nrow = nrow(df_scaled), ncol = 0))
#   
#   for (mi in c(1:length(microbes))) {
#     # microbe <- microbes[1]
#     
#     microbe <- microbes[mi]
#     print(microbe)
#     # 모델 정의
#     # formula_str <- as.formula(
#     #   paste0("incidence ~ ",
#     #          
#     #          "s(",microbe,", k=3) +",
#     #          
#     #          "s(GZ_mean_max, k=3) +",
#     # 
#     #          "s(temp_flower_sampling_10days, k=3) +
#     #           
#     #           s(rhum_flower_sampling_10days, k=3)")
#     # )
#     
    # if(microbe == "Alternaria"){
    #   formula_str <- as.formula(
    #     paste0("incidence ~ ",
    #            "s(",microbe,", k=3) +",
    #            " +
    #           t2(Alternaria,Papiliotrema, k=4) +
    #           t2(Alternaria,Periconia, k=4) +
    #           s(GZ_mean_max, k=3) +
    #           s(temp_flower_sampling_10days, k=3) +
    # 
    #           s(rhum_flower_sampling_10days, k=3)")
    #   )
    # }else if(microbe == "Epicoccum"){
    #   formula_str <- as.formula(
    #     paste0("incidence ~ ",
    #            "s(",microbe,", k=3) +",
    #            " +
    #           t2(Epicoccum, Cladosporium, k=4) +
    #           s(GZ_mean_max, k=3) +
    #           s(temp_flower_sampling_10days, k=3) +
    # 
    #           s(rhum_flower_sampling_10days, k=3)")
    #   )
    # }else if(microbe == "Periconia"){
    #   formula_str <- as.formula(
    #     paste0("incidence ~ ",
    #            "s(",microbe,", k=3) +",
    #            " +
    #           t2(Periconia, Alternaria, k=4) +
    #           s(GZ_mean_max, k=3) +
    #           s(temp_flower_sampling_10days, k=3) +
    # 
    #           s(rhum_flower_sampling_10days, k=3)")
    #   )
    # }else{
    #   formula_str <- as.formula(
    #     paste0("incidence ~ ",
    # 
    #            "s(",microbe,", k=3) +",
    # 
    #            "s(GZ_mean_max, k=3) +",
    # 
    #            "s(temp_flower_sampling_10days, k=3) +
    # 
    #           s(rhum_flower_sampling_10days, k=3)")
    #   )
    # 
    # }
#     
#     
#     
#     
#     model <- brm(
#       formula = formula_str,
#       family = Beta(),
#       data = df_scaled,
#       prior = c(
#         prior(normal(0,1), class="b"),
#         prior(exponential(1), class="sds")
#       ),
#       chains = 4,
#       iter = 6000,
#       warmup = 3000,
#       cores = 4,
#       control = list(adapt_delta=0.995, max_treedepth=15),
#       silent = 2
#     )
#     
#     cat(paste0(microbe, " model finished"))
#     
#     ##### ATE #####
#     
#     new_high <- df_scaled
#     new_low  <- df_scaled
#     
#     new_high[[microbe]] <- quantile(df_scaled[[microbe]],0.9)
#     new_low[[microbe]]  <- quantile(df_scaled[[microbe]],0.1)
#     
#     pred_high <- posterior_epred(model,newdata=new_high)
#     pred_low  <- posterior_epred(model,newdata=new_low)
#     
#     ATE_draws <- pred_high - pred_low
#     ATE_draw_mean <- rowMeans(ATE_draws)
#     ITE_draw_mean <- colMeans(ATE_draws)
#     
#     # save ITE
#     pred_high_df <- as.data.frame(pred_high)
#     pred_low_df <- as.data.frame(pred_low)
#     pred_high_minus_low <- pred_high_df-pred_low_df
#     
#     names(pred_high_minus_low) <- names(pred_low_df) <- names(pred_high_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
#     
#     if(!dir.exists(paste0(imsi_save_dir, "/", ITE_save_version))){
#       dir.create(paste0(imsi_save_dir, "/", ITE_save_version))
#     }
#     
#     writexl::write_xlsx(pred_low_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_low_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(pred_high_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_high_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(pred_high_minus_low, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_results_pred_high_minus_low_",microbe, "_",ITE_save_version2,".xlsx"))
#     
#     print(paste0(microbe, " ITE finished"))
#     
#     # # ATE 계산
#     # ATE <- colMeans(pred_high - pred_low)
#     # 
#     # ATE_results <- cbind(ATE_results, ATE)
#     # names(ATE_results)[mi] <- paste0(microbe, "_ATE")
#     # 
#     # print(paste0(microbe, " ATE finished"))
#     
#     
#     ### CATE -----------------------------------------------------------------------
#     # 데이터 복제
#     new_cate_high_wth_T1 <- df_scaled
#     new_cate_high_wth_T0 <- df_scaled
#     
#     new_cate_low_wth_T1  <- df_scaled
#     new_cate_low_wth_T0  <- df_scaled
#     
#     # 해당 미생물의 90% quantile vs 10% quantile
#     new_cate_high_wth_T1[[microbe]] <- quantile(df_scaled[[microbe]], 0.9, na.rm = TRUE)
#     new_cate_high_wth_T0[[microbe]]  <- quantile(df_scaled[[microbe]], 0.1, na.rm = TRUE)
#     
#     new_cate_low_wth_T1[[microbe]] <- quantile(df_scaled[[microbe]], 0.9, na.rm = TRUE)
#     new_cate_low_wth_T0[[microbe]]  <- quantile(df_scaled[[microbe]], 0.1, na.rm = TRUE)
#     
#     new_cate_high_wth_T1$temp_flower_sampling_10days <-  new_cate_high_wth_T0$temp_flower_sampling_10days <-
#       quantile(df_scaled$temp_flower_sampling_10days, high_weather_quantile)
#     
#     new_cate_high_wth_T1$rhum_flower_sampling_10days <- new_cate_high_wth_T0$rhum_flower_sampling_10days <-
#       quantile(df_scaled$rhum_flower_sampling_10days, high_weather_quantile)
#     
#     new_cate_high_wth_T1$prcp_after_sampling <- new_cate_high_wth_T0$prcp_after_sampling <- 
#       quantile(df_scaled$prcp_after_sampling, high_weather_quantile)
#     
#     new_cate_low_wth_T1$temp_flower_sampling_10days <-  new_cate_low_wth_T0$temp_flower_sampling_10days <-
#       quantile(df_scaled$temp_flower_sampling_10days, low_weather_quantile)
#     
#     new_cate_low_wth_T1$rhum_flower_sampling_10days <- new_cate_low_wth_T0$rhum_flower_sampling_10days <-
#       quantile(df_scaled$rhum_flower_sampling_10days, low_weather_quantile)
#     
#     new_cate_low_wth_T1$prcp_after_sampling <- new_cate_low_wth_T0$prcp_after_sampling <- 
#       quantile(df_scaled$prcp_after_sampling, low_weather_quantile)
#     
#     ## 예측
#     # high pred
#     pred_high_wth_T1 <- posterior_epred(model, newdata = new_cate_high_wth_T1)
#     pred_high_wth_T0 <- posterior_epred(model, newdata = new_cate_high_wth_T0)
#     CATE_high_draws <- pred_high_wth_T1 - pred_high_wth_T0
#     
#     #low pred
#     pred_low_wth_T1 <- posterior_epred(model, newdata = new_cate_low_wth_T1)
#     pred_low_wth_T0 <- posterior_epred(model, newdata = new_cate_low_wth_T0)
#     CATE_low_draws <- pred_low_wth_T1 - pred_low_wth_T0
#     
#     ## save CITE
#     #CATE_high
#     pred_high_wth_T1_df <- as.data.frame(pred_high_wth_T1)
#     pred_high_wth_T0_df <- as.data.frame(pred_high_wth_T0)
#     pred_high_wth_T1_minus_T0 <- pred_high_wth_T1_df - pred_high_wth_T0_df
#     
#     names(pred_high_wth_T1_minus_T0) <- names(pred_high_wth_T1_df) <- names(pred_high_wth_T0_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
#     
#     writexl::write_xlsx(pred_high_wth_T0_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_high_wth_T0_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(pred_high_wth_T1_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_high_wth_T1_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(pred_high_wth_T1_minus_T0, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_high_wth_T1_minus_T0_",microbe, "_",ITE_save_version2,".xlsx"))
#     
#     #CATE_low
#     pred_low_wth_T1_df <- as.data.frame(pred_low_wth_T1)
#     pred_low_wth_T0_df <- as.data.frame(pred_low_wth_T0)
#     pred_low_wth_T1_minus_T0 <- pred_low_wth_T1_df - pred_low_wth_T0_df
#     
#     names(pred_low_wth_T1_minus_T0) <- names(pred_low_wth_T1_df) <- names(pred_low_wth_T0_df) <- paste0(df_for_more_analysis$wth_ID, "_", df_for_more_analysis$year)
#     
#     writexl::write_xlsx(pred_low_wth_T0_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_low_wth_T0_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(pred_low_wth_T1_df, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_low_wth_T1_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(pred_low_wth_T1_minus_T0, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CITE_results_low_wth_T1_minus_T0_",microbe, "_",ITE_save_version2,".xlsx"))
#     
#     
#     # CATE
#     CATE_high <- colMeans(pred_high_wth_T1_minus_T0)
#     CATE_high_results <- cbind(CATE_high_results, CATE_high)
#     names(CATE_high_results)[mi] <- paste0(microbe, "_CATE_high")
#     
#     CATE_low <- colMeans(pred_low_wth_T1_minus_T0)
#     CATE_low_results <- cbind(CATE_low_results, CATE_low)
#     names(CATE_low_results)[mi] <- paste0(microbe, "_CATE_low")
#     
#     ##### observation 평균 → effect per draw #####
#     
#     ATE_draw_mean <- rowMeans(ATE_draws)
#     CATE_high_draw_mean <- rowMeans(CATE_high_draws)
#     CATE_low_draw_mean <- rowMeans(CATE_low_draws)
#     
#     ITE_draw_mean <- colMeans(ATE_draws)
#     CITE_high_draw_mean <- colMeans(CATE_high_draws)
#     CITE_low_draw_mean <- colMeans(CATE_low_draws)
#     ##### Microbe × Weather surface #####
#     # temp
#     microbe_seq <- seq(min(df_scaled[[microbe]]),
#                        max(df_scaled[[microbe]]),
#                        length=5)
#     temp_seq <- seq(min(df_scaled$temp_flower_sampling_10days),
#                     max(df_scaled$temp_flower_sampling_10days),
#                     length=5)
#     
#     CATE_temp_draws_mean <- data.frame(matrix(nrow = 12000, ncol = 0))
#     for(w in c(1:length(temp_seq))){
#       for(m_l in microbe_seq){
#         for(m_h in microbe_seq){
#           #--------------------------
#           if(m_h > m_l){
#             new_imsi_cate_temp_T1 <- df_scaled
#             new_imsi_cate_temp_T0 <- df_scaled
#             
#             new_imsi_cate_temp_T1[[microbe]] <- m_h
#             new_imsi_cate_temp_T0[[microbe]]  <- m_l
#             
#             new_imsi_cate_temp_T1$temp_flower_sampling_10days <-  
#               new_imsi_cate_temp_T0$temp_flower_sampling_10days <- temp_seq[w]
#             
#             pred_imsi_temp_wth_T1 <- posterior_epred(model, newdata = new_imsi_cate_temp_T1)
#             pred_imsi_temp_wth_T0 <- posterior_epred(model, newdata = new_imsi_cate_temp_T0)
#             CATE_imsi_temp_draws <- pred_imsi_temp_wth_T1 - pred_imsi_temp_wth_T0
#             
#             CATE_imsi_temp_draws_mean <- rowMeans(CATE_imsi_temp_draws)  
#             
#             CATE_temp_draws_mean <- cbind(CATE_temp_draws_mean, CATE_imsi_temp_draws_mean)
#             names(CATE_temp_draws_mean)[ncol(CATE_temp_draws_mean)] <- paste0("mh_", which(microbe_seq == m_h), "ml_",which(microbe_seq == m_l))
#             cat("m_l = ", m_l, "m_h = ", m_h)
#           }else{
#             cat("m_l = ", m_l, "m_h = ", m_h)
#           }
#           #--------------------------
#         }
#       }
#     }
#     
#     CATE_rhum_draws_mean<- data.frame(matrix(nrow = 12000, ncol = 0))
#     for(w in c(1:length(rhum_seq))){
#       for(m_l in microbe_seq){
#         for(m_h in microbe_seq){
#           #--------------------------
#           if(m_h > m_l){
#             new_imsi_cate_rhum_T1 <- df_scaled
#             new_imsi_cate_rhum_T0 <- df_scaled
#             
#             new_imsi_cate_rhum_T1[[microbe]] <- m_h
#             new_imsi_cate_rhum_T0[[microbe]]  <- m_l
#             
#             new_imsi_cate_rhum_T1$rhum_flower_sampling_10days <-  
#               new_imsi_cate_rhum_T0$rhum_flower_sampling_10days <- rhum_seq[w]
#             
#             pred_imsi_rhum_wth_T1 <- posterior_epred(model, newdata = new_imsi_cate_rhum_T1)
#             pred_imsi_rhum_wth_T0 <- posterior_epred(model, newdata = new_imsi_cate_rhum_T0)
#             CATE_imsi_rhum_draws <- pred_imsi_rhum_wth_T1 - pred_imsi_rhum_wth_T0
#             
#             CATE_imsi_rhum_draws_mean <- rowMeans(CATE_imsi_rhum_draws)  
#             
#             CATE_rhum_draws_mean <- cbind(CATE_rhum_draws_mean, CATE_imsi_rhum_draws_mean)
#             names(CATE_rhum_draws_mean)[ncol(CATE_rhum_draws_mean)] <- paste0("mh_", which(microbe_seq == m_h), "ml_",which(microbe_seq == m_l))
#             cat("m_l = ", m_l, "m_h = ", m_h)
#           }else{
#             cat("m_l = ", m_l, "m_h = ", m_h)
#           }
#           #--------------------------
#         }
#       }
#     }
#     
#     CATE_prcp_draws_mean <- data.frame(matrix(nrow = 12000, ncol = 0))
#     for(w in c(1:length(prcp_seq))){
#       for(m_l in microbe_seq){
#         for(m_h in microbe_seq){
#           #--------------------------
#           if(m_h > m_l){
#             new_imsi_cate_prcp_T1 <- df_scaled
#             new_imsi_cate_prcp_T0 <- df_scaled
#             
#             new_imsi_cate_prcp_T1[[microbe]] <- m_h
#             new_imsi_cate_prcp_T0[[microbe]]  <- m_l
#             
#             new_imsi_cate_prcp_T1$prcp_after_sampling <-  
#               new_imsi_cate_prcp_T0$prcp_after_sampling <- prcp_seq[w]
#             
#             pred_imsi_prcp_wth_T1 <- posterior_epred(model, newdata = new_imsi_cate_prcp_T1)
#             pred_imsi_prcp_wth_T0 <- posterior_epred(model, newdata = new_imsi_cate_prcp_T0)
#             CATE_imsi_prcp_draws <- pred_imsi_prcp_wth_T1 - pred_imsi_prcp_wth_T0
#             
#             CATE_imsi_prcp_draws_mean <- rowMeans(CATE_imsi_prcp_draws)  
#             
#             CATE_prcp_draws_mean <- cbind(CATE_prcp_draws_mean, CATE_imsi_prcp_draws_mean)
#             names(CATE_prcp_draws_mean)[ncol(CATE_prcp_draws_mean)] <- paste0("mh_", which(microbe_seq == m_h), "ml_",which(microbe_seq == m_l))
#             cat("m_l = ", m_l, "m_h = ", m_h)
#           }else{
#             cat("m_l = ", m_l, "m_h = ", m_h)
#           }
#           #--------------------------
#         }
#       }
#     }
#     
#     
#     
#     # 결과 저장
#     results[[microbe]] <- list(
#       
#       model=model,
#       
#       ATE_draw=ATE_draw_mean,
#       CATE_high_draw=CATE_high_draw_mean,
#       CATE_low_draw=CATE_low_draw_mean, 
#       
#       ITE_draw_mean = ITE_draw_mean, 
#       CITE_high_draw_mean = CITE_high_draw_mean,
#       CITE_low_draw_mean = CITE_low_draw_mean,
#       
#       CATE_temp_draws_mean= CATE_temp_draws_mean,
#       CATE_rhum_draws_mean= CATE_rhum_draws_mean,
#       CATE_prcp_draws_mean= CATE_prcp_draws_mean
#     )
#     
#     cat("// Done:", microbe, "\n")
#     
#     ATE_CATE_draw_mean <- data.frame(ATE_draw_mean = ATE_draw_mean, 
#                                      CATE_high_draw_mean = CATE_high_draw_mean, 
#                                      CATE_low_draw_mean = CATE_low_draw_mean)
#     ITE_CITE_draw_mean <- data.frame(ITE_draw_mean = ITE_draw_mean, 
#                                      CITE_high_draw_mean = CITE_high_draw_mean, 
#                                      CITE_low_draw_mean = CITE_low_draw_mean)
#     
#     writexl::write_xlsx(ATE_CATE_draw_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ATE_CATE_draw_mean_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(ITE_CITE_draw_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ITE_CITE_draw_mean_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(CATE_temp_draws_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_temp_draws_mean_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(CATE_rhum_draws_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_rhum_draws_mean_",microbe, "_",ITE_save_version2,".xlsx"))
#     writexl::write_xlsx(CATE_prcp_draws_mean, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_prcp_draws_mean_",microbe, "_",ITE_save_version2,".xlsx"))
#   }
#   writexl::write_xlsx(ATE_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/ATE_results_",ITE_save_version2,".xlsx"))
#   writexl::write_xlsx(CATE_high_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_high_results_",ITE_save_version2,".xlsx"))
#   writexl::write_xlsx(CATE_low_results, path = paste0(imsi_save_dir, "/", ITE_save_version, "/CATE_low_results_",ITE_save_version2,".xlsx"))
# }
# 
# ATE_results
# CATE_high_results
# CATE_low_results
# 
# colMeans(ATE_results)
# colMeans(CATE_high_results)
# colMeans(CATE_low_results)
# 
# # graph
# ITE_plot_df <- data.frame()
# for(microbe in microbes){
#   ITE <- results[[microbe]]$ITE_draw_mean
#   CITE_high <- results[[microbe]]$CITE_high_draw_mean
#   CITE_low  <- results[[microbe]]$CITE_low_draw_mean
#   
#   ite_tmp <- data.frame(
#     microbe = microbe,
#     effect_type = c("ITE","CITE_high","CITE_low"),
#     mean = c(mean(ITE),
#              mean(CITE_high),
#              mean(CITE_low)),
#     lower = c(quantile(ITE,0.025),
#               quantile(CITE_high,0.025),
#               quantile(CITE_low,0.025)),
#     upper = c(quantile(ITE,0.975),
#               quantile(CITE_high,0.975),
#               quantile(CITE_low,0.975))
#   )
#   
#   ITE_plot_df <- rbind(ITE_plot_df,ite_tmp)
# }
# 
# ITE_g <- ggplot(ITE_plot_df,
#                 aes(x = mean,
#                     y = microbe,
#                     color = effect_type)) +
#   
#   geom_vline(xintercept = 0,
#              linetype = "dashed",
#              color = "grey40") +
#   
#   geom_errorbarh(aes(xmin = lower,
#                      xmax = upper),
#                  height = 0.2,
#                  position = position_dodge(width = 0.5),
#                  size = 1.2) +
#   
#   geom_point(size = 4,
#              shape = "x",
#              position = position_dodge(width = 0.5)) +
#   
#   scale_color_manual(values = c(
#     ITE = "black",
#     CITE_high = "#d73027",
#     CITE_low = "#4575b4"
#   )) +
#   
#   labs(
#     x = "ITE",
#     y = "",
#     color = ""
#   ) +
#   theme_bw(base_size = 15) + 
#   theme(
#     plot.title = element_text(size = 10, face = 'bold'),
#     axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
#     axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
#     axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
#                                size = 15, face = 'bold', color = 'black'),
#     axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
#   )
# ggsave(
#   plot = ITE_g,
#   file = paste0(
#     "./Output/3. obs_pred_graph/2025/",
#     version_pred_obs,
#     "/causal_inference/ITE_ATE_CATE/",
#     ITE_save_version, "/graph/",
#     "ITE_plot_",
#     ITE_save_version2,
#     ".png"
#   ),
#   width = 15,
#   height = 10,
#   units = c("cm")
# )
# 
# ##
# # ATE_plot_df <- data.frame()
# # for(microbe in microbes){
# #   ATE <- results[[microbe]]$ITE_draw_mean
# #   CATE_high <- results[[microbe]]$CATE_high_draw
# #   CATE_low  <- results[[microbe]]$CATE_low_draw
# #   
# #   tmp <- data.frame(
# #     microbe = microbe,
# #     effect_type = c("ATE","CATE_high","CATE_low"),
# #     mean = c(mean(ATE),
# #              mean(CATE_high),
# #              mean(CATE_low)),
# #     lower = c(quantile(ATE,0.025),
# #               quantile(CATE_high,0.025),
# #               quantile(CATE_low,0.025)),
# #     upper = c(quantile(ATE,0.975),
# #               quantile(CATE_high,0.975),
# #               quantile(CATE_low,0.975))
# #   )
# #   
# #   ATE_plot_df <- rbind(ATE_plot_df,tmp)
# # }
# ATE_plot_df <- data.frame()
# for(microbe in microbes){
#   ATE <- results[[microbe]]$ITE_draw_mean
#   # CATE_high <- results[[microbe]]$CATE_high_draw
#   # CATE_low  <- results[[microbe]]$CATE_low_draw
#   
#   tmp <- data.frame(
#     microbe = microbe,
#     effect_type = c("ATE"), #,"CATE_high","CATE_low"
#     mean = c(mean(ATE)#,
#              # mean(CATE_high),
#              # mean(CATE_low)
#     ),
#     lower = c(quantile(ATE,0.025)#,
#               # quantile(CATE_high,0.025),
#               # quantile(CATE_low,0.025)
#     ),
#     upper = c(quantile(ATE,0.975)#,
#               # quantile(CATE_high,0.975),
#               # quantile(CATE_low,0.975)
#     )
#   )
#   
#   ATE_plot_df <- rbind(ATE_plot_df,tmp)
# }
# ATE_plot_df$mean <- ATE_plot_df$mean*100
# ATE_plot_df$lower <- ATE_plot_df$lower*100
# ATE_plot_df$upper <- ATE_plot_df$upper*100
# ATE_plot_df$color_flag <- ifelse(ATE_plot_df$upper > 0, "pos", "default")
# ATE_posterior_g <- ggplot(ATE_plot_df,
#                           aes(x = mean,
#                               y = microbe,
#                               color = color_flag)) +
#   
#   geom_vline(xintercept = 0,
#              linetype = "dashed",
#              color = "grey40") +
#   
#   geom_errorbarh(aes(xmin = lower,
#                      xmax = upper),
#                  height = 0.2,
#                  size = 1.2) +
#   
#   geom_point(size = 4, shape = "x") +
#   
#   scale_color_manual(values = c(
#     default = "#004C99",
#     pos = "grey60"
#   )) +
#   
#   labs(
#     x = "\nEffect on FHB incidence (%)",
#     y = "",
#     color = ""
#   ) +
#   
#   theme_bw(base_size = 15) + 
#   theme(
#     axis.title.x = element_text(size = 15, face = 'bold'),
#     axis.title.y = element_text(size = 15, face = 'bold'),
#     axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
#     axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
#   )
# ggsave(
#   plot = ATE_posterior_g,
#   file = paste0(
#     "./Output/3. obs_pred_graph/2025/",
#     version_pred_obs,
#     "/causal_inference/ITE_ATE_CATE/",
#     ITE_save_version, "/graph/",
#     "ATE_posterior_",
#     ITE_save_version2,
#     ".png"
#   ),
#   width = 15,
#   height = 10,
#   units = c("cm")
# )
# 필요한 패키지 로드
library(tidyverse)
library(WeightIt)
library(cobalt)
library(mgcv)
# ---------------------------------------------------------
# [단계 0] 실제 데이터와 유사한 가상 데이터 생성 (N=98)
# ---------------------------------------------------------
target_m_vec <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia")
causal_infer_save_version <- "2025_v13"
overlap_res_table <- as.data.frame(matrix(nrow = 0, ncol = 5))
CI_res_table <- as.data.frame(matrix(nrow = 0, ncol = 3))
names(CI_res_table) <- c("point_estimate", "CI_lower", "CI_upper")
names(overlap_res_table) <- c("range_min", "range_max", "coef_of_var", "unweighted_ESS", "weighted_ESS")
run_causal_reg <- function(sim_data, indices,target_m) {# 부트스트랩 1회 진행 시 수행할 함수 정의
  # sim_data = sim_data; indices = boot_indices; target_m = target_m
  # 복원 추출로 샘플링된 데이터셋 생성
  boot_data <- sim_data[indices, ]
  
  imsi_fmu <- as.formula(paste(target_m, "~ temp_flower_sampling_10days + rhum_flower_sampling_10days + GZ_mean_max"))
  # 1. 샘플링된 데이터로 가중치 계산
  w_boot <- weightit(
    formula = imsi_fmu,
    data = boot_data,
    method = "ps",
    stabilize = F
  )
  # ---------------------------------------------------------
  # [단계 2] 부트스트랩을 통한 인과적 효과 및 95% 신뢰구간 추정
  # ---------------------------------------------------------
  # 가중치를 데이터셋에 부착
  boot_data$sw <- w_boot$weights
  
  # # 2. 가중 선형 회귀(Weighted OLS) 모델 적합
  # fit <- lm(incidence ~ get(target_m), data = boot_data, weights = sw)
  # 
  # 가중치를 적용한 GAM 모델 적합 (미생물 효과에 스플라인 s() 적용)
  
  imsi_gam_formula <- as.formula(make_formula_inc(microbe = target_m)) #as.formula(paste("incidence ~ s(", target_m, ", k = 3)")) # 데이터가 적으므로 k는 작게 설정
  fit <- gam(formula = imsi_gam_formula, data = boot_data, weights = sw)
  
  # 3. 미생물(Epicoccum)의 계수(인과적 영향력) 반환
  return(coef(fit)[2])
}

make_formula_inc <- function(microbe){
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
  
  return(formula_str)
}

make_formula_m <- function(microbe){
  if(microbe == "Alternaria"){
    formula_str <- as.formula(
      paste0(microbe, " ~ ",
             "t2(Alternaria,Papiliotrema, k=4) +
              t2(Alternaria,Periconia, k=4) +
              s(GZ_mean_max, k=3) +
              s(temp_flower_sampling_10days, k=3) +
              s(rhum_flower_sampling_10days, k=3)")
    )
  }else if(microbe == "Epicoccum"){
    formula_str <- as.formula(
      paste0(microbe, " ~ ",
             "t2(Epicoccum, Cladosporium, k=4) +
              s(GZ_mean_max, k=3) +
              s(temp_flower_sampling_10days, k=3) +

              s(rhum_flower_sampling_10days, k=3)")
    )
  }else if(microbe == "Periconia"){
    formula_str <- as.formula(
      paste0(microbe, " ~ ",
             "t2(Periconia, Alternaria, k=4) +
              s(GZ_mean_max, k=3) +
              s(temp_flower_sampling_10days, k=3) +

              s(rhum_flower_sampling_10days, k=3)")
    )
  }else{
    formula_str <- as.formula(
      paste0(microbe, " ~ ",
             "s(GZ_mean_max, k=3) +",
             "s(temp_flower_sampling_10days, k=3) +
              s(rhum_flower_sampling_10days, k=3)")
    )
    
  }
  
  return(formula_str)
}


for(m in 1:length(target_m_vec)){
  #m=1
  target_m <- target_m_vec[m]#"Epicoccum"
  set.seed(42) # 결과 재현을 위한 시드 설정
  
  sim_data <- df_for_CI_selected #[,c("incidence", target_m, "temp_flower_sampling_10days", "rhum_flower_sampling_10days", "GZ_mean_max")]
  
  # ---------------------------------------------------------
  # [단계 1] 연속형 처치(T)에 대한 Overlap 가정이 충족되는지 확인
  # ---------------------------------------------------------
  # WeightIt 패키지를 사용하여 안정화된 역확률 가중치(Stabilized IPW) 계산
  fmu <- as.formula(paste(target_m, "~ temp_flower_sampling_10days + rhum_flower_sampling_10days + GZ_mean_max"))
  weight_obj <- weightit(
    formula = fmu,
    data = sim_data,
    method = "ps",      # 연속형 변수의 경우 밀도 기반 추정 방식 사용
    stabilize = F    # 가중치 안정화 옵션 (분산 감소 효과)
  )
  weight_obj <- trim(weight_obj, at = .95, lower = TRUE)
  
  # Overlap 및 교란변수 균형도(Balance) 확인
  print("--- 가중치 적용 전/후 교란변수 균형도 확인 ---")
  # print(summary(weight_obj))
  summary_w <- summary(weight_obj)
  overlap_res_table[m,1] <- imsi_range_min <- summary_w$weight.range$all[1]
  overlap_res_table[m,2] <- imsi_range_max <- summary_w$weight.range$all[2]
  overlap_res_table[m,3] <- imsi_coef_of_var <- as.vector(summary_w$coef.of.var)[1]
  overlap_res_table[m,4] <- imsi_unweighted_ESS <- summary_w$effective.sample.size$Total[1]#ESS : Effective sample size
  overlap_res_table[m,5] <- imsi_weighted_ESS <- summary_w$effective.sample.size$Total[2]#ESS : Effective sample size
  
  # 시각적으로 확인하고 싶다면 아래 코드로 가중치 적용 전후의 분포를 볼 수 있습니다.
  imsi_dist_of_weight_graph <- plot(summary(weight_obj))
  imsi_love_plot <- love.plot(weight_obj, stars = "raw", threshold = 0.1)
  ggsave(filename = paste0("./Output/3. obs_pred_graph/2025/",causal_infer_save_version,"/causal_inference/weight_distridution_graph_", target_m, ".png"),
         plot = imsi_dist_of_weight_graph, 
         width = 4, 
         height = 4)
  ggsave(filename = paste0("./Output/3. obs_pred_graph/2025/",causal_infer_save_version,"/causal_inference/love_plot_", target_m, ".png"),
         plot = imsi_love_plot, 
         width = 5, 
         height = 2)
  
  # 부트스트랩 실행 (시간 관계상 1000회 추천, 여기서는 1000회 진행)
  B <- 1000
  boot_results <- numeric(B)
  
  for(i in 1:B) {
    # 98개 중에서 복원 추출할 인덱스 생성
    boot_indices <- sample(1:nrow(sim_data), replace = TRUE)
    # 함수 실행 후 인과 계수 저장
    boot_results[i] <- run_causal_reg(sim_data, boot_indices, target_m = target_m)
  }
  
  # ---------------------------------------------------------
  # [단계 3] 최종 결과 해석 및 신뢰구간 도출
  # ---------------------------------------------------------
  # 점 추정치 (전체 데이터로 구한 값)
  sim_data$sw_final <- weight_obj$weights
  # lm_formula <- as.formula(paste("incidence ~ ", target_m))
  # final_model <- lm(formula = lm_formula, data = sim_data, weights = sw_final)
  gam_formula <- as.formula(make_formula_inc(microbe = target_m)) #as.formula(paste("incidence ~ s(", target_m, ", k = 3)")) # 데이터가 적으므로 k는 작게 설정
  final_model <- gam(formula = gam_formula, data = sim_data, weights = sw_final)
  point_estimate <- as.vector(coef(final_model)[2])
  
  # 95% 신뢰구간 계산 (2.5% 백분위수 ~ 97.5% 백분위수)
  ci_lower <- quantile(boot_results, 0.025)
  ci_upper <- quantile(boot_results, 0.975)
  
  CI_res_table[m,1] <- point_estimate
  CI_res_table[m,2] <- ci_lower
  CI_res_table[m,3] <- ci_upper
  
  cat("\n==================================================\n")
  cat("--- 최종 인과추론 분석 결과 ---\n")
  cat(paste(target_m, "의 점 추정치 (인과 효과):", round(point_estimate, 4), "\n"))
  cat(paste("부트스트랩 기반 95% 신뢰구간: [", round(ci_lower, 4), ",", round(ci_upper, 4), "]\n"))
  
  # 인과적 유의성 판단
  if(ci_lower <= 0 && ci_upper >= 0) {
    cat("결론: 신뢰구간이 0을 포함하므로, 발병률에 유의미한 인과적 효과가 없습니다.\n")
  } else {
    cat("결론: 신뢰구간이 0을 포함하지 않으므로, 발병률에 통계적으로 유의미한 인과적 효과가 있습니다.\n")
  }
  cat("==================================================\n")
  
}
overlap_res_table
write_xlsx(overlap_res_table, path = paste0("./Output/3. obs_pred_graph/2025/",causal_infer_save_version,"/causal_inference/overlap_res_table_", causal_infer_save_version, ".xlsx"))

plot(summary(weight_obj))

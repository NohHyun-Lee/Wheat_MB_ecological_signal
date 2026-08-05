### version 
# v5 : no CLR, Yes interaction, 
# v7 : CLR, Yes interaction, just linear regression
# v8 : CLR, Yes interaction, interaction between MB = with filtered OTUs data, 
# v9 : CLR, Yes interaction, interaction between MB = with filtered OTUs data, genus level
# v10.1 : CLR, Yes interaction, interaction between MB = with filtered OTUs data, genus level (key M = 17 Genus)
# v10.2 : CLR, Yes interaction, interaction between MB = with filtered OTUs data, genus level (key M should be |rho| >= 0.4, total 4 Genus)
# v10.3 : CLR, Yes interaction, interaction between MB = with filtered OTUs data, genus level (key M should be |rho| >= 0.3, total 6 Genus)
# v10.4 : CLR, Yes interaction, interaction between MB = with filtered OTUs data, genus level (key M should be |rho| >= 0.281, total 7 Genus)
# v11 : same as v10.3, but include more analysis (such as compared to the GIBSIM results vs. hybrid model results)
# v12 : same as v11, but using different dataset 
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
library(mgcv)

# conflict_prefer("as.matrix", "base")
### function -------------------------------------------------------------------
#https://yogyui.tistory.com/entry/R-%ED%9A%8C%EA%B7%80%EB%B6%84%EC%84%9D-%EB%AA%A8%EB%8D%B8-%EC%84%B1%EB%8A%A5%ED%8C%90%EB%8B%A8-RMSE-MAE-R-square
rsquare <- function(y1,y2){
  sum((y1- mean(y1))*(y2-mean(y2)))^2 / (sum((y1-mean(y1))^2) * sum((y2-mean(y2))^2))
}


min_max_scaling_for_df <- function(data,Y_value){
  #Y_value <- "FHB_incidence"
  data <- data
  X_data <- data[,which(names(data) != Y_value)]
  max(X_data[,1])
  for(c in c(1:ncol(X_data))){
    #c=1
    X_data[,c] <- min_max_scaling(X_data[,c])
  }
  XY_data <- cbind(X_data, data[,which(names(data) == Y_value)])
  return(XY_data)
}


obs_pred_plot <- function(model, data_contain_X_Y, Y_colname,Title){
  # model = GIBSIM_obs_lm
  # data_contain_X_Y = selected_data_for_lm
  # Y_colname = "FHB_incidence"
  # Title = "FHB_incidence"
  
  predict_tr <- as.data.frame(predict(model, newdata = data_contain_X_Y))
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  predict_tr$pred <- ifelse(predict_tr$pred < 0, 0, predict_tr$pred)
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict(fit2))
  
  # fit <- lm(obs ~ pred, data = predict_tr)
  # r2 <- summary(fit)$r.squared
  # rmse <- Metrics::rmse(predict_tr$obs, predict(fit))
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    # geom_line(aes(y = predict(fit)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(predict_tr$obs), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}

obs_pred_plot_for_gamma3 <- function(model, data_contain_X_Y, Y_colname,Title){
  # model <- glm_gamma_step
  # data_contain_X_Y <- data_gibsim_MB_for_gamma
  # Y_colname = "FHB_incidence"
  # Title <- "t"
  
  predict_tr <- as.data.frame(predict(model, newdata = data_contain_X_Y, type = "response"))
  # predict_tr[,1] <- exp(predict_tr[,1])
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict(fit2))
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    geom_line(aes(y = predict(fit)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(c(predict_tr$obs, predict_tr$pred)), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}

obs_pred_plot_for_GAM <- function(model, data_contain_X_Y, Y_colname,Title){
  # model = gam_model_bf_wth 
  # data_contain_X_Y = data_gibsim_MB_for_beta 
  # Y_colname = "FHB_incidence"
  # Title = "GAM_model_bf_wth"
  
  predict_tr <- as.data.frame(predict(model, type = "response")) #as.data.frame(predict(model, newdata = data_contain_X_Y))
  # predict_tr[,1] <- exp(predict_tr[,1])
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  
  predict_tr$pred <- predict_tr$pred*100
  predict_tr$obs <- predict_tr$obs*100
  
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict_tr$pred)
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    # geom_line(aes(y = predict(fit2)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(predict_tr$obs), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}


obs_pred_plot_for_gamma <- function(model, data_contain_X_Y, Y_colname,Title){
  # model <- glm_stepwise_inc
  # data_contain_X_Y <- wth_intensity_inc_for_lm
  # Y_colname = "FHB_incidence"
  # Title <- "t"
  
  predict_tr <- as.data.frame(predict(model, newdata = data_contain_X_Y))
  predict_tr[,1] <- (1/(100*exp(predict_tr[,1])))
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict(fit2))
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    geom_line(aes(y = predict(fit)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(predict_tr$obs), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}

obs_pred_plot_for_gamma2 <- function(model, data_contain_X_Y, Y_colname,Title){
  # model <- glm_stepwise
  # data_contain_X_Y <- wth_intensity_inc_for_lm
  # Y_colname = "FHB_incidence"
  # Title <- "t"
  
  predict_tr <- as.data.frame(predict(model, newdata = data_contain_X_Y))
  predict_tr[,1] <- exp(predict_tr[,1])
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict(fit2))
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    geom_line(aes(y = predict(fit)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(predict_tr$obs), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}



obs_pred_plot_for_minus_value <- function(model, data_contain_X_Y, Y_colname,Title){
  # model = restricted_model
  # data_contain_X_Y = selected_data_only_MB
  # Y_colname = "FHB_incidence"
  # Title = "FHB_incidence"
  
  predict_tr <- as.data.frame(predict(model, newdata = data_contain_X_Y))
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  
  predict_tr$pred[predict_tr$pred < 0] <- 0
  predict_tr$pred[predict_tr$pred > 100] <- 100
  
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict(fit2))
  
  # fit <- lm(obs ~ pred, data = predict_tr)
  # r2 <- summary(fit)$r.squared
  # rmse <- Metrics::rmse(predict_tr$obs, predict(fit))
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    # geom_line(aes(y = predict(fit)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(predict_tr$obs), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}

obs_pred_plot_for_minus_value_for_restricted_model <- function(model, data_contain_X_Y, Y_colname,Title){
  # model = restricted_model 
  # data_contain_X_Y = data_gibsim_MB_for_model
  # Y_colname = "FHB_incidence"
  # Title = "FHB_incidence"
  
  # restricted_model??? 계수 ?????? 추출
  vars_keep <- names(coef(restricted_model))
  
  # model.matrix ?????? (???체에??? restricted model 변?????? 추출)
  X_restricted <- model.matrix(~ ., data = data_contain_X_Y)[ , vars_keep, drop = FALSE]
  
  # coef(restricted_model)??? named vector??? ????????? 맞춰줘야 ???
  beta <- coef(restricted_model)[vars_keep]
  
  # ???측값 계산
  y_pred <- X_restricted %*% beta
  
  predict_tr <- as.data.frame(y_pred)
  names(predict_tr)[1] <- "pred"
  
  predict_tr$pred[predict_tr$pred < 0] <- 0
  predict_tr$pred[predict_tr$pred > 100] <- 100
  
  #obs
  predict_tr$obs <- as.data.frame(data_contain_X_Y)[,which(names(data_contain_X_Y)==as.character(Y_colname))]
  colnames(predict_tr) <- c("pred", "obs")
  # plot(predict_tr$pred, predict_tr$obs)
  
  fit <- lm(obs ~ pred, data = predict_tr)
  predict_tr$pred <- (predict_tr$pred*summary(fit)[[4]][2,1])+summary(fit)[[4]][1,1] 
  
  predict_tr$pred[predict_tr$pred < 0] <- 0
  predict_tr$pred[predict_tr$pred > 100] <- 100
  
  fit2 <- lm(obs ~ pred, data = predict_tr)
  r2 <- summary(fit2)$r.squared
  rmse <- Metrics::rmse(predict_tr$obs, predict(fit2))
  
  # fit <- lm(obs ~ pred, data = predict_tr)
  # r2 <- summary(fit)$r.squared
  # rmse <- Metrics::rmse(predict_tr$obs, predict(fit))
  # ccc <- CCC(predict_tr$obs, predict(fit))$rho.c[1]
  eq_label <- paste0("R2 = ", round(r2, 2), "\nRMSE = ", round(rmse, 2))
  
  graph  <- ggplot2::ggplot(predict_tr, aes(x=pred, y= obs)) + 
    geom_point(color = "black", alpha = 0.5, shape = 16, size = 4) + 
    theme_classic() + 
    theme(axis.text = element_text(size=27),
          axis.title = element_text(size=27)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5)) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5)) +
    theme(plot.title = element_text(size = 30)) + 
    # scale_y_continuous(limits = c(0, 36)) +
    ggtitle(Title)+
    ylab("Observed\n")+
    xlab("\nPredicted") + #Features\n
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # y=x ??
    # geom_smooth(method = "lm", color = "black", se = F, size = 1, alpha = 0.8) +  # ??????
    guides(fill=FALSE) +
    # geom_line(aes(y = predict(fit)), color = "black", linetype = "solid", linewidth = 1.3, alpha = 0.8) +  # fit ??? 추???
    annotate("text", x = min(predict_tr$pred), y = max(predict_tr$obs), label = eq_label, hjust = 0, vjust = 1, size = 10, color = "darkgray")
  
  return(graph)
}

### data -----------------------------------------------------------------------
# MB_intensity_wth_cor_df <- readxl::read_xlsx("C:/Users/B550/Desktop/work/9. ????????? MB/data/MB_intensity_wth_cor_df.xlsx")
# df_rslt_stat <- read_xlsx("D:/?????????/microbiome_SSU/2024_wheat_data/disease_data/MB_wheat_cor_wth_data_v1.3.xlsx")
df_rslt_stat <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v4.xlsx") 
# df_rslt_stat <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v2.xlsx") #same file = D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/2. merge_meta/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v2.xlsx
# df_rslt_stat <- read_xlsx("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2-1. Run_GIBSIM/Run_GIBSIM/output/for_MB_analysis_v1.3/MB_cor_wth_data/MB_wheat_cor_wth_data_v1.3.xlsx")
# df_rslt_stat <- read_xlsx("F:/?????????/microbiome_SSU/Rcode/FHB_analysis/2-1. Run_GIBSIM/Run_GIBSIM/output/for_MB_analysis_v1.3.3/MB_cor_wth_data/MB_wheat_cor_wth_data_v1.3.3.xlsx")
names(df_rslt_stat)
head(df_rslt_stat)
nrow(df_rslt_stat)


### group by loc -----------------------------------------------------------------------------------------------------
# df_rslt_stat$sample_loc <- NULL
# df_rslt_stat$crop <- NULL
# # df_rslt_stat$loc <- NULL
# df_rslt_stat$Lat <- NULL
# df_rslt_stat$Lon <- NULL
# df_rslt_stat$Alt <- NULL
# 
# table(df_rslt_stat$wth_ID)
# 
# df_rslt_stat <- df_rslt_stat %>% group_by(source, wth_ID, loc, year, heading_date, heading_date_doy) %>% summarise(across(everything(), mean))
# sum(is.na(df_rslt_stat)==T)
# nrow(df_rslt_stat)
# table(df_rslt_stat$loc)
# df_rslt_stat2 <- df_rslt_stat[,c("wth_ID","fhb_incidence_predicted", "FHB_incidence")]
# df_rslt_stat2$minus <- abs(df_rslt_stat2$FHB_incidence*100 - df_rslt_stat2$fhb_incidence_predicted)

### group by loc -----------------------------------------------------------------------------------------------------
nrow(df_rslt_stat)
df_rslt_stat <- df_rslt_stat[complete.cases(df_rslt_stat),]
table(df_rslt_stat$loc)
nrow(df_rslt_stat)

# df_rslt_stat <- df_rslt_stat %>% filter(wth_ID != "ID_043") %>% filter(wth_ID != "ID_124") %>%
#   filter(wth_ID != "ID_022") %>% filter(wth_ID != "ID_021") %>% filter(wth_ID != "ID_027")# %>% filter(wth_ID != "ID_027") %>%
nrow(df_rslt_stat)
# df_rslt_stat$new_wth_ID


# df_rslt_stat <- df_rslt_stat[c(which(df_rslt_stat$wth_ID == "SSU_008"), 
#              which(df_rslt_stat$wth_ID == "ID_043"), 
#              which(df_rslt_stat$wth_ID == "ID_027"), 
#              which(df_rslt_stat$wth_ID == "ID_022"), 
#              which(df_rslt_stat$wth_ID == "ID_021")), ]
# df_rslt_stat$wth_ID[which(df_rslt_stat$FHB_incidence == max(df_rslt_stat$FHB_incidence))]
# df_rslt_stat$FHB_incidence[which(df_rslt_stat$FHB_incidence == max(df_rslt_stat$FHB_incidence))]
# df_rslt_stat$fhb_incidence_predicted[which(df_rslt_stat$FHB_incidence == max(df_rslt_stat$FHB_incidence))]

MB_intensity_wth_cor_df <- df_rslt_stat
colnames(MB_intensity_wth_cor_df)
MB_intensity_wth_cor_df <- MB_intensity_wth_cor_df %>% dplyr::select(-contains("wsd")) # select(!contains("wsd"))
version_pred_obs <- "2025_v12"
# write_xlsx(MB_intensity_wth_cor_df, 
#            path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/MB_intensity_wth_cor_df.xlsx")
sum(is.na(MB_intensity_wth_cor_df == T))
MB_intensity_wth_cor_df <- MB_intensity_wth_cor_df[complete.cases(MB_intensity_wth_cor_df),]
nrow(MB_intensity_wth_cor_df)
# MB_intensity_wth_cor_df$new_wth_ID

# ### only relative abundance >= 0.1 OTUs 
# MB_colnum <- grep("OTU", names(selected_data_for_lm))
# MB_intensity_wth_cor_df_MB <- MB_intensity_wth_cor_df[,MB_colnum]
# MB_intensity_wth_cor_df_not_MB <- MB_intensity_wth_cor_df[,-MB_colnum]
# MB_intensity_wth_cor_df_MB_more_0.1 <- MB_intensity_wth_cor_df_MB[, apply(selected_data_for_lm_MB, 2, function(col) any(col >= 0.1))]
# 
# MB_intensity_wth_cor_df_MB_more_0.1_final <- cbind(MB_intensity_wth_cor_df_not_MB, MB_intensity_wth_cor_df_MB_more_0.1)

###############################################################################################################################################
selected_data_for_lm <- MB_intensity_wth_cor_df
nrow(selected_data_for_lm)
names(selected_data_for_lm)
names(selected_data_for_lm)[which(names(selected_data_for_lm) == "incidence")] <- "FHB_incidence"
selected_data_for_lm$FHB_incidence_predicted <- NULL
names(selected_data_for_lm)[which(names(selected_data_for_lm) == "adj_inc_pred")] <- "fhb_incidence_predicted"
names(selected_data_for_lm)[which(names(selected_data_for_lm) == "adj_FHB_incidence")] <- "fhb_incidence_predicted"
names(selected_data_for_lm)
### GIBSIM vs. obs -------------------------------------------------------------
names(selected_data_for_lm)
selected_data_for_lm$fhb_incidence_predicted

if(max(selected_data_for_lm$FHB_incidence) < 1.1){
  selected_data_for_lm$FHB_incidence <- selected_data_for_lm$FHB_incidence*100
  print("##### Changed #####")
  print(selected_data_for_lm$FHB_incidence)
}else{
  print(selected_data_for_lm$FHB_incidence)
}


GIBSIM_obs_lm <- lm(FHB_incidence ~ fhb_incidence_predicted, data = selected_data_for_lm)
summary(GIBSIM_obs_lm) # R2 = 0.534
plot(selected_data_for_lm$fhb_incidence_predicted, selected_data_for_lm$FHB_incidence)

GIBSIM_obs_lm_summary <- as.data.frame(summary(GIBSIM_obs_lm)$coefficient)
GIBSIM_obs_lm_summary$variables <- rownames(GIBSIM_obs_lm_summary)
GIBSIM_obs_lm_summary <- GIBSIM_obs_lm_summary %>% relocate("variables")
getwd()

# save
save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/"
save_year <- "2025"
version_pred_obs
if(!dir.exists(paste0(save_dir, save_year, "/",version_pred_obs))){
  dir.create(paste0(save_dir, save_year, "/",version_pred_obs))
}
write_xlsx(GIBSIM_obs_lm_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/only_GIBSIM_incidence_",version_pred_obs, ".xlsx"))

### graph --------------------------
names(MB_intensity_wth_cor_df)

only_GIBSIM_inc_pred_obs_plot <- obs_pred_plot(model = GIBSIM_obs_lm, 
                                               data_contain_X_Y = selected_data_for_lm, 
                                               Y_colname = "FHB_incidence", 
                                               Title = "FHB_incidence")

only_GIBSIM_inc_pred_obs_plot


ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/only_GIBSIM_incidence_",version_pred_obs, ".png")),
       plot = only_GIBSIM_inc_pred_obs_plot,
       width = 10, height = 10, bg = "white", dpi = 1500)


selected_data_for_lm$fhb_incidence_predicted_adjust <- selected_data_for_lm$fhb_incidence_predicted * GIBSIM_obs_lm_summary$Estimate[2] + GIBSIM_obs_lm_summary$Estimate[1]
# summary(lm(FHB_incidence ~ fhb_incidence_predicted_adjust, data = selected_data_for_lm))
# ###############################################################################################################################################
# #####################################################################################################################
# #### Genus level ####################################################################################################
# ### ALL - incidence ------------------------------------------------------------
# taxa_data <- readxl::read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/taxa_info_232425_chosen_OTUITSfull_v1.xlsx")
# unique_taxa <- unique(taxa_data$Genus)
# names(taxa_data)
# selected_data_for_lm_no_MB <- selected_data_for_lm[,-grep("OTU", names(selected_data_for_lm))]
# selected_data_for_lm_MB <- selected_data_for_lm[,grep("OTU", names(selected_data_for_lm))]
# imsi_df <- data.frame(matrix(nrow = nrow(selected_data_for_lm), ncol = 0))
# for(t in c(1:length(unique_taxa))){
#   # t=1
#   imsi_same_genus_df <- taxa_data %>% filter(Genus == unique_taxa[t])
#   imsi_same_genus_df_OTU <- imsi_same_genus_df$OTU
#   
#   imsi_same_genus_df_OTU <- imsi_same_genus_df_OTU[imsi_same_genus_df_OTU %in% names(selected_data_for_lm)]
#   
#   imsi_selected_data_for_lm_MB <- selected_data_for_lm_MB %>% dplyr::select(all_of(imsi_same_genus_df_OTU))
#   imsi_df[,t] <- rowSums(imsi_selected_data_for_lm_MB)
#   
#   names(imsi_df)[t] <- stringr::str_replace(unique_taxa[t], "g__", "")
#   
# }
# sum(duplicated(names(imsi_df))) # if more than 0 is error, check again! 
# 
# selected_data_for_lm_genus <- cbind(selected_data_for_lm_no_MB, imsi_df)
# 
# names(selected_data_for_lm_genus)
# 
# # correlation
# 
# 
# # selected_data_for_selected_MB <- min_max_scaling_for_df(selected_data_for_selected_MB, "FHB_incidence")
# 
# 
# selected_data_for_selected_MB <- imsi_df
# selected_data_for_selected_MB$FHB_incidence <- selected_data_for_lm_genus$FHB_incidence
# selected_data_for_selected_MB$fhb_incidence_predicted_adjust <- selected_data_for_lm_genus$fhb_incidence_predicted
# #correlation -------------------------------------------------------------------
# library(ggplot2)
# library(reshape2)
# library("MASS")
# library(corrplot)
# cor_matrix <- cor(selected_data_for_selected_MB, method = "pearson") #pearson, spearman
# p_values <- cor.mtest(selected_data_for_selected_MB)$p
# significant_cor <- cor_matrix * (p_values <= 0.05)
# cor_heatmap <- corrplot(
#   cor_matrix,
#   method = "color",
#   type = "upper",
#   p.mat = p_values,
#   sig.level = 0.05,
#   insig = "pch", #???pch???, ???p-value???, ???blank???, ???n???, ???label_sig???
#   # order = "hclust",
#   outline = "white",
#   tl.col = "black",
#   tl.srt = 45,
#   col = colorRampPalette(c("blue", "white", "red"))(20), 
#   addCoef.col = "white",
#   number.cex = 0.9, tl.cex = 0.9
# )
# 
# # ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/MB_GIB_incidence_",version_pred_obs, ".png")),
# #        plot = cor_heatmap,
# #        width = 10, height = 10, bg = "white", dpi = 1500)
# 
# # library(PerformanceAnalytics)
# # chart.Correlation(selected_data_for_selected_MB, histogram = TRUE, pch = 19)
# 
# ### high correlated Genus
# cor_MB_FHB_df <- as.data.frame(cor_matrix)
# high_cor_MB_df <- cor_MB_FHB_df %>% filter(FHB_incidence >= 0.4 | FHB_incidence <= -0.4)
# 
# selected_colnames_for_high_cor_genus <- rownames(high_cor_MB_df)
# 
# selected_data_for_selected_MB_genus <- selected_data_for_selected_MB %>% dplyr::select(all_of(selected_colnames_for_high_cor_genus))
# ############ -------------------------------------------------------------------
# library(MASS)
# library(car)
# selected_data_for_selected_MB_genus$fhb_incidence_predicted_adjust
# selected_data_for_selected_MB_genus$FHB_incidence
# hist(selected_data_for_selected_MB_genus$FHB_incidence)
# shapiro.test(selected_data_for_selected_MB_genus$FHB_incidence)
# names(selected_data_for_selected_MB_genus)
# 
# selected_data_for_selected_MB_genus$fhb_incidence_predicted_adjust
# MB_GIB_inc_fullmodel <- lm(FHB_incidence ~ . -fhb_incidence_predicted_adjust + offset(fhb_incidence_predicted_adjust), 
#                            data = selected_data_for_selected_MB_genus)
# 
# # MB_GIB_inc_fullmodel <- lm(FHB_incidence ~ ., data = selected_data_for_selected_MB_genus)
# # MB_GIB_inc_fullmodel <- lm(FHB_incidence ~ ., data = selected_data_for_selected_MB_genus)
# MB_GIB_inc_stepwise <- stepAIC(MB_GIB_inc_fullmodel, direction = "both")
# summary(MB_GIB_inc_fullmodel)
# summary(MB_GIB_inc_stepwise)
# # plot(MB_GIB_inc_stepwise)
# vif(MB_GIB_inc_fullmodel)
# vif(MB_GIB_inc_stepwise)
# plot(selected_data_for_selected_MB_genus$Alternaria, selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$Epicoccum, selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$Sporobolomyces, selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$Sphaerulina, selected_data_for_selected_MB_genus$FHB_incidence)
# summary(lm(FHB_incidence ~ Alternaria, data = selected_data_for_selected_MB_genus))
# 
# # MB_GIB_inc_fullmodel <- lm(FHB_incidence ~ OTU0044 , data = selected_data_for_selected_MB_genus)
# # summary(MB_GIB_inc_fullmodel)
# # vif(MB_GIB_inc_fullmodel)
# MB_GIB_inc_stepwise_summary <- as.data.frame(summary(MB_GIB_inc_stepwise)$coefficient)
# MB_GIB_inc_stepwise_summary$variables <- rownames(MB_GIB_inc_stepwise_summary)
# MB_GIB_inc_stepwise_summary <- MB_GIB_inc_stepwise_summary %>% relocate("variables")
# MB_GIB_inc_stepwise_summary[,c(2:5)] <- round(MB_GIB_inc_stepwise_summary[,c(2:5)],2)
# write_xlsx(MB_GIB_inc_stepwise_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/MB_GIB_incidence_",version_pred_obs, ".xlsx"))
# 
# # MB_GIB_inc_inter <- lm(FHB_incidence ~ OTU0007 + OTU0018 + OTU0026 + OTU0034 + OTU0044 + OTU0045 + OTU0053 + fhb_incidence_predicted , data = selected_data_for_selected_MB_genus)
# # MB_GIB_inc_inter <- lm(FHB_incidence ~ OTU0007 + OTU0018 + OTU0026 + OTU0034 + OTU0045 + OTU0053 + fhb_incidence_predicted , data = selected_data_for_selected_MB_genus)
# # summary(MB_GIB_inc_inter)
# 
# stepwise_GIB_inc_pred_obs_plot <- obs_pred_plot(model = MB_GIB_inc_stepwise, 
#                                                 data_contain_X_Y = selected_data_for_selected_MB_genus,
#                                                 Y_colname = "FHB_incidence",
#                                                 Title = "FHB_incidence")
# 
# ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/MB_GIB_incidence_",version_pred_obs, ".png")),
#        plot = stepwise_GIB_inc_pred_obs_plot,
#        width = 10, height = 10, bg = "white", dpi = 1500)
# 
# plot(selected_data_for_selected_MB_genus$OTU0034,selected_data_for_selected_MB_genus$FHB_incidence)
# 
# library(segmented)
# MB_GIB_inc_stepwise
# segmented_MB_GIB_inc_OTU0003 <- segmented(MB_GIB_inc_stepwise, seg.Z = ~ OTU0003)
# segmented_MB_GIB_inc_OTU0034 <- segmented(MB_GIB_inc_stepwise, seg.Z = ~ OTU0034)
# segmented_MB_GIB_inc_OTU0045 <- segmented(MB_GIB_inc_stepwise, seg.Z = ~ OTU0045)
# plot(segmented_MB_GIB_inc_OTU0003)
# plot(segmented_MB_GIB_inc_OTU0034)
# plot(segmented_MB_GIB_inc_OTU0045)
# 
# library(iml)
# rownames(summary(MB_GIB_inc_stepwise)$coefficients)
# x_column_num <- length(rownames(summary(MB_GIB_inc_stepwise)$coefficients))
# X <-rownames(summary(MB_GIB_inc_stepwise)$coefficients)[2:x_column_num]
# X_var <- selected_data_for_selected_MB_genus[,c(X)]
# Y_var <- selected_data_for_selected_MB_genus$FHB_incidence
# Predictor$new(MB_GIB_inc_stepwise, data = X_var, y = Y_var,type = "regression" )
# 
# 
# hist((1/(selected_data_for_selected_MB_genus$FHB_incidence))*100)
# MB_GIB_glm_inc_fullmodel <- glm((1/(FHB_incidence)*100) ~ . , data = selected_data_for_selected_MB_genus, family = Gamma(link = "log"))
# MB_GIB_glm_stepwise <- stepAIC(MB_GIB_glm_inc_fullmodel, direction = "both")
# summary(MB_GIB_glm_stepwise)
# VIF(MB_GIB_glm_stepwise)
# # plot(MB_wth_glm_stepwise)
# summary(MB_GIB_glm_stepwise)
# 
# stepwise_gamma_GIB_inc_pred_obs_plot <- obs_pred_plot_for_gamma(model = MB_GIB_glm_stepwise,
#                                                                 data_contain_X_Y = selected_data_for_selected_MB_genus,
#                                                                 Y_colname = "FHB_incidence",
#                                                                 Title = "FHB_incidence")
# 
# plot(selected_data_for_selected_MB_genus$OTU0003 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0034 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0045 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0053 , selected_data_for_selected_MB_genus$FHB_incidence)
# 
# plot(selected_data_for_selected_MB_genus$OTU0007 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0026 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0034 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0042 , selected_data_for_selected_MB_genus$FHB_incidence)
# plot(selected_data_for_selected_MB_genus$OTU0044 , selected_data_for_selected_MB_genus$FHB_incidence)
# 
# # ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/gamma_MB_GIB_incidence_",version_pred_obs, ".png")),
# #        plot = stepwise_gamma_GIB_inc_pred_obs_plot,
# #        width = 10, height = 10, bg = "white", dpi = 1500)
# # 
# # summary(MB_GIB_glm_stepwise)
# # stepwise_model_GIB_and_MB_summary <- as.data.frame(summary(MB_GIB_glm_stepwise)$coefficients)
# # stepwise_model_GIB_and_MB_summary$variables <- rownames(stepwise_model_GIB_and_MB_summary)
# # stepwise_model_GIB_and_MB_summary <- stepwise_model_GIB_and_MB_summary %>% relocate("variables")
# # write_xlsx(stepwise_model_GIB_and_MB_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_incidence_",version_pred_obs, ".xlsx"))
# # 
# 
# 

##############################################################################################################
##### OTU level ##########################################################################################################################################
MB_wth_cor_result_df <- read_xlsx(path = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/v6.2/MB_wth_cor_result_v6.2.xlsx")
MB_intensity_cor_result_df <- read_xlsx(path = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/v6.2/MB_intensity_cor_result_v6.2.xlsx")
wth_intensity_cor_result_df <- read_xlsx(path = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/v6.2/wth_intensity_cor_result_v6.2.xlsx")
# MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. MB_correlation/v5/MB_filter_df_results_v5.xlsx")
MB_max_filtered_data <- read_xlsx(path = "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/Filtered_genus_more_ITSfull_v3.2.xlsx")
MB_max_filtered_data <- MB_max_filtered_data[-which(MB_max_filtered_data$Taxa == "unclassified"),]
# select high correlated MB (correlated with incidence or severity) -----------------------------------------------------
high_cor_MB_inc_df <- MB_intensity_cor_result_df %>% dplyr::filter(intensity == "incidence") %>% dplyr::filter(rho >= 0.3| rho <= -0.3) %>% dplyr::filter(MB %in% MB_max_filtered_data$Taxa)
# high_cor_MB_inc_df <- MB_intensity_cor_result_df %>% dplyr::filter(intensity == "incidence") %>% dplyr::filter(rho >= 0.0| rho <= -0.0) %>% dplyr::filter(MB %in% MB_max_filtered_data$Taxa)
MB_list_inc <- c(high_cor_MB_inc_df$MB)

high_cor_MB_sev_df <- MB_intensity_cor_result_df %>% dplyr::filter(intensity == "severity") %>% dplyr::filter(rho >= 0.3 | rho <= -0.3) %>% dplyr::filter(MB %in% MB_max_filtered_data$Taxa)
# high_cor_MB_sev_df <- MB_intensity_cor_result_df %>% dplyr::filter(intensity == "severity") %>% dplyr::filter(rho >= 0.0 | rho <= -0.0) %>% dplyr::filter(MB %in% MB_max_filtered_data$Taxa)
MB_list_sev <- c(high_cor_MB_sev_df$MB)

MB_list <- MB_list_inc
# MB_list <- MB_list_inc[MB_list_inc %in% names(selected_data_for_lm)]
# MB_list <- MB_list_origin[MB_list_origin %in% names(MB_intensity_wth_cor_df)]
# MB_list <- MB_list[!MB_list %in% c("OTU0044", "OTU0055", "OTU0058", "OTU0071", "OTU0024", "OTU0042")]

#### only MB ------------------------------------------------------------------- 
### incidence ------------------------------------------------------------------
### MB - incidence ------------------------------------------------------------
# selected_data_for_only_selected_MB <- selected_data_for_lm[,c((MB_list_inc), "FHB_incidence")]
# selected_data_for_only_selected_MB <- min_max_scaling_for_df(selected_data_for_only_selected_MB, "FHB_incidence")

################################################################################
### only MB --------------------------------------------------------------------
library(MASS)
### ALL - incidence ------------------------------------------------------------
selected_data_only_MB <- selected_data_for_lm[,c((MB_list_inc), "FHB_incidence")]
names(selected_data_only_MB)
# selected_data_for_selected_MB_genus <- min_max_scaling_for_df(selected_data_for_selected_MB_genus, "FHB_incidence")

# # CLR for MB data ----
X_mb_raw <- selected_data_only_MB[,c((MB_list_inc))] 
X_mb_raw <- as.matrix(X_mb_raw)
# X_mb_clr <- compositions::clr(X_mb_raw + 1e-6)
# X_mb_clr_df <- as.data.frame(scale(X_mb_clr))

# min(X_mb_clr)

# interaction information : netwrok edges data --- 
# edges_version <- "all_v1" #You should get edges data from that code -->  D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/8. MB_analysis_phyloseq_data/phyloseq_wheat_2025/script/script_for_all_data/6. network_SpiecEasi_v2.R
# edges_df <- readxl::read_xlsx(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/8. MB_analysis_phyloseq_data/phyloseq_wheat_2025/output/v5/network_edges_genus_all_v5.xlsx"))
edges_df <- readxl::read_xlsx(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/8. MB_analysis_phyloseq_data/phyloseq_wheat_2025/output/v5/network_edges_genus_rho_0.3_v5.xlsx"))


# make interaction term ---
MB_interaction_df <- data.frame(matrix(nrow = nrow(X_mb_raw), ncol = 0))
for(k in c(1:nrow(edges_df))){
  # k=1
  i <- edges_df$from[k] 
  j <- edges_df$to[k]
  
  i <- stringr::str_replace(i, "g__", "")
  j <- stringr::str_replace(j, "g__", "")
  
  if(i %in% MB_list_inc && j %in% MB_list_inc){
    imsi_interaction_vector <- X_mb_raw[,i]*X_mb_raw[,j]
    MB_interaction_df$imsi_col <- imsi_interaction_vector
    names(MB_interaction_df)[which(names(MB_interaction_df) == "imsi_col")] <- paste0(i, "_", j)
  }
  
}


# max value of FHB incidence data should be less then 1.0000001 -----------------------
if(max(selected_data_only_MB$FHB_incidence) < 1.001){
  selected_data_only_MB$FHB_incidence <- selected_data_only_MB$FHB_incidence*100
  print("##### change #####")
  print(selected_data_only_MB$FHB_incidence)
  print("##### change #####")
}else{
  print(selected_data_only_MB$FHB_incidence)
  print("##### Not change #####")
}



# data for model -------------------------------------------------------------------------------------------------------------
data_only_MB_for_model <- cbind(as.data.frame(X_mb_raw), MB_interaction_df)
data_only_MB_for_model$FHB_incidence <- selected_data_only_MB$FHB_incidence

# data_only_MB_for_model$FHB_incidence[data_only_MB_for_model$FHB_incidence == 0] <- 0.000001
names(data_only_MB_for_model)
MB_list
# modeling
only_MB_model <- lm(
  FHB_incidence ~ .,
  data = data_only_MB_for_model
)
# nzv_index <- nearZeroVar(data_only_MB_for_model, saveMetrics = TRUE)

summary(only_MB_model)

only_MB_inc_stepwise <- stepAIC(only_MB_model, direction = "both")
summary(only_MB_model)
summary(only_MB_inc_stepwise)
# plot(only_MB_inc_stepwise)
library(car)
vif(only_MB_inc_stepwise)


only_MB_inc_stepwise_summary <- as.data.frame(summary(only_MB_inc_stepwise)$coefficient)
only_MB_inc_stepwise_summary$variables <- rownames(only_MB_inc_stepwise_summary)
only_MB_inc_stepwise_summary <- only_MB_inc_stepwise_summary %>% dplyr::relocate("variables")
only_MB_inc_stepwise_summary[,c(2:5)] <- round(only_MB_inc_stepwise_summary[,c(2:5)],2)
only_MB_inc_stepwise_summary <- only_MB_inc_stepwise_summary %>% dplyr::mutate(`Pr(>|t|)` <= 0.05)

# save
save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/"
save_year <- "2025"
write_xlsx(only_MB_inc_stepwise_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/only_MB_incidence_",version_pred_obs, ".xlsx"))

### graph -----------------
stepwise_only_MB_inc_pred_obs_plot <- obs_pred_plot_for_minus_value(model = only_MB_inc_stepwise, 
                                                                    data_contain_X_Y = data_only_MB_for_model,
                                                                    Y_colname = "FHB_incidence",
                                                                    Title = "FHB_incidence")

ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/only_MB_incidence_",version_pred_obs, ".png")),
       plot = stepwise_only_MB_inc_pred_obs_plot,
       width = 10, height = 10, bg = "white", dpi = 1500)

##################################################################################################################################
### MB + GIBSIM --------------------------------------------------------------------
data_only_MB_for_model$fhb_incidence_predicted <- NULL
data_gibsim_MB_for_model <- data_only_MB_for_model
data_gibsim_MB_for_model$fhb_incidence_predicted <- selected_data_for_lm$fhb_incidence_predicted
data_gibsim_MB_for_model$FHB_incidence <- selected_data_for_lm$FHB_incidence
# selected_data_MB_GIB <- selected_data_for_lm[,c((MB_list_inc),"fhb_incidence_predicted", "FHB_incidence")]
names(data_gibsim_MB_for_model)
# selected_data_for_selected_MB_genus <- min_max_scaling_for_df(selected_data_for_selected_MB_genus, "FHB_incidence")

library(MASS)

if(max(data_gibsim_MB_for_model$FHB_incidence) < 1.1){
  data_gibsim_MB_for_model$FHB_incidence <- data_gibsim_MB_for_model$FHB_incidence*100
  print("##### change #####")
  print(data_gibsim_MB_for_model$FHB_incidence)
  print("##### change #####")
}else{
  print(data_gibsim_MB_for_model$FHB_incidence)
  print("##### Not change #####")
}

writexl::write_xlsx(data_gibsim_MB_for_model, 
                    paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/", 
                           version_pred_obs, 
                           "/data_gibsim_MB_for_model_",version_pred_obs,".xlsx"))
# make model 
GIBSIM_MB_fullmodel <- lm(FHB_incidence ~ ., data = data_gibsim_MB_for_model)
GIBSIM_MB_inc_stepwise <- stepAIC(GIBSIM_MB_fullmodel, direction = "both")

# coefficent of fhb_incidence_predicted should be 1
library(restriktor)
constraints <- 'fhb_incidence_predicted == 1'
restricted_model <- restriktor(GIBSIM_MB_inc_stepwise, constraints = constraints)

summary(GIBSIM_MB_fullmodel)
summary(GIBSIM_MB_inc_stepwise)
summary(restricted_model)
# plot(GIBSIM_MB_inc_stepwise)

GIBSIM_MB_inc_stepwise_summary <- as.data.frame(summary(GIBSIM_MB_inc_stepwise)$coefficient)
GIBSIM_MB_inc_stepwise_summary$variables <- rownames(GIBSIM_MB_inc_stepwise_summary)
GIBSIM_MB_inc_stepwise_summary <- GIBSIM_MB_inc_stepwise_summary %>% relocate("variables")
GIBSIM_MB_inc_stepwise_summary$star <- ifelse(GIBSIM_MB_inc_stepwise_summary$`Pr(>|t|)` < 0.05, "*", "") #GIBSIM_MB_inc_restricted_summary %>% mutate(`Pr(>|t|)` <= 0.05)
GIBSIM_MB_inc_stepwise_summary[,c(2:5)] <- round(GIBSIM_MB_inc_stepwise_summary[,c(2:5)],3)

vif(GIBSIM_MB_inc_stepwise)
# plot(GIBSIM_MB_inc_stepwise)
GIBSIM_MB_inc_restricted_summary <- as.data.frame(summary(restricted_model)$coefficient)
GIBSIM_MB_inc_restricted_summary$variables <- rownames(GIBSIM_MB_inc_restricted_summary)
GIBSIM_MB_inc_restricted_summary <- GIBSIM_MB_inc_restricted_summary %>% relocate("variables")
GIBSIM_MB_inc_restricted_summary$star <- ifelse(GIBSIM_MB_inc_restricted_summary$`Pr(>|t|)` < 0.05, "*", "") #GIBSIM_MB_inc_restricted_summary %>% mutate(`Pr(>|t|)` <= 0.05)
GIBSIM_MB_inc_restricted_summary[,c(2:5)] <- round(GIBSIM_MB_inc_restricted_summary[,c(2:5)],2)
# save
save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/"
save_year <- "2025"
write_xlsx(GIBSIM_MB_inc_stepwise_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_stepwise_incidence_",version_pred_obs, ".xlsx"))
write_xlsx(GIBSIM_MB_inc_restricted_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_restricted_incidence_",version_pred_obs, ".xlsx"))

### graph -----------------
stepwise_GIBSIM_MB_inc_pred_obs_plot <- obs_pred_plot_for_minus_value(model = GIBSIM_MB_inc_stepwise, 
                                                                      data_contain_X_Y = data_gibsim_MB_for_model,
                                                                      Y_colname = "FHB_incidence",
                                                                      Title = "FHB_incidence")


restricted_GIBSIM_MB_inc_pred_obs_plot <- obs_pred_plot_for_minus_value_for_restricted_model(model = restricted_model, 
                                                                                             data_contain_X_Y = data_gibsim_MB_for_model,
                                                                                             Y_colname = "FHB_incidence",
                                                                                             Title = "FHB_incidence")

ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_stepwise_incidence_",version_pred_obs, ".png")),
       plot = stepwise_GIBSIM_MB_inc_pred_obs_plot,
       width = 10, height = 10, bg = "white", dpi = 1500)

ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_restricted_incidence_",version_pred_obs, ".png")),
       plot = restricted_GIBSIM_MB_inc_pred_obs_plot,
       width = 10, height = 10, bg = "white", dpi = 1500)

### GLM (gamma) ========================================================================================================================
library(MASS)        # stepAIC
library(car)         # vif, influencePlot
library(compositions) # ilr / acomp
library(glmnet)      # lasso (??????)
library(betareg)     # (??????: beta regression)
data_gibsim_MB_for_gamma <- data_gibsim_MB_for_model

# 1) ??????변???(???부검???): 0 ??? ?????? ??? 처리 (Gamma??? 0 ?????? ??????)
min(data_gibsim_MB_for_gamma$FHB_incidence)
table_zero <- table(data_gibsim_MB_for_gamma$FHB_incidence == 0)
print(table_zero)
if(any(data_gibsim_MB_for_gamma$FHB_incidence == 0)){
  # 0??? ???물다??? ?????? ????????? ?????? (????????? ??????가 %?????? ?????? ?????? ???)
  minpos <- min(data_gibsim_MB_for_gamma$FHB_incidence[data_gibsim_MB_for_gamma$FHB_incidence > 0], na.rm = TRUE)
  offset <- minpos/5
  message("Zero values detected in FHB_incidence. Replacing zeros with offset = ", offset)
  data_gibsim_MB_for_gamma$FHB_incidence_adj <- ifelse(data_gibsim_MB_for_gamma$FHB_incidence == 0, offset, data_gibsim_MB_for_gamma$FHB_incidence)
} else {
  data_gibsim_MB_for_gamma$FHB_incidence_adj <- data_gibsim_MB_for_gamma$FHB_incidence
}

predictors <- setdiff(names(data_gibsim_MB_for_gamma), c("FHB_incidence","FHB_incidence_adj"))
formula_full <- as.formula(paste("FHB_incidence_adj ~", paste(predictors, collapse = " + ")))
class(formula_full)

#glm (gamma model)
glm_gamma_full <- glm(formula_full, family = Gamma(link = "log"), data = data_gibsim_MB_for_gamma)
summary(glm_gamma_full)

vif_try <- try(vif(glm_gamma_full), silent = TRUE)
if(!inherits(vif_try, "try-error")){
  print(vif_try)
} else {
  message("VIF 계산 ?????? (모델??? ??????공선??? ?????? 기??? 문제 가???).")
}

# stepwise gamma
glm_gamma_step <- stepAIC(glm_gamma_full, direction = "both", trace = TRUE)
summary(glm_gamma_step)

vif_try_step <- try(vif(glm_gamma_step), silent = TRUE)
if(!inherits(vif_try_step, "try-error")){
  print(vif_try_step)
} else {
  message("VIF 계산 ?????? (모델??? ??????공선??? ?????? 기??? 문제 가???).")
}

# 7) ????????? 모델 진단 ??? ??????
df <- data.frame(matrix(ncol = 0, nrow = length(data_gibsim_MB_for_model$FHB_incidence)))
df$FHB_incidence <- data_gibsim_MB_for_model$FHB_incidence
df$pred_gamma_step <- as.vector(predict(glm_gamma_step, type = "response"))
df$pred_gamma_step2 <- ifelse(df$pred_gamma_step > 100, 100, df$pred_gamma_step)
gamma_result <- summary(lm(FHB_incidence ~ pred_gamma_step2, data = df))
gamma_result

# save
GIBSIM_MB_inc_gamma_summary <- as.data.frame(summary(glm_gamma_step)$coefficient)
GIBSIM_MB_inc_gamma_summary$variables <- rownames(GIBSIM_MB_inc_gamma_summary)
GIBSIM_MB_inc_gamma_summary <- GIBSIM_MB_inc_gamma_summary %>% relocate("variables")
GIBSIM_MB_inc_gamma_summary$star <- ifelse(GIBSIM_MB_inc_gamma_summary$`Pr(>|t|)` < 0.05, "*", "") #GIBSIM_MB_inc_gamma_summary %>% mutate(`Pr(>|t|)` <= 0.05)
GIBSIM_MB_inc_gamma_summary[,c(2:5)] <- round(GIBSIM_MB_inc_gamma_summary[,c(2:5)],2)


# save
save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/"
save_year <- "2025"
write_xlsx(GIBSIM_MB_inc_stepwise_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_gamma_stepwise_incidence_",version_pred_obs, ".xlsx"))

gamma_GIBSIM_MB_inc_pred_obs_plot <- obs_pred_plot_for_gamma3(model = glm_gamma_step, 
                                                              data_contain_X_Y = data_gibsim_MB_for_gamma,
                                                              Y_colname = "FHB_incidence",
                                                              Title = "FHB_incidence")
ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_gamma_stepwise_incidence_",version_pred_obs, ".png")),
       plot = gamma_GIBSIM_MB_inc_pred_obs_plot,
       width = 10, height = 10, bg = "white", dpi = 1500)


# ??????/분산(과다분산) 검???
dispersion <- sum(residuals(glm_gamma_step, type = "deviance")^2) / glm_gamma_step$df.residual
message("Dispersion (deviance/df): ", round(dispersion, 3))
plot(glm_gamma_step$fitted.values, residuals(glm_gamma_step, type="deviance"),
     xlab="Fitted", ylab="Deviance residuals"); abline(h=0,lty=2)

try(influencePlot(glm_gamma_step), silent = TRUE)


# #### recheck : select p-val <= 0.05 MB
# recheck_model <- lm(FHB_incidence ~ OTU0004_OTU0038 + OTU0006 + OTU0026_OTU0038 + OTU0038 + fhb_incidence_predicted, data = data_gibsim_MB_for_model)
# recheck_restricted_model <- restriktor(recheck_model, constraints = constraints)
# summary(recheck_restricted_model)
# recheck_model_graph <- obs_pred_plot_for_minus_value_for_restricted_model(recheck_restricted_model, 
#                                                                           data_contain_X_Y = data_gibsim_MB_for_model,
#                                                                           Y_colname = "FHB_incidence",
#                                                                           Title = "FHB_incidence")
# plot(recheck_model)
# 
# plot(data_gibsim_MB_for_model$OTU0026, data_gibsim_MB_for_model$FHB_incidence)

# hist(1/(selected_data_only_MB$FHB_incidence))
# MB_glm_inc_fullmodel <- glm((1/(FHB_incidence)) ~ . , data = selected_data_only_MB, family = Gamma(link = "log"))
# MB_glm_inc_stepwise <- stepAIC(MB_glm_inc_fullmodel, direction = "both")
# summary(MB_glm_inc_stepwise)
# VIF(MB_glm_inc_stepwise)
# # plot(MB_wth_glm_stepwise)
# stepwise_gamma_MB_inc_pred_obs_plot <- obs_pred_plot_for_gamma(model = MB_glm_inc_stepwise, 
#                                                                data_contain_X_Y = selected_data_only_MB, 
#                                                                Y_colname = "FHB_incidence", 
#                                                                Title = "FHB_incidence")
# 
# ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/gamma_only_MB_incidence_",version_pred_obs, ".png")),
#        plot = stepwise_gamma_MB_inc_pred_obs_plot,
#        width = 10, height = 10, bg = "white", dpi = 1500)
# getwd()
# 
# summary(MB_glm_inc_stepwise)
# stepwise_model_only_selected_MB_summary <- as.data.frame(summary(MB_glm_inc_stepwise)$coefficients)
# stepwise_model_only_selected_MB_summary$variables <- rownames(stepwise_model_only_selected_MB_summary)
# stepwise_model_only_selected_MB_summary <- stepwise_model_only_selected_MB_summary %>% relocate("variables")
# write_xlsx(stepwise_model_only_selected_MB_summary, path = paste0(save_dir, save_year, "/",version_pred_obs, "/only_MB_incidence_",version_pred_obs, ".xlsx"))

##########################################################################################################################
##########################################################################################################################
###### scatter plots
names(data_gibsim_MB_for_model)
X_var_df <- data_gibsim_MB_for_model %>% dplyr::select(!(FHB_incidence))
names(X_var_df)
imsi_incidence_list <- list()
for(r in c(1:length(names(X_var_df)))){
  #r=1
  data_for_scatter <- data_gibsim_MB_for_model[, c(names(X_var_df)[r], "FHB_incidence")]
  
  # ggplot_incidence
  G_inci_relation <- ggplot(data=data_for_scatter, aes(x=get(names(X_var_df)[r]), y=FHB_incidence)) + 
    geom_point(size=2, colour="black") + # shape 15: solid square
    labs(title = paste0("Incidence vs. ", names(X_var_df)[r] ,"\n"), 
         x = paste0("\n",names(X_var_df)[r]),
         y = "Incidence\n") +
    #geom_text(aes(label=Freq), hjust=-0.1, color="black", size=2.9)+ 
    theme_classic() + 
    theme(axis.text = element_text(size=18),
          axis.title = element_text(size=30)) + 
    theme(axis.title.x = element_text(size = 30,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 30,hjust = 0.5, face='bold'))
  
  
  # saving incidence plot
  value_type = ""
  assign(paste0("p_", value_type), G_inci_relation)
  
  ggsave(
    plot = get(paste0("p_", value_type)),
    file = paste0(
      "./Output/3. obs_pred_graph/2025/",
      version_pred_obs,
      "/scatter_plot/scatter_Incidence_", 
      names(X_var_df)[r]
      ,".jpg"
    ),
    width = 15,
    height = 15,
    units = c("cm")
  )
  
  
  imsi_incidence_list[[r]] <- G_inci_relation
  
}

# combined_inc_scatter <- grid.arrange(imsi_incidence_list[[1]], imsi_incidence_list[[2]], imsi_incidence_list[[3]], 
#                                   imsi_incidence_list[[4]], imsi_incidence_list[[5]], imsi_incidence_list[[6]],
#                                   imsi_incidence_list[[7]], imsi_incidence_list[[8]], imsi_incidence_list[[9]],
#                                   imsi_incidence_list[[10]], #imsi_incidence_list[[21]], #imsi_graph_list_pdp[[12]],
#                                   #imsi_graph_list_pdp[[13]], imsi_graph_list_pdp[[14]], imsi_graph_list_pdp[[15]],
#                                   #imsi_graph_list_pdp[[16]], imsi_graph_list_pdp[[17]],
#                                   ncol=5)


##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
### ML -------------------------------------------------------------------------
library(caret)
library(randomForest)

### data -----------------------------------------------------------------------
names(data_gibsim_MB_for_model)
selected_data_OTU_level <- data_gibsim_MB_for_model #selected_data_for_lm[,c(6, 8:143, 145)]
names(selected_data_OTU_level)
MB_vec_length <- length(MB_list)
FHB_incidence_num <- which(names(selected_data_OTU_level) == "FHB_incidence")
fhb_incidence_predicted_num <- which(names(selected_data_OTU_level) == "fhb_incidence_predicted")

using_data <- selected_data_OTU_level[,c(1:MB_vec_length, FHB_incidence_num,fhb_incidence_predicted_num)] # no interaction term
names(using_data)
names(using_data)[which(names(using_data) == "fhb_incidence_predicted")] <- "GIBSIM_pred"

### for loop : RF --------------------------------------------------------------
n_folds <- 5   # fold ??? ??????
folds <- cut(seq(1, nrow(using_data)), breaks = n_folds, labels = FALSE)  # ???????????? n_folds??? 분할
set.seed(1)
sample(folds)
RF_result <- data.frame(matrix(ncol = 4))

mtry <- c(2,4,6,8,10)
ntree <- c(200, 300, 400, 500, 750, 1000)

for(mt in c(1:length(mtry))){
  #mt=1
  print("mtry")
  print(mtry[mt])
  for(nt in c(1:length(ntree))){
    #nt=1
    print("ntree")
    print(ntree[nt])
    
    imsi_cv_result <- data.frame()
    imsi_cv_result_avg <- data.frame()
    for (i in 1:n_folds) {
      #i=1
      # train set??? validation set ??????
      train_indices <- which(folds != i)
      train_data <- using_data[train_indices, ]
      validation_data <- using_data[-train_indices, ]
      
      # Random Forest 모델 ??????
      set.seed(i)
      rf_model <- randomForest(FHB_incidence ~ ., mtry = mtry[mt], ntree = ntree[nt],  data = train_data)
      
      # train set ???가 (모델 ????????? ????????? ???????????? ???가)
      train_pred <- predict(rf_model, train_data)
      train_r_squared <- stats::cor(train_data$FHB_incidence, train_pred)^2
      
      # validation set ???가
      validation_pred <- predict(rf_model, validation_data)
      validation_r_squared <- stats::cor(validation_data$FHB_incidence, validation_pred)^2
      
      # 결과 출력
      cat("Fold", i, "- Train set R^2:", train_r_squared, "- Validation set R^2:", validation_r_squared, "\n")
      
      
      imsi_cv_result[i,1] <- train_r_squared
      imsi_cv_result[i,2] <- validation_r_squared
    }
    
    imsi_cv_result_avg[1,1] <- mtry[mt]
    imsi_cv_result_avg[1,2] <- ntree[nt]
    imsi_cv_result_avg[1,3] <- mean(imsi_cv_result$V1) #train_r_squared
    imsi_cv_result_avg[1,4] <- mean(imsi_cv_result$V2) #validation_r_squared
    
    names(imsi_cv_result_avg) <- names(RF_result) <- c("mtry", "ntree", "train_R2", "val_R2")
    RF_result <- rbind(RF_result, imsi_cv_result_avg)
  }
}

RF_result <- RF_result[complete.cases(RF_result), ]
rf_best_parameter_row <- which(RF_result$val_R2 == max(RF_result$val_R2))
rf_best_parameter <- RF_result[rf_best_parameter_row,]
rf_best_parameter


writexl::write_xlsx(RF_result, path = paste0(save_dir, "2025/", version_pred_obs,"/RF_parameter_results_", version_pred_obs, ".xlsx")) #"./output/find_best_parameter/RF_result_without_zero_v1.xlsx")

###find important variable ----------------------------------------------------
RF_important_result_df <- as.data.frame(matrix(ncol = 0,nrow = ncol(using_data)-1)) # ncol(using_data)-1 = number of X variables
RF_important_result_df$X_var <- names(using_data)[-which(names(using_data) == "FHB_incidence")] #X variables
repeat_num <- 1000
for(seed in c(1:repeat_num)){
  #seed=1
  print(seed)
  set.seed(seed)
  best_rf_model <- randomForest(FHB_incidence ~ ., mtry = rf_best_parameter$mtry, ntree = rf_best_parameter$ntree,  data = using_data)
  
  # varImpPlot(best_rf_model)
  imsi_RF_important_result <- data.frame(varImpPlot(best_rf_model)[,1])
  imsi_RF_important_result$X <- rownames(imsi_RF_important_result)
  names(imsi_RF_important_result) <- c(paste0("Round_",seed), "X_var")
  
  ## merge data
  RF_important_result_df <- left_join(RF_important_result_df, imsi_RF_important_result, by = "X_var")
  # writexl::write_xlsx(RF_important_result, path = "./output/find_best_parameter/RF_important_result_with_zero_v1.xlsx")
  
}
RF_important_result_df_t <- as.data.frame(t(RF_important_result_df))
names(RF_important_result_df_t) <- RF_important_result_df_t[1,]
RF_important_result_df_final <- RF_important_result_df_t[-1,]

for(f in c(1:ncol(RF_important_result_df_final))){
  RF_important_result_df_final[,f] <- as.numeric(RF_important_result_df_final[,f])
}
# boxplot(RF_important_result_df_final)

library(ggbeeswarm)
library(tidyr)
df_long <- pivot_longer(RF_important_result_df_final, 
                        cols = everything(), 
                        names_to = "Variable", 
                        values_to = "Importance")

# swarm plot 그리???
swarm_g <- ggplot(df_long, aes(x = fct_reorder(Variable, Importance, .fun = mean),
                               y = Importance)) +
  geom_quasirandom(size = 1.5, bandwidth = 0.1, color = "#353535", alpha=0.2,  shape = 16, stroke = 0) + ##58ACFA ##353535
  stat_summary(fun = mean, geom = "point", color = "orange", size = 10, alpha = 1, stroke = 0, shape = "*") + 
  coord_flip() +
  theme_classic() +
  labs(x = "Variables\n", y = "\nVariable Importance") +
  # theme(axis.text.y = element_text(size = 10)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=15, face='bold',color='black'))+ #vjust=0.4,
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))

swarm_g

ggsave(
  plot = swarm_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/important_swarm_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 20,
  height = 20,
  units = c("cm")
)
#########################################################################################################
### PDP
names(using_data)
set.seed(1)
rf_model_for_analysis <- randomForest(FHB_incidence ~ ., mtry = rf_best_parameter$mtry, ntree = rf_best_parameter$ntree,  data = using_data)

# estimate train set 
RF_pred <- predict(rf_model_for_analysis, using_data)
RF_r_squared <- stats::cor(using_data$FHB_incidence, RF_pred)^2 #0.96
plot(using_data$FHB_incidence, RF_pred)


## goodness-of-fit of RF model
goodness_of_fit_RF <- ggplot(data=using_data, aes(x=RF_pred, y=FHB_incidence)) + 
  geom_point(size=2, colour="black") + # shape 15: solid square
  labs(title = paste0("Goodness-of-fit of RF model\n", "R2 = ", round(RF_r_squared,2)), 
       x = paste0("\nPredicted incidence"),
       y = "Observed incidence\n") +
  #geom_text(aes(label=Freq), hjust=-0.1, color="black", size=2.9)+ 
  theme_classic() + 
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=20)) + 
  theme(axis.title.x = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 20,hjust = 0.5, face='bold'))


# saving incidence plot
ggsave(
  plot = goodness_of_fit_RF,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/goodness_of_fit_RF", 
    version_pred_obs
    ,".png"
  ),
  width = 15,
  height = 15,
  units = c("cm")
)


library(pdp)
names(using_data)
max(using_data$FHB_incidence)
hist(using_data$FHB_incidence)

# pdp_OTU0004 <- pdp::partial(rf_model_for_analysis, pred.var = "OTU0004", train = using_data)
# plot(pdp_OTU0004)


# imp_variables_Top <- c("fhb_incidence_predicted","OTU0004_OTU0038", "OTU0004", "OTU0004_OTU0026", "OTU0027", "OTU0006", 
#                        "OTU0038", "OTU0056", "OTU0075")

imp_variables_Top <- names(using_data)[-which(names(using_data) == "FHB_incidence")]
length(imp_variables_Top)
imsi_graph_list_pdp <- list()
for(p in c(1:length(imp_variables_Top))){
  #p=6
  imp_var <- imp_variables_Top[p]
  imsi_pdp <- pdp::partial(rf_model_for_analysis, pred.var = imp_var, train = using_data)
  
  title_lab <- paste0("PDP : ",imp_var,"\n") 
  xlab <- paste0("\n",imp_var,"\n") 
  
  imsi_graph <- ggplot(imsi_pdp, aes_string(x = imp_var, y = "yhat")) +
    geom_line(color = "#303030", linewidth = 1) +
    labs(
      title = title_lab,
      x = xlab,
      y = "\nFHB incidence\n" #PDP of incidence # y hat
    ) +
    theme_classic() #+
  # theme(
  #   axis.text = element_text(size = 15),
  #   axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
  #   axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
  #   plot.title = element_text(size = 17, face = 'bold')
  # )
  
  
  imsi_graph_list_pdp[[p]] <- imsi_graph
}

# combind graph
# combined_plot_PDP <- grid.arrange(imsi_graph_list_pdp[[1]], imsi_graph_list_pdp[[2]], imsi_graph_list_pdp[[3]], 
#                                   imsi_graph_list_pdp[[4]], imsi_graph_list_pdp[[5]], imsi_graph_list_pdp[[6]],
#                                   ncol=3)

combined_plot_PDP <- grid.arrange(imsi_graph_list_pdp[[1]], imsi_graph_list_pdp[[2]], imsi_graph_list_pdp[[3]], 
                                  imsi_graph_list_pdp[[4]], imsi_graph_list_pdp[[5]], imsi_graph_list_pdp[[6]],
                                  imsi_graph_list_pdp[[7]], #imsi_graph_list_pdp[[8]], imsi_graph_list_pdp[[9]],
                                  # imsi_graph_list_pdp[[10]], imsi_graph_list_pdp[[11]], imsi_graph_list_pdp[[12]],
                                  # imsi_graph_list_pdp[[13]], imsi_graph_list_pdp[[14]], imsi_graph_list_pdp[[15]],
                                  # imsi_graph_list_pdp[[16]], imsi_graph_list_pdp[[17]],imsi_graph_list_pdp[[18]],
                                  ncol=4)

ggsave(
  plot = combined_plot_PDP,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/combined_PDP_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 25,
  height = 12,
  units = c("cm")
)
### ================================================================================================
### ICE
imp_variables_Top <- names(using_data)[-which(names(using_data) == "FHB_incidence")]
length(imp_variables_Top)
imsi_graph_list_ice <- list()
for(p in c(1:length(imp_variables_Top))){
  #p=6
  imp_var <- imp_variables_Top[p]
  imsi_ice <- pdp::partial(rf_model_for_analysis, pred.var = imp_var, train = using_data, ice = T)
  # plotPartial(ice, alpha = 0.1)
  
  title_lab <- paste0("ICE : ",imp_var,"\n") 
  xlab <- paste0("\n",imp_var,"\n") 
  
  imsi_graph <- autoplot(imsi_ice, center = TRUE, alpha = 0.2, rug = F, train = using_data) + 
    # geom_line(color = "#303030", linewidth = 1) +
    labs(
      title = title_lab,
      x = xlab,
      y = "\nFHB incidence\n" #ice of incidence # y hat
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title.x = element_text(size = 12, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 12, hjust = 0.5, face = 'bold'),
      plot.title = element_text(size = 14, face = 'bold')
    )
  
  
  imsi_graph_list_ice[[p]] <- imsi_graph
}

# combind graph
# combined_plot_ice <- grid.arrange(imsi_graph_list_ice[[1]], imsi_graph_list_ice[[2]], imsi_graph_list_ice[[3]], 
#                                   imsi_graph_list_ice[[4]], imsi_graph_list_ice[[5]], imsi_graph_list_ice[[6]],
#                                   ncol=3)

combined_plot_ice <- grid.arrange(imsi_graph_list_ice[[1]], imsi_graph_list_ice[[2]], imsi_graph_list_ice[[3]], 
                                  imsi_graph_list_ice[[4]], imsi_graph_list_ice[[5]], imsi_graph_list_ice[[6]],
                                  imsi_graph_list_ice[[7]], #imsi_graph_list_ice[[8]], imsi_graph_list_ice[[9]],
                                  # imsi_graph_list_ice[[10]], #imsi_graph_list_ice[[11]], #imsi_graph_list_ice[[12]],
                                  # imsi_graph_list_ice[[13]], imsi_graph_list_ice[[14]], imsi_graph_list_ice[[15]],
                                  # imsi_graph_list_ice[[16]], imsi_graph_list_ice[[17]],imsi_graph_list_ice[[18]],
                                  ncol=4)


ggsave(
  plot = combined_plot_ice,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/combined_ICE_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 28,
  height = 13,
  units = c("cm")
)
### ============================================================================================================
### ============================================================================================================
### Section 2 : ALE ############################################################
### data : with_zero_df // model : rf_model_with_zero --------------------------
data_for_iml <- using_data
model_for_iml <- rf_model_for_analysis # best_brt_model | rf_model_with_zero
plot(using_data$Alternaria, using_data$Cercospora)

### X variable
names(data_for_iml)
X_df <- data_for_iml[which(names(data_for_iml) != "FHB_incidence")]
X_df <- data.frame(X_df)

names(data_for_iml)[which(names(data_for_iml) == "FHB_incidence")] <- "Y"

predict_function <- function(model, newdata) {
  preds <- predict(model, newdata)
  as.numeric(preds)
  return(preds)
}

library(iml)
predictor <- Predictor$new(model_for_iml, data = X_df, y = data_for_iml$Y,  predict.fun = predict_function)


### ALE
# what is the most important variable?
#imp_variables_Top5 <- c("rhum_after_flower_10day", "prcp_before_hd_10day", "prcp_hd_10day", "rhum_after_flower_10day", "rhum_flower_10day")
# imp_variables_Top <- c("fhb_incidence_predicted","OTU0004_OTU0038", "OTU0004", "OTU0004_OTU0026", "OTU0027", "OTU0006")
imp_variables_Top <- names(using_data)[-which(names(using_data) == "FHB_incidence")]

imsi_graph_list <- list()
for(p in c(1:length(imp_variables_Top))){
  #p=6
  imp_var <- imp_variables_Top[p]
  ale <- FeatureEffect$new(predictor, feature = imp_var, method = "ale")
  
  title_lab <- paste0(imp_var,"\n") 
  xlab <- paste0("\n",imp_var,"\n") 
  
  imsi_graph <- ale$plot() +
    labs(title = title_lab, 
         x = xlab,
         y = "\nALE of incidence\n") +
    geom_line(col = "#303030", linewidth = 1) + 
    theme_classic()+
    theme(axis.text = element_text(size=11)) + #, color = "black"
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold')) +
    theme(plot.title = element_text(size = 13,face='bold'))
  
  imsi_graph_list[[p]] <- imsi_graph
}

# combind graph
# combined_plot_ALE <- grid.arrange(imsi_graph_list[[1]], imsi_graph_list[[2]], imsi_graph_list[[3]], 
#                                   imsi_graph_list[[4]], imsi_graph_list[[5]], imsi_graph_list[[6]],
#                                   ncol=3)

combined_plot_ALE <- grid.arrange(imsi_graph_list[[1]], imsi_graph_list[[2]], imsi_graph_list[[3]], 
                                  imsi_graph_list[[4]], imsi_graph_list[[5]], imsi_graph_list[[6]],
                                  imsi_graph_list[[7]], #imsi_graph_list[[8]], imsi_graph_list[[9]],
                                  # imsi_graph_list[[10]], imsi_graph_list[[11]], imsi_graph_list[[12]],
                                  # imsi_graph_list[[13]], imsi_graph_list[[14]], imsi_graph_list[[15]],
                                  # imsi_graph_list[[16]], imsi_graph_list[[17]],imsi_graph_list[[18]],
                                  ncol=4)
ggsave(
  plot = combined_plot_ALE,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/combined_ALE_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 25,
  height = 12,
  units = c("cm")
)
### ============================================================================================================
### ============================================================================================================

###SHAP
library(iml)

# Predictor 객체 ??????
X <- using_data[, -which(names(using_data) == "FHB_incidence")]  
names(X)[which(names(X) == "fhb_incidence_predicted" )] <- "GIBSIM_pred"
predictor <- Predictor$new(rf_model_for_analysis, data = X, y = using_data$FHB_incidence)

# SHAP ??? 계산
shap <- Shapley$new(predictor, x.interest = X[1, ])
plot(shap)


### interaction ------------------
interaction <- Interaction$new(predictor)
plot(interaction)
interaction_Hannaella <- Interaction$new(predictor, feature = "Hannaella")
# plot(interaction_Cercospora)
interaction_Alternaria <- Interaction$new(predictor, feature = "Alternaria")
plot(interaction_Alternaria)
interaction_GIBSIM <- Interaction$new(predictor, feature = "GIBSIM_pred")
plot(interaction_GIBSIM)


total_interaction_graph <- plot(interaction) +
  labs(title = "Total interaction") +
  theme_classic() +
  theme(
    axis.text = element_text(size=20),
    axis.title.x = element_text(size = 20,hjust = 0.5, face='bold'),
    axis.title.y = element_text(size = 20,hjust = 0.5, face='bold'),
    plot.title = element_text(size = 20, face='bold')
  )

Hannaella_interaction_graph <- plot(interaction_Hannaella) +
  labs(title = "Interaction with Hannaella") +
  theme_classic() +
  theme(
    axis.text = element_text(size=20),
    axis.title.x = element_text(size = 20,hjust = 0.5, face='bold'),
    axis.title.y = element_text(size = 20,hjust = 0.5, face='bold'),
    plot.title = element_text(size = 20, face='bold')
  )

Alternaria_interaction_graph <- plot(interaction_Alternaria) +
  labs(title = "Interaction with Alternaria") +
  theme_classic() +
  theme(
    axis.text = element_text(size=20),
    axis.title.x = element_text(size = 20,hjust = 0.5, face='bold'),
    axis.title.y = element_text(size = 20,hjust = 0.5, face='bold'),
    plot.title = element_text(size = 20, face='bold')
  )

GIBSIM_interaction_graph <- plot(interaction_GIBSIM) +
  labs(title = "Interaction with GIBSIM_pred") +
  theme_classic() +
  theme(
    axis.text = element_text(size=20),
    axis.title.x = element_text(size = 20,hjust = 0.5, face='bold'),
    axis.title.y = element_text(size = 20,hjust = 0.5, face='bold'),
    plot.title = element_text(size = 20, face='bold')
  )



ggsave(
  plot = total_interaction_graph,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/interaction_total_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 20,
  height = 20,
  units = c("cm")
)

ggsave(
  plot = Hannaella_interaction_graph,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/interaction_Hannaella_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 20,
  height = 20,
  units = c("cm")
)

ggsave(
  plot = Alternaria_interaction_graph,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/interaction_Alternaria_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 20,
  height = 20,
  units = c("cm")
)
ggsave(
  plot = GIBSIM_interaction_graph,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/interaction_GIBSIM_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 20,
  height = 20,
  units = c("cm")
)
#####################################################################################################################################
#####################################################################################################################################
### 2D PDP
pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c("Hannaella", "GIBSIM_pred"), train = using_data, chull = TRUE)
plotPartial(pdp_2d, levelplot = FALSE, contour = TRUE)
# 
# pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0026", "OTU0006"), train = using_data, chull = TRUE)
# plotPartial(pdp_2d, levelplot = FALSE, contour = TRUE)
# 
# pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0004", "OTU0038"), train = using_data, chull = TRUE)
# pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0004", "OTU0056"), train = using_data, chull = TRUE)
# pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0001", "OTU0075"), train = using_data, chull = TRUE)

# # Heatmap ????????? 2D PDP
# plotPartial(
#   pdp_2d,
#   levelplot   = TRUE, 
#   # contour     = TRUE,
#   # contour.color = "#a0522d",
#   col.regions = viridis::rocket(100, direction = -1) #colorRampPalette(c("yellow","red"))(100)
# )


# for loop
imp_variables_Top <- names(using_data)[-which(names(using_data) == "FHB_incidence")]
imp_variables_Top <- imp_variables_Top[-which(imp_variables_Top == "GIBSIM_pred")]

pdp_list <- lapply(imp_variables_Top, function(var) {
  pdp::partial(rf_model_for_analysis, pred.var = c(var, "GIBSIM_pred"),
               train = using_data, chull = TRUE)
})

all_values <- unlist(lapply(pdp_list, function(df) df$yhat))
global_min <- min(all_values, na.rm = TRUE)
global_max <- max(all_values, na.rm = TRUE)


imsi_graph_list_2d_pdp <- list()
for(p in c(1:length(imp_variables_Top))){
  #p=1
  imp_var <- imp_variables_Top[p]
  pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c(imp_var, "GIBSIM_pred"), train = using_data, chull = TRUE)
  
  title_lab <- paste0("")
  # xlab <- paste0("\n",imp_var,"\n") 
  
  imsi_2d_pdp_graph <- plotPartial(
    pdp_2d,
    levelplot   = TRUE, 
    # contour     = TRUE,
    # contour.color = "#a0522d",
    col.regions = viridis::rocket(100, direction = -1),
    # at          = seq(floor(global_min/10)*10 , ceiling(global_max / 10) * 10, length.out = 100),
    at          = seq(global_min, global_max, length.out = 100),
    xlab        = list(cex=1.2, font=2),
    ylab        = list(label="GIBSIM", cex=1.2, font=2),
    main        = list(label=title_lab, cex=1.3, font=2),
    par.settings = list(axis.text = list(cex=1.1, col="black"))
  )
  
  imsi_graph_list_2d_pdp[[p]] <- imsi_2d_pdp_graph
  print(p)
}
combined_plot_2d_pdp <- do.call(grid.arrange, 
                                c(imsi_graph_list_2d_pdp, ncol = 3))
ggsave(
  plot = combined_plot_2d_pdp,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/2D_PDP_with_GIBSIM_pred_plot_fixed_legend_", 
    version_pred_obs,
    ".png"
  ),
  width = 30,
  height = 20,
  units = c("cm")
)

# ------
imsi_graph_list_2d_pdp_2 <- list()
for(p in c(1:length(imp_variables_Top))){
  #p=1
  imp_var <- imp_variables_Top[p]
  pdp_2d <- pdp::partial(rf_model_for_analysis, pred.var = c(imp_var, "GIBSIM_pred"), train = using_data, chull = TRUE)
  
  title_lab <- paste0("")
  # xlab <- paste0("\n",imp_var,"\n") 
  
  imsi_2d_pdp_graph <- plotPartial(
    pdp_2d,
    levelplot   = TRUE, 
    # contour     = TRUE,
    # contour.color = "#a0522d",
    col.regions = viridis::rocket(100, direction = -1),
    # at          = seq(floor(global_min/10)*10 , ceiling(global_max / 10) * 10, length.out = 100),
    xlab        = list(cex=1.2, font=2),
    ylab        = list(label="GIBSIM", cex=1.2, font=2),
    main        = list(label=title_lab, cex=1.3, font=2),
    par.settings = list(axis.text = list(cex=1.1, col="black"))
  )
  
  imsi_graph_list_2d_pdp_2[[p]] <- imsi_2d_pdp_graph
  print(p)
}
combined_plot_2d_pdp_not_fixed_legend <- do.call(grid.arrange, 
                                                 c(imsi_graph_list_2d_pdp_2, ncol = 3))
ggsave(
  plot = combined_plot_2d_pdp_not_fixed_legend,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/2D_PDP_with_GIBSIM_pred_plot_not_fixed_legend", 
    version_pred_obs,
    ".png"
  ),
  width = 30,
  height = 20,
  units = c("cm")
)

# =================================================================================================================================
# combined_plot_2d_pdp <- grid.arrange(
#   imsi_graph_list_2d_pdp[[1]], imsi_graph_list_2d_pdp[[2]], imsi_graph_list_2d_pdp[[3]], 
#   imsi_graph_list_2d_pdp[[4]], imsi_graph_list_2d_pdp[[5]], imsi_graph_list_2d_pdp[[6]],
#   imsi_graph_list_2d_pdp[[7]], imsi_graph_list_2d_pdp[[8]], imsi_graph_list_2d_pdp[[9]],
#   imsi_graph_list_2d_pdp[[10]],
#   ncol = 5
# )


# ### local 2d pdp (specific OTU - specific OTU )
# # OTU0004_OTU0056
# pdp_2d_OTU0004_OTU0056 <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0004", "OTU0056"), train = using_data, chull = TRUE)
# 
# # Heatmap ????????? 2D PDP
# plotPartial(
#   pdp_2d_OTU0004_OTU0056,
#   levelplot   = TRUE, 
#   # contour     = TRUE,
#   # contour.color = "#a0522d",
#   col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#   xlab        = list(cex=1.2, font=2),
#   ylab        = list(cex=1.2, font=2),
#   main        = list(label=title_lab, cex=1.3, font=2),
#   par.settings = list(axis.text = list(cex=1.1, col="black"))
# )
# 
# 
# #OTU0004_OTU0038
# pdp_2d_OTU0004_OTU0038 <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0004", "OTU0038"), train = using_data, chull = TRUE)
# 
# # Heatmap ????????? 2D PDP
# plotPartial(
#   pdp_2d_OTU0004_OTU0038,
#   levelplot   = TRUE, 
#   # contour     = TRUE,
#   # contour.color = "#a0522d",
#   col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#   xlab        = list(cex=1.2, font=2),
#   ylab        = list(cex=1.2, font=2),
#   main        = list(label=title_lab, cex=1.3, font=2),
#   par.settings = list(axis.text = list(cex=1.1, col="black"))
# )
# 
# #OTU0006_OTU0026
# pdp_2d_OTU0006_OTU0026 <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0006", "OTU0026"), train = using_data, chull = TRUE)
# 
# # Heatmap ????????? 2D PDP
# plotPartial(
#   pdp_2d_OTU0006_OTU0026,
#   levelplot   = TRUE, 
#   # contour     = TRUE,
#   # contour.color = "#a0522d",
#   col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#   xlab        = list(cex=1.2, font=2),
#   ylab        = list(cex=1.2, font=2),
#   main        = list(label=title_lab, cex=1.3, font=2),
#   par.settings = list(axis.text = list(cex=1.1, col="black"))
# )
# 
# #OTU0004_OTU0026
# pdp_2d_OTU0004_OTU0006 <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0004", "OTU0006"), train = using_data, chull = TRUE)
# 
# # Heatmap ????????? 2D PDP
# plotPartial(
#   pdp_2d_OTU0004_OTU0006,
#   levelplot   = TRUE, 
#   # contour     = TRUE,
#   # contour.color = "#a0522d",
#   col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#   xlab        = list(cex=1.2, font=2),
#   ylab        = list(cex=1.2, font=2),
#   main        = list(label=title_lab, cex=1.3, font=2),
#   par.settings = list(axis.text = list(cex=1.1, col="black"))
# )

### ----------------------------------------------------------------------------------------------------------------------------------------
interaction_vec <- names(data_only_MB_for_model)[grep("_", names(data_only_MB_for_model))]
interaction_vec <- interaction_vec[-which(interaction_vec == "FHB_incidence")]
for_loop_df_for_2d_PDP <- data.frame(interaction = interaction_vec) %>%
  tidyr::separate(interaction, into = c("interaction_v1", "interaction_v2"), sep = "_")

# interaction_v1 <- interaction_df$V1 #c("OTU0001", "OTU0001", "OTU0004", "OTU0004", "OTU0004", "OTU0004", "OTU0006", "OTU0026", "OTU0027", "OTU0039")
# interaction_v2 <- interaction_df$V2 #c("OTU0004", "OTU0039", "OTU0006", "OTU0026", "OTU0038", "OTU0056", "OTU0026", "OTU0038", "OTU0075", "OTU0056")
# for_loop_df_for_2d_PDP <- data.frame(interaction_v1, interaction_v2)

all_vals <- c()
for (x in 1:nrow(for_loop_df_for_2d_PDP)) {
  imsi_pdp <- pdp::partial(rf_model_for_analysis, pred.var = c(for_loop_df_for_2d_PDP$interaction_v1[x], for_loop_df_for_2d_PDP$interaction_v2[x]), 
                           train = using_data, chull = TRUE)
  
  all_vals <- c(all_vals,imsi_pdp$yhat)
  print(x)
}

# (2) ?????? global_min / global_max 계산 + floor/ceiling
MB_MB_global_min <- floor(min(all_vals))
MB_MB_global_max <- ceiling(max(all_vals))

imsi_2d_pdp_combind_not_fixed <- list()
imsi_2d_pdp_combind_fixed <- list()

for(x in c(1:nrow(for_loop_df_for_2d_PDP))){
  imsi_pdp <- pdp::partial(rf_model_for_analysis, pred.var = c(for_loop_df_for_2d_PDP$interaction_v1[x], for_loop_df_for_2d_PDP$interaction_v2[x]), 
                           train = using_data, chull = TRUE)
  
  imsi_2d_pdp <- plotPartial(
    imsi_pdp,
    levelplot   = TRUE, 
    # contour     = TRUE,
    # contour.color = "#a0522d",
    col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
    xlab        = list(cex=1.2, font=2),
    ylab        = list(cex=1.2, font=2),
    main        = list(label=title_lab, cex=1.3, font=2),
    par.settings = list(axis.text = list(cex=1.1, col="black"))
  )
  
  imsi_2d_pdp_combind_not_fixed[[x]] <- imsi_2d_pdp
  
  imsi_2d_pdp_fixed <- plotPartial(
    imsi_pdp,
    levelplot   = TRUE, 
    # contour     = TRUE,
    # contour.color = "#a0522d",
    col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
    at          = seq(MB_MB_global_min , MB_MB_global_max, length.out = 100),
    xlab        = list(cex=1.2, font=2),
    ylab        = list(cex=1.2, font=2),
    main        = list(label=title_lab, cex=1.3, font=2),
    par.settings = list(axis.text = list(cex=1.1, col="black"))
  )
  
  imsi_2d_pdp_combind_fixed[[x]] <- imsi_2d_pdp_fixed
  
  
  print(x)
}

combined_plot_MB_MB_2d_pdp_not_fixed <- grid.arrange(
  imsi_2d_pdp_combind_not_fixed[[1]], imsi_2d_pdp_combind_not_fixed[[2]], 
  imsi_2d_pdp_combind_not_fixed[[3]], imsi_2d_pdp_combind_not_fixed[[4]],
  imsi_2d_pdp_combind_not_fixed[[5]], #imsi_2d_pdp_combind_not_fixed[[6]],
  # imsi_2d_pdp_combind_not_fixed[[7]], imsi_2d_pdp_combind_not_fixed[[8]],
  # imsi_2d_pdp_combind_not_fixed[[9]], imsi_2d_pdp_combind_not_fixed[[10]],
  # imsi_2d_pdp_combind_not_fixed[[11]], imsi_2d_pdp_combind_not_fixed[[12]],
  # imsi_2d_pdp_combind_not_fixed[[13]], imsi_2d_pdp_combind_not_fixed[[14]],
  # imsi_2d_pdp_combind_not_fixed[[15]], imsi_2d_pdp_combind_not_fixed[[16]],
  ncol = 3
)

combined_plot_MB_MB_2d_pdp_fixed <- grid.arrange(
  imsi_2d_pdp_combind_fixed[[1]], imsi_2d_pdp_combind_fixed[[2]], 
  imsi_2d_pdp_combind_fixed[[3]],
  imsi_2d_pdp_combind_fixed[[4]],
  imsi_2d_pdp_combind_fixed[[5]], #imsi_2d_pdp_combind_fixed[[6]],
  # imsi_2d_pdp_combind_fixed[[7]], imsi_2d_pdp_combind_fixed[[8]],
  # imsi_2d_pdp_combind_fixed[[9]], imsi_2d_pdp_combind_fixed[[10]],
  # imsi_2d_pdp_combind_fixed[[11]], imsi_2d_pdp_combind_fixed[[12]],
  # imsi_2d_pdp_combind_fixed[[13]], imsi_2d_pdp_combind_fixed[[14]],
  # imsi_2d_pdp_combind_fixed[[15]], imsi_2d_pdp_combind_fixed[[16]],
  ncol = 3
)

ggsave(
  plot = combined_plot_MB_MB_2d_pdp_not_fixed,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/2D_PDP_MB_MB_plot_not_fixed_legend", 
    version_pred_obs,
    ".png"
  ),
  width = 30,
  height = 20,
  units = c("cm")
)


ggsave(
  plot = combined_plot_MB_MB_2d_pdp_fixed,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/2D_PDP_MB_MB_plot_fixed_legend", 
    version_pred_obs,
    ".png"
  ),
  width = 30,
  height = 20,
  units = c("cm")
)


# 
# pdp_2d_OTU0004_OTU0001 <- pdp::partial(rf_model_for_analysis, pred.var = c("OTU0004", "OTU0001"), train = using_data, chull = TRUE)
# 
# # Heatmap ????????? 2D PDP
# plotPartial(
#   pdp_2d_OTU0004_OTU0001,
#   levelplot   = TRUE, 
#   # contour     = TRUE,
#   # contour.color = "#a0522d",
#   col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#   xlab        = list(cex=1.2, font=2),
#   ylab        = list(cex=1.2, font=2),
#   main        = list(label=title_lab, cex=1.3, font=2),
#   par.settings = list(axis.text = list(cex=1.1, col="black"))
# )

# 
# 
# # # LOOCV setting
# # n <- nrow(selected_data_OTU_level) #selected_data_OTU_level, selected_data_for_selected_MB_genus
# # preds <- rep(NA, n)  # ???측값 ??????
# # actuals <- selected_data_OTU_level$FHB_incidence  #obs value
# # for (i in 1:n) {
# #   train_data <- selected_data_OTU_level[-i, ]  # i번째 ?????? ??????
# #   test_data <- selected_data_OTU_level[i, , drop = FALSE]  # i번째 ??????
# #   
# #   # ?????? ???????????? 모델 ??????
# #   rf_model <- randomForest(FHB_incidence ~ ., data = train_data, ntree = 100)
# #   
# #   # ????????? ????????? ??????
# #   preds[i] <- predict(rf_model, test_data)
# # }
# # 
# # # ????????? 계산
# # LOOCV_result <- data.frame(
# #   "pred" <- preds,
# #   "obs" <- actuals
# # )
# # names(LOOCV_result) <- c("pred", "obs")
# # LOOCV_lm_result <- lm(obs ~ pred, data = LOOCV_result)
# # summary(LOOCV_lm_result)
# 
# ###
# set.seed(1)
# # random_numbers <- sample(1:nrow(selected_data_OTU_level), 20)
# random_numbers <- sample(1:nrow(selected_data_OTU_level), nrow(selected_data_OTU_level))
# RF_tr <- selected_data_OTU_level[random_numbers,]
# RF_te <- selected_data_OTU_level[-random_numbers, ]
# nrow(RF_tr)
# nrow(RF_te)
# #
# set.seed(5)
# ctrl <- trainControl(method = "cv", number = 5)
# rf_model_v2 <- randomForest(FHB_incidence ~ ., data = RF_tr, ntree = 100, mtry = 4, cv.folds = 3)
# varImpPlot(rf_model_v2)
# # train 
# predict_tr <- as.data.frame(predict(rf_model_v2, newdata = RF_tr))
# predict_tr$obs <- RF_tr$FHB_incidence
# colnames(predict_tr) <- c("pred", "obs")
# plot(predict_tr$pred, predict_tr$obs)
# summary(lm(pred ~ obs, data = predict_tr))
# 
# rsquare(predict_tr$pred, predict_tr$obs) # 0.72
# RMSE(predict_tr$pred, predict_tr$obs) # 4.658
# 
# # test 
# predict_te <- as.data.frame(predict(rf_model, newdata = RF_te))
# predict_te$obs <- RF_te$FHB_incidence
# colnames(predict_te) <- c("pred", "obs")
# plot(predict_te$obs, predict_te$pred)
# summary(lm(pred ~ obs, data = predict_te))
# rsquare(predict_te$pred, predict_te$obs) #0.50
# RMSE(predict_te$pred, predict_te$obs) # 8.037
# 
# # for loop----
# for_loop_results_RF_R2 <- data.frame()
# for(x in c(1:100)){
#   set.seed(x)
#   random_numbers <- sample(1:nrow(selected_data_OTU_level), round(nrow(selected_data_OTU_level)/10*3))
#   RF_tr <- selected_data_OTU_level[-random_numbers,]
#   RF_te <- selected_data_OTU_level[random_numbers, ]
#   nrow(RF_tr)
#   nrow(RF_te)
#   #
#   ctrl <- trainControl(method = "cv", number = 3)
#   rf_model_v2 <- randomForest(FHB_incidence ~ ., data = RF_tr, ntree = 100, mtry = 4, cv.folds = 3)
#   varImpPlot(rf_model_v2)
#   # train 
#   predict_tr <- as.data.frame(predict(rf_model_v2, newdata = RF_tr))
#   predict_tr$obs <- RF_tr$FHB_incidence
#   colnames(predict_tr) <- c("pred", "obs")
#   plot(predict_tr$pred, predict_tr$obs)
#   summary(lm(pred ~ obs, data = predict_tr))
#   
#   rsquare(predict_tr$pred, predict_tr$obs) # 0.72
#   RMSE(predict_tr$pred, predict_tr$obs) # 4.658
#   
#   # test 
#   predict_te <- as.data.frame(predict(rf_model, newdata = RF_te))
#   predict_te$obs <- RF_te$FHB_incidence
#   colnames(predict_te) <- c("pred", "obs")
#   plot(predict_te$obs, predict_te$pred)
#   summary(lm(pred ~ obs, data = predict_te))
#   rsquare(predict_te$pred, predict_te$obs) #0.50
#   RMSE(predict_te$pred, predict_te$obs) # 8.037
#   
#   #results
#   for_loop_results_RF_R2[x,1] <- rsquare(predict_tr$pred, predict_tr$obs)
#   for_loop_results_RF_R2[x,2] <- rsquare(predict_te$pred, predict_te$obs) #0.50
# }
# 
# names(for_loop_results_RF_R2) <- c("RF_train_R2", "RF_test_R2")
# hist(for_loop_results_RF_R2$RF_train_R2)
# hist(for_loop_results_RF_R2$RF_test_R2)
### beta model (GAM) --------------------------------------------------------------------------------------------------------------
library(mgcv)
data_gibsim_MB_for_beta <- data_gibsim_MB_for_model

if(max(data_gibsim_MB_for_beta$FHB_incidence) > 1.1){
  data_gibsim_MB_for_beta$FHB_incidence <- data_gibsim_MB_for_beta$FHB_incidence/100
  print("### change ###")
}else{
  print("### not change ###")
}

if(max(data_gibsim_MB_for_beta$fhb_incidence_predicted) > 1.1){
  data_gibsim_MB_for_beta$fhb_incidence_predicted <- data_gibsim_MB_for_beta$fhb_incidence_predicted/100
  print("### change ###")
}else{
  print("### not change ###")
}

data_gibsim_MB_for_beta$fhb_incidence_predicted <- ifelse(data_gibsim_MB_for_beta$fhb_incidence_predicted < 0, 0, data_gibsim_MB_for_beta$fhb_incidence_predicted) 
names(data_gibsim_MB_for_beta)
gam_model <- gam(
  FHB_incidence ~ 
    s(Alternaria, k=3) + 
    s(Cladosporium, k=3) + 
    s(Epicoccum, k=3) + 
    s(Hannaella, k=3) + 
    s(Papiliotrema, k=3) + 
    s(Periconia, k=3) + 
    # s(Vishniacozyma, k=3) +
    s(fhb_incidence_predicted, k=3) +
    ti(Alternaria, Cladosporium, k=c(3,3)) +
    ti(Alternaria, Papiliotrema, k=c(3,3)) +
    ti(Alternaria, Periconia, k=c(3,3)) +
    ti(Cladosporium, Epicoccum, k=c(3,3)) +
    ti(Epicoccum, Periconia, k=c(3,3)),
  family = betar(link="logit"),
  data = data_gibsim_MB_for_beta, 
  select=TRUE
)

summary_gam <- summary(gam_model)
summary_gam_res <- as.data.frame(summary_gam$s.table)
summary_gam_res$X <- rownames(summary_gam_res)
summary_gam_res$star <- ifelse(summary_gam_res$`p-value` <= 0.01, "**", 
                               ifelse(summary_gam_res$`p-value` <= 0.05, "*", 
                                      ifelse(summary_gam_res$`p-value` <= 0.1, ".", "")))
summary_gam_res
writexl::write_xlsx(summary_gam_res, path = paste0(
  "./Output/3. obs_pred_graph/2025/",
  version_pred_obs,
  "/GAM_model_coef_", 
  version_pred_obs,
  ".xlsx"
))
# plot(gam_model, pages = 5, shade = TRUE)
GAM_GIBSIM_MB_pred_obs_plot <- obs_pred_plot_for_GAM(model = gam_model, 
                                                     data_contain_X_Y = data_gibsim_MB_for_beta, 
                                                     Y_colname = "FHB_incidence", 
                                                     Title = "GAM_model"
)
GAM_GIBSIM_MB_pred_obs_plot
ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/GIBSIM_MB_GAM_beta_",version_pred_obs, ".png")),
       plot = GAM_GIBSIM_MB_pred_obs_plot,
       width = 10, height = 10, bg = "white", dpi = 1500)

### GIBSIM vs. GAM vs. obs
GAM_predict_tr <- as.data.frame(predict(gam_model, type = "response"))*100 #as.data.frame(predict(model, newdata = data_contain_X_Y))
GAM_predict_tr$GIBSIM <- data_gibsim_MB_for_beta$fhb_incidence_predicted*100
GAM_predict_tr$obs <- data_gibsim_MB_for_beta$FHB_incidence*100
names(GAM_predict_tr) <- c("GAM", "GIBSIM", "obs")

GAM_predict_tr$obs_minus_GAM <- GAM_predict_tr$obs - GAM_predict_tr$GAM
GAM_predict_tr$obs_minus_GIBSIM <- GAM_predict_tr$obs - GAM_predict_tr$GIBSIM

GAM_predict_tr$abs_obs_minus_GAM <- abs(GAM_predict_tr$obs_minus_GAM)
GAM_predict_tr$abs_obs_minus_GIBSIM <- abs(GAM_predict_tr$obs_minus_GIBSIM)

GAM_predict_tr$which_is_good <- ifelse(GAM_predict_tr$abs_obs_minus_GAM < GAM_predict_tr$abs_obs_minus_GIBSIM, "GAM_is_good", "GIBSIM_is_good")
head(GAM_predict_tr)
GAM_predict_tr$sample_id <- rownames(GAM_predict_tr)

#Graph
GAM_predict_tr2 <- GAM_predict_tr %>%
  mutate(
    diff_abs = abs(obs_minus_GAM - obs_minus_GIBSIM)
  ) %>%
  dplyr::arrange(
    which_is_good,        # block 분리 (GIBSIM ??? GAM)
    desc(diff_abs)        # block ???부 ??????
  ) %>%
  mutate(
    sample_id_ord = factor(sample_id, levels = sample_id)
  )

plot_df <- GAM_predict_tr2 %>%
  pivot_longer(
    cols = c(obs_minus_GAM, obs_minus_GIBSIM),
    names_to = "model",
    values_to = "residual"
  )

table(plot_df$which_is_good)
plot_df_GAM <- plot_df %>% filter(which_is_good == "GAM_is_good")
plot_df_GIBSIM <- plot_df %>% filter(which_is_good == "GIBSIM_is_good")

library(ggplot2)

G_GAM_is_good <- ggplot(
  plot_df_GAM,
  aes(
    x = sample_id_ord,
    y = residual,
    color = model
  )
) +
  geom_line(
    aes(group = sample_id_ord),
    color = "grey60",
    linewidth = 0.8
  ) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.0)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_color_manual(
    values = c(
      obs_minus_GAM    = "orange",
      obs_minus_GIBSIM = "#708090"
    ),
    labels = c(
      obs_minus_GAM    = "Obs ??? GAM",
      obs_minus_GIBSIM = "Obs ??? GIBSIM"
    )
  ) +
  labs(
    title = "GAM is good",
    x = "Sample ID\n",
    y = "\nResidual (obs ??? prediction)",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x  = element_text(
      size = 15,
      # face = "bold",
      color = "black",
      angle = 0,
      hjust = 0.5,
      vjust = 0.4
    ),
    axis.text.y  = element_text(
      size = 10,
      # face = "bold",
      color = "black"
    ),
    panel.grid.major.x = element_blank()
  ) + 
  coord_flip()
G_GAM_is_good
G_GIBSIM_is_good <- ggplot(
  plot_df_GIBSIM,
  aes(x = sample_id_ord, y = residual,color = model)) +
  geom_line(
    aes(group = sample_id_ord),
    color = "grey60",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3
  ) +
  geom_hline(yintercept = 0,linetype = "dashed",color = "grey40") +
  scale_color_manual(
    values = c(obs_minus_GAM    = "darkorange", obs_minus_GIBSIM = "#708090"),
    labels = c(obs_minus_GAM    = "Obs ??? GAM", obs_minus_GIBSIM = "Obs ??? GIBSIM")
  ) +
  labs(
    title = "GIBSIM is good",
    x = "Sample ID\n",
    y = "\nResidual (obs ??? prediction)",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold"),
    axis.text.x  = element_text(
      size = 15,
      # face = "bold",
      color = "black",
      angle = 0,
      hjust = 0.5,
      vjust = 0.4
    ),
    axis.text.y  = element_text(
      size = 10,
      # face = "bold",
      color = "black"
    ),
    panel.grid.major.x = element_blank()
  ) + 
  coord_flip()

ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/G_GAM_is_good_",version_pred_obs, ".png")),
       plot = G_GAM_is_good,
       width = 7, height = 10, bg = "white", dpi = 1500)
ggsave(filename = file.path(paste0(save_dir, save_year, "/",version_pred_obs, "/G_GIBSIM_is_good_",version_pred_obs, ".png")),
       plot = G_GIBSIM_is_good,
       width = 7, height = 6, bg = "white", dpi = 1500)
### ICE
GAM_imp_variables_Top <- names(data_gibsim_MB_for_beta)[-which(names(data_gibsim_MB_for_beta) == "FHB_incidence")]
length(GAM_imp_variables_Top)
GAM_imsi_graph_list_ice <- list()
for(p in c(1:length(GAM_imp_variables_Top))){
  #p=1
  GAM_imp_var <- GAM_imp_variables_Top[p]
  GAM_imsi_ice <- pdp::partial(gam_model, pred.var = GAM_imp_var, train = data_gibsim_MB_for_beta,   
                               grid.resolution = 30,
                               pred.fun = function(object, newdata) {
                                 predict(object, newdata = newdata, type = "response")
                               },
                               ice = T)
  # plotPartial(ice, alpha = 0.1)
  
  title_lab <- paste0("ICE : ",GAM_imp_var,"\n") 
  xlab <- paste0("\n",GAM_imp_var,"\n") 
  
  if(GAM_imp_var != "fhb_incidence_predicted"){
    imsi_graph <- autoplot(GAM_imsi_ice, center = TRUE, alpha = 0.2, rug = F, train = data_gibsim_MB_for_beta) + 
      # geom_line(color = "#303030", linewidth = 1) +
      labs(
        title = title_lab,
        x = xlab,
        y = "\nFHB incidence (%)\n" #ice of incidence # y hat
      ) +
      scale_y_continuous(
        labels = function(x) x * 100
      ) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 12, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 14, face = 'bold')
      )
  }else if(GAM_imp_var == "fhb_incidence_predicted"){
    imsi_graph <- autoplot(GAM_imsi_ice, center = TRUE, alpha = 0.2, rug = F, train = data_gibsim_MB_for_beta) + 
      # geom_line(color = "#303030", linewidth = 1) +
      labs(
        title = title_lab,
        x = xlab,
        y = "\nFHB incidence (%)\n" #ice of incidence # y hat
      ) +
      scale_y_continuous(
        labels = function(x) x * 100
      ) +
      scale_x_continuous(
        labels = function(x) x * 100
      ) + 
      theme_classic() +
      theme(
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 12, hjust = 0.5, face = 'bold'),
        axis.title.y = element_text(size = 12, hjust = 0.5, face = 'bold'),
        plot.title = element_text(size = 14, face = 'bold')
      )
  }
  
  
  
  GAM_imsi_graph_list_ice[[p]] <- imsi_graph
}
GAM_combined_plot_ice <- grid.arrange(GAM_imsi_graph_list_ice[[1]], GAM_imsi_graph_list_ice[[2]], GAM_imsi_graph_list_ice[[3]], 
                                      GAM_imsi_graph_list_ice[[4]], GAM_imsi_graph_list_ice[[5]], GAM_imsi_graph_list_ice[[6]],
                                      # GAM_imsi_graph_list_ice[[7]], GAM_imsi_graph_list_ice[[8]], GAM_imsi_graph_list_ice[[9]],
                                      # GAM_imsi_graph_list_ice[[10]], GAM_imsi_graph_list_ice[[11]], 
                                      GAM_imsi_graph_list_ice[[12]],
                                      # GAM_imsi_graph_list_ice[[13]], GAM_imsi_graph_list_ice[[14]], GAM_imsi_graph_list_ice[[15]],
                                      # GAM_imsi_graph_list_ice[[16]], GAM_imsi_graph_list_ice[[17]],GAM_imsi_graph_list_ice[[18]],
                                      ncol=4)


ggsave(
  plot = GAM_combined_plot_ice,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/combined_GAM_ICE_plot_", 
    version_pred_obs,
    ".png"
  ),
  width = 28,
  height = 13,
  units = c("cm")
)


### why GAM is better than GIBSIM
df_for_more_analysis <- cbind(df_rslt_stat, GAM_predict_tr)
names(df_for_more_analysis)
df_for_more_analysis$abs_obs_minus_GIBSIM
plot(df_for_more_analysis$obs,df_for_more_analysis$adj_FHB_incidence)
hist(df_for_more_analysis$abs_obs_minus_GIBSIM)
hist(df_for_more_analysis$abs_obs_minus_GAM)

df_for_more_analysis <- df_for_more_analysis %>%
  mutate(
    cluster = case_when(
      abs(abs_obs_minus_GAM - abs_obs_minus_GIBSIM) > 10 &
        which_is_good == "GAM_is_good" ~ "GAM_good",
      
      abs(abs_obs_minus_GAM - abs_obs_minus_GIBSIM) < 10 ~ "no_diff",
      
      abs(abs_obs_minus_GAM - abs_obs_minus_GIBSIM) > 10 &
        which_is_good == "GIBSIM_is_good" ~ "GIBSIM_good",
      TRUE ~ NA_character_
    )
  ) %>% 
  mutate(
    Check_GIBSIM_res = case_when(
      
      abs_obs_minus_GIBSIM < 20 ~ "GIBSIM_good",
      abs_obs_minus_GIBSIM < 110 ~ "GIBSIM_bad",
      abs_obs_minus_GIBSIM < 110 ~ "GIBSIM_very_bad",
      
      TRUE ~ NA_character_
    )
  )%>% 
  mutate(
    Check_GAM_res = case_when(
      
      abs_obs_minus_GAM < 20 ~ "GAM_good",
      abs_obs_minus_GAM < 100 ~ "GAM_bad",
      abs_obs_minus_GAM < 110 ~ "GAM_very_bad",
      
      TRUE ~ NA_character_
    )
  )
df_for_more_analysis$cluster <- as.factor(df_for_more_analysis$cluster)
df_for_more_analysis$Check_GIBSIM_res <- factor(df_for_more_analysis$Check_GIBSIM_res, level = c("GIBSIM_good", "GIBSIM_bad", "GIBSIM_very_bad"))
df_for_more_analysis$Check_GAM_res <- factor(df_for_more_analysis$Check_GAM_res, level = c("GAM_good", "GAM_bad", "GAM_very_bad"))
table(df_for_more_analysis$Check_GIBSIM_res)
table(df_for_more_analysis$Check_GAM_res)
names(df_for_more_analysis)

writexl::write_xlsx(df_for_more_analysis, 
                    paste0(
                      "./Output/3. obs_pred_graph/2025/",
                      version_pred_obs,
                      "/df_for_more_analysis_good_20_bad_100",
                      ".xlsx"
                    ))


GIBSIM_obs_pred_residual_g <- ggplot(data = df_for_more_analysis, mapping = aes(x = GIBSIM, y = obs, color = Check_GIBSIM_res)) + 
  geom_point(size = 3, alpha = 0.5)+ 
  scale_color_manual(values = c("GIBSIM_good" = "forestgreen",
                                "GIBSIM_bad" = "orange",
                                "GIBSIM_very_bad" = "red")) + 
  labs(title = "", 
       y = "FHB incidence\n", 
       x = "\nGIBSIM predicted") + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  )+ 
  geom_hline(yintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)+
  geom_vline(xintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)  

GAM_obs_pred_residual_g <- ggplot(data = df_for_more_analysis, mapping = aes(x = GAM, y = obs, color = Check_GAM_res)) + 
  geom_point(size = 3, alpha = 0.5)+ 
  # # geom_abline(slope = 1, intercept = 0,  linetype = "dashed",linewidth = 1,
  # #             color = "forestgreen", alpha = 0.2) +
  # geom_abline(slope = 1, intercept = 10, linetype = "dashed",linewidth = 1,
  #             color = "orange", alpha = 0.2) +
  # geom_abline(slope = 1, intercept = -10, linetype = "dashed",linewidth = 1,
  #             color = "orange", alpha = 0.2) +
  # geom_abline(slope = 1, intercept = 20, linetype = "dashed",linewidth = 1,
  #             color = "red", alpha = 0.2) +
  # geom_abline(slope = 1, intercept = -20, linetype = "dashed",linewidth = 1,
  #             color = "red", alpha = 0.2) +
  scale_color_manual(values = c("GAM_good" = "forestgreen",
                                "GAM_bad" = "orange",
                                "GAM_very_bad" = "red")) + 
  labs(title = "", 
       y = "FHB incidence\n", 
       x = "\nGAM predicted") + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  )+ 
  geom_hline(yintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)+
  geom_vline(xintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)  

GIBSIM_obs_pred_GAM_label_residual_g <- ggplot(data = df_for_more_analysis, mapping = aes(x = GIBSIM, y = obs, color = Check_GAM_res)) + 
  geom_point(size = 3, alpha = 0.5)+ 
  scale_color_manual(values = c("GAM_good" = "forestgreen",
                                "GAM_bad" = "orange",
                                "GAM_very_bad" = "red")) + 
  labs(title = "", 
       y = "FHB incidence\n", 
       x = "\nGIBSIM predicted") + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  )+ 
  geom_hline(yintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)+
  geom_vline(xintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)  

GAM_obs_pred_GIBSIM_label_residual_g <- ggplot(data = df_for_more_analysis, mapping = aes(x = GAM, y = obs, color = Check_GIBSIM_res)) + 
  geom_point(size = 3, alpha = 0.5)+ 
  scale_color_manual(values = c("GIBSIM_good" = "forestgreen",
                                "GIBSIM_bad" = "orange",
                                "GIBSIM_very_bad" = "red")) + 
  labs(title = "", 
       y = "FHB incidence\n", 
       x = "\nGAM predicted") + 
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  ) + 
  geom_hline(yintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)+
  geom_vline(xintercept=20, linetype='dashed', color='tomato', linewidth = 1, alpha = 0.8)  


GIBSIM_obs_pred_residual_g
GAM_obs_pred_residual_g
GIBSIM_obs_pred_GAM_label_residual_g
GAM_obs_pred_GIBSIM_label_residual_g


ggsave(
  plot = GIBSIM_obs_pred_residual_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/GIBSIM_obs_pred_residual_g_",
    version_pred_obs,
    ".png"
  ),
  width = 14,
  height = 11,
  units = c("cm")
)

ggsave(
  plot = GAM_obs_pred_residual_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/GAM_obs_pred_residual_g_",
    version_pred_obs,"(bad_20, no_verybad)",
    ".png"
  ),
  width = 14,
  height = 11,
  units = c("cm")
)

ggsave(
  plot = GIBSIM_obs_pred_GAM_label_residual_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/GIBSIM_obs_pred_GAM_label_residual_g_",
    version_pred_obs,"(bad_20, no_verybad)",
    ".png"
  ),
  width = 14,
  height = 11,
  units = c("cm")
)

ggsave(
  plot = GAM_obs_pred_GIBSIM_label_residual_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/GAM_obs_pred_GIBSIM_label_residual_g_",
    version_pred_obs,"(bad_20, no_verybad)",
    ".png"
  ),
  width = 14,
  height = 11,
  units = c("cm")
)

GIBSIM_very_bad_df <- df_for_more_analysis %>% filter(Check_GIBSIM_res == "GIBSIM_very_bad")
GIBSIM_good_df <- df_for_more_analysis %>% filter(Check_GIBSIM_res == "GIBSIM_good")

names(GIBSIM_very_bad_df)
table(df_for_more_analysis$year)
table(GIBSIM_very_bad_df$year)
GIBSIM_very_bad_df$under_over <- ifelse(GIBSIM_very_bad_df$obs_minus_GIBSIM < 0, "over", "under")
GIBSIM_good_df_lowinc <- GIBSIM_good_df %>% filter(obs < 20)
table(GIBSIM_very_bad_df$under_over)

GIBSIM_very_bad_df_over <- GIBSIM_very_bad_df %>% filter(under_over == "over")

wilcox.test(GIBSIM_very_bad_df_over$rhum_hd,GIBSIM_good_df_lowinc$rhum_hd)
wilcox.test(GIBSIM_very_bad_df_over$rhum_after_sampling_div2_1st,GIBSIM_good_df_lowinc$rhum_after_sampling_div2_1st)
wilcox.test(GIBSIM_very_bad_df_over$temp_before_hd,GIBSIM_good_df_lowinc$temp_before_hd)
boxplot(GIBSIM_very_bad_df_over$rhum_after_sampling_div2_1st,GIBSIM_good_df_lowinc$rhum_after_sampling_div2_1st)
boxplot(GIBSIM_very_bad_df_over$rhum_after_sampling,GIBSIM_good_df_lowinc$rhum_after_sampling)

boxplot(GIBSIM_very_bad_df_over$temp_before_hd,GIBSIM_good_df_lowinc$temp_before_hd)
boxplot(GIBSIM_very_bad_df_over$temp_hd,GIBSIM_good_df_lowinc$temp_hd)
boxplot(GIBSIM_very_bad_df_over$temp_flower,GIBSIM_good_df_lowinc$temp_flower)
boxplot(GIBSIM_very_bad_df_over$temp_after_sampling_div2_1st,GIBSIM_good_df_lowinc$temp_after_sampling_div2_1st)
boxplot(GIBSIM_very_bad_df_over$temp_after_sampling_div2_2nd,GIBSIM_good_df_lowinc$temp_after_sampling_div2_2nd)


GIBSIM_very_bad_df_under <- GIBSIM_very_bad_df %>% filter(under_over == "under")
GIBSIM_very_bad_df_under[,c(425, 424, 396:422)]
GIBSIM_very_bad_df_over[,c(425, 424, 396:422)]
names(GIBSIM_very_bad_df_under)








GIBSIM_good_df <- df_for_more_analysis %>% filter(cluster == "GIBSIM_good")
no_diff_df <- df_for_more_analysis %>% filter(cluster == "no_diff")
GAM_good_df <- df_for_more_analysis %>% filter(cluster == "GAM_good")

boxplot(GIBSIM_good_df$temp_before_hd, no_diff_df$temp_before_hd, GAM_good_df$temp_before_hd)
aov(GIBSIM_good_df$temp_after_sampling, GAM_normal_df$temp_after_sampling, GAM_good_df$temp_after_sampling)
boxplot(GIBSIM_good_df$temp_after_flower, GAM_normal_df$temp_after_flower, GAM_good_df$temp_after_flower)


boxplot(df_for_more_analysis$obs_minus_GAM, df_for_more_analysis$obs_minus_GIBSIM)

### additional analysis based on ICE results (boxplot)
hist(df_for_more_analysis$Alternaria)
more_analysis_ICE_res <- df_for_more_analysis %>% dplyr::select("Alternaria", "Cladosporium", "Epicoccum", "Hannaella", 
                                                                "Papiliotrema","Periconia", "incidence")
if(max(more_analysis_ICE_res$incidence) < 1.1 ){
  more_analysis_ICE_res$incidence <- more_analysis_ICE_res$incidence*100  
}

hist(more_analysis_ICE_res$Periconia)
more_analysis_ICE_res$Alternaria_cat <- ifelse(more_analysis_ICE_res$Alternaria < 7, "low_abu", 
                                               ifelse(more_analysis_ICE_res$Alternaria > 7, "high_abu", "mid_abu"))
more_analysis_ICE_res$Epicoccum_cat <- ifelse(more_analysis_ICE_res$Epicoccum < 2.5, "low_abu", 
                                              ifelse(more_analysis_ICE_res$Epicoccum > 7.5, "high_abu", "mid_abu"))

more_analysis_ICE_res$Hannaella_cat <- ifelse(more_analysis_ICE_res$Hannaella < 3, "low_abu", 
                                              ifelse(more_analysis_ICE_res$Hannaella > 5, "high_abu", "mid_abu"))

more_analysis_ICE_res$Periconia_cat <- ifelse(more_analysis_ICE_res$Periconia < 4, "low_abu", 
                                              ifelse(more_analysis_ICE_res$Periconia > 4, "high_abu", "mid_abu"))

more_analysis_ICE_res$Cladosporium_cat <- ifelse(more_analysis_ICE_res$Cladosporium < 9, "low_abu", 
                                                 ifelse(more_analysis_ICE_res$Cladosporium > 9, "high_abu", "mid_abu"))
more_analysis_ICE_res$Papiliotrema_cat <- ifelse(more_analysis_ICE_res$Papiliotrema < 7, "low_abu", 
                                                 ifelse(more_analysis_ICE_res$Papiliotrema > 7, "high_abu", "mid_abu"))

Alt_boxplot <- ggplot(more_analysis_ICE_res, aes(x = factor(Alternaria_cat), y = incidence)) + 
  geom_boxplot()

Epi_boxplot <- ggplot(more_analysis_ICE_res, aes(x = factor(Epicoccum_cat), y = incidence)) + 
  geom_boxplot()

Han_boxplot <- ggplot(more_analysis_ICE_res, aes(x = factor(Hannaella_cat), y = incidence)) + 
  geom_boxplot()

Per_boxplot <- ggplot(more_analysis_ICE_res, aes(x = factor(Periconia_cat), y = incidence)) + 
  geom_boxplot()

Cla_boxplot <- ggplot(more_analysis_ICE_res, aes(x = factor(Cladosporium_cat), y = incidence)) + 
  geom_boxplot()
Pap_boxplot <- ggplot(more_analysis_ICE_res, aes(x = factor(Papiliotrema_cat), y = incidence)) + 
  geom_boxplot()


### correlation
if(!dir.exists(paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve"))){
  dir.create(paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve"))
}
names(df_rslt_stat)
df_for_more_analysis3 <- df_rslt_stat %>% dplyr::select("Alternaria", "Cladosporium", "Epicoccum", "Hannaella", 
                                                        "Papiliotrema","Periconia", "incidence", 
                                                        "temp_before_hd", "rhum_before_hd", "prcp_before_hd", 
                                                        "temp_hd", "rhum_hd", "prcp_hd", 
                                                        "temp_before_hd_sampling_10days", "rhum_before_hd_sampling_10days", "prcp_before_hd_sampling_10days", 
                                                        "temp_hd_sampling_10days", "rhum_hd_sampling_10days", "prcp_hd_sampling_10days")
# df_for_more_analysis3[,c(1:6)] <- -df_for_more_analysis3[,c(1:6)]
gam_microbe_DVR <- function(response, data) {
  gam(
    as.formula(
      paste0(response, " ~ ",
             "s(temp_before_hd, k=5) + ",
             "s(rhum_before_hd, k=5) + ",
             "s(prcp_before_hd, k=5) + ",
             
             "s(temp_hd, k=5) + ",
             "s(rhum_hd, k=5) + ",
             "s(prcp_hd, k=5)")
    ),
    data = data,
    method = "REML",
    family = gaussian()
  )
}

gam_microbe_sampling <- function(response, data) {
  gam(
    as.formula(
      paste0(response, " ~ ",
             "s(temp_before_hd_sampling_10days, k=5) + ",
             "s(rhum_before_hd_sampling_10days, k=5) + ",
             "s(prcp_before_hd_sampling_10days, k=5) + ",
             
             "s(temp_hd_sampling_10days, k=5) + ",
             "s(rhum_hd_sampling_10days, k=5) + ",
             "s(prcp_hd_sampling_10days, k=5)")
    ),
    data = data,
    method = "REML",
    family = gaussian()
  )
}

head(as.data.frame(df_for_more_analysis3))
m_alt  <- gam_microbe_DVR("Alternaria", df_for_more_analysis3)
m_epi  <- gam_microbe_DVR("Epicoccum", df_for_more_analysis3)
m_han  <- gam_microbe_DVR("Hannaella", df_for_more_analysis3)
m_per  <- gam_microbe_DVR("Periconia", df_for_more_analysis3)

m_alt_sampling  <- gam_microbe_sampling("Alternaria", df_for_more_analysis3)
m_epi_sampling  <- gam_microbe_sampling("Epicoccum", df_for_more_analysis3)
m_han_sampling  <- gam_microbe_sampling("Hannaella", df_for_more_analysis3)
m_per_sampling  <- gam_microbe_sampling("Periconia", df_for_more_analysis3)

#
df_for_more_analysis3$incidence <- ifelse(df_for_more_analysis3$incidence == 0, 0.001, df_for_more_analysis3$incidence)
df_for_more_analysis3$incidence <- ifelse(df_for_more_analysis3$incidence == 1, 0.999, df_for_more_analysis3$incidence)

m_fhb <- gam(
  incidence ~ s(temp_before_hd, k=5) + s(rhum_before_hd, k=5) + s(prcp_before_hd, k=5) + s(temp_hd, k=5) + s(rhum_hd, k=5) + s(prcp_hd, k=5),
  data = df_for_more_analysis3,
  method = "REML",
  family = betar(link = "logit")   # 비율?????? beta regression 권장
)

m_fhb_sampling <- gam(
  incidence ~ s(temp_before_hd_sampling_10days, k=5) + s(rhum_before_hd_sampling_10days, k=5) + s(prcp_before_hd_sampling_10days, k=5) + 
    s(temp_hd_sampling_10days, k=5) + s(rhum_hd_sampling_10days, k=5) + s(prcp_hd_sampling_10days, k=5),
  data = df_for_more_analysis3,
  method = "REML",
  family = betar(link = "logit")   # 비율?????? beta regression 권장
)

#save
summary_m_alt <- as.data.frame(summary(m_alt)$s.table)
summary_m_alt$names <- rownames(summary_m_alt)
write_xlsx(summary_m_alt, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/DVR_GAM_X_is_wth_Y_is_", 
                                        "Alternaria_",version_pred_obs, ".xlsx"))
summary_m_epi <- as.data.frame(summary(m_epi)$s.table)
summary_m_epi$names <- rownames(summary_m_epi)
write_xlsx(summary_m_epi, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/DVR_GAM_X_is_wth_Y_is_", 
                                        "Epicoccum_",version_pred_obs, ".xlsx"))
summary_m_han <- as.data.frame(summary(m_han)$s.table)
summary_m_han$names <- rownames(summary_m_han)
write_xlsx(summary_m_han, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/DVR_GAM_X_is_wth_Y_is_", 
                                        "Hannaella_",version_pred_obs, ".xlsx"))
summary_m_per <- as.data.frame(summary(m_per)$s.table)
summary_m_per$names <- rownames(summary_m_per)
write_xlsx(summary_m_per, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/DVR_GAM_X_is_wth_Y_is_", 
                                        "Periconia_",version_pred_obs, ".xlsx"))
summary_m_fhb <- as.data.frame(summary(m_fhb)$s.table)
summary_m_fhb$names <- rownames(summary_m_fhb)
write_xlsx(summary_m_fhb, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/DVR_GAM_X_is_wth_Y_is_", 
                                        "FHB_incidence_",version_pred_obs, ".xlsx"))

summary_m_alt_sampling <- as.data.frame(summary(m_alt_sampling)$s.table)
summary_m_alt_sampling$names <- rownames(summary_m_alt_sampling)
write_xlsx(summary_m_alt_sampling, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/sampling_GAM_X_is_wth_Y_is_", 
                                        "Alternaria_",version_pred_obs, ".xlsx"))
summary_m_epi_sampling <- as.data.frame(summary(m_epi_sampling)$s.table)
summary_m_epi_sampling$names <- rownames(summary_m_epi_sampling)
write_xlsx(summary_m_epi_sampling, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/sampling_GAM_X_is_wth_Y_is_", 
                                        "Epicoccum_",version_pred_obs, ".xlsx"))
summary_m_han_sampling <- as.data.frame(summary(m_han_sampling)$s.table)
summary_m_han_sampling$names <- rownames(summary_m_han_sampling)
write_xlsx(summary_m_han_sampling, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/sampling_GAM_X_is_wth_Y_is_", 
                                        "Hannaella_",version_pred_obs, ".xlsx"))
summary_m_per_sampling <- as.data.frame(summary(m_per_sampling)$s.table)
summary_m_per_sampling$names <- rownames(summary_m_per_sampling)
write_xlsx(summary_m_per_sampling, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/sampling_GAM_X_is_wth_Y_is_", 
                                        "Periconia_",version_pred_obs, ".xlsx"))
summary_m_fhb_sampling <- as.data.frame(summary(m_fhb_sampling)$s.table)
summary_m_fhb_sampling$names <- rownames(summary_m_fhb_sampling)
write_xlsx(summary_m_fhb_sampling, path = paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve/sampling_GAM_X_is_wth_Y_is_", 
                                        "FHB_incidence_",version_pred_obs, ".xlsx"))




extract_smooth_var_DVR <- function(model, varname, data) {
  # ???균값?????? 고정??? 변?????? ??????
  newdata <- data.frame(
    temp_before_hd  = mean(data$temp_before_hd),
    rhum_before_hd  = mean(data$rhum_before_hd),
    prcp_before_hd  = mean(data$prcp_before_hd),
    temp_hd         = mean(data$temp_hd),
    rhum_hd         = mean(data$rhum_hd),
    prcp_hd         = mean(data$prcp_hd)
  )
  
  # ?????? 변?????? sequence??? 변???
  newdata <- newdata[rep(1, 100), ]
  newdata[[varname]] <- seq(
    min(data[[varname]]),
    max(data[[varname]]),
    length = 100
  )
  
  pred <- predict(model, newdata, type = "terms", se.fit = TRUE)
  
  smooth_index <- grep(paste0("s\\(", varname), colnames(pred$fit))
  
  data.frame(
    x = newdata[[varname]],
    est = pred$fit[, smooth_index],
    se  = pred$se.fit[, smooth_index]
  )
}

extract_smooth_var_sampling <- function(model, varname, data) {
  
  # ???균값?????? 고정??? 변?????? ??????
  newdata <- data.frame(
    temp_before_hd_sampling_10days  = mean(data$temp_before_hd_sampling_10days),
    rhum_before_hd_sampling_10days  = mean(data$rhum_before_hd_sampling_10days),
    prcp_before_hd_sampling_10days  = mean(data$prcp_before_hd_sampling_10days),
    temp_hd_sampling_10days         = mean(data$temp_hd_sampling_10days),
    rhum_hd_sampling_10days         = mean(data$rhum_hd_sampling_10days),
    prcp_hd_sampling_10days         = mean(data$prcp_hd_sampling_10days)
  )
  
  # ?????? 변?????? sequence??? 변???
  newdata <- newdata[rep(1, 100), ]
  newdata[[varname]] <- seq(
    min(data[[varname]]),
    max(data[[varname]]),
    length = 100
  )
  
  pred <- predict(model, newdata, type = "terms", se.fit = TRUE)
  
  smooth_index <- grep(paste0("s\\(", varname), colnames(pred$fit))
  
  data.frame(
    x = newdata[[varname]],
    est = pred$fit[, smooth_index],
    se  = pred$se.fit[, smooth_index]
  )
}

# for loop
print("for loop")
mb_v <- c("Alternaria", "Cladosporium", "Epicoccum", "Hannaella", "Papiliotrema","Periconia")
w_v <- c("temp_before_hd", "rhum_before_hd", "prcp_before_hd", "temp_hd", "rhum_hd", "prcp_hd")
for(wv in c(1:length(w_v))){
  #wv <- 1
  target_w <- w_v[wv]
  imsi_alt  <- extract_smooth_var_DVR(m_alt, target_w, df_for_more_analysis3)
  imsi_epi  <- extract_smooth_var_DVR(m_epi, target_w, df_for_more_analysis3)
  imsi_han  <- extract_smooth_var_DVR(m_han, target_w, df_for_more_analysis3)
  imsi_per  <- extract_smooth_var_DVR(m_per, target_w, df_for_more_analysis3)
  imsi_fhb  <- extract_smooth_var_DVR(m_fhb, target_w, df_for_more_analysis3)
  
  imsi_alt$taxon <- "Alternaria"
  imsi_epi$taxon <- "Epicoccum"
  imsi_han$taxon <- "Hannaella"
  imsi_per$taxon <- "Periconia"
  imsi_fhb$taxon <- "FHB"
  
  imsi_all <- rbind(imsi_alt, imsi_epi, imsi_han, imsi_per, imsi_fhb)
  imsi_all$taxon <- factor(imsi_all$taxon, levels = c("FHB", "Alternaria", "Epicoccum", "Hannaella", "Periconia"))
  
  imsi_g <- ggplot(imsi_all, aes(x = x, y = est, color = taxon, fill = taxon)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = est - se, ymax = est + se),
                alpha = 0.15, color = NA) +
    scale_color_manual(values = c("FHB" = "black", "Alternaria" = "#ff6347",
                                  "Epicoccum" = "#daa520", "Hannaella" = "#006400",
                                  "Periconia" = "#4682b4")) +
    scale_fill_manual(values = c("FHB" = "black", "Alternaria" = "#ff6347",
                                 "Epicoccum" = "#daa520", "Hannaella" = "#006400",
                                 "Periconia" = "#4682b4")) +
    labs(
      x = paste0("\n",target_w),
      y = paste0("Smooth effect","\n"),
      title = paste0("Shared niche test: Response to ", target_w)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                                 size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    )
  # imsi_g  
  # geom_hline(yintercept=0, linetype='dashed', color='red', size = 0.75)
  if(!dir.exists(paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve"))){
    dir.create(paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve"))
  }
  ggsave(
    plot = imsi_g,
    file = paste0(
      "./Output/3. obs_pred_graph/2025/",
      version_pred_obs,
      "/response_curve/DVR_response_curve_",target_w,"_",
      version_pred_obs,
      ".png"
    ),
    width = 12,
    height = 9,
    units = c("cm")
  )
}

w_v2 <- c("temp_before_hd_sampling_10days", "rhum_before_hd_sampling_10days", "prcp_before_hd_sampling_10days", 
          "temp_hd_sampling_10days", "rhum_hd_sampling_10days", "prcp_hd_sampling_10days")
for(wv2 in c(1:length(w_v2))){
  #wv2 <- 1
  target_w <- w_v2[wv2]
  imsi_alt  <- extract_smooth_var_sampling(m_alt_sampling, target_w, df_for_more_analysis3)
  imsi_epi  <- extract_smooth_var_sampling(m_epi_sampling, target_w, df_for_more_analysis3)
  imsi_han  <- extract_smooth_var_sampling(m_han_sampling, target_w, df_for_more_analysis3)
  imsi_per  <- extract_smooth_var_sampling(m_per_sampling, target_w, df_for_more_analysis3)
  imsi_fhb  <- extract_smooth_var_sampling(m_fhb_sampling, target_w, df_for_more_analysis3)
  
  imsi_alt$taxon <- "Alternaria"
  imsi_epi$taxon <- "Epicoccum"
  imsi_han$taxon <- "Hannaella"
  imsi_per$taxon <- "Periconia"
  imsi_fhb$taxon <- "FHB"
  
  imsi_all <- rbind(imsi_alt, imsi_epi, imsi_han, imsi_per, imsi_fhb)
  imsi_all$taxon <- factor(imsi_all$taxon, levels = c("FHB", "Alternaria", "Epicoccum", "Hannaella", "Periconia"))
  
  imsi_g <- ggplot(imsi_all, aes(x = x, y = est, color = taxon, fill = taxon)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = est - se, ymax = est + se),
                alpha = 0.15, color = NA) +
    scale_color_manual(values = c("FHB" = "black", "Alternaria" = "#ff6347",
                                  "Epicoccum" = "#daa520", "Hannaella" = "#006400",
                                  "Periconia" = "#4682b4")) +
    scale_fill_manual(values = c("FHB" = "black", "Alternaria" = "#ff6347",
                                 "Epicoccum" = "#daa520", "Hannaella" = "#006400",
                                 "Periconia" = "#4682b4")) +
    labs(
      x = paste0("\n",target_w),
      y = paste0("Smooth effect","\n"),
      title = paste0("Shared niche test: Response to ", target_w)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                                 size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    )
  # imsi_g  
  # geom_hline(yintercept=0, linetype='dashed', color='red', size = 0.75)
  if(!dir.exists(paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve"))){
    dir.create(paste0("./Output/3. obs_pred_graph/2025/",version_pred_obs,"/response_curve"))
  }
  ggsave(
    plot = imsi_g,
    file = paste0(
      "./Output/3. obs_pred_graph/2025/",
      version_pred_obs,
      "/response_curve/sampling_response_curve_",target_w,"_",
      version_pred_obs,
      ".png"
    ),
    width = 12,
    height = 9,
    units = c("cm")
  )
}

# correlation
cor(rh_hd_alt$est, rh_hd_fhb$est, method = "spearman")
cor(rh_hd_epi$est, rh_hd_fhb$est, method = "spearman")
cor(rh_hd_han$est, rh_hd_fhb$est, method = "spearman")
cor(rh_hd_per$est, rh_hd_fhb$est, method = "spearman")

cor(rh_hd_alt$est, rh_hd_fhb$est, method = "pearson")
cor(rh_hd_epi$est, rh_hd_fhb$est, method = "pearson")
cor(rh_hd_han$est, rh_hd_fhb$est, method = "pearson")
cor(rh_hd_per$est, rh_hd_fhb$est, method = "pearson")

weather_vars <- c("temp_before_hd", "rhum_before_hd", "prcp_before_hd",
                  "temp_hd", "rhum_hd", "prcp_hd")
weather_vars_sampling <- c("temp_before_hd_sampling_10days", "rhum_before_hd_sampling_10days", "prcp_before_hd_sampling_10days", 
                  "temp_hd_sampling_10days", "rhum_hd_sampling_10days", "prcp_hd_sampling_10days")

compare_to_fhb_DVR <- function(model_microbe, model_fhb, varname, data) {
  s1 <- extract_smooth_var_DVR(model_microbe, varname, data)
  s2 <- extract_smooth_var_DVR(model_fhb, varname, data)
  stats::cor(s1$est, s2$est, method = "pearson")
}

compare_to_fhb_sampling <- function(model_microbe, model_fhb, varname, data) {
  s1 <- extract_smooth_var_sampling(model_microbe, varname, data)
  s2 <- extract_smooth_var_sampling(model_fhb, varname, data)
  stats::cor(s1$est, s2$est, method = "pearson")
}

alt_result_matrix <- sapply(weather_vars, function(v)compare_to_fhb_DVR(m_alt, m_fhb, v, df_for_more_analysis3))
epi_result_matrix <- sapply(weather_vars, function(v)compare_to_fhb_DVR(m_epi, m_fhb, v, df_for_more_analysis3))
han_result_matrix <- sapply(weather_vars, function(v)compare_to_fhb_DVR(m_han, m_fhb, v, df_for_more_analysis3))
per_result_matrix <- sapply(weather_vars, function(v)compare_to_fhb_DVR(m_per, m_fhb, v, df_for_more_analysis3))

result_df_DVR <- rbind(
  Alternaria = alt_result_matrix,
  Epicoccum = epi_result_matrix,
  Hannaella = han_result_matrix,
  Periconia = per_result_matrix
)

alt_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_alt_sampling, m_fhb_sampling, v, df_for_more_analysis3))
epi_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_epi_sampling, m_fhb_sampling, v, df_for_more_analysis3))
han_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_han_sampling, m_fhb_sampling, v, df_for_more_analysis3))
per_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_per_sampling, m_fhb_sampling, v, df_for_more_analysis3))

result_df_sampling <- rbind(
  Alternaria = alt_result_matrix_sampling,
  Epicoccum = epi_result_matrix_sampling,
  Hannaella = han_result_matrix_sampling,
  Periconia = per_result_matrix_sampling
)

# Similarity in the effects of weather on FHB and microorganisms
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. effect matrix??? ??????
heat_df <- as.data.frame(result_df_DVR) %>%
  tibble::rownames_to_column("Microbe") %>%
  tidyr::pivot_longer(-Microbe,
                      names_to = "Weather",
                      values_to = "Correlation")

# 2. p-value + edf 추출 ??????
extract_smooth_info <- function(model, name){
  s_tab <- summary(model)$s.table
  data.frame(
    Microbe = name,
    Weather = weather_vars, #rownames(s_tab)
    edf = s_tab[,"edf"],
    pval = s_tab[,"p-value"]
  )
}

smooth_df <- bind_rows(
  extract_smooth_info(m_alt, "Alternaria"),
  extract_smooth_info(m_epi, "Epicoccum"),
  extract_smooth_info(m_han, "Hannaella"),
  extract_smooth_info(m_per, "Periconia"),
  extract_smooth_info(m_fhb, "FHB")
)

smooth_df <- smooth_df %>%
  group_by(Weather) %>%
  mutate(edf2 = edf / edf[Microbe == "FHB"]) %>%
  ungroup()

# 3. 병합
plot_df <- left_join(heat_df, smooth_df,
                     by=c("Microbe","Weather"))
as.data.frame(plot_df)
# p > 0.05 ??? 경우 ?????? 처리??? 변???
plot_df <- plot_df %>%
  mutate(sig = ifelse(pval <= 0.05, "sig", "ns"))

plot_df$Weather <- factor(plot_df$Weather, levels = c("temp_before_hd", "rhum_before_hd", "prcp_before_hd", 
                                                      "temp_hd", "rhum_hd", "prcp_hd"))

FHB_M_sim_g <- ggplot(plot_df, aes(x=Weather, y=Microbe)) +
  geom_point(aes(size=edf2,
                 fill=ifelse(sig=="sig", Correlation, NA)),
             shape=21, color= "black") +
  geom_text(aes(label = round(ifelse(sig=="sig", Correlation, NA), 2)), 
            size = 4, fontface = "bold", color = "black") +
  labs(
    title = "Similarity of weather-response patterns \nbetween FHB incidence and microbes", 
    x = "\nWeather variables",
    y = "FHB incidence and microbial responses\n"
  )+ 
  # xlab("\nWeather variables") + 
  # ylab("Microbes\n") + 
  scale_fill_gradient2(
    name = "Weather-response \nsimilarity (Rho)",
    low = "#08979D",#"#2166AC",
    mid = "white",#"white",
    high = "#8474A1",#"#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value="#d3d3d3"
  ) + 
  scale_size(name   = "Similarity of \nresponse complexity",range=c(10,20),breaks = c(0,1,2,3),limits = c(0,3)) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,face='bold')) + 
  theme(axis.title.y = element_text(size = 15, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, size=15, face='bold',color='black', hjust = 1))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black')) + 
  theme(legend.title = element_text(size = 15, face="bold"),
        legend.text  = element_text(size = 14))
# theme(panel.grid.major = element_blank())
FHB_M_sim_g
ggsave(
  plot = FHB_M_sim_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/response_curve/Similarity_wth_effects_on_FHB_and_mb_",
    version_pred_obs,
    ".png"
  ),
  width = 22,
  height = 18,
  units = c("cm")
)

####################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
weather_vars_sampling <- c("temp_before_hd_sampling_10days", "rhum_before_hd_sampling_10days", "prcp_before_hd_sampling_10days", 
                           "temp_hd_sampling_10days", "rhum_hd_sampling_10days", "prcp_hd_sampling_10days")

compare_to_fhb_sampling <- function(model_microbe, model_fhb, varname, data) {
  s1 <- extract_smooth_var_sampling(model_microbe, varname, data)
  s2 <- extract_smooth_var_sampling(model_fhb, varname, data)
  stats::cor(s1$est, s2$est, method = "pearson")
}

alt_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_alt_sampling, m_fhb_sampling, v, df_for_more_analysis3))
epi_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_epi_sampling, m_fhb_sampling, v, df_for_more_analysis3))
han_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_han_sampling, m_fhb_sampling, v, df_for_more_analysis3))
per_result_matrix_sampling <- sapply(weather_vars_sampling, function(v)compare_to_fhb_sampling(m_per_sampling, m_fhb_sampling, v, df_for_more_analysis3))

result_df_sampling <- rbind(
  Alternaria = alt_result_matrix_sampling,
  Epicoccum = epi_result_matrix_sampling,
  Hannaella = han_result_matrix_sampling,
  Periconia = per_result_matrix_sampling
)

# Similarity in the effects of weather on FHB and microorganisms
# 1. effect matrix??? ??????
heat_df <- as.data.frame(result_df_sampling) %>%
  tibble::rownames_to_column("Microbe") %>%
  tidyr::pivot_longer(-Microbe,
                      names_to = "Weather",
                      values_to = "Correlation")

# 2. p-value + edf 추출 ??????
extract_smooth_info <- function(model, name){
  s_tab <- summary(model)$s.table
  data.frame(
    Microbe = name,
    Weather = weather_vars_sampling, #rownames(s_tab)
    edf = s_tab[,"edf"],
    pval = s_tab[,"p-value"]
  )
}

smooth_df <- bind_rows(
  extract_smooth_info(m_alt_sampling, "Alternaria"),
  extract_smooth_info(m_epi_sampling, "Epicoccum"),
  extract_smooth_info(m_han_sampling, "Hannaella"),
  extract_smooth_info(m_per_sampling, "Periconia"),
  extract_smooth_info(m_fhb_sampling, "FHB")
)

smooth_df <- smooth_df %>%
  group_by(Weather) %>%
  mutate(edf2 = edf / edf[Microbe == "FHB"]) %>%
  ungroup()

# 3. 병합
plot_df <- left_join(heat_df, smooth_df,
                     by=c("Microbe","Weather"))
as.data.frame(plot_df)
# p > 0.05 ??? 경우 ?????? 처리??? 변???
plot_df <- plot_df %>%
  mutate(sig = ifelse(pval <= 0.05, "sig", "ns"))
plot_df$Weather <- str_replace(plot_df$Weather, "_sampling_10days", "")
plot_df$Weather <- factor(plot_df$Weather, levels = c("temp_before_hd", "rhum_before_hd", "prcp_before_hd", 
                                                      "temp_hd", "rhum_hd", "prcp_hd"))
plot_df$edf2
FHB_M_sim_g <- ggplot(plot_df, aes(x=Weather, y=Microbe)) +
  geom_point(aes(size=edf2,
                 fill=ifelse(sig=="sig", Correlation, NA)),
             shape=21, color= "black") +
  geom_text(aes(label = round(ifelse(sig=="sig", Correlation, NA), 2)), 
            size = 4, fontface = "bold", color = "black") +
  labs(
    title = "Similarity of weather-response patterns \nbetween FHB incidence and microbes", 
    x = "\nWeather variables",
    y = "FHB incidence and microbial responses\n"
  )+ 
  # xlab("\nWeather variables") + 
  # ylab("Microbes\n") + 
  scale_fill_gradient2(
    name = "Weather-response \nsimilarity (Rho)",
    low = "#08979D",#"#2166AC",
    mid = "white",#"white",
    high = "#8474A1",#"#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value="#d3d3d3"
  ) + 
  scale_size(name   = "Similarity of \nresponse complexity",range=c(10,20),breaks = c(0,1,2,3,4),limits = c(0,4)) +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,face='bold')) + 
  theme(axis.title.y = element_text(size = 15, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, size=15, face='bold',color='black', hjust = 1))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black')) + 
  theme(legend.title = element_text(size = 15, face="bold"),
        legend.text  = element_text(size = 14))
# theme(panel.grid.major = element_blank())
FHB_M_sim_g
ggsave(
  plot = FHB_M_sim_g,
  file = paste0(
    "./Output/3. obs_pred_graph/2025/",
    version_pred_obs,
    "/response_curve/sampling_Similarity_wth_effects_on_FHB_and_mb_",
    version_pred_obs,
    ".png"
  ),
  width = 22,
  height = 18,
  units = c("cm")
)

# correlation MB&FHB vs.wth (bipartite correlation network)
head(as.data.frame(df_for_more_analysis3))
stats::cor(df_for_more_analysis3)

print("Run `4.correlation_plot_v12.R` code for graph")

# library(Hmisc)
# bio_vars <- c("Alternaria","Cladosporium","Epicoccum",
#               "Hannaella","Papiliotrema","Periconia","incidence")
# 
# weather_vars <- c("temp_before_hd","rhum_before_hd","prcp_before_hd",
#                   "temp_hd","rhum_hd","prcp_hd")
# sub_df <- df_for_more_analysis3[, c(bio_vars, weather_vars)]
# cor_res <- rcorr(as.matrix(sub_df), type="pearson")
# cor_mat <- cor_res$r
# p_mat   <- cor_res$P
# 
# edges <- expand.grid(from = bio_vars,
#                      to   = weather_vars)
# 
# edges$r <- mapply(function(x,y) cor_mat[x,y],
#                   edges$from, edges$to)
# 
# edges$p <- mapply(function(x,y) p_mat[x,y],
#                   edges$from, edges$to)
# 
# # ????????? 것만
# edges <- edges[edges$p < 0.05, ]
# 
# edges$weight <- abs(edges$r)
# edges$color  <- ifelse(edges$r > 0, "positive", "negative")
# 
# nodes <- data.frame(
#   name = c(bio_vars, weather_vars),
#   type = c(rep("Biological", length(bio_vars)),
#            rep("Weather", length(weather_vars)))
# )
# library(tidygraph)
# edges_com <- edges[complete.cases(edges),]
# net <- tbl_graph(nodes = nodes,
#                  edges = edges_com,
#                  directed = FALSE)
# 
# library(dplyr)
# 
# bio_n <- length(bio_vars)
# wea_n <- length(weather_vars)
# 
# nodes <- nodes %>%
#   mutate(
#     angle = c(seq(pi*0.75, pi*1.25, length.out = bio_n),
#               seq(-pi*0.25, pi*0.25, length.out = wea_n)),
#     x = cos(angle),
#     y = sin(angle)
#   )
# 
# net <- tbl_graph(nodes = nodes,
#                  edges = edges_com,
#                  directed = FALSE)
# library(ggraph)
# 
# ggraph(net, layout = "manual",
#        x = nodes$x,
#        y = nodes$y) +
#   
#   geom_edge_link(aes(width = weight,
#                      color = color),
#                  alpha = 0.8) +
#   
#   scale_edge_width(range = c(0.5, 3)) +
#   
#   scale_edge_color_manual(values = c(
#     positive = "#B2182B",
#     negative = "#2166AC"
#   )) +
#   
#   geom_node_point(aes(fill = type),
#                   shape = 21,
#                   size = 8,
#                   color = "black") +
#   
#   geom_node_text(aes(label = name),
#                  repel = TRUE,
#                  size = 4) +
#   
#   scale_fill_manual(values = c(
#     Biological = "#E69F00",
#     Weather    = "#56B4E9"
#   )) +
#   
#   theme_void() +
#   theme(legend.position = "right")

########################################################################################




# 
# ggplot(heat_df, aes(x = Weather, y = Microbe, fill = Correlation)) +
#   geom_tile(color = "white") +
#   geom_text(aes(label = round(Correlation, 2)), size = 4) +
# scale_fill_gradient2(
#   low = "#2166AC",
#   mid = "white",
#   high = "#B2182B",
#   midpoint = 0,
#   limits = c(-1, 1)
# ) +
#   theme_minimal(base_size = 14) +
#   labs(fill = "Correlation",
#        x = "Weather variable",
#        y = "Microbe") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 
# 
# ### Does the distribution of key_m really differ depending on the FHB incidence?
# df_for_more_analysis2 <- df_rslt_stat
# names(df_for_more_analysis2)
# 
# hist(df_for_more_analysis2$incidence)
# df_for_more_analysis2 <- df_for_more_analysis2 %>%
#   mutate(
#     cluster = case_when(
#       incidence < 0.2 ~ "low",
#       incidence < 0.5 ~ "mid",
#       incidence > 0.5 ~ "high",
#       TRUE ~ "NA"
#     )
#   )
# df_for_more_analysis2$cluster <- factor(df_for_more_analysis2$cluster, 
#                                         levels = c("low", "mid", "high"))
# table(df_for_more_analysis2$cluster)
# 
# summary(gam_model)
# df_for_more_analysis2_for_boxplot <- df_for_more_analysis2 %>% dplyr::select("Alternaria", "Cladosporium", "Epicoccum",
#                                                                              "Hannaella", "Papiliotrema", "Periconia", 
#                                                                              "incidence", "cluster")
# # boxplot
# library(agricolae)
# boxplot_mb_vs_inc <- function(target_microbe){
#   target_microbe <- target_microbe #"Alternaria"
#   
#   kw_res <- kruskal(
#     y = df_for_more_analysis2_for_boxplot[,target_microbe],
#     trt = df_for_more_analysis2_for_boxplot$cluster,
#     alpha = 0.05,
#     group = TRUE
#   )
#   
#   letter_df <- kw_res$groups %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("cluster")
#   
#   df_for_more_analysis2_for_boxplot$cluster <- as.factor(df_for_more_analysis2_for_boxplot$cluster)
#   y_pos <- as.data.frame(df_for_more_analysis2_for_boxplot) %>%
#     dplyr::group_by(cluster) %>% 
#     dplyr::summarise(y = max(get(target_microbe)) * 1.05)
#   
#   letter_df <- left_join(letter_df, y_pos, by = "cluster")
#   
#   g <- ggplot(df_for_more_analysis2_for_boxplot, aes(x = cluster, y = get(target_microbe))) +
#     geom_boxplot(outlier.shape = NA, fill = "grey90") +
#     geom_jitter(width = 0.15, alpha = 0.6) +
#     geom_text(
#       data = letter_df,
#       aes(x = cluster, y = y, label = groups),
#       size = 5
#     ) +
#     theme_bw() +
#     labs(
#       title = paste0(target_microbe, " abundance by cluster"), 
#       y = target_microbe
#     ) + 
#     theme_bw() +
#     theme(
#       axis.title = element_text(size = 14, face = "bold"),
#       axis.text  = element_text(size = 12, face = "bold"),
#       plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
#     )
#   return(g)
# }
# 
# boxplot_mb_vs_inc(names(df_for_more_analysis2_for_boxplot)[1])
# boxplot_mb_vs_inc(names(df_for_more_analysis2_for_boxplot)[2])
# boxplot_mb_vs_inc(names(df_for_more_analysis2_for_boxplot)[3])
# boxplot_mb_vs_inc(names(df_for_more_analysis2_for_boxplot)[4])
# boxplot_mb_vs_inc(names(df_for_more_analysis2_for_boxplot)[5])
# boxplot_mb_vs_inc(names(df_for_more_analysis2_for_boxplot)[6])
# 
# plot(df_for_more_analysis2_for_boxplot$Alternaria, df_for_more_analysis2_for_boxplot$incidence)


# ### 2D PDP for GAM
# library(pdp)
# 
# cms_AC <- pdp::partial(
#   object = gam_model,
#   pred.var = c("Alternaria", "Cladosporium"),
#   train = data_gibsim_MB_for_beta,
#   grid.resolution = 30,
#   pred.fun = function(object, newdata) {
#     predict(object, newdata = newdata, type = "response")
#   }
# )
# 
# ggplot(cms_AC,
#        aes(x = Alternaria,
#            y = Cladosporium,
#            fill = yhat)) +
#   geom_raster(interpolate = TRUE) +
#   geom_contour(aes(z = yhat),
#                color = "black",
#                alpha = 0.5) +
#   scale_fill_viridis_c(
#     name = "Conditional mean\nFHB incidence"
#   ) +
#   geom_point(
#     data = data_gibsim_MB_for_beta,
#     aes(x = Alternaria, y = Cladosporium),
#     inherit.aes = FALSE,
#     size = 0.7,
#     alpha = 0.4
#   )+
#   theme_bw() +
#   labs(
#     title = "Conditional Mean Surface",
#     subtitle = "Alternaria × Cladosporium (GAM)",
#     x = "Alternaria",
#     y = "Cladosporium"
#   )
# 
# cms_A  <- pdp::partial(gam_model, "Alternaria", train = data_gibsim_MB_for_beta,
#                        grid.resolution = 30,
#                        pred.fun = function(object, newdata) {
#                          predict(object, newdata = newdata, type = "response")
#                        })
# cms_C  <- pdp::partial(gam_model, "Cladosporium", train = data_gibsim_MB_for_beta,
#                        grid.resolution = 30,
#                        pred.fun = function(object, newdata) {
#                          predict(object, newdata = newdata, type = "response")
#                        })
# cms_2d <- pdp::partial(
#   gam_model,
#   pred.var = c("Alternaria", "Cladosporium"),
#   train = data_gibsim_MB_for_beta,
#   grid.resolution = 30,
#   pred.fun = function(object, newdata) {
#     predict(object, newdata, type = "response")
#   }
# )
# 
# names(cms_2d)[3] <- "yhat_2d"
# names(cms_A)[2] <- "yhat_A"
# names(cms_C)[2] <- "yhat_C"
# names(cms_A)[3] <- names(cms_C)[3] <- "yhat.id"
# additive_cms <- cms_2d %>%
#   left_join(cms_A,
#             by = c("Alternaria", "yhat.id")) %>%
#   left_join(cms_C,
#             by = c("Cladosporium", "yhat.id")) %>%
#   mutate(
#     yhat_add = yhat_A + yhat_C - mean(cms_2d$yhat_2d)
#   )
# interaction_cms <- additive_cms %>%
#   mutate(
#     yhat_int = yhat_2d - yhat_add
#   )
# 
# library(ggplot2)
# 
# ggplot(interaction_cms,
#        aes(x = Alternaria,
#            y = Cladosporium,
#            fill = yhat_int)) +
#   geom_raster(interpolate = F) +
#   geom_contour(aes(z = yhat_int),
#                color = "black",
#                alpha = 0.6) +
#   scale_fill_gradient2(
#     low = "blue",
#     mid = "white",
#     high = "red",
#     midpoint = 0,
#     name = "Interaction\n(surface)"
#   ) +
#   theme_bw() +
#   labs(
#     title = "GAM interaction surface",
#     subtitle = "2D CMS ??? additive CMS",
#     x = "Alternaria",
#     y = "Cladosporium"
#   )
# 
# 
# 
# # GAM??? ti()가 ?????? ?????? 가???
# vis.gam(
#   gam_model,
#   view = c("Alternaria", "Cladosporium"),
#   plot.type = "contour",
#   too.far = 0.1
# )
# 
# 
# 
# head(data_gibsim_MB_for_beta)
# 
# interaction_vec_GAM <- GAM_imp_variables_Top[grep("_", GAM_imp_variables_Top)]
# interaction_vec <- interaction_vec[-which(interaction_vec == "FHB_incidence")]
# for_loop_df_for_2d_PDP <- data.frame(interaction = interaction_vec) %>%
#   tidyr::separate(interaction, into = c("interaction_v1", "interaction_v2"), sep = "_")
# 
# # interaction_v1 <- interaction_df$V1 #c("OTU0001", "OTU0001", "OTU0004", "OTU0004", "OTU0004", "OTU0004", "OTU0006", "OTU0026", "OTU0027", "OTU0039")
# # interaction_v2 <- interaction_df$V2 #c("OTU0004", "OTU0039", "OTU0006", "OTU0026", "OTU0038", "OTU0056", "OTU0026", "OTU0038", "OTU0075", "OTU0056")
# # for_loop_df_for_2d_PDP <- data.frame(interaction_v1, interaction_v2)
# 
# all_vals <- c()
# for (x in 1:nrow(for_loop_df_for_2d_PDP)) {
#   imsi_pdp <- pdp::partial(rf_model_for_analysis, 
#                            pred.var = c(for_loop_df_for_2d_PDP$interaction_v1[x], for_loop_df_for_2d_PDP$interaction_v2[x]), 
#                            train = using_data, chull = TRUE)
#   
#   gam_imsi_pdp_2 <- pdp::partial(gam_model, 
#                                  pred.var = c("Epicoccum", "Cladosporium"),
#                                  train = data_gibsim_MB_for_beta)
#   pdp_2d <- pdp::partial(
#     object = gam_model,
#     pred.var = c("Epicoccum", "Cladosporium"),
#     train = data_gibsim_MB_for_beta,
#     grid.resolution = 30,
#     pred.fun = function(object, newdata) {
#       predict(object, newdata = newdata, type = "response")
#     }
#   )
#   
#   table(pdp_2d$yhat.id)
#   hist(pdp_2d$yhat)
#   
#   # Enhanced plotting with ggplot2
#   ggplot(pdp_2d, aes(x = Epicoccum, y = Cladosporium, fill = yhat)) +
#     geom_tile() +
#     scale_fill_viridis_c(option = "rocket") + # Use a colorblind-friendly color map
#     labs(title = "2D Partial Dependence Plot for Median House Value",
#          x = "Epicoccum",
#          y = "Cladosporium",
#          fill = "FHB incidence") +
#     theme_minimal()
#   
#   pdp::partial(gam_model, pred.var = c("Alternaira", "Cladosporium"), plot = TRUE, chull = TRUE)
#   
#   GAM_imsi_2d_pdp <- plotPartial(
#     pdp_2d,
#     levelplot   = TRUE, 
#     # contour     = TRUE,
#     # contour.color = "#a0522d",
#     col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#     xlab        = list(cex=1.2, font=2),
#     ylab        = list(cex=1.2, font=2),
#     main        = list(label=title_lab, cex=1.3, font=2),
#     par.settings = list(axis.text = list(cex=1.1, col="black"))
#   )
#   
#   ggplot(pdp_2d, aes(
#     x = Epicoccum,
#     y = Cladosporium,
#     fill = yhat
#   )) +
#     geom_raster(interpolate = TRUE) +
#     geom_contour(aes(z = yhat), color = "black", alpha = 0.5) +
#     scale_fill_viridis_c(name = "Predicted\nFHB incidence") +
#     theme_minimal() +
#     labs(
#       title = "2D Partial Dependence Plot",
#       subtitle = "Alternaria × Cladosporium",
#       x = "Alternaria (CLR / transformed)",
#       y = "Cladosporium (CLR / transformed)"
#     )
#   
#   ggplot(pdp_2d, aes(
#     x = Epicoccum,
#     y = Cladosporium,
#     fill = yhat
#   )) +
#     geom_raster(interpolate = TRUE) +
#     geom_contour(aes(z = yhat), color = "black", alpha = 0.4) +
#     geom_point(
#       data = data_gibsim_MB_for_beta,
#       aes(x = Alternaria, y = Cladosporium),
#       inherit.aes = FALSE,
#       size = 0.8,
#       alpha = 0.4
#     ) +
#     scale_fill_viridis_c() +
#     theme_bw()
#   
#   
#   
#   
#   
#   
#   all_vals <- c(all_vals,imsi_pdp$yhat)
#   print(x)
# }
# 
# # (2) ?????? global_min / global_max 계산 + floor/ceiling
# MB_MB_global_min <- floor(min(all_vals))
# MB_MB_global_max <- ceiling(max(all_vals))
# 
# imsi_2d_pdp_combind_not_fixed <- list()
# imsi_2d_pdp_combind_fixed <- list()
# 
# for(x in c(1:nrow(for_loop_df_for_2d_PDP))){
#   imsi_pdp <- pdp::partial(rf_model_for_analysis, pred.var = c(for_loop_df_for_2d_PDP$interaction_v1[x], for_loop_df_for_2d_PDP$interaction_v2[x]), 
#                            train = using_data, chull = TRUE)
#   
#   imsi_2d_pdp <- plotPartial(
#     imsi_pdp,
#     levelplot   = TRUE, 
#     # contour     = TRUE,
#     # contour.color = "#a0522d",
#     col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#     xlab        = list(cex=1.2, font=2),
#     ylab        = list(cex=1.2, font=2),
#     main        = list(label=title_lab, cex=1.3, font=2),
#     par.settings = list(axis.text = list(cex=1.1, col="black"))
#   )
#   
#   imsi_2d_pdp_combind_not_fixed[[x]] <- imsi_2d_pdp
#   
#   imsi_2d_pdp_fixed <- plotPartial(
#     imsi_pdp,
#     levelplot   = TRUE, 
#     # contour     = TRUE,
#     # contour.color = "#a0522d",
#     col.regions = viridis::rocket(100, direction = -1), #colorRampPalette(c("yellow","red"))(100)
#     at          = seq(MB_MB_global_min , MB_MB_global_max, length.out = 100),
#     xlab        = list(cex=1.2, font=2),
#     ylab        = list(cex=1.2, font=2),
#     main        = list(label=title_lab, cex=1.3, font=2),
#     par.settings = list(axis.text = list(cex=1.1, col="black"))
#   )
#   
#   imsi_2d_pdp_combind_fixed[[x]] <- imsi_2d_pdp_fixed
#   
#   
#   print(x)
# }
# 
# combined_plot_MB_MB_2d_pdp_not_fixed <- grid.arrange(
#   imsi_2d_pdp_combind_not_fixed[[1]], imsi_2d_pdp_combind_not_fixed[[2]], 
#   imsi_2d_pdp_combind_not_fixed[[3]], imsi_2d_pdp_combind_not_fixed[[4]],
#   imsi_2d_pdp_combind_not_fixed[[5]], #imsi_2d_pdp_combind_not_fixed[[6]],
#   # imsi_2d_pdp_combind_not_fixed[[7]], imsi_2d_pdp_combind_not_fixed[[8]],
#   # imsi_2d_pdp_combind_not_fixed[[9]], imsi_2d_pdp_combind_not_fixed[[10]],
#   # imsi_2d_pdp_combind_not_fixed[[11]], imsi_2d_pdp_combind_not_fixed[[12]],
#   # imsi_2d_pdp_combind_not_fixed[[13]], imsi_2d_pdp_combind_not_fixed[[14]],
#   # imsi_2d_pdp_combind_not_fixed[[15]], imsi_2d_pdp_combind_not_fixed[[16]],
#   ncol = 3
# )
# 
# combined_plot_MB_MB_2d_pdp_fixed <- grid.arrange(
#   imsi_2d_pdp_combind_fixed[[1]], imsi_2d_pdp_combind_fixed[[2]], 
#   imsi_2d_pdp_combind_fixed[[3]],
#   imsi_2d_pdp_combind_fixed[[4]],
#   imsi_2d_pdp_combind_fixed[[5]], #imsi_2d_pdp_combind_fixed[[6]],
#   # imsi_2d_pdp_combind_fixed[[7]], imsi_2d_pdp_combind_fixed[[8]],
#   # imsi_2d_pdp_combind_fixed[[9]], imsi_2d_pdp_combind_fixed[[10]],
#   # imsi_2d_pdp_combind_fixed[[11]], imsi_2d_pdp_combind_fixed[[12]],
#   # imsi_2d_pdp_combind_fixed[[13]], imsi_2d_pdp_combind_fixed[[14]],
#   # imsi_2d_pdp_combind_fixed[[15]], imsi_2d_pdp_combind_fixed[[16]],
#   ncol = 3
# )
# 
# ggsave(
#   plot = combined_plot_MB_MB_2d_pdp_not_fixed,
#   file = paste0(
#     "./Output/3. obs_pred_graph/2025/",
#     version_pred_obs,
#     "/2D_PDP_MB_MB_plot_not_fixed_legend",
#     version_pred_obs,
#     ".png"
#   ),
#   width = 30,
#   height = 20,
#   units = c("cm")
# )
# 
# 
# ggsave(
#   plot = combined_plot_MB_MB_2d_pdp_fixed,
#   file = paste0(
#     "./Output/3. obs_pred_graph/2025/",
#     version_pred_obs,
#     "/2D_PDP_MB_MB_plot_fixed_legend", 
#     version_pred_obs,
#     ".png"
#   ),
#   width = 30,
#   height = 20,
#   units = c("cm")
# )

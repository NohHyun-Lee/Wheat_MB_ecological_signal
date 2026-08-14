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
df_rslt_stat <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v4.xlsx") 
names(df_rslt_stat)
head(df_rslt_stat)
nrow(df_rslt_stat)

### group by loc -----------------------------------------------------------------------------------------------------
nrow(df_rslt_stat)
table(df_rslt_stat$loc)

MB_intensity_wth_cor_df <- df_rslt_stat
colnames(MB_intensity_wth_cor_df)
MB_intensity_wth_cor_df <- MB_intensity_wth_cor_df %>% dplyr::select(-contains("wsd")) # select(!contains("wsd"))
version_pred_obs <- "2025_v12"

sum(is.na(MB_intensity_wth_cor_df == T))
MB_intensity_wth_cor_df <- MB_intensity_wth_cor_df[complete.cases(MB_intensity_wth_cor_df),]
nrow(MB_intensity_wth_cor_df)

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

##############################################################################################################
##### Genus level ##########################################################################################################################################
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

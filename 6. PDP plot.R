# library(gratia)
library(ggplot2)

# 모델에 포함된 모든 smooth 이름 가져오기
smooth_names <- smooths(gam_model)

# 저장 경로 지정
save_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/gratia/"

for (var in vars) {
  
  x_seq <- seq(min(data_gibsim_MB_for_beta[[var]], na.rm = TRUE),
               max(data_gibsim_MB_for_beta[[var]], na.rm = TRUE),
               length.out = 50)
  
  ice_list <- list()
  
  for (i in 1:nrow(data_gibsim_MB_for_beta)) {
    
    row_i <- data_gibsim_MB_for_beta[i, ]
    
    newdata <- row_i[rep(1, length(x_seq)), ]
    newdata[[var]] <- x_seq
    
    # baseline
    baseline <- newdata
    baseline[[var]] <- min(data_gibsim_MB_for_beta[[var]], na.rm = TRUE)
    
    pred <- predict(gam_model, newdata = newdata, type = "response")
    pred_base <- predict(gam_model, newdata = baseline, type = "response")
    
    delta <- (pred - pred_base) * 100
    
    ice_list[[i]] <- data.frame(
      x = x_seq,
      delta = delta,
      id = i
    )
  }
  
  ice_df <- bind_rows(ice_list)
  
  # 평균 (PDP)
  pdp_df <- ice_df %>%
    group_by(x) %>%
    summarise(mean_delta = mean(delta), .groups = "drop")
  
  p <- ggplot() +
    geom_line(data = ice_df,
              aes(x = x, y = delta, group = id),
              size = 0.5, alpha = 0.4, color = "black") +
    geom_line(data = pdp_df,
              aes(x = x, y = mean_delta),
              size = 1.5, color = "red") +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = var,
      y = "Δ FHB incidence (%)",
      title = paste("ICE plot :", var)
    ) +
    theme_classic()+
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                                 size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    )
  
  print(p)
  
  ggsave( filename = paste0(save_dir, "GAM_delta_ICE_plot_", var, ".png"), plot = p, width = 4, height = 3 )
  
  
}

vars1 <- c("Alternaria", "Alternaria", "Alternaria", "Epicoccum",  "Epicoccum" )
vars2 <- c("Cladosporium", "Papiliotrema", "Periconia" ,   "Cladosporium", "Periconia")
vars1 <- "Alternaria"
vars2 <- "Periconia"
for(k in c(1:length(vars1))){
  
  var1 <- vars1[k]
  var2 <- vars2[k]
  
  x_seq <- seq(min(data_gibsim_MB_for_beta[[var1]], na.rm = TRUE),
               max(data_gibsim_MB_for_beta[[var1]], na.rm = TRUE),
               length.out = 50)
  
  y_seq <- seq(min(data_gibsim_MB_for_beta[[var2]], na.rm = TRUE),
               max(data_gibsim_MB_for_beta[[var2]], na.rm = TRUE),
               length.out = 50)
  
  grid <- expand.grid(x_seq, y_seq)
  colnames(grid) <- c(var1, var2)
  
  base <- data.frame(lapply(data_gibsim_MB_for_beta, function(x) {
    if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1]
  }))
  
  newdata <- base[rep(1, nrow(grid)), ]
  newdata[[var1]] <- grid[[var1]]
  newdata[[var2]] <- grid[[var2]]
  
  # baseline (둘 다 최소)
  baseline <- newdata
  baseline[[var1]] <- min(data_gibsim_MB_for_beta[[var1]], na.rm = TRUE)
  baseline[[var2]] <- min(data_gibsim_MB_for_beta[[var2]], na.rm = TRUE)
  
  pred <- predict(gam_model, newdata = newdata, type = "response")
  pred_base <- predict(gam_model, newdata = baseline, type = "response")
  
  newdata$delta <- (pred - pred_base)*100
  
  print(min(newdata$delta))
  print(max(newdata$delta))
  
  p <- ggplot(newdata, aes_string(x = var1, y = var2, fill = "delta")) +
    # geom_tile() +
    # geom_contour(aes(z = delta), color = "white") +
    geom_tile(data = newdata,
              aes_string(x = var1, y = var2, fill = "delta")) +
    geom_contour(data = newdata,
                 aes_string(x = var1, y = var2, z = "delta"),
                 color = "white") +
    labs(fill = "Δ FHB incidence (%)", 
         x = paste0("\n", var1), 
         y = paste0(var2, "\n")) +
    theme_minimal()+
    scale_fill_gradient2(
      low = "#3B4CC0",
      mid = "white",
      high = "#B40426",
      midpoint = 0,
      limits = c(-100, 100),
      breaks = seq(-100, 100, by = 50),
      oob = scales::squish
    ) +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                                 size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    )
  
  print(p)
  
  ggsave(
    filename = paste0(save_dir, "GAM_delta_plot_", var1, "_", var2, ".png"),
    plot = p,
    width = 5.2,
    height = 4
  )
}
# # 반복문으로 각 smooth term을 그려서 저장
# for (s in smooth_names) {
#   p <- draw(gam_model, select = s, type = "response") +
#     theme_minimal() +
#     theme(
#       text = element_text(size = 15),
#       axis.title = element_text(size = 16),
#       axis.text = element_text(size = 15)
#     )
# 
#   save_name <- gsub("[() ,]", "_", s)  # 파일명에 안전하게 변환
#   ggsave(
#     filename = paste0(save_dir, "GAM_smooth_plot_", save_name, ".png"),
#     plot = p,
#     width = 5,
#     height = 5
#   )
# }
# 
# 
# #####################################################################################
# ############################################################################################
# library(ggplot2)
# 
# vars <- c("Alternaria", "Cladosporium", "Epicoccum",
#           "Hannaella", "Papiliotrema", "Periconia",
#           "fhb_incidence_predicted")
# 
# for (var in vars) {
#   # var <- vars[1]
#   x_seq <- seq(min(data_gibsim_MB_for_beta[[var]], na.rm = TRUE),
#                max(data_gibsim_MB_for_beta[[var]], na.rm = TRUE),
#                length.out = 100)
#   
#   # 나머지 변수는 평균으로 고정
#   newdata <- data.frame(lapply(data_gibsim_MB_for_beta, function(x) {
#     if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1]
#   }))
#   
#   newdata <- newdata[rep(1, 100), ]
#   newdata[[var]] <- x_seq
#   
#   pred <- predict(gam_model, newdata = newdata,
#                   type = "response", se.fit = TRUE)
#   
#   newdata$fit <- pred$fit
#   newdata$se <- pred$se.fit
#   
#   p <- ggplot(newdata, aes_string(x = var, y = "fit")) +
#     geom_line(size = 1.2) +
#     geom_ribbon(aes(ymin = fit - 2*se,
#                     ymax = fit + 2*se),
#                 alpha = 0.2) +
#     labs(y = "FHB incidence (predicted)") +
#     theme_minimal()
#   
#   print(p)
#   save_name <- var
#   ggsave(
#     filename = paste0(save_dir, "GAM_smooth_plot_", save_name, ".png"),
#     plot = p,
#     width = 5,
#     height = 5
#   )
# }
# 
# ###=======
# vars1 <- c("Alternaria", "Alternaria", "Alternaria", "Epicoccum", "Epicoccum")
# vars2 <- c("Cladosporium","Papiliotrema", "Periconia", "Cladosporium", "Periconia")
# 
# for(k in c(1:length(vars1))){
#   var1 <- vars1[k]
#   var2 <- vars2[k]
#     
#   x_seq <- seq(min(data_gibsim_MB_for_beta[[var1]], na.rm = TRUE),
#                max(data_gibsim_MB_for_beta[[var1]], na.rm = TRUE),
#                length.out = 50)
#   
#   y_seq <- seq(min(data_gibsim_MB_for_beta[[var2]], na.rm = TRUE),
#                max(data_gibsim_MB_for_beta[[var2]], na.rm = TRUE),
#                length.out = 50)
#   
#   grid <- expand.grid(x_seq, y_seq)
#   colnames(grid) <- c(var1, var2)
#   
#   # 나머지 변수 평균 고정
#   base <- data.frame(lapply(data_gibsim_MB_for_beta, function(x) {
#     if (is.numeric(x)) mean(x, na.rm = TRUE) else x[1]
#   }))
#   
#   newdata <- base[rep(1, nrow(grid)), ]
#   newdata[[var1]] <- grid[[var1]]
#   newdata[[var2]] <- grid[[var2]]
#   
#   pred <- predict(gam_model, newdata = newdata, type = "response")
#   
#   newdata$fit <- pred
#   
#   p <- ggplot(newdata, aes_string(x = var1, y = var2, fill = "fit")) +
#     geom_tile() +
#     geom_contour(aes(z = fit), color = "white") +
#     scale_fill_viridis_c() +
#     labs(fill = "FHB incidence") +
#     theme_minimal()
#   
#   print(p)
#   save_name <- paste0(var1, "_", var2)
#   ggsave(
#     filename = paste0(save_dir, "GAM_smooth_plot_", save_name, ".png"),
#     plot = p,
#     width = 6,
#     height = 5
#   )
#   
# }
# 
# # ####
# # # 사용할 데이터 (모델에 사용한 데이터)
# # dat <- model.frame(gam_model)
# # library(dplyr)
# # # 고정값 (median 추천)
# # fixed_vals <- dat %>%
# #   summarise(across(everything(), median))
# # 
# # # 변수 선택
# # var <- "Alternaria"
# # 
# # # x grid 생성
# # x_seq <- seq(min(dat[[var]]), max(dat[[var]]), length.out = 100)
# # 
# # # newdata 생성
# # newdata <- fixed_vals[rep(1, 100), ]
# # newdata[[var]] <- x_seq
# # 
# # # baseline (x 최소값)
# # baseline <- newdata
# # baseline[[var]] <- min(dat[[var]])
# # 
# # # 예측 (response scale!)
# # pred <- predict(gam_model, newdata = newdata, type = "response")
# # pred_base <- predict(gam_model, newdata = baseline, type = "response")
# # 
# # # Δy 계산
# # delta_y <- pred - pred_base
# # 
# # # plot
# # df_plot <- data.frame(
# #   x = x_seq,
# #   delta_y = delta_y
# # )
# # 
# # ggplot(df_plot, aes(x = x, y = delta_y)) +
# #   geom_line(size = 1.2) +
# #   geom_hline(yintercept = 0, linetype = "dashed") +
# #   labs(
# #     x = var,
# #     y = "Change in FHB incidence (Δ)",
# #     title = paste("PDP (Δy) -", var)
# #   ) +
# #   theme_classic()
# # 
# # 
# # #
# # var1 <- "Epicoccum"
# # var2 <- "Cladosporium"
# # 
# # x_seq <- seq(min(dat[[var1]]), max(dat[[var1]]), length.out = 50)
# # y_seq <- seq(min(dat[[var2]]), max(dat[[var2]]), length.out = 50)
# # 
# # grid <- expand.grid(
# #   Alternaria = x_seq,
# #   Cladosporium = y_seq
# # )
# # 
# # # 나머지 변수 고정
# # for (v in setdiff(names(dat), c(var1, var2))) {
# #   grid[[v]] <- fixed_vals[[v]]
# # }
# # 
# # # baseline (둘 다 최소)
# # baseline <- grid
# # baseline[[var1]] <- min(dat[[var1]])
# # baseline[[var2]] <- min(dat[[var2]])
# # 
# # # 예측
# # pred <- predict(gam_model, newdata = grid, type = "response")
# # pred_base <- predict(gam_model, newdata = baseline, type = "response")
# # 
# # grid$delta_y <- pred - pred_base
# # 
# # # plot
# # ggplot(grid, aes(x = Alternaria, y = Cladosporium, fill = delta_y)) +
# #   geom_tile() +
# #   scale_fill_gradient2(midpoint = 0) +
# #   labs(
# #     fill = "Δ FHB incidence",
# #     title = paste("2D PDP (Δy):", var1, "×", var2)
# #   ) +
# #   theme_classic()

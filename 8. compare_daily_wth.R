### library --------------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)

### Function -------------------------------------------------------------------
make_linegraph_incidence <- function(data_for_linegraph_join, y_var){
  line_g <- ggplot(data_for_linegraph_join, aes(x = day_num, 
                                                y = get(y_var), 
                                                group = id_year, 
                                                color = incidence,
                                                linetype = GIBSIM_group)) +
    geom_line(alpha = 0.5, linewidth = 1.2) +
    labs(x = "Day Number", 
         y = y_var,
         color = "Obs Incidence",
         linetype = "GIBSIM Group",
         title = paste0("Daily ",y_var," by Weather id and years"))+
    scale_colour_gradient2(high = "red", mid = "yellow", low = "forestgreen", 
                           midpoint = 0.2, limits = c(0,1)) + 
    theme_minimal() #+
  # theme(legend.position = "none")
  
  return(line_g)
}

make_linegraph_GIBSIM <- function(data_for_linegraph_join, y_var){
  line_g2 <- ggplot(data_for_linegraph_join, aes(x = day_num, 
                                                 y =get(y_var), 
                                                 group = id_year, 
                                                 color = GIBSIM,
                                                 linetype = GIBSIM_group)) + #pred value
    geom_line(alpha = 0.5, linewidth = 1.2) +
    labs(x = "Day Number", 
         y = y_var,
         color = "GIBSIM Result",
         title = paste0("Daily ",y_var," by Weather id and years"))+
    scale_colour_gradient2(high = "red", mid = "yellow", low = "forestgreen", midpoint = 20, limits = c(0,100)) + 
    theme_minimal()# +
  theme(legend.position = "none")
  
  return(line_g2)
}

make_linegraph_compare <- function(data_for_g, y_var){
  # data_for_g = over_compare_df
  # y_var = "tavg"
  line_g2 <- ggplot(data_for_g, aes(x = day_num, 
                                    y =get(y_var), 
                                    group = id_year, 
                                    color = GIBSIM_group)) + #pred value
    geom_line(alpha = 0.5, linewidth = 1.2) +
    labs(x = "Day Number", 
         y = y_var,
         color = "GIBSIM_group",
         title = paste0("Daily ",y_var," by Weather id and years"))+
    # scale_colour_gradient2(high = "red", mid = "yellow", low = "forestgreen", midpoint = 20, limits = c(0,100)) + 
    theme_minimal()
  
  return(line_g2)
}

make_linegraph_compare_vs_inc <- function(data_for_g, y_var, save_path,prefix){
  # data_for_g = over_heading_df_for_bad_G4_vs_inc
  # y_var = "tavg"
  # bad = "bad_G4"
  
  line_g2 <- ggplot(data_for_g, aes(x = day_num, 
                                    y = get(y_var), 
                                    group = id_year, 
                                    color = GIBSIM_group)) + 
    geom_line(alpha = 0.5, linewidth = 0.8) +
    labs(x = "Day Number", 
         y = y_var,
         color = "GIBSIM_group",
         title = paste0("Daily ", y_var, " by Weather id and years")) +
    scale_color_manual(values = c(
      "bad_G4" = "#A0A0A0",
      "low_inc" = "#90C590",
      "mid_inc" = "#FFD27F",
      "high_inc" = "#FF7F7F"
    )) + 
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    ) + 
    geom_vline(xintercept= 15, linetype='dashed', color='#8080C0', size = 0.75, alpha = 0.5)+ 
    geom_vline(xintercept= 25, linetype='dashed', color='#8080C0', size = 0.75, alpha = 0.5) + 
    geom_vline(xintercept= 35, linetype='dashed', color='#8080C0', size = 0.75, alpha = 0.5)
  
  ggsave(
    plot = line_g2,
    file = paste0(save_path, "/", prefix, "_", y_var, ".png"),
    width = 17,
    height = 10,
    units = "cm"
  )
  
  return(line_g2)
}
make_kw_boxplot_for_bad_G4 <- function(data, y_var, save_path, prefix){
  # data = over_flowering_df_0_to_15
  # y_var = v
  # save_path = save_dir
  # prefix = "box_over_bad_G4_0_15"
  
  df <- data %>% dplyr::select(GIBSIM_group, all_of(y_var)) %>% na.omit()
  
  # Kruskal-Wallis
  kw_res <- kruskal_test(df, formula(paste(y_var, "~ GIBSIM_group")))
  
  # Dunn test
  dunn_res <- df %>%
    dunn_test(formula(paste(y_var, "~ GIBSIM_group")), p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "GIBSIM_group")
  
  dunn_res_sig <- dunn_res %>% filter(p.adj < 0.05)
  
  # plot
  p <- ggplot(df, aes(x = GIBSIM_group, y = .data[[y_var]], fill = GIBSIM_group)) +
    geom_boxplot() +
    geom_jitter(color = "black", size = 1, alpha = 0.1) +
    stat_pvalue_manual(
      dunn_res_sig,
      label = "p.adj.signif",
      tip.length = 0.01,
      inherit.aes = FALSE
    ) +
    # stat_compare_means(
    #   method = "kruskal.test"
    #   # label.y = max(df[[y_var]], na.rm = TRUE) * 1.2
    # ) +
    labs(
      # title = y_var,
      x = "",
      y = paste0(y_var, "\n")
    ) +
    scale_fill_manual(values = c(
      "bad_G4" = "#A0A0A0",
      "low_inc" = "#90C590",
      "mid_inc" = "#FFD27F",
      "high_inc" = "#FF7F7F"
    )) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.text.x = element_text(size = 12, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 12, face = 'bold', color = 'black')
    ) + theme(legend.position = "none")
  
  ggsave(
    plot = p,
    file = paste0(save_path, "/", prefix, "_", y_var, ".png"),
    width = 10,
    height = 9,
    units = "cm"
  )
  
  return(p)
}

make_kw_boxplot_for_bad_G1_for_m <- function(data, target_microbe, save_path, prefix){
  # 그룹 순서 지정
  data$GIBSIM_group <- factor(data$GIBSIM_group, 
                              levels = c("bad_G1","low_inc", "mid_inc", "high_inc"))
  
  # 분석용 데이터 준비
  df <- data %>% dplyr::select(GIBSIM_group, all_of(target_microbe)) %>% na.omit()
  
  # Kruskal-Wallis test
  kw_res <- kruskal.test(df[[target_microbe]] ~ df$GIBSIM_group)
  print(kw_res)
  
  # Dunn test 실행
  dunn_res <- dunn.test::dunn.test(df[[target_microbe]], df$GIBSIM_group, method = "bonferroni")
  
  # Dunn test 결과를 tibble로 변환
  dunn_df <- tibble::tibble(
    group1 = sapply(strsplit(dunn_res$comparisons, " - "), `[`, 1),
    group2 = sapply(strsplit(dunn_res$comparisons, " - "), `[`, 2),
    p.adj = dunn_res$P.adjusted,
    p.adj.signif = dplyr::case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    # y.position을 그룹별 최대값 위에 조금 띄워서 표시
    y.position = seq(from = max(df[[target_microbe]], na.rm = TRUE) * 1.05,
                     by = 0.3,
                     length.out = length(dunn_res$comparisons))
  )
  dunn_df_sig <- dunn_df %>% dplyr::filter(p.adj.signif != "ns")
  # plot
  p <- ggplot(df, aes(x = GIBSIM_group, y = .data[[target_microbe]], fill = GIBSIM_group)) +
    geom_boxplot() +
    geom_jitter(color = "black", size = 1, alpha = 0.1) +
    stat_pvalue_manual(
      dunn_df_sig,
      label = "p.adj.signif",
      tip.length = 0.01,
      inherit.aes = FALSE
    ) +
    stat_compare_means(
      method = "kruskal.test"
      # label.y = max(df[[target_microbe]], na.rm = TRUE) * 1.2
    ) +
    labs(
      title = target_microbe,
      x = "\nGIBSIM_group",
      y = paste0(target_microbe, "\n")
    ) +
    scale_fill_manual(values = c(
      "bad_G1" = "#A0A0A0",
      "low_inc" = "#90C590",
      "mid_inc" = "#FFD27F",
      "high_inc" = "#FF7F7F"
    )) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, face = 'bold'),
      axis.title.y = element_text(size = 15, face = 'bold'),
      axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    )
  
  # 저장
  ggsave(
    plot = p,
    file = paste0(save_path, "/", prefix, "_", target_microbe, ".png"),
    width = 17,
    height = 15,
    units = "cm"
  )
  
  return(p)
}


make_kw_boxplot_for_bad_G4_for_m <- function(data, target_microbe, save_path, prefix){
  # 그룹 순서 지정
  data$GIBSIM_group <- factor(data$GIBSIM_group, 
                              levels = c("bad_G4","low_inc", "mid_inc", "high_inc"))
  
  # 분석용 데이터 준비
  df <- data %>% dplyr::select(GIBSIM_group, all_of(target_microbe)) %>% na.omit()
  
  # Kruskal-Wallis test
  kw_res <- kruskal.test(df[[target_microbe]] ~ df$GIBSIM_group)
  print(kw_res)
  
  # Dunn test 실행
  dunn_res <- dunn.test::dunn.test(df[[target_microbe]], df$GIBSIM_group, method = "bonferroni")
  
  # Dunn test 결과를 tibble로 변환
  dunn_df <- tibble::tibble(
    group1 = sapply(strsplit(dunn_res$comparisons, " - "), `[`, 1),
    group2 = sapply(strsplit(dunn_res$comparisons, " - "), `[`, 2),
    p.adj = dunn_res$P.adjusted,
    p.adj.signif = dplyr::case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    # y.position을 그룹별 최대값 위에 조금 띄워서 표시
    y.position = seq(from = max(df[[target_microbe]], na.rm = TRUE) * 1.05,
                     by = 0.3,
                     length.out = length(dunn_res$comparisons))
  )
  dunn_df_sig <- dunn_df %>% dplyr::filter(p.adj.signif != "ns")
  # plot
  p <- ggplot(df, aes(x = GIBSIM_group, y = .data[[target_microbe]], fill = GIBSIM_group)) +
    geom_boxplot() +
    geom_jitter(color = "black", size = 1, alpha = 0.1) +
    stat_pvalue_manual(
      dunn_df_sig,
      label = "p.adj.signif",
      tip.length = 0.01,
      inherit.aes = FALSE
    ) +
    stat_compare_means(
      method = "kruskal.test",
      label.y = max(df[[target_microbe]], na.rm = TRUE) * 1.2
    ) +
    labs(
      title = target_microbe,
      x = "\nGIBSIM_group",
      y = paste0(target_microbe, "\n")
    ) +
    scale_fill_manual(values = c(
      "bad_G4" = "#A0A0A0",
      "low_inc" = "#90C590",
      "mid_inc" = "#FFD27F",
      "high_inc" = "#FF7F7F"
    )) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, face = 'bold'),
      axis.title.y = element_text(size = 15, face = 'bold'),
      axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    )
  
  # 저장
  ggsave(
    plot = p,
    file = paste0(save_path, "/", prefix, "_", target_microbe, ".png"),
    width = 17,
    height = 15,
    units = "cm"
  )
  
  return(p)
}



make_linegraph_compare_vs_inc_bad_G1 <- function(data_for_g, y_var, save_path,prefix){
  # data_for_g = over_heading_df_for_bad_G1_vs_inc
  # y_var = "tavg"
  # bad = "bad_G4"
  
  line_g2 <- ggplot(data_for_g, aes(x = day_num, 
                                    y = get(y_var), 
                                    group = id_year, 
                                    color = GIBSIM_group)) + 
    geom_line(alpha = 0.5, linewidth = 0.8) +
    labs(x = "Day Number", 
         y = y_var,
         color = "GIBSIM_group",
         title = paste0("Daily ", y_var, " by Weather id and years")) +
    scale_color_manual(values = c(
      "bad_G1" = "#A0A0A0",
      "low_inc" = "#90C590",
      "mid_inc" = "#FFD27F",
      "high_inc" = "#FF7F7F"
    )) + 
    theme_classic() + 
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
      axis.text.x = element_text(size = 15, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
    ) + 
    geom_vline(xintercept= 15, linetype='dashed', color='#8080C0', size = 0.75, alpha = 0.5)+ 
    geom_vline(xintercept= 25, linetype='dashed', color='#8080C0', size = 0.75, alpha = 0.5) + 
    geom_vline(xintercept= 35, linetype='dashed', color='#8080C0', size = 0.75, alpha = 0.5)
  
  ggsave(
    plot = line_g2,
    file = paste0(save_path, "/", prefix, "_", y_var, ".png"),
    width = 17,
    height = 10,
    units = "cm"
  )
  
  return(line_g2)
}

make_kw_boxplot_for_bad_G1 <- function(data, y_var, save_path, prefix){
  df <- data %>% dplyr::select(GIBSIM_group, all_of(y_var)) %>% na.omit()
  
  # Kruskal-Wallis
  kw_res <- kruskal_test(df, formula(paste(y_var, "~ GIBSIM_group")))
  
  # Dunn test
  dunn_res <- df %>%
    dunn_test(formula(paste(y_var, "~ GIBSIM_group")), p.adjust.method = "bonferroni") %>%
    add_xy_position(x = "GIBSIM_group")
  
  dunn_res_sig <- dunn_res %>% filter(p.adj < 0.05)
  
  # plot
  p <- ggplot(df, aes(x = GIBSIM_group, y = .data[[y_var]], fill = GIBSIM_group)) +
    geom_boxplot() +
    geom_jitter(color = "black", size = 1, alpha = 0.1) +
    stat_pvalue_manual(
      dunn_res_sig,
      label = "p.adj.signif",
      tip.length = 0.01,
      inherit.aes = FALSE
    ) +
    # stat_compare_means(
    #   method = "kruskal.test"
    #   # label.y = max(df[[y_var]], na.rm = TRUE) * 1.2
    # ) +
    labs(
      # title = y_var,
      x = "",
      y = paste0(y_var, "\n")
    ) +
    scale_fill_manual(values = c(
      "bad_G1" = "#A0A0A0",
      "low_inc" = "#90C590",
      "mid_inc" = "#FFD27F",
      "high_inc" = "#FF7F7F"
    )) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.text.x = element_text(size = 12, face = 'bold', color = 'black'),
      axis.text.y = element_text(size = 12, face = 'bold', color = 'black')
    ) + 
    theme(legend.position = "none")
  
  ggsave(
    plot = p,
    file = paste0(save_path, "/", prefix, "_", y_var, ".png"),
    width = 10,
    height = 9,
    units = "cm"
  )
  
  return(p)
}
### data dir--------------------------------------------------------------------
GIBSIM_daily_res_dir <- "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/1.Run_GIBSIM/GIBSIM_daily_results/MB_v3.3/"
GIBSIM_daily_res_filelist <- list.files(GIBSIM_daily_res_dir)
df_for_more_analysis <- read_xlsx("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v12/df_for_more_analysis_good_20_bad_100.xlsx")
unique(df_for_more_analysis$flower_sampling_date)
df_for_more_analysis <- df_for_more_analysis %>%
  mutate(GIBSIM_group = case_when(
    Check_GIBSIM_res == "GIBSIM_bad"  & GIBSIM < 20 & obs >= 20 ~ "bad_G1",
    Check_GIBSIM_res == "GIBSIM_bad"  & GIBSIM >= 20 & obs >= 20 ~ "bad_G2",
    Check_GIBSIM_res == "GIBSIM_bad"  & GIBSIM < 20 & obs < 20 ~ "bad_G3",
    Check_GIBSIM_res == "GIBSIM_bad"  & GIBSIM >= 20 & obs < 20 ~ "bad_G4",
    Check_GIBSIM_res == "GIBSIM_good" & GIBSIM < 20 & obs >= 20 ~ "good_G1",
    Check_GIBSIM_res == "GIBSIM_good" & GIBSIM >= 20 & obs >= 20 ~ "good_G2",
    Check_GIBSIM_res == "GIBSIM_good" & GIBSIM < 20 & obs < 20 ~ "good_G3",
    Check_GIBSIM_res == "GIBSIM_good" & GIBSIM >= 20 & obs < 20 ~ "good_G4",
    TRUE ~ "error"
  ))

table(df_for_more_analysis$GIBSIM_group)

Check_GIBSIM_res_df <- df_for_more_analysis %>% 
  dplyr::select(wth_ID, year, Check_GIBSIM_res, incidence, GIBSIM, GIBSIM_group, 
                Alternaria, Epicoccum, Hannaella, Periconia) %>% 
  dplyr::mutate(id_year = paste0(wth_ID, "_", year))
duplicated_id_year <- unique(Check_GIBSIM_res_df[duplicated(Check_GIBSIM_res_df$id_year),]$id_year)

Check_GIBSIM_res_df2 <- Check_GIBSIM_res_df[!duplicated(Check_GIBSIM_res_df$id_year),]

# Check_GIBSIM_res_df2 <- Check_GIBSIM_res_df
# Check_GIBSIM_res_df2[Check_GIBSIM_res_df2$wth_ID == "wth_057",]
# 
# Check_GIBSIM_res_df[duplicated(Check_GIBSIM_res_df$id_year),]
# 
# Check_GIBSIM_res_df %>% filter(id_year == duplicated_id_year[1])
# Check_GIBSIM_res_df %>% filter(id_year == duplicated_id_year[2])
# Check_GIBSIM_res_df %>% filter(id_year == duplicated_id_year[3])
# Check_GIBSIM_res_df %>% filter(id_year == duplicated_id_year[4])
# Check_GIBSIM_res_df %>% filter(id_year == duplicated_id_year[5])
names(df_for_more_analysis)
heading # 4/16

hour_dir <- "P:/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/weather_data/KMA_500m_grid_2022_2025/hourly_wth"
### merge all GIBSIM_daily data ------------------------------------------------
data_for_linegraph <- data.frame()
for(i in c(1:length(GIBSIM_daily_res_filelist))){
  #i=1
  imsi_GIBSIM_daily_res_xlsx <- read_xlsx(paste0(GIBSIM_daily_res_dir, GIBSIM_daily_res_filelist[i]))
  
  ##
  imsi_target_wth <- stringr::str_sub(GIBSIM_daily_res_filelist[i], 18,24) 
  imsi_target_year <- stringr::str_sub(GIBSIM_daily_res_filelist[i], 26,29) 
  
  imsi_hour_raw_df <- read_xlsx(paste0(hour_dir,"/", imsi_target_wth, ".xlsx"))
  imsi_hour_raw_df$Year <- as.character(stringr::str_sub(imsi_hour_raw_df$data_dt_str, 1,4)) 
  imsi_hour_raw_df$Mon <- as.character(stringr::str_sub(imsi_hour_raw_df$data_dt_str, 5,6)) 
  imsi_hour_raw_df$Day <- as.character(stringr::str_sub(imsi_hour_raw_df$data_dt_str, 7,8)) 
  imsi_hour_raw_df$Date <- as.Date(paste0(imsi_hour_raw_df$Year, "-", imsi_hour_raw_df$Mon, "-", imsi_hour_raw_df$Day)) 
  imsi_hour_raw_df$RH80 <- ifelse(imsi_hour_raw_df$hm > 80, 1, 0)
  
  imsi_hour_to_daily <- imsi_hour_raw_df %>% group_by(Date, Year, Mon, Day) %>% summarize(RH80 = sum(RH80))
  imsi_hour_to_daily_filtered <- imsi_hour_to_daily %>% dplyr::filter(Year == imsi_target_year) %>%
    dplyr::filter(Mon == "04" | Mon == "05")
  ##
  
  Apr_to_May_df <- imsi_GIBSIM_daily_res_xlsx %>% 
    dplyr::filter(month == 4 | month == 5) %>% 
    dplyr::select(id, year, date, doy, heading_date, disease_sampling_date, 
                  tavg, tmin, tmax, prcp, rhum, rsds, 
                  tavg48, tavg72, rhum48, rhum72, rainy_day, gz_mean, inf)
  Apr_to_May_df <- Apr_to_May_df %>%
    dplyr::mutate(id_year = paste0(id, "_", year), 
                  day_num = c(1:nrow(Apr_to_May_df))) %>%
    dplyr::relocate(id_year, .after = year) %>%
    dplyr::relocate(day_num, .before = doy)
  
  # v3
  Apr_to_May_df$RH80 <- imsi_hour_to_daily_filtered$RH80
  Apr_to_May_df <- Apr_to_May_df %>% dplyr::mutate(Teq = tavg/35, 
                                                   Teq_Rossi = (tavg-5)/30) %>% 
    dplyr::mutate(Ino_Manstretta = (5.317*(Teq^1.501)*(1-Teq))^4.983,  #Manstretta, V., & Rossi, V. (2016). Ascospore discharge by Fusarium graminearum as affected by temperature and relative humidity. European journal of plant pathology, 146(1), 191-197.
                  Ino_Rossi = (25.98*(Teq_Rossi^8.59)*(1-Teq_Rossi))^0.24/(1+exp(5.52-0.51*0.80)), #Rossi, V., Pattori, E., Ravanetti, A., & Giosue, S. (2002). Effect of constant and fluctuating temperature regimes on sporulation of four fungi causing head blight of wheat. Journal of Plant Pathology, 95-105.
                  Ino_Dufault = 2.34+0.0732*tavg-0.0192*tavg^2-0.000740*tavg^3, #Dufault, N. S., De Wolf, E. D., Lipps, P. E., & Madden, L. V. (2006). Role of temperature and moisture in the production and maturation of Gibberella zeae perithecia. Plant Disease, 90(5), 637-644.
                  Ino_Rossi2 = log10(ifelse(-682.3+45.68*tavg+21.5*RH80 + 107.0*prcp < 0, 0, -682.3+45.68*tavg+21.5*RH80 + 107.0*prcp)+1), #Rossi, V., Languasco, L., Pattori, E., & Giosue, S. (2002). Dynamics of airborne Fusarium macroconidia in wheat fields naturally affected by head blight. Journal of Plant Pathology, 53-64.
                  Ino_GIBSIM = gz_mean, 
    )
  
  
  data_for_linegraph <- rbind(data_for_linegraph, Apr_to_May_df)
  print(i)
}


# d_before_hd_sampling = date %in% seq(flower_sampling_date-24,flower_sampling_date-15,1), 
# d_hd_sampling = date %in% seq(flower_sampling_date-14,flower_sampling_date-5,1), 
# d_flower_sampling = date %in% seq(flower_sampling_date-4,flower_sampling_date+5,1), 
# d_after_flower_sampling = date %in% seq(flower_sampling_date+6,flower_sampling_date+15,1) 


# data_for_linegraph_start_DVR_heading <- data.frame()
# for(i in c(1:length(GIBSIM_daily_res_filelist))){
#   #i = 1
#   imsi_GIBSIM_daily_res_xlsx <- read_xlsx(paste0(GIBSIM_daily_res_dir, GIBSIM_daily_res_filelist[i]))
#   imsi_heading_date <- unique(imsi_GIBSIM_daily_res_xlsx$heading_date)
#   imsi_disease_date <- unique(imsi_GIBSIM_daily_res_xlsx$disease_sampling_date)
#   imsi_flowering_date <- unique(imsi_GIBSIM_daily_res_xlsx$flower_sampling_date)
#   
#   heading_raw <- which(imsi_GIBSIM_daily_res_xlsx$date == imsi_heading_date)
#   flowering_raw <- which(imsi_GIBSIM_daily_res_xlsx$date == imsi_flowering_date)
#   disease_raw <- which(imsi_GIBSIM_daily_res_xlsx$date == imsi_disease_date)
#   
#   # heading_minus_20_raw <- heading_raw-20
#   # heading_plus_40_raw <- heading_raw+40
#   flowering_minus_24_raw <- flowering_raw-24
#   
#   imsi_GIBSIM_daily_res_xlsx_filtered_by_heading <- imsi_GIBSIM_daily_res_xlsx[c(flowering_minus_24_raw:disease_raw),]
#   
#   ##
#   imsi_target_wth_heading <- stringr::str_sub(GIBSIM_daily_res_filelist[i], 18,24) 
#   imsi_target_year_heading <- stringr::str_sub(GIBSIM_daily_res_filelist[i], 26,29) 
#   
#   imsi_hour_raw_df_heading <- read_xlsx(paste0(hour_dir,"/", imsi_target_wth, ".xlsx"))
#   imsi_hour_raw_df_heading$Year <- as.character(stringr::str_sub(imsi_hour_raw_df_heading$data_dt_str, 1,4)) 
#   imsi_hour_raw_df_heading$Mon <- as.character(stringr::str_sub(imsi_hour_raw_df_heading$data_dt_str, 5,6)) 
#   imsi_hour_raw_df_heading$Day <- as.character(stringr::str_sub(imsi_hour_raw_df_heading$data_dt_str, 7,8)) 
#   imsi_hour_raw_df_heading$Date <- as.Date(paste0(imsi_hour_raw_df_heading$Year, "-", imsi_hour_raw_df_heading$Mon, "-", imsi_hour_raw_df_heading$Day)) 
#   imsi_hour_raw_df_heading$RH80 <- ifelse(imsi_hour_raw_df_heading$hm > 80, 1, 0)
#   
#   imsi_hour_to_daily_heading <- imsi_hour_raw_df_heading %>% dplyr::group_by(Date, Year, Mon, Day) %>% dplyr::summarize(RH80 = sum(RH80))
#   
#   start_day <- imsi_GIBSIM_daily_res_xlsx_filtered_by_heading$date[1]
#   end_day <- imsi_GIBSIM_daily_res_xlsx_filtered_by_heading$date[nrow(imsi_GIBSIM_daily_res_xlsx_filtered_by_heading)]
#   start_day_rawnum <- which(imsi_hour_to_daily_heading$Date == start_day)
#   end_day_rawnum <- which(imsi_hour_to_daily_heading$Date == end_day)
#   
#   imsi_hour_to_daily_filtered_heading <- imsi_hour_to_daily_heading[c(start_day_rawnum:end_day_rawnum),]
#   ##
#   before_heading_to_disease_df <- imsi_GIBSIM_daily_res_xlsx_filtered_by_heading %>% 
#     dplyr::select(id, year, date, doy, heading_date, disease_sampling_date, 
#                   tavg, tmin, tmax, prcp, rhum, rsds, 
#                   tavg48, tavg72, rhum48, rhum72, rainy_day, gz_mean, inf) %>%
#     mutate(id_year = paste0(id, "_", year), 
#            day_num = c(1:nrow(imsi_GIBSIM_daily_res_xlsx_filtered_by_heading))) %>%
#     relocate(id_year, .after = year) %>%
#     relocate(day_num, .before = doy) 
#   
#   before_heading_to_disease_df$RH80 <- imsi_hour_to_daily_filtered_heading$RH80
#   
#   before_heading_to_disease_df <- before_heading_to_disease_df %>% dplyr::mutate(Teq = tavg/35, 
#                                                  Teq_Rossi = (tavg-5)/30) %>% 
#     dplyr::mutate(Ino_Manstretta = ifelse(tavg < 0, 0, ((5.317*(Teq^1.501)*(1-Teq))^4.983)) ,  #Manstretta, V., & Rossi, V. (2016). Ascospore discharge by Fusarium graminearum as affected by temperature and relative humidity. European journal of plant pathology, 146(1), 191-197.
#                   Ino_Rossi = ifelse(tavg < 5, 0, (25.98*(Teq_Rossi^8.59)*(1-Teq_Rossi))^0.24/(1+exp(5.52-0.51*0.80))) , #Rossi, V., Pattori, E., Ravanetti, A., & Giosue, S. (2002). Effect of constant and fluctuating temperature regimes on sporulation of four fungi causing head blight of wheat. Journal of Plant Pathology, 95-105.
#                   Ino_Dufault = 2.34+0.0732*tavg-0.0192*tavg^2-0.000740*tavg^3, #Dufault, N. S., De Wolf, E. D., Lipps, P. E., & Madden, L. V. (2006). Role of temperature and moisture in the production and maturation of Gibberella zeae perithecia. Plant Disease, 90(5), 637-644.
#                   Ino_Rossi2 = log10(ifelse(-682.3+45.68*tavg+21.5*RH80 + 107.0*prcp < 0, 0, -682.3+45.68*tavg+21.5*RH80 + 107.0*prcp)+1), #Rossi, V., Languasco, L., Pattori, E., & Giosue, S. (2002). Dynamics of airborne Fusarium macroconidia in wheat fields naturally affected by head blight. Journal of Plant Pathology, 53-64.
#                   Ino_GIBSIM = gz_mean, 
#     )
#   
#   data_for_linegraph_start_DVR_heading <- rbind(data_for_linegraph_start_DVR_heading, before_heading_to_disease_df)
#   print(i)
# }

data_for_linegraph_start_DVR_flowering <- data.frame()
for(i in c(1:length(GIBSIM_daily_res_filelist))){
  #i = 1
  imsi_GIBSIM_daily_res_xlsx <- read_xlsx(paste0(GIBSIM_daily_res_dir, GIBSIM_daily_res_filelist[i]))
  imsi_disease_date <- unique(imsi_GIBSIM_daily_res_xlsx$disease_sampling_date)
  imsi_flowering_date <- unique(imsi_GIBSIM_daily_res_xlsx$flower_sampling_date)
  
  flowering_raw <- which(imsi_GIBSIM_daily_res_xlsx$date == imsi_flowering_date)
  disease_raw <- which(imsi_GIBSIM_daily_res_xlsx$date == imsi_disease_date)
  
  # flowering_minus_20_raw <- flowering_raw-20
  # flowering_plus_40_raw <- flowering_raw+40
  flowering_minus_24_raw <- flowering_raw-24
  
  imsi_GIBSIM_daily_res_xlsx_filtered_by_flowering <- imsi_GIBSIM_daily_res_xlsx[c(flowering_minus_24_raw:disease_raw),]
  
  ##
  imsi_target_wth_flowering <- stringr::str_sub(GIBSIM_daily_res_filelist[i], 18,24) 
  imsi_target_year_flowering <- stringr::str_sub(GIBSIM_daily_res_filelist[i], 26,29) 
  
  imsi_hour_raw_df_flowering <- read_xlsx(paste0(hour_dir,"/", imsi_target_wth, ".xlsx"))
  imsi_hour_raw_df_flowering$Year <- as.character(stringr::str_sub(imsi_hour_raw_df_flowering$data_dt_str, 1,4)) 
  imsi_hour_raw_df_flowering$Mon <- as.character(stringr::str_sub(imsi_hour_raw_df_flowering$data_dt_str, 5,6)) 
  imsi_hour_raw_df_flowering$Day <- as.character(stringr::str_sub(imsi_hour_raw_df_flowering$data_dt_str, 7,8)) 
  imsi_hour_raw_df_flowering$Date <- as.Date(paste0(imsi_hour_raw_df_flowering$Year, "-", imsi_hour_raw_df_flowering$Mon, "-", imsi_hour_raw_df_flowering$Day)) 
  imsi_hour_raw_df_flowering$RH80 <- ifelse(imsi_hour_raw_df_flowering$hm > 80, 1, 0)
  
  imsi_hour_to_daily_flowering <- imsi_hour_raw_df_flowering %>% dplyr::group_by(Date, Year, Mon, Day) %>% dplyr::summarize(RH80 = sum(RH80))
  
  start_day <- imsi_GIBSIM_daily_res_xlsx_filtered_by_flowering$date[1]
  end_day <- imsi_GIBSIM_daily_res_xlsx_filtered_by_flowering$date[nrow(imsi_GIBSIM_daily_res_xlsx_filtered_by_flowering)]
  start_day_rawnum <- which(imsi_hour_to_daily_flowering$Date == start_day)
  end_day_rawnum <- which(imsi_hour_to_daily_flowering$Date == end_day)
  
  imsi_hour_to_daily_filtered_flowering <- imsi_hour_to_daily_flowering[c(start_day_rawnum:end_day_rawnum),]
  ##
  before_flowering_to_disease_df <- imsi_GIBSIM_daily_res_xlsx_filtered_by_flowering %>% 
    dplyr::select(id, year, date, doy, flower_sampling_date, disease_sampling_date, 
                  tavg, tmin, tmax, prcp, rhum, rsds, 
                  tavg48, tavg72, rhum48, rhum72, rainy_day, gz_mean, inf) %>%
    mutate(id_year = paste0(id, "_", year), 
           day_num = c(1:nrow(imsi_GIBSIM_daily_res_xlsx_filtered_by_flowering))) %>%
    relocate(id_year, .after = year) %>%
    relocate(day_num, .before = doy) 
  
  before_flowering_to_disease_df$RH80 <- imsi_hour_to_daily_filtered_flowering$RH80
  #Matengu, T. T., Bullock, P. R., Mkhabela, M. S., Zvomuya, F., Henriquez, M. A., Ojo, E. R., & Fernando, W. D. (2024). Weather???based models for forecasting Fusarium head blight risks in wheat and barley: A review. Plant Pathology, 73(3), 492-505.
  before_flowering_to_disease_df <- before_flowering_to_disease_df %>% dplyr::mutate(Teq = tavg/35, 
                                                                                     Teq_Rossi = (tavg-5)/30) %>% 
    dplyr::mutate(Ino_Manstretta = ifelse(tavg < 0, 0, ((5.317*(Teq^1.501)*(1-Teq))^4.983)) ,  #Manstretta, V., & Rossi, V. (2016). Ascospore discharge by Fusarium graminearum as affected by temperature and relative humidity. European journal of plant pathology, 146(1), 191-197.
                  Ino_Rossi = ifelse(tavg < 5, 0, (25.98*(Teq_Rossi^8.59)*(1-Teq_Rossi))^0.24/(1+exp(5.52-0.51*0.80))) , #Rossi, V., Pattori, E., Ravanetti, A., & Giosue, S. (2002). Effect of constant and fluctuating temperature regimes on sporulation of four fungi causing head blight of wheat. Journal of Plant Pathology, 95-105.
                  # Ino_Dufault = 2.34+0.0732*tavg-0.0192*tavg^2-0.000740*tavg^3, #Dufault, N. S., De Wolf, E. D., Lipps, P. E., & Madden, L. V. (2006). Role of temperature and moisture in the production and maturation of Gibberella zeae perithecia. Plant Disease, 90(5), 637-644. 
                  # Ino_Rossi2 = log10(ifelse(-682.3+45.68*tavg+21.5*RH80 + 107.0*prcp < 0, 0, -682.3+45.68*tavg+21.5*RH80 + 107.0*prcp)+1), #Rossi, V., Languasco, L., Pattori, E., & Giosue, S. (2002). Dynamics of airborne Fusarium macroconidia in wheat fields naturally affected by head blight. Journal of Plant Pathology, 53-64.
                  Ino_GIBSIM = gz_mean, 
    )
  
  data_for_linegraph_start_DVR_flowering <- rbind(data_for_linegraph_start_DVR_flowering, before_flowering_to_disease_df)
  print(i)
}

nrow(data_for_linegraph)
nrow(data_for_linegraph_start_DVR_flowering)

### merge Check_GIBSIM_res_df2 & data_for_linegraph
# Check_GIBSIM_res_df3 <- Check_GIBSIM_res_df2 %>% 
#   dplyr::select(id_year, incidence, Check_GIBSIM_res, GIBSIM, GIBSIM_group, 
#                 Alternaria, Epicoccum, Hannaella, Periconia)
# 
# data_for_linegraph_join <- left_join(data_for_linegraph, Check_GIBSIM_res_df3, by = "id_year")
# data_for_linegraph_flowering_join <- left_join(data_for_linegraph_start_DVR_flowering, Check_GIBSIM_res_df3, by = "id_year")


Check_GIBSIM_res_df4 <- Check_GIBSIM_res_df %>% 
  dplyr::select(id_year, incidence, Check_GIBSIM_res, GIBSIM, GIBSIM_group, 
                Alternaria, Epicoccum, Hannaella, Periconia)
data_for_linegraph_join <- left_join(Check_GIBSIM_res_df4,data_for_linegraph, by = "id_year")
data_for_linegraph_flowering_join <- left_join(Check_GIBSIM_res_df4,data_for_linegraph_start_DVR_flowering,  by = "id_year")



### graph
library(ggplot2)
# NA_df <- data_for_linegraph_join[which(is.na(data_for_linegraph_join$incidence) == T),]
# NA_df$id_year
# 
# sum(is.na(data_for_linegraph_join$incidence))
# hist(data_for_linegraph_join$incidence)

### prepare the dataset for line graph------------------------------------------
G1_good_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "good_G1")
G2_good_df_for_over <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "good_G2") %>% dplyr::filter(GIBSIM < 50)
G2_good_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "good_G2")
G3_good_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "good_G3")
G4_good_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "good_G4")

G1_bad_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "bad_G1")
G2_bad_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "bad_G2")
G3_bad_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "bad_G3")
G4_bad_df <- data_for_linegraph_join %>% dplyr::filter(GIBSIM_group == "bad_G4")


G1_good_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "good_G1")
G2_good_flowering_df_for_over <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "good_G2") %>% dplyr::filter(GIBSIM < 50)
G2_good_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "good_G2")
G3_good_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "good_G3")
G4_good_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "good_G4")

G1_bad_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "bad_G1")
G2_bad_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "bad_G2")
G3_bad_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "bad_G3")
G4_bad_flowering_df <- data_for_linegraph_flowering_join %>% dplyr::filter(GIBSIM_group == "bad_G4")



low_inc_df_for_badG4 <- data_for_linegraph_flowering_join %>% dplyr::filter(incidence <= 0.2) %>% dplyr::filter(GIBSIM_group != "bad_G4")
mid_inc_df_for_badG4 <- data_for_linegraph_flowering_join %>% dplyr::filter(incidence > 0.2) %>% dplyr::filter(incidence <= 0.5)
high_inc_df_for_badG4 <- data_for_linegraph_flowering_join %>% dplyr::filter(incidence > 0.5)
low_inc_df_for_badG4$GIBSIM_group <- "low_inc"
mid_inc_df_for_badG4$GIBSIM_group <- "mid_inc"
high_inc_df_for_badG4$GIBSIM_group <- "high_inc"

low_inc_df_for_badG1 <- data_for_linegraph_flowering_join %>% dplyr::filter(incidence <= 0.2)
mid_inc_df_for_badG1 <- data_for_linegraph_flowering_join %>% dplyr::filter(incidence > 0.2) %>% dplyr::filter(incidence <= 0.5)  %>% dplyr::filter(GIBSIM_group != "bad_G1") 
high_inc_df_for_badG1 <- data_for_linegraph_flowering_join %>% dplyr::filter(incidence > 0.5)  %>% dplyr::filter(incidence != "bad_G1")
low_inc_df_for_badG1$GIBSIM_group <- "low_inc"
mid_inc_df_for_badG1$GIBSIM_group <- "mid_inc"
high_inc_df_for_badG1$GIBSIM_group <- "high_inc"


### make subset -----------------------------------------------------------------
under_flowering_df_for_bad_G1 <- rbind(G1_bad_flowering_df, G2_good_flowering_df, G3_good_flowering_df)
over_flowering_df_for_bad_G4 <- rbind(G4_bad_flowering_df, G2_good_flowering_df_for_over, G3_good_flowering_df)

under_flowering_df_for_bad_G1_vs_inc <- rbind(G1_bad_flowering_df, low_inc_df_for_badG1, mid_inc_df_for_badG1, high_inc_df_for_badG1)
over_flowering_df_for_bad_G4_vs_inc <- rbind(G4_bad_flowering_df,  low_inc_df_for_badG4, mid_inc_df_for_badG4, high_inc_df_for_badG4)
names(under_flowering_df_for_bad_G1_vs_inc)

###########################################################################################################################################
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggplot2)
### bad_G4 vs. good_G2 & good_G3 ======================================================================================================
# select wth vars
ino_vars <- c("Ino_Manstretta", "Ino_Rossi", "Ino_Dufault", "Ino_Rossi2", "Ino_GIBSIM")
vars <- c(
  "gz_mean",
  "tavg",
  "tmin",
  "rhum",
  "rhum48",
  "prcp_log_plus_1",
  "rainy_day" , 
  ino_vars
)


save_dir <- paste0(
  "./Output/3. obs_pred_graph/2025/",
  version_pred_obs,
  "/compare_daily_wth"
)

# for loops
over_flowering_df_for_bad_G4_vs_inc$GIBSIM_group <- factor(over_flowering_df_for_bad_G4_vs_inc$GIBSIM_group,
                                                           levels = c("bad_G4", "low_inc", "mid_inc", "high_inc"))
over_flowering_df_for_bad_G4_vs_inc$prcp_log_plus_1 <- log(over_flowering_df_for_bad_G4_vs_inc$prcp + 1)
over_flowering_df_0_to_15 <- over_flowering_df_for_bad_G4_vs_inc %>% filter(day_num <= 10)
over_flowering_df_15_to_50 <- over_flowering_df_for_bad_G4_vs_inc %>% filter(day_num >= 10) %>% filter(day_num < 50)

for(v in vars){
  print(v)
  make_kw_boxplot_for_bad_G4(
    data = over_flowering_df_0_to_15,
    y_var = v,
    save_path = save_dir,
    prefix = "box_over_bad_G4_0_15"
  )
  
  make_kw_boxplot_for_bad_G4(
    data = over_flowering_df_15_to_50,
    y_var = v,
    save_path = save_dir,
    prefix = "box_over_bad_G4_15_50"
  )
  
  make_linegraph_compare_vs_inc(over_flowering_df_for_bad_G4_vs_inc, 
                                y_var = v, 
                                save_path = save_dir,
                                prefix = "line_over_bad_G4")
}

### bad_G1 vs. good_G2 & good_G3 ======================================================================================================
# select wth vars
vars <- c(
  "gz_mean",
  "tavg",
  "tmin",
  "rhum",
  "rhum48",
  "prcp_log_plus_1",
  "rainy_day", ino_vars
)

save_dir <- paste0(
  "./Output/3. obs_pred_graph/2025/",
  version_pred_obs,
  "/compare_daily_wth"
)

# for loops
under_flowering_df_for_bad_G1_vs_inc$GIBSIM_group <- factor(under_flowering_df_for_bad_G1_vs_inc$GIBSIM_group,
                                                            levels = c("bad_G1", "low_inc", "mid_inc", "high_inc"))
under_flowering_df_for_bad_G1_vs_inc$prcp_log_plus_1 <- log(under_flowering_df_for_bad_G1_vs_inc$prcp + 1)
under_flowering_df_0_to_15 <- under_flowering_df_for_bad_G1_vs_inc %>% filter(day_num <= 10)
under_flowering_df_15_to_50 <- under_flowering_df_for_bad_G1_vs_inc %>% filter(day_num >= 10) %>% filter(day_num < 50)

for(v in vars){
  print(v)
  make_kw_boxplot_for_bad_G1(
    data = under_flowering_df_0_to_15,
    y_var = v,
    save_path = save_dir,
    prefix = "box_under_bad_G1_0_15"
  )
  
  make_kw_boxplot_for_bad_G1(
    data = under_flowering_df_15_to_50,
    y_var = v,
    save_path = save_dir,
    prefix = "box_under_bad_G1_15_50"
  )
  
  make_linegraph_compare_vs_inc_bad_G1(under_flowering_df_for_bad_G1_vs_inc, 
                                       y_var = v, 
                                       save_path = save_dir,
                                       prefix = "line_under_bad_G1")
}

### Microbes ----------------------------------------------------------------------------------------------------
under_flowering_df_for_bad_G1_vs_inc_groupby <- under_flowering_df_for_bad_G1_vs_inc %>% dplyr::group_by(id, year, id_year) %>% 
  dplyr::summarise(incidence = mean(incidence), 
                   GIBSIM = mean(GIBSIM), 
                   Check_GIBSIM_res = unique(Check_GIBSIM_res), 
                   GIBSIM_group = unique(GIBSIM_group), 
                   Alternaria = mean(Alternaria), 
                   Epicoccum = mean(Epicoccum), 
                   Hannaella = mean(Hannaella), 
                   Periconia = mean(Periconia)
  )

over_flowering_df_for_bad_G4_vs_inc_groupby <- over_flowering_df_for_bad_G4_vs_inc %>% dplyr::group_by(id, year, id_year) %>% 
  dplyr::summarise(incidence = mean(incidence), 
                   GIBSIM = mean(GIBSIM), 
                   Check_GIBSIM_res = unique(Check_GIBSIM_res), 
                   GIBSIM_group = unique(GIBSIM_group), 
                   Alternaria = mean(Alternaria), 
                   Epicoccum = mean(Epicoccum), 
                   Hannaella = mean(Hannaella), 
                   Periconia = mean(Periconia)
  )

mvars <- c("Alternaria", "Epicoccum", "Hannaella", "Periconia")
for(mv in mvars){
  make_kw_boxplot_for_bad_G1_for_m(data = under_flowering_df_for_bad_G1_vs_inc_groupby, 
                                   target_microbe = mv, 
                                   save_path = save_dir, 
                                   prefix = "box_under_bad_G1")
  make_kw_boxplot_for_bad_G4_for_m(data = over_flowering_df_for_bad_G4_vs_inc_groupby, 
                                   target_microbe = mv, 
                                   save_path = save_dir, 
                                   prefix = "box_over_bad_G4")
  
}

cor.test(bf_heading_for_cor$bf_heading_prcp, bf_heading_for_cor$bf_heading_rhum)
plot(bf_heading_for_cor$bf_heading_prcp, bf_heading_for_cor$incidence)
plot(bf_heading_for_cor$bf_heading_rhum, bf_heading_for_cor$incidence)
plot(bf_heading_for_cor$bf_heading_tavg, bf_heading_for_cor$incidence)
### correlation ---------------------------------------------------------------------------------------------------------------------
bf_heading_for_cor_raw <- data_for_linegraph_flowering_join %>% dplyr::filter(day_num <= 10)
bf_heading_for_cor <- bf_heading_for_cor_raw %>% dplyr::group_by(id, year, id_year) %>% dplyr::summarize(bf_heading_tavg = mean(tavg), 
                                                                                                         bf_heading_tmin = mean(tmin), 
                                                                                                         bf_heading_tmax = mean(tmax), 
                                                                                                         bf_heading_rhum = mean(rhum), 
                                                                                                         bf_heading_prcp = mean(log(prcp+1)),
                                                                                                         Alternaria = mean(Alternaria, na.rm = T), 
                                                                                                         Epicoccum = mean(Epicoccum, na.rm = T), 
                                                                                                         Hannaella = mean(Hannaella, na.rm = T), 
                                                                                                         Periconia = mean(Periconia, na.rm = T), 
                                                                                                         bf_heading_Ino_Manstretta_mean = mean(Ino_Manstretta), 
                                                                                                         bf_heading_Ino_Rossi_mean = mean(Ino_Rossi), 
                                                                                                         bf_heading_Ino_GIBSIM_mean = mean(Ino_GIBSIM),
                                                                                                         bf_heading_Ino_Manstretta_max = max(Ino_Manstretta), 
                                                                                                         bf_heading_Ino_Rossi_max = max(Ino_Rossi), 
                                                                                                         bf_heading_Ino_GIBSIM_max = max(Ino_GIBSIM),
                                                                                                         incidence = mean(incidence)
)
as.data.frame(bf_heading_for_cor[which(bf_heading_for_cor$id_year == "wth_067_2025"),])

names(bf_heading_for_cor2)
head(bf_heading_for_cor2)
bf_heading_for_cor2 <- bf_heading_for_cor[,c(4:ncol(bf_heading_for_cor))]
round(cor(bf_heading_for_cor2),2)

library(dplyr)
library(tidyr)
library(purrr)

# 모든 변수 이름
vars <- names(bf_heading_for_cor2)

# 모든 조합 생성 (중복 제거: combn 사용)
var_pairs <- t(combn(vars, 2)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  setNames(c("Variable1", "Variable2"))

# correlation + p-value 계산 함수
calc_cor <- function(v1, v2, data) {
  x <- data[[v1]]
  y <- data[[v2]]
  
  # NA 제거
  df <- na.omit(data.frame(x, y))
  
  # Pearson
  p_test <- cor.test(df$x, df$y, method = "pearson")
  
  # Spearman
  s_test <- cor.test(df$x, df$y, method = "spearman")
  
  tibble(
    pearson = p_test$estimate,
    pearson_pval = p_test$p.value,
    spearman = s_test$estimate,
    spearman_pval = s_test$p.value
  )
}

# 모든 조합에 대해 계산
cor_table <- var_pairs %>%
  mutate(results = map2(Variable1, Variable2, ~calc_cor(.x, .y, bf_heading_for_cor2))) %>%
  unnest(results)

# 결과 확인
cor_table
cor_table$spearman_pval <- round(cor_table$spearman_pval, 3)
cor_table %>% dplyr::filter(Variable1 == "Alternaria" | Variable2 == "Alternaria")
cor_table %>% dplyr::filter(Variable1 == "Epicoccum" | Variable2 == "Epicoccum")
cor_table %>% dplyr::filter(Variable1 == "Hannaella" | Variable2 == "Hannaella")
cor_table %>% dplyr::filter(Variable1 == "Periconia" | Variable2 == "Periconia")

cor_table %>% dplyr::filter(Variable1 == "incidence" | Variable2 == "incidence")


writexl::write_xlsx(cor_table, 
                    path = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/" , 
                                  version_pred_obs,"/Correlation_GZ_and_M_and_wth_",
                                  version_pred_obs,".xlsx") 
)
writexl::write_xlsx(bf_heading_for_cor_raw, 
                    path = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/" , 
                                  version_pred_obs,"/bf_heading_for_cor_raw_",
                                  version_pred_obs,".xlsx") 
)
writexl::write_xlsx(bf_heading_for_cor, 
                    path = paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/" , 
                                  version_pred_obs,"/bf_heading_for_cor_",
                                  version_pred_obs,".xlsx") 
)

plot(bf_heading_for_cor$incidence, bf_heading_for_cor$bf_heading_Ino_Manstretta_mean)
plot(bf_heading_for_cor$incidence, bf_heading_for_cor$bf_heading_Ino_Rossi_mean)
plot(bf_heading_for_cor$incidence, bf_heading_for_cor$bf_heading_Ino_GIBSIM_mean)

plot(bf_heading_for_cor$Alternaria, bf_heading_for_cor$bf_heading_Ino_Manstretta_mean)
plot(bf_heading_for_cor$Alternaria, bf_heading_for_cor$bf_heading_Ino_Rossi_mean)
plot(bf_heading_for_cor$Alternaria, bf_heading_for_cor$bf_heading_Ino_GIBSIM_mean)

plot(bf_heading_for_cor$Epicoccum, bf_heading_for_cor$bf_heading_Ino_Manstretta_mean)
plot(bf_heading_for_cor$Epicoccum, bf_heading_for_cor$bf_heading_Ino_Rossi_mean)
plot(bf_heading_for_cor$Epicoccum, bf_heading_for_cor$bf_heading_Ino_GIBSIM_mean)

plot(bf_heading_for_cor$Hannaella, bf_heading_for_cor$bf_heading_Ino_Manstretta_mean)
plot(bf_heading_for_cor$Hannaella, bf_heading_for_cor$bf_heading_Ino_Rossi_mean)
plot(bf_heading_for_cor$Hannaella, bf_heading_for_cor$bf_heading_Ino_GIBSIM_mean)


##########################################################################################################################################


# df_long_for_corGraph <- bf_heading_for_cor %>%
#   select(
#     incidence,
#     Alternaria, Epicoccum, Hannaella, Periconia,
#     bf_heading_Ino_Manstretta_mean,
#     bf_heading_Ino_Rossi_mean,
#     bf_heading_Ino_GIBSIM_mean
#   ) %>%
#   pivot_longer(
#     cols = c(incidence, Alternaria, Epicoccum, Hannaella, Periconia),
#     names_to = "bio_var",
#     values_to = "bio_value"
#   ) %>%
#   pivot_longer(
#     cols = starts_with("bf_heading_Ino"),
#     names_to = "ino_var",
#     values_to = "ino_value"
#   )
# 
# df_long_for_corGraph$ino_var <- recode(df_long_for_corGraph$ino_var,
#                           bf_heading_Ino_Manstretta_mean = "Manstretta",
#                           bf_heading_Ino_Rossi_mean = "Rossi",
#                           bf_heading_Ino_GIBSIM_mean = "GIBSIM"
# )
# df_long_for_corGraph$bio_var <- factor(df_long_for_corGraph$bio_var,
#                           levels = c("incidence","Alternaria","Epicoccum","Hannaella", "Periconia")
# )
# p <- ggplot(df_long_for_corGraph, aes(x = bio_value, y = ino_value)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = T, linewidth = 0.8, col = "red") +
#   
#   facet_grid(ino_var ~ bio_var, scales = "free") +
#   
#   theme_bw() +
#   theme(
#     strip.text = element_text(size = 11, face = "bold"),
#     axis.title = element_text(size = 12),
#     axis.text  = element_text(size = 10)
#   ) +
#   
#   labs(
#     x = "Biological / Incidence",
#     y = "Inoculum level"
#   ) + 
#   
#   stat_cor( method = "spearman", 
#             aes(label = paste(..r.label.., ..p.label.., sep = "~,~")), 
#             label.x.npc = 0.02, 
#             label.y.npc = 1, 
#             size = 3, 
#             color = "grey50" )
# 
# p

plot(bf_heading_for_cor$Periconia, bf_heading_for_cor$bf_heading_Ino_Manstretta_mean)
plot(bf_heading_for_cor$Periconia, bf_heading_for_cor$bf_heading_Ino_Rossi_mean)
plot(bf_heading_for_cor$Periconia, bf_heading_for_cor$bf_heading_Ino_GIBSIM_mean)

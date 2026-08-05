### version3.3 : remove IS3_2024, IS4_2024

rm(list=ls()) # remove all variables
getwd()
source("./Script/Functions/Function_GIBSIM_for_MB_v2.R") #D:\Microbiome\000.codes\FHB_microbiome_code\FHB_analysis\2025\MB_code_2025\Script\Functions
library(dplyr)
library(conflicted)
detach("package:plyr", unload=TRUE)
conflict_prefer("summarise", "dplyr")
conflict_prefer("group_by", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("year", "lubridate")
conflict_prefer("month", "lubridate")
conflict_prefer("day", "lubridate")
conflict_prefer("yday", "lubridate")

# source("E:/Microbiome/000.codes/FHB_pb_model/Script/Functions/Functions_integrated_fhb_model_v2.0_NH_v0.2.R")
library("stringr")

##### VERSION INFORMATION #####
version <- "MB_v3.3"

# make folder for saving results
pre_path_folder <- "./Output/1.Run_GIBSIM"
if(!file.exists(pre_path_folder)){
  dir.create(pre_path_folder)
}

path_folder <- sprintf("%s%s%s%s", pre_path_folder,"/", version, "/")
if(!file.exists(path_folder)){
  dir.create(path_folder)
}


##### SET PARAMETER VALUES #####

##### _PARAMETER FOR PRELIMINARY TEST #####
NO_DISEASE_MASKING_THERSHOLD = 1

##### _ORIGINAL VALUES: #####
HNG_1 = -0.0127
HNG_2 = 2.4352
FIRST_ANTHER = 3
ANTHER_RSDS_THRESHOLD = 10
MINIMUM_ANTHER_LONGEVITY = 3
MAXIMUM_ANTHER_LONGEVITY = 5
# INFECTION_RHUM_THRESHOLD_1 = 66.01787#72.233
# INFECTION_RHUM_THRESHOLD_2 = 78.01564#65.177
# # INFECTION_RHUM_THRESHOLD_3 = 69.81#69.90
# INF_COEF_1 = 9.99e-03#1.625433e-03
# INF_COEF_2 = 0.269#0.21
INFECTION_RHUM_THRESHOLD_1 = 67.5473
INFECTION_RHUM_THRESHOLD_2 = 75.416
INFECTION_RHUM_THRESHOLD_3 = 79.031
INF_COEF_1 = 0.011
INF_COEF_2 = 0.0001
PRCP_COEF = 0.90



##### _NEW VALUES: #####
# HNG_1 = -0.0447
# HNG_2 = 2.5041
# FIRST_ANTHER = 3
# ANTHER_RSDS_THRESHOLD = 10
# MINIMUM_ANTHER_LONGEVITY = 3
# MAXIMUM_ANTHER_LONGEVITY = 5
# INFECTION_RHUM_THRESHOLD_1 = 75
# INFECTION_RHUM_THRESHOLD_2 = 75
# INFECTION_RHUM_THRESHOLD_3 = 80
# INF_COEF_1 = 0.001029
# INF_COEF_2 = 0.1957

# ##### _NEW VALUES_v1.2: #####
# HNG_1 = -0.0127
# HNG_2 = 2.4352
# FIRST_ANTHER = 3
# ANTHER_RSDS_THRESHOLD = 10
# MINIMUM_ANTHER_LONGEVITY = 3
# MAXIMUM_ANTHER_LONGEVITY = 5
# INFECTION_RHUM_THRESHOLD_1 = 70
# INFECTION_RHUM_THRESHOLD_2 = 70
# INFECTION_RHUM_THRESHOLD_3 = 75
# INF_COEF_1 = 0.001029
# INF_COEF_2 = 0.1957


##### ___________________ #####
# ORIGINAL:
# Gv = 43.10
# Thv = 6.40
# Lc = 9.35
# Av = 1.75
# B = 0.42
# HNG_1 = -0.0127
# HNG_2 = 2.4352
FIRST_ANTHER = 3
ANTHER_RSDS_THRESHOLD = 10
MINIMUM_ANTHER_LONGEVITY = 2
MAXIMUM_ANTHER_LONGEVITY = 5
# INFECTION_RHUM_THRESHOLD_1 <- 80
# INFECTION_RHUM_THRESHOLD_2 <- 80
# INFECTION_RHUM_THRESHOLD_3 <- 85
ST_PERIOD_1 = 7
ST_PERIOD_2 = 7
GZ_INTERCEPT = -0.6306
GZ_RHUM = 0.0152
GZ_CRD = 0.1076
# INF_RHUM_THRESHOLD <- 80
# INF_COEF_1 = 0.001029
# INF_COEF_2 = 0.1957
# INF_INTERCEPT = 0

# host parameter for winter wheat
HNG_1 <- -0.0447
HNG_2 <- 2.5041

# SANN parameters in R: R2 = 0.72, RMSE = 3.471 <= BEST SHAPE
Gv = 43.08508
Thv = 6.212119
Lc = 9.314734
Av = 0.11566
B = 0.553594
GIB_COEF = 4.7
# GIB_COEF = 1
# 
# INF_TAVG = 0.263955031
# INF_RHUM = 0.013343999
# INF_RHUM_THRESHOLD = 80
# INF_PRCP = 0.011793784
# INF_INTERCEPT = -4.441917626
# GIB_COEF = 0.356062754


##### FHB OBSERVATION DATA #####

# data processing
# v0.71: by_NH
# v0.72: by_SSU
#new part by NH (time : 2024/08/30) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
wth_dir <- "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/weather_data/KMA_500m_grid_2022_2025/daily_wth_rsds_for_MB"
path_FHB_obs <- "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/FHB_intensity/23_to_25_FHB_intensity_data_summary_for_MB_v3.xlsx"#"./data/FHB_obs/23_to_25_FHB_intensity_data_summary.xlsx"

df_obs_raw <- read_xlsx(path_FHB_obs)
df_obs_raw <- df_obs_raw[complete.cases(df_obs_raw),]
# filter only wheat and barley
df_obs <- df_obs_raw

df_obs <- df_obs #%>% filter(year == 2024)
nrow(df_obs)

# ### histogram --------------------------------------------------------------------------------------
# df_obs_graph <- df_obs %>% filter(year == 2025)
# df_obs_graph$FHB_incidence <- df_obs_graph$FHB_incidence * 100
# df_obs_graph$FHB_severity <- df_obs_graph$FHB_severity * 100
# hist(df_obs$FHB_incidence)
# ggplot(df_obs_graph, aes(x=FHB_severity)) + 
#   geom_histogram(binwidth = 0.1, color="black", fill="darkgrey") + 
#   theme_classic()+
#   xlab("\nFHB severity (%)") + 
#   ylab("Count\n") + 
#   ylim(c(0,13)) +
#   theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
#   theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
#   theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
#   theme(axis.text.y = element_text(size=15, face='bold',color='black'))
# ###-------------------------------------------------------------------------------------------------

# df_obs[which(df_obs$loc == "?Î¿?1"),]$source <- "epi_2024"
# df_obs[which(df_obs$loc == "?Ø³?2"),]$source <- "epi_2024"
# 
# df_obs[which(df_obs$loc == "?Î¿?1"),]$wth_ID <- "ID_124"
# df_obs[which(df_obs$loc == "?Ø³?2"),]$wth_ID <- "ID_123"
names(df_obs) 
df_obs <- df_obs %>% dplyr::select("source", "new_wth_ID", "year", "FHB_incidence", "FHB_severity")
names(df_obs) <- c("source", "ID", "year", "FHB_incidence" ,"FHB_severity")

df_rslt <- df_obs %>% 
  group_by(source, ID, year) %>% 
  summarise(FHB_incidence = mean(FHB_incidence), 
            FHB_severity = mean(FHB_severity)
  ) %>%
  # ?????? ?? ?ß°? ?? À§Ä¡ Á¶Á¤
  dplyr::mutate(
    FHB_index_predicted = NA,
    FHB_incidence_predicted = NA,
    FHB_severity_predicted = NA
  ) %>%
  relocate(
    FHB_index_predicted, .before = FHB_incidence
  ) %>%
  relocate(
    FHB_incidence_predicted, .after = FHB_incidence
  ) %>%
  relocate(
    FHB_severity_predicted, .after = FHB_severity
  )

df_rslt
nrow(df_rslt)
table(df_rslt$year)
# df_rslt <- df_rslt%>% filter(ID != "ID_020") %>% filter(ID != "ID_049") %>% filter(ID != "ID_025") %>% filter(ID != "ID_047") %>% filter(ID != "ID_044") %>% filter(ID != "ID_046") 
# df_rslt <- df_rslt # %>% filter(ID != "SSU_014") %>% filter(ID != "ID_041") %>% filter(ID != "ID_047") %>% filter(ID != "ID_020")    #%>% filter(ID != "ID_033") %>% filter(ID != "SSU_016") 
nrow(df_rslt)

# we will use incidence data!!! ---------------------------------please check this !!
df_rslt$FHB_severity <- df_rslt$FHB_incidence

# incidence & severity?? percentage?? ??È¯
df_rslt
if(max(df_rslt$FHB_incidence) < 1.1){
  df_rslt$FHB_incidence <- df_rslt$FHB_incidence*100
}
if(max(df_rslt$FHB_severity) < 1.1){
  df_rslt$FHB_severity <- df_rslt$FHB_severity*100
}



### wth info file
wth_info <- readxl::read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/2025_wheat_data/MB_wth_info.xlsx")
df_rslt2 <- df_rslt

df_rslt <- df_rslt2#[c(1:26, 28:97),] #there is no 27 wth file.. so delect that raw

nrow(df_rslt)
head(df_rslt)


##### RUN INTEGRATED MODEL: DVR + GIBSIM #####

# rm(list=ls())
# source("./Script/Functions/Functions_integrated_FHB_model_v0.1.R")

tic("total code run time") #-------------------------------------------------------------------------------------------------------------------------------------------
for(target_weather_data in c('specific_wth_epinet')){
  # target_weather_data='specific_wth_epinet'
  
  tic(paste0(target_weather_data))
  
  ##### LOAD OBS DATA AND WEATEHR DATA INFO #####
  # Read obs data and weather data information ('specific_wth_epinet': site-specific high-resolution weather data from EPINET)
  if(target_weather_data == 'specific_wth_epinet'){
    # obs data with EPINET weather data
    # path_folder_wth <-"./Data/obs_FHB_NICS_SSU/epinet_weather_data/merged_data" 
    path_folder_wth <- wth_dir
  } else {
    stop("Call wrong weather data")
  }
  
  # Add empty columns to assign predicted values
  #head(df_rslt)
  df_rslt <- df_rslt %>% 
    mutate(
      Lat = NA,
      Lon = NA,
      Alt = NA,
      # FHB_index_predicted = NA,
      FHB_incidence_predicted = NA,
      FHB_severity_predicted = NA,
      heading_date = NA,
      heading_date_doy = NA,
      flower_sampling_date = NA,
      disease_sampling_date = NA, 
      
      GZ_mean_max = NA,
      GZ_mean_mean = NA,
      GZ_mean_median = NA
      # DVR_heading_date = NA,
      # DVR_heading_date_doy = NA,
      # Fixed_heading_date = NA,
      # Fixed_heading_date_doy = NA
      
    ) %>% 
    relocate(
      c(Lat, Lon, Alt, heading_date, heading_date_doy,flower_sampling_date,disease_sampling_date,  FHB_incidence_predicted, FHB_severity_predicted, ), .before = FHB_incidence
    )
  # %>% relocate(
  #   FHB_index_predicted, .before = FHB_incidence
  # ) %>% 
  # relocate(
  #   FHB_incidence_predicted, .after = FHB_incidence
  # ) %>% 
  # relocate(
  #   FHB_severity_predicted, .after = FHB_severity
  # )
  
  # empty data frame for saving weather variables used in statistical analysis
  wth_var_stat <- data.frame()
  # for loop
  for(m in 1:nrow(df_rslt)){
    # m=1
    # m=2
    # m=51
    names(df_rslt)[which(names(df_rslt) == "ID")] <- "wth_ID"
    
    target_weather_data = 'specific_wth_epinet'
    target_source <- df_rslt$source[m]
    target_ID <- df_rslt$wth_ID[m]
    target_year <- df_rslt$year[m]
    
    tic(paste0(m,". ",target_ID, " / ", target_year))
    
    ##### LOAD WEATHER DATA #####
    if(target_weather_data == 'specific_wth_epinet'){
      # file path of EPINET weather data
      path_folder_wth <- "D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/weather_data/KMA_500m_grid_2022_2025/daily_wth_rsds_for_MB"
      path_wth <-  paste0(path_folder_wth, "/", target_ID, ".csv")
      df_wth <- read.csv(path_wth) %>%
        mutate(date = paste0(Year, sprintf("%02d", Mon), sprintf("%02d", Day))) %>%
        mutate(id = target_ID) %>% 
        dplyr::select(id, date, tavg, tmin, tmax, prcp, rhum, rsds)
      
      # if(target_source == "epi_2023"){
      #   path_wth <-  paste0(path_folder_wth, "/epinet_for_2023/merged_data/", target_ID, ".xlsx")
      #   
      # df_wth <- read_xlsx(path_wth) %>%
      #   rename(
      #     date = datetime,
      #     rsds = ins,
      #     wsd = wsa
      #   )
      #   
      #   
      # }else if(target_source == "epi_2024"){
      #   path_wth <-  paste0(path_folder_wth, "/epinet_2020_2024/epi_wth_2024_spilt_by_ID/epi_wth_", target_ID, "_2024.xlsx")
      #   
      #   df_wth <- read_xlsx(path_wth)
      #   df_wth$alt <- 28
      #   df_wth$sunshine <- -99
      #   df_wth <- df_wth%>% select("ID",  "Longitude","Latitude", "alt","date", "tavg", "tmin", "tmax","ins", "sunshine","wsa", "rain", "hm") %>%
      #     rename(
      #       id = ID, 
      #       lon = Longitude, 
      #       lat = Latitude, 
      #       rsds = ins, 
      #       wsd = wsa, 
      #       prcp = rain, 
      #       rhum = hm
      #     )
      # }else{
      #   print("error : line 208") 
      # }
    }
    #### wth info file (for Lat, Lon info)
    imsi_rownum <- which(wth_info$Wth_ID == target_ID)[1]
    
    
    #####################################################################################################################
    #####################################################################################################################
    ##### RUN GIBSIM MODEL #####
    
    
    ##### PARAMETERS #####
    QTY_PLANTS <- 300
    QTY_HEADS_BY_PLANTS <- 2
    QTY_ANTHERS_BY_HEADS <- 80
    
    wth <- df_wth
    yr <- target_year
    lon <- wth_info$Longitude[imsi_rownum]#df_wth$lon[1]
    lat <- wth_info$Latitude[imsi_rownum]#df_wth$lat[1]
    alt <- 30 #df_wth$alt[1]
    
    valid_year_range <- sort(unique(lubridate::year(as.Date(wth$date, format = "%Y%m%d"))))[-1]
    if(!yr %in% valid_year_range){
      next
    }
    
    # run FHB simulation: GIBSIM
    # model_result <- integrated_FHB_model(wth, yr, lon, lat, alt, QTY_PLANTS, QTY_HEADS_BY_PLANTS, QTY_ANTHERS_BY_HEADS,
    #                                      HNG_1,
    #                                      HNG_2,
    #                                      FIRST_ANTHER,
    #                                      ANTHER_RSDS_THRESHOLD,
    #                                      MINIMUM_ANTHER_LONGEVITY,
    #                                      MAXIMUM_ANTHER_LONGEVITY,
    #                                      INFECTION_RHUM_THRESHOLD_1,
    #                                      INFECTION_RHUM_THRESHOLD_2,
    #                                      INFECTION_RHUM_THRESHOLD_3,
    #                                      INF_COEF_1,
    #                                      INF_COEF_2
    # )
    model_result <- integrated_fhb_model(wth, yr, lon, lat, alt, QTY_PLANTS, QTY_HEADS_BY_PLANTS, QTY_ANTHERS_BY_HEADS,
                                         Gv,
                                         Thv,
                                         Lc,
                                         Av,
                                         B,
                                         HNG_1,
                                         HNG_2,
                                         FIRST_ANTHER,
                                         ANTHER_RSDS_THRESHOLD,
                                         MINIMUM_ANTHER_LONGEVITY,
                                         MAXIMUM_ANTHER_LONGEVITY,
                                         ST_PERIOD_1,
                                         ST_PERIOD_2,
                                         GZ_INTERCEPT,
                                         GZ_RHUM,
                                         GZ_CRD,
                                         INFECTION_RHUM_THRESHOLD_1,
                                         INFECTION_RHUM_THRESHOLD_2,
                                         INFECTION_RHUM_THRESHOLD_3,
                                         INF_COEF_1,
                                         INF_COEF_2,
                                         PRCP_COEF,
                                         # INF_TAVG,
                                         # INF_RHUM,
                                         # INF_RHUM_THRESHOLD,
                                         # INF_PRCP,
                                         # INF_INTERCEPT,
                                         GIB_COEF
    )
    #####################################################################################################################
    #####################################################################################################################
    if(!dir.exists(paste0("./Output/1.Run_GIBSIM/GIBSIM_daily_results/", version))){
      dir.create(paste0("./Output/1.Run_GIBSIM/GIBSIM_daily_results/", version))
    }
    writexl::write_xlsx(model_result, 
                        paste0("./Output/1.Run_GIBSIM/GIBSIM_daily_results/", version, "/GIBSIM_daily_res_",  target_ID,"_",target_year, ".xlsx"))
    
    heading_date <- max(model_result$heading_date , na.rm = T)
    heading_date_doy <- max(model_result$heading_date_doy , na.rm = T)
    # DVR_heading_date <- max(model_result$DVR_heading_date , na.rm = T)
    # DVR_heading_date_doy <- max(model_result$DVR_heading_date_doy , na.rm = T)
    # Fixed_heading_date <- max(model_result$Fixed_heading_date, na.rm = T)
    # Fixed_heading_date_doy <- max(model_result$Fixed_heading_date_doy, na.rm = T)
    disease_index <- max(model_result$gib4_inc, na.rm = T)
    incidence <- max(model_result$gib4_final, na.rm = T)
    flower_sampling_date<- unique(model_result$flower_sampling_date)[1]
    disease_sampling_date <- unique(model_result$disease_sampling_date)[1]
    
    GZ_mean_max <- max(model_result$gz_mean[which(model_result$st!=0)])
    GZ_mean_mean <- mean(model_result$gz_mean[which(model_result$st!=0)])
    GZ_mean_median <- median(model_result$gz_mean[which(model_result$st!=0)])
    
    # inf_day <- max(model_result$day_inf)
    # tavg <- max(model_result$tavg_inf)
    # prcp <- max(model_result$prcp_inf)
    # rainy_day <- max(model_result$rainy_day_inf)
    # rhum <- max(model_result$rhum_inf)
    
    # plot(model_result$ant[100:150])
    # plot(model_result$st[100:150])
    
    df_rslt[m,]$Lat <- lat
    df_rslt[m,]$Lon <- lon
    df_rslt[m,]$Alt <- alt
    df_rslt[m,]$FHB_severity_predicted <- disease_index
    df_rslt[m,]$FHB_incidence_predicted <- incidence
    df_rslt[m,]$heading_date <- heading_date
    df_rslt[m,]$heading_date_doy <- heading_date_doy
    df_rslt[m,]$flower_sampling_date <- as.Date(flower_sampling_date)  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    df_rslt[m,]$disease_sampling_date <- as.Date(disease_sampling_date)
    
    df_rslt[m,]$GZ_mean_max <- GZ_mean_max
    df_rslt[m,]$GZ_mean_mean <- GZ_mean_mean
    df_rslt[m,]$GZ_mean_median <- GZ_mean_median
    
    
    day_diff_1st_2nd_sample <- as.numeric(disease_sampling_date - flower_sampling_date)
    div_2 <- round(day_diff_1st_2nd_sample/2) 
    div_3 <- round(day_diff_1st_2nd_sample/3) 
    
    # df_rslt[i,]$DVR_heading_date <- DVR_heading_date
    # df_rslt[i,]$DVR_heading_date_doy <- DVR_heading_date_doy
    # df_rslt[i,]$Fixed_heading_date <- Fixed_heading_date
    # df_rslt[i,]$Fixed_heading_date_doy <- Fixed_heading_date_doy
    
    ##### _WEATHER VARIABLES FOR STATISTICAL ANALYSIS ####
    wth_final <- wth
    
    # wth_final$GZ_mean <- model_result$gz_mean
    # wth_final$inf <- model_result$inf
    
    wth_var_stat_row <- wth_final %>%
      mutate(
        date = ymd(date),
        year = lubridate::year(date),
        month = month(date),
        day = day(date),
        doy = yday(date)
      ) %>%
      relocate(
        c(year, month, day, doy), .after = date
      ) %>%
      dplyr::filter(year == yr) %>%
      mutate(
        # d_around_hd = date %in% seq(heading_date-5,heading_date+25,1),
        
        # d_before_hd_10day = date %in% seq(Fixed_heading_date-16,Fixed_heading_date-6,1),
        # d_hd_10day = date %in% seq(Fixed_heading_date-5,Fixed_heading_date+5,1),
        # d_flower_10day = date %in% seq(Fixed_heading_date+6,Fixed_heading_date+16,1),
        # d_after_flower_10day = date %in% seq(Fixed_heading_date+17,Fixed_heading_date+27,1),
        
        d_before_hd_10day = date %in% seq(heading_date-16,heading_date-6,1),
        d_hd_10day = date %in% seq(heading_date-5,heading_date+5,1),
        d_flower_10day = date %in% seq(heading_date+6,heading_date+16,1),
        d_after_flower_10day = date %in% seq(heading_date+17,heading_date+27,1) , #!!!!!!!!!!!!!!!!!!!!!!!!
        
        d_flower_to_disease = date %in% seq(flower_sampling_date,disease_sampling_date,1), 
        
        d_flower_to_disease_div2_1st = date %in% seq(flower_sampling_date,flower_sampling_date+div_2,1), 
        d_flower_to_disease_div2_2nd = date %in% seq(flower_sampling_date+div_2+1,disease_sampling_date,1), 
        
        d_flower_to_disease_div3_1st = date %in% seq(flower_sampling_date,flower_sampling_date+div_3,1), 
        d_flower_to_disease_div3_2nd = date %in% seq(flower_sampling_date+div_3+1,flower_sampling_date+div_3+div_3,1), 
        d_flower_to_disease_div3_3rd = date %in% seq(flower_sampling_date+div_3+div_3+1,disease_sampling_date,1), 
        
        # version 3.2
        d_before_hd_sampling = date %in% seq(flower_sampling_date-24,flower_sampling_date-15,1), 
        d_hd_sampling = date %in% seq(flower_sampling_date-14,flower_sampling_date-5,1), 
        d_flower_sampling = date %in% seq(flower_sampling_date-4,flower_sampling_date+5,1), 
        d_after_flower_sampling = date %in% seq(flower_sampling_date+6,flower_sampling_date+15,1) 
        
        
        # # d_01_31 = day %in% 1:31,
        # d_01_10 = day %in% 1:10,
        # d_11_20 = day %in% 11:20,
        # d_21_31 = day %in% 21:31,
        # 
        # # d_01_31 = day %in% 1:31,
        # d_01_10 = day %in% 1:10,
        # d_11_20 = day %in% 11:20,
        # d_21_31 = day %in% 21:31,
        
      ) %>%
      rename(ID = id) %>%
      group_by(ID, year) %>%
      summarise(
        ## use pure weather data (not GIBSIM parameter)
        
        # temp_around_hd = mean(ifelse(d_around_hd,tavg,NA), na.rm = T),
        # rhum_around_hd = mean(ifelse(d_around_hd,rhum,NA), na.rm = T),
        # prcp_around_hd = sum(ifelse(d_around_hd,prcp,NA), na.rm = T),
        # rsds_around_hd = mean(ifelse(d_around_hd,rsds,NA), na.rm = T),
        # hm_day_around_hd = sum(ifelse(d_around_hd,humid_day,NA), na.rm = T),
        # dew_day_around_hd = sum(ifelse(d_around_hd,dew_day,NA), na.rm = T),
        
        temp_before_hd = mean(ifelse(d_before_hd_10day,tavg,NA), na.rm = T),
        rhum_before_hd = mean(ifelse(d_before_hd_10day,rhum,NA), na.rm = T),
        prcp_before_hd = sum(ifelse(d_before_hd_10day,prcp,NA), na.rm = T),
        # wsd_before_hd = sum(ifelse(d_before_hd_10day,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_before_hd = mean(ifelse(d_before_hd_10day,Tmax_Tmin_diff,NA), na.rm = T),
        
        # rsds_before_hd_10day = mean(ifelse(d_before_hd_10day,rsds,NA), na.rm = T),
        # hm_day_before_hd_10day = sum(ifelse(d_before_hd_10day,humid_day,NA), na.rm = T),
        # dew_day_before_hd_10day = sum(ifelse(d_before_hd_10day,dew_day,NA), na.rm = T),
        
        temp_hd = mean(ifelse(d_hd_10day,tavg,NA), na.rm = T),
        rhum_hd = mean(ifelse(d_hd_10day,rhum,NA), na.rm = T),
        prcp_hd = sum(ifelse(d_hd_10day,prcp,NA), na.rm = T),
        # wsd_hd_ = sum(ifelse(d_hd_10day,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_hd = mean(ifelse(d_hd_10day,Tmax_Tmin_diff,NA), na.rm = T),
        # rsds_hd_10day = mean(ifelse(d_hd_10day,rsds,NA), na.rm = T),
        # hm_day_hd_10day = sum(ifelse(d_hd_10day,humid_day,NA), na.rm = T),
        # dew_day_hd_10day = sum(ifelse(d_hd_10day,dew_day,NA), na.rm = T),
        
        temp_flower = mean(ifelse(d_flower_10day,tavg,NA), na.rm = T),
        rhum_flower = mean(ifelse(d_flower_10day,rhum,NA), na.rm = T),
        prcp_flower = sum(ifelse(d_flower_10day,prcp,NA), na.rm = T),
        # wsd_flower = sum(ifelse(d_flower_10day,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_flower = mean(ifelse(d_flower_10day,Tmax_Tmin_diff,NA), na.rm = T),
        # rsds_flower_10day = mean(ifelse(d_flower_10day,rsds,NA), na.rm = T),
        # hm_day_flower_10day = sum(ifelse(d_flower_10day,humid_day,NA), na.rm = T),
        # dew_day_flower_10day = sum(ifelse(d_flower_10day,dew_day,NA), na.rm = T),
        
        temp_after_flower = mean(ifelse(d_after_flower_10day,tavg,NA), na.rm = T),
        rhum_after_flower = mean(ifelse(d_after_flower_10day,rhum,NA), na.rm = T),
        prcp_after_flower = sum(ifelse(d_after_flower_10day,prcp,NA), na.rm = T) , #!!!!!!!!!!!!!!!!!!!!!
        # wsd_after_flower = sum(ifelse(d_after_flower_10day,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_after_flower = mean(ifelse(d_after_flower_10day,Tmax_Tmin_diff,NA), na.rm = T),
        # rsds_after_flower_10day = mean(ifelse(d_after_flower_10day,rsds,NA), na.rm = T),
        # hm_day_after_flower_10day = sum(ifelse(d_after_flower_10day,humid_day,NA), na.rm = T),
        # dew_day_after_flower_10day = sum(ifelse(d_after_flower_10day,dew_day,NA), na.rm = T),
        
        temp_after_sampling = mean(ifelse(d_flower_to_disease,tavg,NA), na.rm = T),
        rhum_after_sampling = mean(ifelse(d_flower_to_disease,rhum,NA), na.rm = T),
        prcp_after_sampling = sum(ifelse(d_flower_to_disease,prcp,NA), na.rm = T), 
        
        temp_after_sampling_div2_1st = mean(ifelse(d_flower_to_disease_div2_1st,tavg,NA), na.rm = T),
        rhum_after_sampling_div2_1st = mean(ifelse(d_flower_to_disease_div2_1st,rhum,NA), na.rm = T),
        prcp_after_sampling_div2_1st = sum(ifelse(d_flower_to_disease_div2_1st,prcp,NA), na.rm = T),
        
        temp_after_sampling_div2_2nd = mean(ifelse(d_flower_to_disease_div2_2nd,tavg,NA), na.rm = T),
        rhum_after_sampling_div2_2nd = mean(ifelse(d_flower_to_disease_div2_2nd,rhum,NA), na.rm = T),
        prcp_after_sampling_div2_2nd = sum(ifelse(d_flower_to_disease_div2_2nd,prcp,NA), na.rm = T),
        
        temp_after_sampling_div3_1st = mean(ifelse(d_flower_to_disease_div3_1st,tavg,NA), na.rm = T),
        rhum_after_sampling_div3_1st = mean(ifelse(d_flower_to_disease_div3_1st,rhum,NA), na.rm = T),
        prcp_after_sampling_div3_1st = sum(ifelse(d_flower_to_disease_div3_1st,prcp,NA), na.rm = T),
        
        temp_after_sampling_div3_2nd = mean(ifelse(d_flower_to_disease_div3_2nd,tavg,NA), na.rm = T),
        rhum_after_sampling_div3_2nd = mean(ifelse(d_flower_to_disease_div3_2nd,rhum,NA), na.rm = T),
        prcp_after_sampling_div3_2nd = sum(ifelse(d_flower_to_disease_div3_2nd,prcp,NA), na.rm = T),
        
        temp_after_sampling_div3_1st = mean(ifelse(d_flower_to_disease_div3_3rd,tavg,NA), na.rm = T),
        rhum_after_sampling_div3_1st = mean(ifelse(d_flower_to_disease_div3_3rd,rhum,NA), na.rm = T),
        prcp_after_sampling_div3_1st = sum(ifelse(d_flower_to_disease_div3_3rd,prcp,NA), na.rm = T),
        
        
        # version 3.2
        temp_before_hd_sampling_10days = mean(ifelse(d_before_hd_sampling,tavg,NA), na.rm = T),
        rhum_before_hd_sampling_10days = mean(ifelse(d_before_hd_sampling,rhum,NA), na.rm = T),
        prcp_before_hd_sampling_10days = sum(ifelse(d_before_hd_sampling,prcp,NA), na.rm = T),
        
        temp_hd_sampling_10days = mean(ifelse(d_hd_sampling,tavg,NA), na.rm = T),
        rhum_hd_sampling_10days = mean(ifelse(d_hd_sampling,rhum,NA), na.rm = T),
        prcp_hd_sampling_10days = sum(ifelse(d_hd_sampling,prcp,NA), na.rm = T),
        
        temp_flower_sampling_10days = mean(ifelse(d_flower_sampling,tavg,NA), na.rm = T),
        rhum_flower_sampling_10days = mean(ifelse(d_flower_sampling,rhum,NA), na.rm = T),
        prcp_flower_sampling_10days = sum(ifelse(d_flower_sampling,prcp,NA), na.rm = T),
        
        temp_after_flower_sampling_10days = mean(ifelse(d_after_flower_sampling,tavg,NA), na.rm = T),
        rhum_after_flower_sampling_10days = mean(ifelse(d_after_flower_sampling,rhum,NA), na.rm = T),
        prcp_after_flower_sampling_10days = sum(ifelse(d_after_flower_sampling,prcp,NA), na.rm = T)
        
        # temp_april = mean(ifelse(d_01_31,tavg,NA), na.rm = T),
        # rhum_april = mean(ifelse(d_01_31,rhum,NA), na.rm = T),
        # prcp_april = sum(ifelse(d_01_31,prcp,NA), na.rm = T),
        # rsds_april = mean(ifelse(d_01_31,rsds,NA), na.rm = T),
        #
        # temp_may = mean(ifelse(d_01_31 & month == 5,tavg,NA), na.rm = T),
        # rhum_may = mean(ifelse(d_01_31 & month == 5,rhum,NA), na.rm = T),
        # prcp_may = sum(ifelse(d_01_31 & month == 5,prcp,NA), na.rm = T),
        # rsds_may = mean(ifelse(d_01_31 & month == 5,rsds,NA), na.rm = T),
        #
        # temp_04_01_04_10 = mean(ifelse(d_01_10 & month == 4,tavg,NA), na.rm = T),
        # rhum_04_01_04_10 = mean(ifelse(d_01_10 & month == 4,rhum,NA), na.rm = T),
        # prcp_04_01_04_10 = sum(ifelse(d_01_10 & month == 4,prcp,NA), na.rm = T),
        # wsd_04_01_04_10 = mean(ifelse(d_01_10 & month == 4,wsd,NA), na.rm = T),
        # hm_day_04_01_04_10 = sum(ifelse(d_01_10 & month == 4,humid_day,NA), na.rm = T),
        # dew_day_04_01_04_10 = sum(ifelse(d_01_10 & month == 4,dew_day,NA), na.rm = T),
        
        # temp_04_11_04_20 = mean(ifelse(d_11_20 & month == 4,tavg,NA), na.rm = T),
        # rhum_04_11_04_20 = mean(ifelse(d_11_20 & month == 4,rhum,NA), na.rm = T),
        # prcp_04_11_04_20 = sum(ifelse(d_11_20 & month == 4,prcp,NA), na.rm = T),
        # wsd_04_11_04_20 = mean(ifelse(d_11_20 & month == 4,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_04_11_04_20 = mean(ifelse(d_11_20 & month == 4,Tmax_Tmin_diff,NA), na.rm = T),
        # # hm_day_04_11_04_20 = sum(ifelse(d_11_20 & month == 4,humid_day,NA), na.rm = T),
        # # dew_day_04_11_04_20 = sum(ifelse(d_11_20 & month == 4,dew_day,NA), na.rm = T),
        # 
        # temp_04_21_04_30 = mean(ifelse(d_21_31 & month == 4,tavg,NA), na.rm = T),
        # rhum_04_21_04_30 = mean(ifelse(d_21_31 & month == 4,rhum,NA), na.rm = T),
        # prcp_04_21_04_30 = sum(ifelse(d_21_31 & month == 4,prcp,NA), na.rm = T),
        # wsd_04_21_04_30 = mean(ifelse(d_21_31 & month == 4,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_04_21_04_30 = mean(ifelse(d_21_31 & month == 4,Tmax_Tmin_diff,NA), na.rm = T),
        # # hm_day_04_21_04_30 = sum(ifelse(d_21_31 & month == 4,humid_day,NA), na.rm = T),
        # # dew_day_04_21_04_30 = sum(ifelse(d_21_31 & month == 4,dew_day,NA), na.rm = T),
        # 
        # temp_05_01_05_10 = mean(ifelse(d_01_10 & month == 5,tavg,NA), na.rm = T),
        # rhum_05_01_05_10 = mean(ifelse(d_01_10 & month == 5,rhum,NA), na.rm = T),
        # prcp_05_01_05_10 = sum(ifelse(d_01_10 & month == 5,prcp,NA), na.rm = T),
        # wsd_05_01_05_10 = mean(ifelse(d_01_10 & month == 5,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_05_01_05_10 = mean(ifelse(d_01_10 & month == 5,Tmax_Tmin_diff,NA), na.rm = T),
        # # hm_day_05_01_05_10 = sum(ifelse(d_01_10 & month == 5,humid_day,NA), na.rm = T),
        # # dew_day_05_01_05_10 = sum(ifelse(d_01_10 & month == 5,dew_day,NA), na.rm = T),
        # #
        # temp_05_11_05_20 = mean(ifelse(d_11_20 & month == 5,tavg,NA), na.rm = T),
        # rhum_05_11_05_20 = mean(ifelse(d_11_20 & month == 5,rhum,NA), na.rm = T),
        # prcp_05_11_05_20 = sum(ifelse(d_11_20 & month == 5,prcp,NA), na.rm = T),
        # wsd_05_11_05_20 = mean(ifelse(d_11_20 & month == 5,wsd,NA), na.rm = T),
        # Tmax_Tmin_diff_05_11_05_20 = mean(ifelse(d_11_20 & month == 5,Tmax_Tmin_diff,NA), na.rm = T),
        # # hm_day_05_11_05_20 = sum(ifelse(d_11_20 & month == 5,humid_day,NA), na.rm = T),
        # # dew_day_05_11_05_20 = sum(ifelse(d_11_20 & month == 5,dew_day,NA), na.rm = T),
        # #
        # # temp_05_21_05_31 = mean(ifelse(d_21_31 & month == 5,tavg,NA), na.rm = T),
        # # rhum_05_21_05_31 = mean(ifelse(d_21_31 & month == 5,rhum,NA), na.rm = T),
        # # prcp_05_21_05_31 = sum(ifelse(d_21_31 & month == 5,prcp,NA), na.rm = T),
        # # rsds_05_21_05_31 = mean(ifelse(d_21_31 & month == 5,rsds,NA), na.rm = T),
        # # hm_day_05_21_05_31 = sum(ifelse(d_21_31 & month == 5,humid_day,NA), na.rm = T),
        # # dew_day_05_21_05_31 = sum(ifelse(d_21_31 & month == 5,dew_day,NA), na.rm = T)
        
        
        # ## use GIBSIM parameter
        # # GIBSIM_temp_before_hd_10day = mean(ifelse(d_before_hd_10day,GIBSIM_temp,NA), na.rm = T),
        # # GIBSIM_rhum_before_hd_10day = mean(ifelse(d_before_hd_10day,GIBSIM_rhum,NA), na.rm = T),
        # # GIBSIM_prcp_before_hd_10day = sum(ifelse(d_before_hd_10day,GIBSIM_prcp,NA), na.rm = T),
        # GIBSIM_INF_before_hd = mean(ifelse(d_before_hd_10day,inf,NA), na.rm = T),
        # GIBSIM_GZ_sum_before_hd = sum(ifelse(d_before_hd_10day,GZ_mean,NA), na.rm = T),
        # 
        # # GIBSIM_temp_hd_10day = mean(ifelse(d_hd_10day,GIBSIM_temp,NA), na.rm = T),
        # # GIBSIM_rhum_hd_10day = mean(ifelse(d_hd_10day,GIBSIM_rhum,NA), na.rm = T),
        # # GIBSIM_prcp_hd_10day = sum(ifelse(d_hd_10day,GIBSIM_prcp,NA), na.rm = T),
        # GIBSIM_INF_hd = mean(ifelse(d_hd_10day,inf,NA), na.rm = T),
        # GIBSIM_GZ_sum_hd = sum(ifelse(d_hd_10day,GZ_mean,NA), na.rm = T),
        # 
        # 
        # # GIBSIM_temp_flower_10day = mean(ifelse(d_flower_10day,GIBSIM_temp,NA), na.rm = T),
        # # GIBSIM_rhum_flower_10day = mean(ifelse(d_flower_10day,GIBSIM_rhum,NA), na.rm = T),
        # # GIBSIM_prcp_flower_10day = sum(ifelse(d_flower_10day,GIBSIM_prcp,NA), na.rm = T),
        # GIBSIM_INF_flower = mean(ifelse(d_flower_10day,inf,NA), na.rm = T),
        # GIBSIM_GZ_sum_flower = sum(ifelse(d_flower_10day,GZ_mean,NA), na.rm = T),
        # 
        # 
        # # GIBSIM_temp_after_flower_10day = mean(ifelse(d_after_flower_10day,GIBSIM_temp,NA), na.rm = T),
        # # GIBSIM_rhum_after_flower_10day = mean(ifelse(d_after_flower_10day,GIBSIM_rhum,NA), na.rm = T),
        # # GIBSIM_prcp_after_flower_10day = sum(ifelse(d_after_flower_10day,GIBSIM_prcp,NA), na.rm = T),
        # GIBSIM_INF_after_flower = mean(ifelse(d_after_flower_10day,inf,NA), na.rm = T),
        # GIBSIM_GZ_sum_after_flower = sum(ifelse(d_after_flower_10day,GZ_mean,NA), na.rm = T),
        
        
        # INF_04_11_04_20 = mean(ifelse(d_11_20 & month == 4,INF,NA), na.rm = T),
        # INF_04_21_04_30 = mean(ifelse(d_21_31 & month == 4,INF,NA), na.rm = T),
        # INF_05_01_05_10 = mean(ifelse(d_01_10 & month == 5,INF,NA), na.rm = T),
        # INF_05_11_05_20 = mean(ifelse(d_11_20 & month == 5,INF,NA), na.rm = T),
        # 
        # GZ_04_11_04_20 = sum(ifelse(d_11_20 & month == 4,GIBSIM_GZ,NA), na.rm = T),
        # GZ_04_21_04_30 = sum(ifelse(d_21_31 & month == 4,GIBSIM_GZ,NA), na.rm = T),
        # GZ_05_01_05_10 = sum(ifelse(d_01_10 & month == 5,GIBSIM_GZ,NA), na.rm = T),
        # GZ_05_11_05_20 = sum(ifelse(d_11_20 & month == 5,GIBSIM_GZ,NA), na.rm = T),
        # 
        # GZ_mean_04_11_04_20 = sum(ifelse(d_11_20 & month == 4,GIBSIM_GZ_mean,NA), na.rm = T),
        # GZ_mean_04_21_04_30 = sum(ifelse(d_21_31 & month == 4,GIBSIM_GZ_mean,NA), na.rm = T),
        # GZ_mean_05_01_05_10 = sum(ifelse(d_01_10 & month == 5,GIBSIM_GZ_mean,NA), na.rm = T),
        # GZ_mean_05_11_05_20 = sum(ifelse(d_11_20 & month == 5,GIBSIM_GZ_mean,NA), na.rm = T),
      )
    
    wth_var_stat <- rbind(wth_var_stat, wth_var_stat_row)
    
    toc()
  }
  
  # merge both results
  wth_var_stat <- wth_var_stat[!duplicated(wth_var_stat),]
  #is it true? : length(wth_var_stat$ID) == length(unique(wth_var_stat$ID))
  names(wth_var_stat)[which(names(wth_var_stat)=="ID")] <- "wth_ID"
  df_rslt_stat <- left_join(df_rslt,wth_var_stat, by = c("wth_ID", "year"))
  
  names(df_rslt_stat)
  df_rslt_stat$year.y <- NULL
  names(df_rslt_stat)[which(names(df_rslt_stat) == "year.x")] <- "year"
  
  # save result
  if(nrow(df_rslt) == nrow(df_rslt_stat)){
    path_rslt <-  paste0(path_folder,"/","MB_wheat_cor_wth_data",
                         "_", version,
                         ".xlsx")
    writexl::write_xlsx(df_rslt_stat, path_rslt)
  }else{
    print("ERROR : nrow(df_rslt) == nrow(df_rslt_stat) is not true")
  }
  
}

# df_rslt_stat$FHB_severity_predicted <- NULL
df_rslt_stat
nrow(df_rslt_stat)
nrow(df_rslt)
nrow(wth_var_stat)
table(wth_var_stat$wth_ID)
plot(df_rslt_stat$FHB_incidence_predicted, df_rslt_stat$FHB_incidence)
summary(lm(FHB_incidence ~ FHB_incidence_predicted, data = df_rslt_stat))


df_rslt_stat_file_path <- paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/1.Run_GIBSIM/",
                                 version,"/MB_wheat_cor_wth_data_",
                                 version ,".xlsx")

df_rslt_stat <- readxl::read_xlsx(df_rslt_stat_file_path)
# df_rslt_stat <- readxl::read_xlsx("E:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/1.Run_GIBSIM/MB_v0.2/MB_wheat_cor_wth_data_MB_v0.2.xlsx")
plot(df_rslt_stat$FHB_incidence_predicted, df_rslt_stat$FHB_incidence)
complete_df <- df_rslt_stat

# filter outlier
complete_df <- complete_df[-c(which(df_rslt_stat$wth_ID == "wth_063" & df_rslt_stat$year == 2024),
                              which(df_rslt_stat$wth_ID == "wth_061" & df_rslt_stat$year == 2024)
),]

# complete_df <- complete_df[-c(which(df_rslt_stat$wth_ID == "wth_066" & df_rslt_stat$year == 2024),
#                               which(df_rslt_stat$wth_ID == "wth_064" & df_rslt_stat$year == 2024),
#                               which(df_rslt_stat$wth_ID == "wth_062" & df_rslt_stat$year == 2024)
# ),]
# 
# complete_df <- complete_df[-c(which(df_rslt_stat$wth_ID == "wth_070" & df_rslt_stat$year == 2025),
#                               which(df_rslt_stat$wth_ID == "wth_019" & df_rslt_stat$year == 2025),
#                               which(df_rslt_stat$wth_ID == "wth_067" & df_rslt_stat$year == 2025)
# ),]

# complete_df <- complete_df[-which(complete_df$wth_ID == "ID_041"),]
nrow(complete_df)
plot(complete_df$FHB_incidence_predicted, complete_df$FHB_incidence)
obs_pred_inc <- lm(FHB_incidence_predicted ~FHB_incidence, data = complete_df)
summary(obs_pred_inc)
abline(obs_pred_inc)
# plot(complete_df$FHB_severity, complete_df$FHB_index_predicted)
# obs_pred_sev <- lm(FHB_severity_predicted ~FHB_severity, data = complete_df)
# summary(obs_pred_sev)






##### PLOTTING #####
# df_plot <- objfun_result(new_parameter, df_rslt)$sim
# original model performance
# df_plot <- objfun_result(c(4.7), df_rslt)$sim

# section 1: inputting values (x = predicted values, y = observed values)
df_plot <- complete_df %>% dplyr::select(wth_ID, source, FHB_incidence, FHB_incidence_predicted)
df_plot2 <- df_plot %>% group_by(wth_ID, source) %>% summarise(FHB_incidence = mean(FHB_incidence), 
                                                               FHB_incidence_predicted = mean(FHB_incidence_predicted))
nrow(df_plot2)


# df_plot2 <- df_plot2[-c(91,88,85,93),]
# df_plot2 <- df_plot2[-c(91,88,93),]
if(df_plot2$FHB_incidence[1] < 1){
  df_plot2$FHB_incidence <- df_plot2$FHB_incidence*100
}else{
  df_plot2$FHB_incidence <- df_plot2$FHB_incidence
}

df_plot2$FHB_incidence_predicted[which(df_plot2$FHB_incidence_predicted < 0)] <- 0


x_element = df_plot2$FHB_incidence_predicted
y_element = df_plot2$FHB_incidence
x_label = "\nPredicted FHB incidence"
y_label = "Observed FHB incidence\n"
## pred_value_type: If you don't want to show the type of your predicted value, then pred_value_type = ""
pred_value_type = target_value = "incidence"
# pred_value_type = ""


# section 2: plotting
res = lm(y_element ~ x_element)
summary_res <- summary(res)
coef_df <- as.data.frame(summary_res$coefficients)
coef_df$Estimate
x_element <- df_plot2$FHB_incidence_predicted <- df_plot2$FHB_incidence_predicted*coef_df$Estimate[2] + coef_df$Estimate[1]
x_element[which(x_element < 0)] <- 0
x_element[which(x_element > 100)] <- 100
rmse <- sqrt(mean((df_plot2$FHB_incidence - df_plot2$FHB_incidence_predicted)^2))
res2 = lm(y_element ~ x_element)
summary_res2 <- summary(res2)
p1 <- ggplot(mapping = aes(x = x_element))+
  theme_classic() +
  geom_point(mapping = aes(y = y_element), shape = 1) +
  stat_smooth(mapping = aes(y = x_element, color = "1 : 1"), formula = y ~ x, method = 'lm', se = F) +
  stat_smooth(mapping = aes(y = y_element, color = "observed ~ predicted"), formula = y ~ x, method = 'lm', se = F) +
  scale_color_manual(name = "",
                     breaks=c('1 : 1', 'observed ~ predicted'),
                     values=c('1 : 1'= adjustcolor("red",alpha=0.8), 'observed ~ predicted' =adjustcolor("black",alpha=1))) +
  # scale_x_continuous(limits = c(0,max(x_element))) +
  # scale_y_continuous(limits = c(0,max(y_element))) +
  # plot.margin = unit(c(0, 0, 0, 0), "cm") +
  labs(title = paste0("Observed VS Predicted values", if(!pred_value_type == "") paste0(" (", pred_value_type, ")") ),
       subtitle = bquote(y == .(res2$coefficients[2]) * x + .(res2$coefficients[1]) ~ " / " ~ n == .(length(x_element)) ~ " / " ~ R^2 == .(round(summary_res2$r.squared, 2)) ~ " / " ~ RMSE == .(round(rmse, 2))),
       # caption = "caption",
       x= x_label,
       y= y_label)+
  theme(
    # legend.background = element_rect(fill="white", linewidth = 0.5, linetype="solid", colour = "black"),
    legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  # theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=13,color='black'))+ #, face='bold'
  theme(axis.text.y = element_text(size=13,color='black')) #, face='bold'

print(p1)
assign(paste0("p_", pred_value_type), p1)

# section 3: saving plot
ggsave(
  plot = get(paste0("p_", pred_value_type)),
  file = paste0(
    path_folder,"/",
    "GIBSIM", "_",
    "matching_obs_pred",
    "_", version,
    ".jpg"
  ),
  width = 15,
  height = 15,
  units = c("cm")
)
hist(x_element)



complete_df$adj_FHB_incidence <- complete_df$FHB_incidence_predicted*coef_df$Estimate[2] + coef_df$Estimate[1]
complete_df_save_file_path <- paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/1.Run_GIBSIM/",
                                     version,"/MB_wheat_wth_data_adj_inc_",
                                     version,".xlsx")
writexl::write_xlsx(complete_df, path = complete_df_save_file_path)

plot(complete_df$GZ_mean_median, complete_df$FHB_incidence)
complete_df2 <- complete_df
complete_df2$year <- as.factor(complete_df2$year)

GZ_mean_median_g <- ggplot2::ggplot(complete_df2, aes(x = GZ_mean_median, y = FHB_incidence, color = year)) +
  geom_point() + 
  labs(#title = "GZ_mean_median vs. FHB incidence", 
    x = "\nGZ_mean_median",
    y = "FHB incidence\n") +
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))

GZ_mean_mean_g <- ggplot2::ggplot(complete_df2, aes(x = GZ_mean_mean, y = FHB_incidence, color = year)) +
  geom_point() + 
  labs(#title = "GZ_mean_median vs. FHB incidence", 
    x = "\nGZ_mean_mean",
    y = "FHB incidence\n") +
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))

GZ_mean_max_g <- ggplot2::ggplot(complete_df2, aes(x = GZ_mean_max, y = FHB_incidence, color = year)) +
  geom_point() + 
  labs(#title = "GZ_mean_median vs. FHB incidence", 
    x = "\nGZ_mean_max",
    y = "FHB incidence\n") +
  theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))

ggsave(
  plot = GZ_mean_median_g,
  file = paste0(
    "./Output/1.Run_GIBSIM/",
    version,
    "/GZ_mean_median_", 
    version,
    ".png"
  ),
  width = 15.5,
  height = 12,
  units = c("cm")
)
ggsave(
  plot = GZ_mean_mean_g,
  file = paste0(
    "./Output/1.Run_GIBSIM/",
    version,
    "/GZ_mean_mean_", 
    version,
    ".png"
  ),
  width = 15.5,
  height = 12,
  units = c("cm")
)
ggsave(
  plot = GZ_mean_max_g,
  file = paste0(
    "./Output/1.Run_GIBSIM/",
    version,
    "/GZ_mean_max_", 
    version,
    ".png"
  ),
  width = 15.5,
  height = 12,
  units = c("cm")
)
# 
# ##### PLOTTING OPTIMIZATION TRACE #####
# 
# target_N <- representative_row$Nth_CV
# target_fold <- representative_row$fold
# id <- paste0(target_N ,"_", target_fold)
# df_trace_plot <- list_trace_matrix[[id]] %>%
#   as_tibble() %>% 
#   mutate(index = index(list_trace_matrix[[id]])) %>% 
#   relocate(index)
# # section 1: inputting values (x = predicted values, y = observed values)
# x_element = df_trace_plot$index
# y_element = df_trace_plot$current.minimum
# x_label = "Index?"
# y_label = "RMSE"
# ## pred_value_type: If you don't want to show the type of your predicted value, then pred_value_type = ""
# pred_value_type = target_value = paste("RMSE")
# # pred_value_type = ""
# plot_counts <- df_total_result[df_total_result$Nth_CV == target_N & df_total_result$fold == target_fold, "counts"]
# plot_train_RMSE <- df_total_result[df_total_result$Nth_CV == target_N & df_total_result$fold == target_fold, "train_RMSE"] %>% round(1)
# plot_train_R2 <- df_total_result[df_total_result$Nth_CV == target_N & df_total_result$fold == target_fold, "train_R2"] %>% round(2)
# plot_test_RMSE <- df_total_result[df_total_result$Nth_CV == target_N & df_total_result$fold == target_fold, "test_RMSE"] %>% round(1)
# plot_test_R2 <- df_total_result[df_total_result$Nth_CV == target_N & df_total_result$fold == target_fold, "test_R2"] %>% round(2)
# 
# # section 2: plotting
# p1 <- ggplot(mapping = aes(x = x_element))+
#   theme_classic() +
#   geom_line(mapping = aes(y = y_element), linetype = 1, linewidth = 1, color = "Blue") +
#   # scale_x_continuous(limits = c(0,max(x_element))) +
#   # scale_y_continuous(limits = c(0,max(y_element))) +
#   # plot.margin = unit(c(0, 0, 0, 0), "cm") +
#   labs(title = paste0("Optimization Process Trace of ", target_N, "-", target_fold, " fold", if(!pred_value_type == "") paste0(" (", pred_value_type, ")") ),
#        subtitle = bquote(count == .(plot_counts) ~ "/" ~ train ~ RMSE == .(plot_train_RMSE) ~ "," ~ R^2 == .(plot_train_R2) ~ "/" ~ test ~ RMSE == .(plot_test_RMSE) ~ "," ~ R^2 == .(plot_test_R2)),
#        # caption = "caption",
#        x= x_label,
#        y= y_label)+
#   theme(
#     # legend.background = element_rect(fill="white", linewidth = 0.5, linetype="solid", colour = "black"),
#     legend.position = "bottom") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# print(p1)
# assign(paste0("p_", pred_value_type), p1)
# 
# # section 3: saving plot
# ggsave(
#   plot = get(paste0("p_", pred_value_type)),
#   file = paste0(
#     path_folder,"/",
#     "GIBSIM", "_",
#     "trace",
#     "_", version,
#     ".jpg"
#   ),
#   width = 15,
#   height = 15,
#   units = c("cm")
# )
# 

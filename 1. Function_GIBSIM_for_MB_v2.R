##### VERSION DISCRIPTION #####
##### MB_v1 : fix 2025` heading date, results contain sampling date (flowering, disease sampling)
##### MB_v2 : using DVR to calculate heading date

# additional option: decadal interval option / past 10 years interval option
# incidence => intensity: intensity ranges from 0 to 100


# rm(list=ls())

# library("data.table")
# #install.packages("Metrics")
# library("Metrics")
# # Fetch weather for year 2000 wet season for a 120 day rice variety at the IRRI
# # Zeigler Experiment Station
# library(readxl)
# library(dplyr)
# library(tidyverse)
# library(magrittr)
# library(lubridate)
# install.packages("geosphere")
# install.packages("imputeTS")
# install.packages("nloptr")
# install.packages('parrallel')
# install.packages("cvTools")

library(tidyverse)
library(lubridate)
library(zoo)
library(tictoc)
library(data.table)
library(patchwork)
library(ggtext)
library(geosphere)
library(imputeTS)
library(car)
library(readxl)
library(writexl)
library(openxlsx)
library(nloptr)
library(parallel)
library(doParallel)
library(foreach)
library(cvTools)
library(dplyr)

##### _______________ #####

##### FUNCTION: READ WEATHER DATA OF MODEL #####
read_model <- function(model, path){
  tic(model)
  file_read <- list.files(paste(path, model, sep = "/"), pattern = ".csv")
  model_read <- paste(path,  model, file_read, sep = "/") %>% lapply(read.csv)
  names(model_read) <- file_read
  toc()
  return(model_read)
}

##### _______________ #####

##### FUNCTION: READ SECOND SHEET OF EXCEL #####
read_xlsx_2sheet <- function(path){
  return(read_xlsx(path, sheet = 2))
}

##### _______________ #####



##### FUNCTION: INTEGRATED FHB MODEL: DVR + GIBSIM #####

# # ORIGINAL:
# hd_diff = 0
# Gv = 43.10
# Thv = 6.40
# Lc = 9.35
# Av = 1.75
# B = 0.42
# HNG_1 = -0.0127
# HNG_2 = 2.4352
# FIRST_ANTHER = 3
# ANTHER_RSDS_THRESHOLD = 10
# MINIMUM_ANTHER_LONGEVITY = 2
# MAXIMUM_ANTHER_LONGEVITY = 5
# ST_PERIOD_1 = 7
# ST_PERIOD_2 = 7
# GZ_INTERCEPT = -0.6306
# GZ_RHUM = 0.0152
# GZ_CRD = 0.1076
# INFECTION_RHUM_THRESHOLD_1 = 80
# INFECTION_RHUM_THRESHOLD_2 = 80
# INFECTION_RHUM_THRESHOLD_3 = 85
# INF_COEF_1 = 0.001029
# INF_COEF_2 = 0.1957
# INF_INTERCEPT = 0
# GIB_COEF = 1


integrated_fhb_model <- function(wth, yr, lon, lat, alt, QTY_PLANTS, QTY_HEADS_BY_PLANTS, QTY_ANTHERS_BY_HEADS,
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
){
  
  ##### WEATHER DATA PROCESSING #####
  wth <- wth %>% 
    mutate(
      date = ymd(date),
      year = year(date),
      month = month(date),
      day = day(date),
      doy = yday(date)
    ) %>% 
    relocate(
      c(year, month, day, doy), .after = date
    )
  
  # If there is no tavg column, make one.
  if(!"tavg" %in% colnames(wth)){
    wth$tavg <- (wth$tmax + wth$tmin) / 2
  }
  
  ##### DVR MODEL #####
  # x = daily minimum temperature in January averaged over 10 years at original version
  yr_range <- unique(wth$year)
  
  # if(yr %in% 1981:1990){
  #   yr_range <- 1981:1990
  # } else if(yr %in% 1991:2000){
  #   yr_range <- 1991:2000
  # } else if(yr %in% 2001:2010){
  #   yr_range <- 2001:2010
  # } else if(yr %in% 2011:2022){
  #   yr_range <- 2011:2020
  # }
  
  # if(target_stn_name == "asos256"){
  #   if(yr %in% 2011:2012){
  #     yr_range <- 2001:2010
  #   } else if(yr > 2012){
  #     next
  #   }
  # }
  
  
  target_tmin_JAN <- x <- wth %>% 
    filter(year %in% yr_range) %>% 
    group_by(month) %>%
    summarise(tmin_JAN = mean(tmin)) %>% 
    filter(month == 1) %>% pull()
  
  target_tmin_JAN
  
  # if(x < -9 & x > 10){
  #   next
  # }
  
  ##### _Calculate optimal sowing date #####
  if(alt < 100){
    # If the altitude is below 100m
    y = (-0.1017 * x^2) + (2.2899 * x) + 305.98
  } else if(alt >= 100){
    # If the altitude is more than 100m
    y = (-0.081 * x^2) + (2.2603 * x) + 299.35
  }
  
  sowing_doy <- y
  
  ##### _Filter the weather data and designate sowing date #####
  target_sim <- wth %>% 
    filter(year %in% (yr-1):yr) %>% 
    mutate(
      sowing_date_logical = FALSE,
      heading_date_logical = FALSE
    ) %>% 
    relocate(doy, .after = day)
  
  ##### _predicted sowing date: DVI = 0 #####
  if(!(yr - 1) %in% target_sim$year) {
    next
  }
  target_sim[(target_sim$year == yr-1) & (target_sim$doy == ceiling(y)),]$sowing_date_logical <- TRUE
  target_sim <- target_sim[which(target_sim$sowing_date_logical):nrow(target_sim),]
  
  
  ##### _Run model #####
  target_sim <- target_sim %>% 
    mutate(
      # target_sim[1,]$
      
      # parameters
      Gv = Gv,
      Thv = Thv,
      Lc = Lc,
      Av = Av,
      B = B,
      
      # # new parameters: R2 = 0.74
      # Gv = 49.3119,
      # Thv = 6.96108,
      # Lc = 9.66437,
      # Av = 0.070191,
      # B = 1.36713,
      
      ##### _L: daily photoperiod #####
      # phi: latitude
      phi = pi/180 * lat,
      # delta: the declination of the sun
      # delta = 23.45 * sin(pi/180 * 360/365 * (284 + doy)),
      delta = 23.45 * cos(pi/180 * 360/365 * (doy - 172)),
      # L: daily photoperiod
      L = 24/pi * (acos(-1 * tan(phi) * tan(pi/180 * delta))),
      
      ##### _DVR #####
      DVR = ifelse(B*(L-Lc) >= 0, 1/Gv * (1-exp(-B*(L-Lc))) / (1+exp(-Av*(tavg-Thv))),0),
      
      ##### _DVI: accumulated DVR #####
      DVI = cumsum(DVR)
    )
  
  
  ##### _predicted heading date: DVI = 1 #####
  target_sim[min(which(target_sim$DVI >= 1)), ]$heading_date_logical <- TRUE
  
  
  ##### _SOWING AND HEADING DATE BY DVR MODEL #####
  
  sowing_date_logical <- target_sim[min(which(target_sim$sowing_date_logical)), ]$date %>% ymd()
  heading_date_40 <- target_sim[min(which(target_sim$heading_date_logical)), ]$date %>% ymd()
  
  temp <- seq(1,25)
  hng_acc <- 1 - exp(HNG_1 * (temp ^ HNG_2))
  hd_diff <- min(which(hng_acc > 0.4))
  
  #version 0.2.5
  # hd_diff_change = -10
  heading_starting_date <- heading_date_40 - hd_diff + 1
  
  
  
  ##### heading date ##### 
  # version : Functions_GIBSIM_for_MB_v1
  if(yr == 2023){
    heading_starting_date <- heading_starting_date
    flower_sampling_date <- as.Date("2023-05-01")
    disease_sampling_date <- as.Date("2023-05-22")
  }else if(yr == 2024){
    heading_starting_date <- heading_starting_date
    flower_sampling_date <- as.Date("2024-04-27")
    disease_sampling_date <- as.Date("2024-05-14")
  }else if(yr == 2025){
    heading_starting_date <- heading_starting_date #as.Date("2025-04-26")
    flower_sampling_date <- as.Date("2025-04-29")
    disease_sampling_date <- as.Date("2025-05-19")
  }else{
    heading_starting_date <- heading_starting_date
    flower_sampling_date <- heading_starting_date + 10
    disease_sampling_date <- flower_sampling_date + 21
  }
  
  
  ##### GIBSIM MODEL #####
  
  ##### _WEATHER DATA PROCESSING #####
  
  # Filter wth by yr
  wth <- target_sim %>% 
    filter(year %in% (yr-1):(yr))
  
  # days after heading date column: should be revised
  # v0.9
  wth <- wth %>% 
    mutate(day_hd = doy - yday(heading_starting_date) + 1)
  wth[1:which(wth$date == heading_starting_date - 1),]$day_hd <- 0
  
  
  
  # Add tavg48, rhum48 columns: 48-hour average value of temperature, relative humidity
  wth$tavg48 <- rollmean(wth$tavg, k = 2, fill = wth$tavg[1], align = "right")
  wth$tavg72 <- rollmean(wth$tavg, k = 3, fill = wth$tavg[1], align = "right")
  wth$rhum48 <- rollmean(wth$rhum, k = 2, fill = wth$rhum[1], align = "right")
  wth$rhum72 <- rollmean(wth$rhum, k = 3, fill = wth$rhum[1], align = "right")
  
  
  
  ##### _MODEL START #####
  
  # Set variables
  hng <- hng_daily_diff <- head_daily_emerged <- a <- b <- anther <- ant <- st <- rainy_day <- crd <- gz <- gz_mean <- inf <- inf_mean <- gib1 <- gib1_acc_percent <- gib2 <- gib2_acc_percent <- gib3 <- gib3_acc_percent <- gib4 <- gib4_acc_percent <- 
    rep(0, times = length(wth$date))
  
  
  for(i in 1:length(wth$date)) {
    # i = which(wth$date == heading_starting_date)
    # i = 110
    
    ##### _HEAD MODULE #####
    
    # HNG: the daily cumulative proportion of heads emerged
    hng[i] <- 1 - exp(-0.0127 * (wth$day_hd[i] ^ 2.4352))
    
    # HNG daily diff
    if(i == 1){
      hng_daily_diff[i] <- hng[i]
    } else if (hng[i] - hng[i-1] > 0) {
      hng_daily_diff[i] <- hng[i] - hng[i-1]
    } else {
      hng_daily_diff[i] <- 0
    }
    
    # heads daily emerged
    head_daily_emerged[i] <- hng_daily_diff[i] * QTY_HEADS_BY_PLANTS * QTY_PLANTS
    
    
    
    ##### _ANTHER MODULE #####
    
    # a
    if(wth$day_hd[i] > 0){
      a[i] <- 0.255 + (-0.029 * wth$tavg[i]) + (0.0009 * (wth$tavg[i]) ^ 2)
    }
    
    # b
    if(wth$day_hd[i] > 0){
      b[i] <- -5.773 + (0.966 * wth$tavg[i]) + (-0.0278 * (wth$tavg[i]) ^ 2)
    }
    
    ##### _GZ(Gibberella zeae) MODULE: inoculum factor #####
    
    # consecutive rainy days: select threshold
    if(i == 1){
      if(wth$prcp[i] > 0.3){
        rainy_day[i] <- 1
      } else{
        rainy_day[i] <- 0
      }
    } else if(wth$prcp[i] > 0.3){
      rainy_day[i] <- rainy_day[i-1] + 1
    } else{
      rainy_day[i] <- 0
    }
    
    # CRD: dummy variable for a position of a rainy( >0.3mm) day in a consecutive period of rainy days.
    if(rainy_day[i] == 0){
      crd[i] <- 0
    } else if(rainy_day[i] == 1){
      crd[i] <- 1
    } else if(rainy_day[i] == 2){
      crd[i] <- 2
    } else if(rainy_day[i] == 3){
      crd[i] <- 2.5
    } else if(rainy_day[i] >= 4){
      crd[i] <- 0.3
    } else{
      crd[i] <- 0
    }
    
    # GZ: the daily relative density of a GZ spore cloud
    gz[i] <- (GZ_INTERCEPT + (GZ_RHUM * wth$rhum[i]) + (GZ_CRD * crd[i])) ^ 2
    
    # GZ mean: 48-hour average value of GZ
    # This code can be replaced with rollmean (ex. tavg48)
    if(i == 1){
      gz_mean[i] <- gz[i]
    } else {
      gz_mean[i] <- mean(c(gz[i], gz[i-1]))
    }
    
    
    
    
    
    ##### _INFECTION MODULE: environmental factor #####
    
    # INF: the proportion of susceptible tissue likely to be infected at any time.
    # v0.7: revision for counting infection event
    
    if(i == 1){
      inf[i] <- 0
    } else if((rainy_day[i] >= 2 & wth$rhum48[i] >= INFECTION_RHUM_THRESHOLD_1) |
              ((wth$prcp[i-1] >= PRCP_COEF & wth$rhum[i-1] >= INFECTION_RHUM_THRESHOLD_2) & (wth$prcp[i] < PRCP_COEF & wth$rhum[i] >= INFECTION_RHUM_THRESHOLD_3)) |
              ((wth$prcp[i] >= PRCP_COEF & wth$rhum[i] >= INFECTION_RHUM_THRESHOLD_2) & (wth$prcp[i-1] < PRCP_COEF & wth$rhum[i-1] >= INFECTION_RHUM_THRESHOLD_3))){
      inf[i] <- INF_COEF_1 * exp(INF_COEF_2 * wth$tavg48[i])
    }else{
      inf[i] <- 0
    }
    # ##### _INFECTION MODULE: environmental factor #####
    # 
    # # INF: the proportion of susceptible tissue likely to be infected at any time.
    # # v0.7: revision for counting infection event
    # 
    # if(i == 1){
    #   inf[i] <- 0
    # } else {
    #   inf[i] <- INF_TAVG * ifelse(wth$tavg[i] > 30 | wth$tavg[i] < 10, 0, wth$tavg[i]) + INF_RHUM * ifelse(wth$rhum[i] < INF_RHUM_THRESHOLD, 0, wth$rhum[i]) + INF_PRCP * ifelse(wth$prcp[i] < 0.3, 0, wth$prcp[i]) + INF_INTERCEPT
    # }
    # if(i == 1){
    #   inf[i] <- 0
    # } else if(rainy_day[i] >= 3 | wth$rhum72[i] >= INFECTION_RHUM_THRESHOLD){
    #   inc_t=72
    #   inf[i] <- - 0.099 - 0.363*inc_t + 0.07808*wth$tavg72[i]*inc_t - 0.00591*wth$tavg72[i]^2*inc_t + 0.000199*wth$tavg72[i]^3*inc_t - 0.0000024*wth$tavg72[i]^4*inc_t
    # } else if(rainy_day[i] == 2 | wth$rhum48[i] >= INFECTION_RHUM_THRESHOLD){
    #   inc_t=48
    #   inf[i] <- - 0.099 - 0.363*inc_t + 0.07808*wth$tavg48[i]*inc_t - 0.00591*wth$tavg48[i]^2*inc_t + 0.000199*wth$tavg48[i]^3*inc_t - 0.0000024*wth$tavg48[i]^4*inc_t
    # } else if(rainy_day[i] == 1 | wth$rhum[i] >= INFECTION_RHUM_THRESHOLD){
    #   inc_t=24
    #   inf[i] <- - 0.099 - 0.363*inc_t + 0.07808*wth$tavg[i]*inc_t - 0.00591*wth$tavg[i]^2*inc_t + 0.000199*wth$tavg[i]^3*inc_t - 0.0000024*wth$tavg[i]^4*inc_t
    # } else {
    #   inf[i] <- 0
    # }
    
    # inf[i] <- (INF_TAVG_1 * exp(INF_TAVG_2 * wth$tavg48[i])) +
    #   (INF_RHUM_1 * (wth$rhum48[i] - INF_RHUM_2)^2 + INF_RHUM_3) +
    #   (INF_PRCP * (wth$prcp[i] >= 0.3))
    # inf[i] <- INF_TAVG * wth$tavg[i] + INF_RHUM * wth$rhum[i] + INF_PRCP * wth$prcp[i]
    # inf[i] <- TRH_COEF * ifelse(wth$tavg[i] >= 15 & wth$tavg[i] <= 30 & wth$rhum[i] >= 80, 1, 0)
    # } else if((rainy_day[i] >= 2 & wth$rhum48[i] >= INFECTION_RHUM_THRESHOLD_1) |
    #           ((wth$prcp[i-1] >= 0.3 & wth$rhum[i-1] >= INFECTION_RHUM_THRESHOLD_2) & (wth$prcp[i] < 0.3 & wth$rhum[i] >= INFECTION_RHUM_THRESHOLD_3)) |
    #           ((wth$prcp[i] >= 0.3 & wth$rhum[i] >= INFECTION_RHUM_THRESHOLD_2) & (wth$prcp[i-1] < 0.3 & wth$rhum[i-1] >= INFECTION_RHUM_THRESHOLD_3))){
    #   inf[i] <- INF_COEF_1 * exp(INF_COEF_2 * wth$tavg48[i])
    # }
  }
  
  
  # Many code above in for loop can be replaced with more effective form!!!!!!!!!!
  
  
  ##### _DATAFRAME PREPARATION #####
  
  head_present_doy <- which(head_daily_emerged != 0)
  anther_duration <- length(head_present_doy) + 10
  
  
  # Useful format for making empty dataframe!
  
  # Making empty dataframes
  list_empty <- list()
  for(i in head_present_doy){
    list_empty <- c(list_empty, list(rep(0,anther_duration)))
  }
  df_empty <- (do.call(cbind,list_empty)) %>%
    as.data.frame()
  
  df_head <- df_antext <- df_antext_diff <- df_anther <- df_empty
  
  head_name <- antext_name <- antext_diff_name <- anther_name <- empty_name <- c()
  for(i in head_present_doy){
    empty_name <- c(empty_name, paste0("empty_",i))
    head_name <- c(head_name, paste0("head_",i))
    antext_name <- c(antext_name, paste0("antext_",i))
    antext_diff_name <- c(antext_diff_name, paste0("antext_diff_",i))
    anther_name <- c(anther_name, paste0("anther_",i))
  }
  colnames(df_head) <- head_name
  colnames(df_antext) <- antext_name
  colnames(df_antext_diff) <- antext_diff_name
  colnames(df_anther) <- anther_name
  colnames(df_empty) <- empty_name
  
  
  ### do.call() and reduce() can be used for empty dataframe!!!
  # df_antext <- list_antext %>%
  #   reduce(cbind) %>%
  #   as.data.frame()
  # colnames(df_antext) <- antext_name
  
  
  
  
  
  ##### _DAILY EMERGED HEADS #####
  
  # Fill df_head
  for(i in 1:ncol(df_head)){
    # i=1
    for(j in i:nrow(df_head)){
      # j=1
      df_head[j, i] <- head_daily_emerged[i + head_present_doy[1] - 1]
    }
  }
  
  
  
  
  ##### _ANTHERS EXTRUSION #####
  
  # ANText: the daily rate of cumulative proportion of extruded anthers in a cohort of heads
  # ANText daily diff: the daily rate of proportion of extruded anthers
  
  # v0.9: The first anther extrusion timing was revised
  # Fill df_antext and df_antext_diff
  for(i in 1:ncol(df_antext)){
    # i=1
    # i=2
    for(j in (i + FIRST_ANTHER):nrow(df_antext)){
      # j=4
      # j=5
      df_antext[j, i] <- 1 - exp(-a[j + head_present_doy[1] - 1] * ((j - i - FIRST_ANTHER + 1) ^ b[j + head_present_doy[1] - 1]))
      
      # 3 days for first anther extrusion
      if(j < i + FIRST_ANTHER){
        df_antext_diff[j, i] <- 0
      } else if(j == i + FIRST_ANTHER){
        df_antext_diff[j, i] <- df_antext[j, i]
        antext_temp <- df_antext[j, i]
      } else if(df_antext[j, i] - antext_temp >= 0){
        df_antext_diff[j, i] <- df_antext[j, i] - antext_temp
        antext_temp <- df_antext[j, i]
      } else {
        df_antext_diff[j, i] <- 0
      }
      
      # antext > 0.99, antext = 1
      if(df_antext[j, i] > 0.99) {
        df_antext[j:nrow(df_antext), i] <- 1
        break
      }
    }
  }
  
  
  
  
  
  ##### _DAILY EXTRUDED ANTHERS #####
  
  
  # calculate df_anther
  df_anther <- df_head * df_antext_diff * QTY_ANTHERS_BY_HEADS
  
  
  
  
  ##### _LIVE ANTHERS #####
  
  anther_present_doy <- head_present_doy[1]:(head_present_doy[1] + anther_duration - 1)
  TOTAL_ANTHER <- QTY_PLANTS * QTY_HEADS_BY_PLANTS * QTY_ANTHERS_BY_HEADS 
  
  # anther column
  anther[anther_present_doy] <- rowSums(df_anther)
  
  # accumulative anther column
  anther_acc <- cumsum(anther)
  
  # Empty live anther column
  anther_live <- rep(0, length(anther))
  
  # v0.9: NEW calculation code containing specific value, minimum 2 ~ maximum 5
  # Fill live anther column: Anther's longevity
  for(i in anther_present_doy) {
    # i = anther_present_doy[5]
    # i = anther_present_doy[7]
    # i = anther_present_doy[26]
    if(i %in% anther_present_doy[1:(3 + MINIMUM_ANTHER_LONGEVITY)]){
      anther_live[i] <- anther_acc[i]
    } else if(wth$rsds[i-1] < ANTHER_RSDS_THRESHOLD){
      anther_live[i] <- anther[i] + anther_live[i-1]
      if(identical(wth$rsds[(i - MAXIMUM_ANTHER_LONGEVITY + 1):(i - 1)] < ANTHER_RSDS_THRESHOLD, rep(TRUE, MAXIMUM_ANTHER_LONGEVITY - 1))) {
        anther_live[i] <- anther_live[i] - anther[i - MAXIMUM_ANTHER_LONGEVITY]
      }
    } else if(wth$rsds[i - 1] >= ANTHER_RSDS_THRESHOLD) {
      anther_live[i] <- sum(anther[(i - MINIMUM_ANTHER_LONGEVITY + 1):i]) 
    }
  }
  
  # ant column
  ant <- anther_live / TOTAL_ANTHER
  
  # st column
  st <- ant
  st_max <- max(ant)
  
  # coefficients for post peak flowering infections
  for(i in anther_present_doy){
    # i = 130
    if (i > max(which(st == st_max))) {
      if(st[i] >= 0.01 & st[i] < 0.25) {
        # After peak flowering & ANT < 0.25
        st[i] <- 0.25
      } else if(st[i] < st_max & st[i] < 0.01){
        # After flowering
        # ST = 0.25 for next seven days
        st[i:(i+(ST_PERIOD_1-1))] <- 0.25
        # ST = 0.01 from eight to 14 days
        st[(i+ST_PERIOD_1):(i+(ST_PERIOD_1+ST_PERIOD_2-1))] <- 0.10
        break
      }
    }
  }
  
  # ant[which(ant >= 0.01 & ant < 0.25)]
  # st[which(ant >= 0.01 & ant < 0.25)]
  # sum((st == 0.25))
  # plot(ant[anther_present_doy])
  # plot(st[anther_present_doy])
  
  ant_mean <- rollmean(ant, k = 2, fill = ant[1], align = "right")
  # plot(ant_mean[anther_present_doy])
  
  st_mean <- rollmean(st, k = 2, fill = st[1], align = "right")
  # plot(st_mean[anther_present_doy])
  
  
  ##### _GIBBERELLA RISK MODULE #####
  
  inf <- ifelse(inf <=0, 0, inf)
  
  # GIB1 = ANT * INF
  gib1 <- ant * ifelse(inf <=0, 0, inf)
  
  # GIB2 = ANT * INF * GZ
  gib2 <- ant * inf * gz_mean
  
  # GIB3 = ST * INF
  gib3 <- st * inf
  
  # GIB4 = ST * INF * GZ
  gib4 <- st * inf * gz_mean * GIB_COEF
  
  
  # GIB % = SUM(GIB * 100)
  gib1_acc_percent <- cumsum(gib1) * 100
  gib2_acc_percent <- cumsum(gib2) * 100
  gib3_acc_percent <- cumsum(gib3) * 100
  gib4_acc_percent <- cumsum(gib4) * 100
  
  gib1_final <- gib1_acc_percent[length(wth$date)]
  gib2_final <- gib2_acc_percent[length(wth$date)]
  gib3_final <- gib3_acc_percent[length(wth$date)]
  gib4_final <- gib4_acc_percent[length(wth$date)]
  # 
  # gib1_final <- ifelse(gib1_final < 0, 0, ifelse(gib1_final > 100, 100, gib1_final))
  # gib2_final <- ifelse(gib2_final < 0, 0, ifelse(gib2_final > 100, 100, gib2_final))
  # gib3_final <- ifelse(gib3_final < 0, 0, ifelse(gib3_final > 100, 100, gib3_final))
  # gib4_final <- ifelse(gib4_final < 0, 0, ifelse(gib4_final > 100, 100, gib4_final))
  
  
  eq_calibration <- function(x) {GIB_COEF * x}
  
  gib1_inc <- eq_calibration(gib1_final)
  gib2_inc <- eq_calibration(gib2_final)
  gib3_inc <- eq_calibration(gib3_final)
  gib4_inc <- eq_calibration(gib4_final)
  
  # gib1_sev <- 0.7442 + 1.46 * gib1_final
  # gib2_sev <- 0.7442 + 1.46 * gib2_final
  # gib3_sev <- 0.7442 + 1.46 * gib3_final
  # gib4_sev <- 0.7442 + 1.46 * gib4_final
  
  
  # v0.9: anther and anther_live columns were added
  # Merge all the columns
  res <-
    cbind(
      wth,
      rainy_day,
      hng,
      anther,
      anther_live,
      ant,
      st,
      gz_mean,
      inf,
      gib1_acc_percent,
      gib1_final,
      gib1_inc,
      gib2_acc_percent,
      gib2_final,
      gib2_inc,
      gib3_acc_percent,
      gib3_final,
      gib3_inc,
      gib4_acc_percent,
      gib4_final,
      gib4_inc
    )
  
  # v0.5: check doy and weather variables of infection period 
  # v0.6: revised version, there was error about counting infection days. How to calculate data when incidence == 0?
  doy_inf <- yday(heading_starting_date) : max(yday(heading_starting_date), min(which(round(res$gib4_acc_percent, 2) == round(res$gib4_final, 2))))
  
  gib1_final <- ifelse(gib1_final < 0, 0, ifelse(gib1_final > 100, 100, gib1_final))
  gib2_final <- ifelse(gib2_final < 0, 0, ifelse(gib2_final > 100, 100, gib2_final))
  gib3_final <- ifelse(gib3_final < 0, 0, ifelse(gib3_final > 100, 100, gib3_final))
  gib4_final <- ifelse(gib4_final < 0, 0, ifelse(gib4_final > 100, 100, gib4_final))
  
  gib1_inc <- ifelse(gib1_inc < 0, 0, ifelse(gib1_inc > 100, 100, gib1_inc))
  gib2_inc <- ifelse(gib2_inc < 0, 0, ifelse(gib2_inc > 100, 100, gib2_inc))
  gib3_inc <- ifelse(gib3_inc < 0, 0, ifelse(gib3_inc > 100, 100, gib3_inc))
  gib4_inc <- ifelse(gib4_inc < 0, 0, ifelse(gib4_inc > 100, 100, gib4_inc))
  
  day_inf <- length(doy_inf)
  # day_inf <- ifelse(length(doy_inf) == 1, 0, length(doy_inf))
  tavg_inf <- mean(res[doy_inf,]$tavg, na.rm = T)
  prcp_inf <- sum(res[doy_inf,]$prcp, na.rm = T)
  rainy_day_inf <- sum(res[doy_inf,]$rainy_day != 0, na.rm = T)
  rhum_inf <- mean(res[doy_inf,]$rhum, na.rm = T)
  
  res <-
    cbind(
      wth,
      rainy_day,
      hng,
      anther,
      anther_live,
      ant,
      st,
      gz_mean,
      inf,
      gib1_acc_percent,
      gib1_final,
      gib1_inc,
      gib2_acc_percent,
      gib2_final,
      gib2_inc,
      gib3_acc_percent,
      gib3_final,
      gib3_inc,
      gib4_acc_percent,
      gib4_final,
      gib4_inc
    )
  
  
  
  
  
  # v0.9
  res <- res %>% 
    mutate(
      "heading_date" = heading_date_40,
      "heading_date_doy" = yday(heading_date_40),
      # "heading_starting_date" = heading_starting_date, 
      "flower_sampling_date" = flower_sampling_date, 
      "disease_sampling_date" = disease_sampling_date,
      "day_inf" = day_inf,
      "tavg_inf" = tavg_inf,
      "prcp_inf" = prcp_inf,
      "rainy_day_inf" = rainy_day_inf,
      "rhum_inf" = rhum_inf,
      'Gv' = Gv,
      'Thv' = Thv,
      'Lc' = Lc,
      'Av' = Av,
      'B' = B,
      'HNG_1' = HNG_1,
      'HNG_2' = HNG_2,
      'FIRST_ANTHER' = FIRST_ANTHER,
      'ANTHER_RSDS_THRESHOLD' = ANTHER_RSDS_THRESHOLD,
      'MINIMUM_ANTHER_LONGEVITY' = MINIMUM_ANTHER_LONGEVITY,
      'MAXIMUM_ANTHER_LONGEVITY' = MAXIMUM_ANTHER_LONGEVITY,
      'ST_PERIOD_1' = ST_PERIOD_1,
      'ST_PERIOD_2' = ST_PERIOD_2,
      'GZ_INTERCEPT' = GZ_INTERCEPT,
      'GZ_RHUM' = GZ_RHUM,
      'GZ_CRD' = GZ_CRD,
      'INFECTION_RHUM_THRESHOLD_1' = INFECTION_RHUM_THRESHOLD_1,
      'INFECTION_RHUM_THRESHOLD_2' = INFECTION_RHUM_THRESHOLD_2,
      'INFECTION_RHUM_THRESHOLD_3' = INFECTION_RHUM_THRESHOLD_3,
      # "INF_COEF_1" = INF_COEF_1,
      # "INF_COEF_2" = INF_COEF_2,
      'GIB_COEF' = GIB_COEF
    )
  
  return(res)
}




##### _______________ #####

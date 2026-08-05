### library 
# remotes::install_version("piecewiseSEM", version = "2.1.2")
library(piecewiseSEM)
library(readxl)
library(betareg)
### data
df_rslt_stat <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v4.xlsx") 
nrow(df_rslt_stat)

### group by loc -----------------------------------------------------------------------------------------------------
nrow(df_rslt_stat)
df_rslt_stat <- df_rslt_stat[complete.cases(df_rslt_stat),]
names(df_rslt_stat)

df_SEM <- df_rslt_stat %>% dplyr::select("incidence","Alternaria", "Epicoccum", "Hannaella", "Periconia", 
                                         #"Cladosporium", "Papiliotrema", 
                                         "temp_before_hd_sampling_10days",
                                         "rhum_before_hd_sampling_10days",
                                         # "prcp_before_hd_sampling_10days",
                                         "prcp_after_sampling",
                                         "temp_flower_sampling_10days",
                                         "rhum_flower_sampling_10days")
head(as.data.frame(df_SEM))
nrow(df_SEM)
### analysis
mb_pca <- prcomp(df_SEM[, c("Alternaria","Epicoccum","Hannaella",
                            "Periconia")],
                 scale. = TRUE)

summary(mb_pca)
df_SEM$MB_axis1 <- mb_pca$x[,1]
df_SEM$MB_axis2 <- mb_pca$x[,2]

###
df_SEM_std <- as.data.frame(scale(df_SEM))
head(df_SEM_std)
df_SEM_std$incidence <- df_rslt_stat$incidence
df_SEM_std$incidence <- ifelse(df_SEM_std$incidence == 0, 0.00001, df_SEM_std$incidence)
df_SEM_std$incidence <- ifelse(df_SEM_std$incidence == 1, 0.9999, df_SEM_std$incidence)
nrow(df_SEM_std)
###########################################################################################################
library(piecewiseSEM)
library(betareg)

# 기존 모델 그대로 사용
Before_heading_m_MB1 <- lm(MB_axis1 ~ temp_before_hd_sampling_10days + rhum_before_hd_sampling_10days,
                           data = df_SEM_std)

Before_heading_m_MB2 <- lm(MB_axis2 ~ temp_before_hd_sampling_10days + rhum_before_hd_sampling_10days,
                           data = df_SEM_std)

Before_heading_m_inc <- glm(
  incidence ~ MB_axis1 + MB_axis2 +
    temp_before_hd_sampling_10days +
    rhum_before_hd_sampling_10days,
  data = df_SEM_std,
  family = quasibinomial(link = "logit")
)

#???? 핵심: piecewise SEM 구성
sem_before <- psem(
  Before_heading_m_MB1,
  Before_heading_m_MB2,
  Before_heading_m_inc
)

# ???? 1) 전체 모델 적합도 (가장 중요)
summary(sem_before)

# ???? 2) d-separation test (independence claims)
dSep(sem_before)

# ???? 3) 경로계수 (표준화 포함)
coefs(sem_before, standardize = "scale")

# ???? 4) R2
rsquared(sem_before)

# pred_Before_heading_m_inc <- predict(Before_heading_m_inc, type = "response")
# 1 - (Before_heading_m_inc$deviance / Before_heading_m_inc$null.deviance)

####################
Flowering_m_MB1 <- lm(MB_axis1 ~ temp_flower_sampling_10days + 
                        rhum_flower_sampling_10days,
                      data = df_SEM_std)

Flowering_m_MB2 <- lm(
  MB_axis2 ~ temp_flower_sampling_10days +
    rhum_flower_sampling_10days, 
    
  data = df_SEM_std
)

Flowering_m_inc <- glm(incidence ~ MB_axis1 + MB_axis2 +
                             temp_flower_sampling_10days +
                             rhum_flower_sampling_10days,
                       data = df_SEM_std,
                       family = quasibinomial(link = "logit")
                       )
sem_flowering <- psem(
  Flowering_m_MB1,
  Flowering_m_MB2,
  Flowering_m_inc
)

summary(sem_flowering)
dSep(sem_flowering)
coefs(sem_flowering, standardize = "scale")
rsquared(sem_flowering)
dSep(sem_flowering)

### Graph
library(DiagrammeR)

grViz("
digraph SEM_Before {

  graph [layout = dot, rankdir = LR]

  node [shape = rectangle, style = filled, fillcolor = lightgrey]

  temp [label = 'Temperature']
  rhum [label = 'Humidity']
  MB1 [label = 'MB_axis1']
  MB2 [label = 'MB_axis2']
  inc [label = 'Disease incidence']

  # edges
  temp -> MB1 [label = '-0.64***']
  rhum -> MB1 [label = '-0.06']

  temp -> MB2 [label = '0.06']
  rhum -> MB2 [label = '-0.23']

  MB1 -> inc [label = '-0.39**']
  MB2 -> inc [label = '-0.33**']

  temp -> inc [label = '0.83***']
  rhum -> inc [label = '0.24']

}
")


grViz("
digraph SEM_Flowering {

  graph [layout = dot, rankdir = LR]

  node [shape = rectangle, style = filled, fillcolor = lightgrey]

  temp [label = 'Temperature']
  rhum [label = 'Humidity']
  MB1 [label = 'MB_axis1']
  MB2 [label = 'MB_axis2']
  inc [label = 'Disease incidence']

  # edges
  temp -> MB1 [label = '-0.23*']
  rhum -> MB1 [label = '-0.46***']

  temp -> MB2 [label = '-0.27**']
  rhum -> MB2 [label = '-0.13']

  MB1 -> inc [label = '-0.21*']
  MB2 -> inc [label = '0.07']

  temp -> inc [label = '0.62***']
  rhum -> inc [label = '0.99***']

}
")

## indirect effect table
library(dplyr)

indirect_before <- data.frame(
  Path = c(
    "Temp → MB1 → Incidence",
    "Temp → MB2 → Incidence",
    "Rhum → MB1 → Incidence",
    "Rhum → MB2 → Incidence"
  ),
  Effect = c(
    -0.6355 * -0.3853,
    0.0626 * -0.3330,
    -0.0648 * -0.3853,
    -0.2256 * -0.3330
  )
)

indirect_before

indirect_flowering <- data.frame(
  Path = c(
    "Temp → MB1 → Incidence",
    "Temp → MB1 → MB2 → Incidence",
    "Temp → MB2 → Incidence",
    "Rhum → MB1 → Incidence",
    "Rhum → MB1 → MB2 → Incidence",
    "Rhum → MB2 → Incidence"
  ),
  Effect = c(
    -0.2270 * -0.2087,
    -0.2270 * -0.3669 * 0.0683,
    -0.3520 * 0.0683,
    -0.4638 * -0.2087,
    -0.4638 * -0.3669 * 0.0683,
    -0.3079 * 0.0683
  )
)

indirect_flowering

indirect_before %>%
  mutate(Stage = "Before") %>%
  bind_rows(indirect_flowering %>% mutate(Stage = "Flowering"))

total_effect <- data.frame(
  Variable = c("Temperature", "Humidity"),
  Direct_before = c(0.8282, 0.2350),
  Indirect_before = c(
    -0.6355 * -0.3853 + 0.0626 * -0.3330,
    -0.0648 * -0.3853 + -0.2256 * -0.3330
  ),
  Direct_flowering = c(0.6217, 0.9914),
  Indirect_flowering = c(
    -0.2270 * -0.2087 + -0.2270 * -0.3669 * 0.0683 + -0.3520 * 0.0683,
    -0.4638 * -0.2087 + -0.4638 * -0.3669 * 0.0683 + -0.3079 * 0.0683
  )
)

total_effect

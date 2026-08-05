library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(emmeans)
library(multcomp)
library(multcompView)

vivo_df <- readxl::read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/길항test_결과_사진/wheat_head/wheat_head_260702_NH.xlsx")
vivo_df_sel <- vivo_df %>% dplyr::select("treatment", "bio_rep", "total_spikelet", "disease_spikelet", "severity")
names(vivo_df_sel)[which(names(vivo_df_sel) == "bio_rep")] <- "rep"

# treatment 순서 지정 (원하는 순서로 변경 가능)
vivo_df_sel$treatment <- factor(
  vivo_df_sel$treatment,
  levels = c("PC", "Alternaria_F.g", "Epicoccum_F.g", 
             "NC","only_Alternaria","only_Epicoccum")
)

# 평균과 SD 계산
summary_df <- vivo_df_sel %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(severity),
    sd = sd(severity),
    n = n(),
    .groups = "drop"
  )

# one-way ANOVA
fit_aov <- aov(severity ~ treatment, data = vivo_df_sel)
anova(fit_aov)
emm_aov <- emmeans(fit_aov, "treatment")
# 
letters_aov <- cld(
  emm_aov,
  adjust = "tukey",
  Letters = letters
)
emm_vivo <- emmeans(fit_aov, "treatment")
vivo_res <- contrast(
  emm_vivo,
  method = "trt.vs.ctrl",
  ref = "PC",
  adjust = "dunnett"
)

summary_df$letter <-
  letters_aov$.group[
    match(summary_df$treatment,
          letters_aov$treatment)
  ]

summary_df <- summary_df %>% mutate(letter = trimws(letter))
summary_df$mean <- summary_df$mean*100
summary_df$sd <- summary_df$sd*100
wheat_head_graph <- ggplot(summary_df, aes(x = treatment, y = mean)) +
  
  geom_col(width = 0.7,
           fill = "#707070",
           color = "black") +
  
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.15,
                linewidth = 0.8) +
  
  geom_quasirandom(
    data = vivo_df_sel,
    aes(x = treatment,
        y = severity*100),
    inherit.aes = FALSE,
    width = 0.15,
    size = 2.5,
    shape = 21,
    fill = alpha("white", 0.6),
    color = "#707070",
    stroke = 0.8
  ) +
  
  coord_cartesian(ylim = c(0, 120)) +
  
  labs(
    x = NULL,
    y = "FHB severity\n"
  ) +
  
  theme_classic(base_size = 14)+
  theme(
    axis.text.x = element_text(
      angle = 90,            # 90도 회전
      vjust = 0.5,
      hjust = 1
    )
  ) 
wheat_head_graph
summary_df

###########
vitro_df <- readxl::read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/길항test_결과_사진/test_results/in_vitro_test_260703.xlsx")
vitro_df_sel <- vitro_df %>% dplyr::select("treatment", "teq_rep","bio_rep", "Fusarium_growth", "control_growth", "inhibition_rate", "inhibition_zone")
names(vitro_df_sel)[which(names(vitro_df_sel) == "bio_rep")] <- "rep"

# treatment 순서 지정 (원하는 순서로 변경 가능)
vitro_df_sel$treatment <- factor(
  vitro_df_sel$treatment,
  levels = c("control", "KJ5_1", "KJ1_E_inh1")
)

vitro_plot_df <- vitro_df_sel %>%
  filter(treatment != "control") %>%
  mutate(
    treatment = factor(treatment,
                       levels = c("KJ5_1", "KJ1_E_inh1"))
  )

# 평균과 SD 계산
vitro_summary_rate <- vitro_plot_df %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(inhibition_rate),
    sd = sd(inhibition_rate),
    .groups = "drop"
  )

inh_rate_graph <- ggplot(vitro_summary_rate,
       aes(x = treatment, y = mean)) +
  
  geom_col(
    fill = "#707070",
    color = "black",
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(ymin = mean - sd,
        ymax = mean + sd),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_quasirandom(
    data = vitro_plot_df,
    aes(y = inhibition_rate,
    # inherit.aes = FALSE,
    x = treatment),
    ## ★ 수정
    size = 2.5,
    shape = 21,
    fill = alpha("white", 0.6),
    color = "#707070",
    stroke = 0.8
  ) +
  
  labs(
    x = NULL,
    y = "Inhibition rate (%)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )
inh_rate_graph
t.test(
  Fusarium_growth ~ treatment,
  data = subset(vitro_df_sel,
                treatment %in% c("control","KJ5_1"))
)

t.test(
  Fusarium_growth ~ treatment,
  data = subset(vitro_df_sel,
                treatment %in% c("control","KJ1_E_inh1"))
)

###------
vitro_summary_zone <- vitro_plot_df %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(inhibition_zone),
    sd = sd(inhibition_zone),
    .groups = "drop"
  )

inh_zone_graph <- ggplot(vitro_summary_zone,
       aes(x = treatment, y = mean)) +
  
  geom_col(
    fill = "#707070",
    color = "black",
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(ymin = mean - sd,
        ymax = mean + sd),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_quasirandom(
    data = vitro_plot_df,
    aes(y = inhibition_zone,x = treatment),
    ## ★ 수정
    size = 2.5,
    shape = 21,
    fill = alpha("white", 0.6),
    color = "#707070",
    stroke = 0.8
  ) +
  
  labs(
    x = NULL,
    y = "Inhibition zone (mm)"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )
inh_zone_graph
#
library(dplyr)
library(ggplot2)
library(ggbeeswarm)

## Control 포함
growth_df <- vitro_df_sel %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c("control", "KJ5_1", "KJ1_E_inh1")
    )
  )

## 평균 ± SD
growth_summary <- growth_df %>%
  group_by(treatment) %>%
  summarise(
    mean = mean(Fusarium_growth),
    sd = sd(Fusarium_growth),
    .groups = "drop"
  )

growth_graph <- ggplot(growth_summary,
                       aes(x = treatment, y = mean)) +
  
  geom_col(
    fill = "#707070",
    color = "black",
    width = 0.7
  ) +
  
  geom_errorbar(
    aes(ymin = mean - sd,
        ymax = mean + sd),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_quasirandom(
    data = growth_df,
    aes(
      x = treatment,
      y = Fusarium_growth
    ),
    inherit.aes = FALSE,
    width = 0.12,
    ## ★ 수정
    size = 2.5,
    shape = 21,
    fill = alpha("white", 0.6),
    color = "#707070",
    stroke = 0.8
  ) +
  
  labs(
    x = NULL,
    y = "Fusarium growth (mm)"
  ) +
  coord_cartesian(ylim = c(37, 55)) +
  theme_classic(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

growth_graph

library(emmeans)

fit <- aov(Fusarium_growth ~ treatment, data = growth_df)

## ANOVA
anova(fit)

## Dunnett (control과만 비교)
emm <- emmeans(fit, "treatment")

growth_res <- contrast(
  emm,
  method = "trt.vs.ctrl",
  ref = "control",
  adjust = "dunnett"
)
growth_res


### save ----------------------------------------------------------------------
wheat_head_graph
growth_graph
inh_rate_graph
inh_zone_graph

ggsave(filename = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/6. antifungal_test/wheat_head_graph.png", 
       plot = wheat_head_graph, 
       width = 10/2.5, 
       height = 10.5/2.5)

ggsave(filename = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/6. antifungal_test/growth_graph.png", 
       plot = growth_graph, 
       width = 4.5/2.5, 
       height = 7.6/2.5)

ggsave(filename = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/6. antifungal_test/inh_rate_graph.png", 
       plot = inh_rate_graph, 
       width = 4.5/2.5, 
       height = 7.6/2.5)

ggsave(filename = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/6. antifungal_test/inh_zone_graph.png", 
       plot = inh_zone_graph, 
       width = 4.5/2.5, 
       height = 7.6/2.5)


# stat
vivo_res # in vivo result
summary_df #in vivo result
growth_res #growth result
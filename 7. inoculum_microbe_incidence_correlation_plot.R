library(Hmisc)
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)
df_rslt_stat2 <- read_xlsx("D:/Microbiome/000.data/FHB_microbiome_data/FHB_intensity_and_wth_data/For_phyloseq/2025/merge_meta_MB_GIBSIM_group_by_loc_only_flowering_ITSfull_v4.xlsx") 
df_rslt_stat2 <- df_rslt_stat2[complete.cases(df_rslt_stat2),]
bf_heading_for_cor <- read_xlsx(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/" , 
                                           "2025_v12","/bf_heading_for_cor_",
                                           "2025_v12",".xlsx"))
bf_heading_for_cor_ino <- bf_heading_for_cor %>% dplyr::select(id, year, bf_heading_Ino_Manstretta_mean, bf_heading_Ino_Rossi_mean, bf_heading_Ino_GIBSIM_mean)
names(bf_heading_for_cor_ino)[which(names(bf_heading_for_cor_ino) == "id")] <- "wth_ID"
df_rslt_stat3 <- left_join(df_rslt_stat2, bf_heading_for_cor_ino, by = c("wth_ID", "year"))


df_for_more_analysis3 <- df_rslt_stat3 %>% dplyr::select("Alternaria", "Cladosporium", "Epicoccum", "Hannaella", 
                                                        "Papiliotrema","Periconia", "incidence", 
                                                        "temp_before_hd", "rhum_before_hd", "prcp_before_hd", 
                                                        "temp_hd", "rhum_hd", "prcp_hd", 
                                                        "temp_before_hd_sampling_10days", "rhum_before_hd_sampling_10days", "prcp_before_hd_sampling_10days", 
                                                        "temp_hd_sampling_10days", "rhum_hd_sampling_10days", "prcp_hd_sampling_10days", 
                                                        "bf_heading_Ino_Manstretta_mean", "bf_heading_Ino_Rossi_mean", "bf_heading_Ino_GIBSIM_mean")


# bio_vars <- c("Alternaria","Cladosporium","Epicoccum",
#               "Hannaella","Papiliotrema","Periconia","incidence")

incidence_Var <- "incidence"

bio_vars <- c("Alternaria","Epicoccum",
              "Hannaella","Periconia")

weather_vars <- c("temp_before_hd_sampling_10days","rhum_before_hd_sampling_10days","prcp_before_hd_sampling_10days") #weather_vars <- c("temp_before_hd","rhum_before_hd","prcp_before_hd","temp_hd","rhum_hd","prcp_hd")

ino_vars <- c("bf_heading_Ino_Manstretta_mean", "bf_heading_Ino_Rossi_mean", "bf_heading_Ino_GIBSIM_mean")

###
cor_res <- rcorr(as.matrix(sub_df), type = "spearman")
cor_mat <- cor_res$r
p_mat   <- cor_res$P

all_vars <- c(incidence_Var, bio_vars, ino_vars)

edges <- expand.grid(
  from = all_vars,
  to   = all_vars,
  stringsAsFactors = FALSE
)

# 자기 자신 제거
edges <- edges[edges$from != edges$to, ]

# 중복 제거 (A-B, B-A 중 하나만 유지)
edges <- edges %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(from, to)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

edges$r <- mapply(function(x, y) cor_mat[x, y],
                  edges$from, edges$to)

edges$p <- mapply(function(x, y) p_mat[x, y],
                  edges$from, edges$to)

edges <- edges %>%
  filter(!is.na(r), !is.na(p)) %>%
  filter(p < 0.05) %>%        # 유의한 것만
  filter(r <= -0.3 | r >= 0.3) %>%        # 유의한 것만
  mutate(
    weight = abs(r),
    sign   = ifelse(r > 0, "positive", "negative")
  )

nodes <- data.frame(
  name = all_vars,
  type = c(
    "Incidence",
    rep("Biological", length(bio_vars)),
    rep("Inoculum", length(ino_vars))
  ),
  stringsAsFactors = FALSE
)

inc_n <- 1
bio_n <- length(bio_vars)
ino_n <- length(ino_vars)

nodes$x <- c(
  rep(-0.5, inc_n),
  rep(0, bio_n),
  rep( 0.5, ino_n)
)

nodes$y <- c(
  0,
  c(-0.5, -0.3, 0.3, 0.5),
  seq(-0.3, 0.3, length.out = ino_n)
)

nodes$x <- c(-0.3, -0.5, -0.45, -0.3, 0, 0.3, 0.5, 0.3 ) #Incidence, Alt, Epi, Han, Per, Man, Rossi, GIBSIM
nodes$y <- c(0.3, 0, -0.25, -0.43, -0.5, 0.3, 0, -0.43)

nodes$hjust <- 0.5 #ifelse(nodes$x < 0, 1, 0)

nodes$name <- str_replace(nodes$name, "bf_heading_Ino_", "")
nodes$name <- str_replace(nodes$name, "_mean", "")

edges$from <- str_replace(edges$from, "bf_heading_Ino_", "")
edges$from <- str_replace(edges$from, "_mean", "")
edges$to <- str_replace(edges$to, "bf_heading_Ino_", "")
edges$to <- str_replace(edges$to, "_mean", "")

net <- tbl_graph(nodes = nodes,
                 edges = edges,
                 directed = FALSE)

p <- ggraph(net, layout = "manual",
            x = nodes$x,
            y = nodes$y) +
  
  geom_edge_link(aes(width = weight,
                     color = sign,
                     alpha = weight)) +
  
  scale_edge_width(range = c(0.5, 5)) +
  scale_edge_alpha(range = c(0.5, 0.9), guide = "none") +
  
  scale_edge_color_manual(
    values = c(
      positive = "#B2182B",
      negative = "#2166AC"
    )
  ) +
  
  geom_node_point(aes(fill = type),
                  size = 10,
                  shape = 21,
                  color = "white") +
  
  scale_fill_manual(
    values = c(
      Incidence  = "#000000",
      Biological = "#2e8b57",
      Inoculum   = "orange"
    )
  ) +
  
  geom_node_text(aes(label = name, hjust = 0.5, vjust = -1.1),
                 size = 4.5,
                 fontface = "bold") +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

p

      # ### choose the data for analysis
      # # chosen_wth_var <- weather_vars #weather_vars_sampling
      # names(df_for_more_analysis3)
      # sub_df <- df_for_more_analysis3[, c(incidence_Var,bio_vars, chosen_wth_var, ino_vars)]
      # head(as.data.frame(sub_df))
      # # as.data.frame(bf_heading_for_cor[which(bf_heading_for_cor$id_year == "wth_067_2025"),])
      # # sub_df[which(sub_df$Alternaria == bf_heading_for_cor$Alternaria[which(bf_heading_for_cor$id_year == "wth_067_2025")]) ,]
      # 
      # cor_res <- rcorr(as.matrix(sub_df), type="spearman")
      # cor_mat <- cor_res$r
      # p_mat   <- cor_res$P
      # 
      # edges <- expand.grid(from = bio_vars,
      #                      to   = chosen_wth_var,
      #                      stringsAsFactors = FALSE)
      # 
      # edges$r <- mapply(function(x,y) cor_mat[x,y],
      #                   edges$from, edges$to)
      # 
      # edges$p <- mapply(function(x,y) p_mat[x,y],
      #                   edges$from, edges$to)
      # 
      # edges <- edges %>%
      #   filter(p < 0.05) %>%
      #   mutate(weight = abs(r),
      #          sign = ifelse(r > 0, "positive", "negative"))
      # 
      # nodes <- data.frame(
      #   name = c(bio_vars, chosen_wth_var),
      #   type = c(rep("Biological", length(bio_vars)),
      #            rep("Weather", length(chosen_wth_var))),
      #   stringsAsFactors = FALSE
      # )
      # 
      # bio_n <- length(bio_vars)
      # wea_n <- length(chosen_wth_var)
      # 
      # nodes$x <- c(rep(-1, bio_n),
      #              rep( 1, wea_n))
      # 
      # nodes$y <- c(seq(-1, 1, length.out = bio_n),
      #              seq(-1, 1, length.out = wea_n))
      # nodes$hjust <- ifelse(nodes$type %in% c("Biological","Incidence"), 1, 0)
      # nodes$type <- ifelse(nodes$name == "incidence", "Incidence", nodes$type)
      # 
      # net <- tbl_graph(nodes = nodes,
      #                  edges = edges,
      #                  directed = FALSE) %>%
      #   activate(nodes) %>%
      #   mutate(degree = centrality_degree())
      # 
      # 
      # 
      # p <- ggraph(net, layout = "manual",
      #             x = nodes$x,
      #             y = nodes$y) +
      #   
      #   # edges
      #   geom_edge_link(aes(width = weight,
      #                      color = sign,
      #                      alpha = weight),
      #                  show.legend = TRUE) +
      #   
      #   scale_edge_width(range = c(0.5, 2.5),
      #                    name = "|Correlation|") +
      #   
      #   scale_edge_alpha(range = c(0.4, 0.9),
      #                    guide = "none") +
      #   
      #   scale_edge_color_manual(
      #     values = c(positive = "#B2182B",
      #                negative = "#2166AC"),
      #     name = "Correlation sign"
      #   ) +
      #   
      #   # nodes
      #   # geom_node_point(aes(size = degree,
      #   #                     fill = type),
      #   #                 shape = 21,
      #   #                 color = "black",
      #   #                 stroke = 1) +
      #   # 
      #   # scale_size(range = c(5,10),
      #   #            name = "Node degree") +
      #   
      # geom_node_point(aes(size = 10, fill = type),
      #                 shape = 21,
      #                 color = "white",
      #                 stroke = 1) +
      #   
      #   scale_fill_manual(
      #     values = c(Biological = "#ffa500",
      #                Weather = "#2e8b57", 
      #                Incidence = "#404040"),
      #     name = "Node type"
      #   ) +
      #   
      #   # labels
      #   # geom_node_text(aes(label = name),
      #   #                repel = TRUE,
      #   #                size = 4,
      #   #                fontface = "bold") +
      #   
      #   geom_node_text(aes(label = name,
      #                      hjust = hjust),
      #                  size = 4.5,
      #                  fontface = "bold") +
      #   
      #   theme_void() +
      #   theme(
      #     legend.position = "right",
      #     legend.title = element_text(size=12, face="bold"),
      #     legend.text  = element_text(size=11)
      #   )
      # p
# # PDF (vector)
# ggsave("FHB_weather_microbe_bipartite.pdf",
#        plot = p,
#        width = 8,
#        height = 6)

# High-res PNG
# ggsave("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/2025_v11/response_curve/FHB_weather_microbe_bipartite_no_name.png",
#        plot = p,
#        width = 5,
#        height = 5,
#        dpi = 400)


ggsave(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/",
              "2025_v12",
              "/response_curve/FHB_weather_microbe_bipartite_name (spearman).png"),
       plot = p,
       width = 8,
       height = 6,
       dpi = 400)

# no name
np <- ggraph(net, layout = "manual",
            x = nodes$x,
            y = nodes$y) +
  
  geom_edge_link(aes(width = weight,
                     color = sign,
                     alpha = weight)) +
  
  scale_edge_width(range = c(0.5, 5)) +
  scale_edge_alpha(range = c(0.5, 0.9), guide = "none") +
  
  scale_edge_color_manual(
    values = c(
      positive = "#B2182B",
      negative = "#2166AC"
    )
  ) +
  
  geom_node_point(aes(fill = type),
                  size = 10,
                  shape = 21,
                  color = "white") +
  
  scale_fill_manual(
    values = c(
      Incidence  = "#000000",
      Biological = "#2e8b57",
      Inoculum   = "orange"
    )
  ) +
  
  # geom_node_text(aes(label = name, hjust = 0.5, vjust = -1.1),
  #                size = 4.5,
  #                fontface = "bold") +
  
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  )

np
ggsave(paste0("D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/3. obs_pred_graph/2025/",
              "2025_v12",
              "/response_curve/FHB_weather_microbe_bipartite_no_name (spearman).png"),
       plot = np,
       width = 7,
       height = 5,
       dpi = 400)


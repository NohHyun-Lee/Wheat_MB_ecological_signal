isolated_Fungi_df = data.frame("Genus" = c("Alternaria","Epicoccum","Cladosporium",
                                "Sporobolomyces", "Papiliotrema","Penicillium", 
                                "Chaetomium", "Colletotrichum", "Eutypella", 
                                "Moesziomyces", "Phaeobotryon","Vishniacozyma"), 
                    "numbers" = c(43/109,29/109,15/109,
                                  4/109,
                                  3/109,3/109,
                                  1/109,1/109,1/109,1/109,1/109,1/109))

isolated_Fungi_df$count <- c(43,29,15,
                  4,
                  3,3,
                  1,1,1,1,1,1)

isolated_Fungi_graph <- ggplot(
  isolated_Fungi_df,
  aes(
    x = reorder(Genus, count),
    y = count
  )
) +  labs(title = "", 
          x = "Genus\n", 
          y = "\nNumber of isolates") + 
  
  geom_col(
    width = 0.8,
    fill = "#003366",
    color = "white"
  ) +
  coord_flip() +
  geom_text(
    aes(label = count),
    hjust = -0.11,
    size = 5
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_blank()
  ) + 
  theme(
    plot.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.title.y = element_text(size = 15, hjust = 0.5, face = 'bold'),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.4,
                               size = 15, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 15, face = 'bold', color = 'black')
  )
isolated_Fungi_graph


ggsave(filename = "D:/Microbiome/000.codes/FHB_microbiome_code/FHB_analysis/2025/MB_code_2025/Output/4. isolated_fungi/isolated_Fungi_graph.png",
       plot = isolated_Fungi_graph, 
       width = 7.5, 
       height = 8)

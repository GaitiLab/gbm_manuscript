# Code to reproduce Figure 2b

library(dplyr)
library(ggplot2)
library(ggrepel)

output_dir <- ""

# Creating directories needed for outputs 
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plot_dir <- paste0(output_dir, "/path_to_output_directory")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

interactions <- read.csv("interactions_summary.csv")
interactions <- interactions %>%  
  filter(Region == "PT") %>% 
  filter(source_target == "Glutamatergic__Invasive-high OPC/NPC1" | source_target == "Invasive-high OPC/NPC1__Glutamatergic") 

interactions <- interactions %>% 
  filter(Interaction_type != "NA") %>% 
  group_by(Interaction_type) %>% 
  summarise(n = n())

interactions$rank <- rank(-interactions$n, ties.method = "first")
interactions <- interactions %>% arrange(rank)
interactions$pathway_label <- ifelse(interactions$rank <= 4, interactions$Interaction_type, NA)

ggplot(interactions, aes(x=rank, y=n)) +
  geom_point(color = "black", stroke = 1, shape = 21, size = 5) +
  geom_text_repel(aes(label = pathway_label), nudge_x = 5, nudge_y = 2, na.rm = TRUE) +
  theme_classic()
ggsave("cci_pathway_rank_new.pdf", path = plot_dir)

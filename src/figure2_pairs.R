# Below is the R script for figures 2c and 2d. 


library(tidyverse)
library(readxl)
library(here)

# Import data from (Supp Table 6) Landscape of digenic dependencies in uveal melanoma 

Supplementary_Table_6 <- read_excel(here("data/Supplementary Table 6.xlsx"), 
                                    skip =2)

hit_1 <- Supplementary_Table_6 |> 
         filter(is_bassik_hit == 1)

hit_2 <- hit_1 |> 
        filter(!targetA__is_single_depleted == 1, 
               !targetB__is_single_depleted == 1)

filtered_data <- hit_2 |>
  group_by(gene_pair_id) |>
  filter(n() >= 6) |>
  ungroup()

# 160 pairs in "filtered_data" and also in 230820_PairedScreen_DotPlot_Rinput.csv

# Filtered_data is identical to 230820_PairedScreen_DotPlot_Rinput.csv with 
# the exception that in 230820_PairedScreen_DotPlot_Rinput.csv
# the previous hits are indicated by a colour 
# Also column N in 230820_PairedScreen_DotPlot_Rinput.csv.

PairedResults_6lines_SinglesEssentialsRemoved <- read_csv("data/230820_PairedScreen_DotPlot_Rinput.csv")

jpeg(here("results/plots/2024_PairedDotPlot_final.jpeg"), units = "in", width = 4.8, height = 7, res = 600)
DotPlotFinal <- ggplot(PairedResults_6lines_SinglesEssentialsRemoved) +
  geom_point(size = 2, aes(x = mean_norm_gi,
                           y = reorder(sorted_gene_pair, 
                           mean_GI, descending = FALSE), color = colour, 
                           dotsize = 4))
DotPlotFinal <- DotPlotFinal + 
                theme_classic() + 
                theme(axis.line.x = 
                element_line(linewidth = 0.2, color = "black")) + 
                theme(axis.line.y = 
                element_line(linewidth = 0.2, color = "black")) + 
                scale_color_manual(values = c("#217CA3", "grey"), labels = c("Previous SL screen hit", "Previous SL screen non-hit")) + 
                labs(x = "GI score", y = "") + 
                theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) + 
                theme(legend.title = element_blank()) + 
                guides(colour = guide_legend(nrow = 2)) + 
                theme(axis.text.y = element_text(face = "italic")) + 
                scale_y_discrete(limits = rev) + 
                theme(text = element_text(family = "Helvetica", size = 14))
DotPlotFinal
dev.off()

######### 25th June 24
SinglesResults_ParalogPairDerived <- read.table(here("data/common_SL_hit_gi_scores.csv"),
                     header = TRUE, sep = ",") 
jpeg(here("results/plots","2024_SinglesDotPlotParalogs25th_6_24(AHCYCL1_fix).jpeg"), 
    units = "in", width = 5.1, height = 7, res = 600)

DotPlot_SINGLESparalogs <- ggplot(SinglesResults_ParalogPairDerived) +
  geom_point(size = 2.2, aes(x = X921, 
  y = (reorder(GENE, average_LFC_all, decreasing = TRUE)), 
  color = PARALOG_BINARY), size = 2) + 
  theme_classic() + 
  theme(axis.line.x = element_line(linewidth = 0.2, color = "black")) + 
  theme(axis.line.y = element_line(linewidth = 0.2, color = "black")) + 
  scale_color_manual(values = c("grey", "#C52930"), 
  labels = c("Member of previous SL screen non-hit", 
             "Member of previous SL screen hit")) + 
  labs(x = bquote("Normalised log2 fold change"), y = "") + 
  theme(legend.position = "bottom", 
        legend.box = "vertical", 
        legend.margin = margin()) + 
  theme(legend.title = element_blank()) + 
  guides(colour = guide_legend(nrow = 2, reverse = TRUE)) + 
  theme(axis.text.y = element_text(face = "italic")) + 
  theme(text = element_text(family = "Helvetica", size = 14))
print(DotPlot_SINGLESparalogs)
dev.off()

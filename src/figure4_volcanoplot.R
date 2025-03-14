# 6th_March_24. redraw volcano plot to switch axes and change colours. The figure is for Figure 4 and the data is from Supp table 9.
# Updated July 2024

library(tidyverse)
library(ggtext)
library(plotly)
library(ggrepel)
library(here)

# Import data - this is supp data table 9 from the originally submitted paper to NG. I trimmed the heading to help the import and converted to a .csv
Supp_Table_9 <- read_csv(here("data/Supplementary Table 9.csv")) 

# Selected required columns and add gene status column (upregulated, down or non significant)
Supp_Table_9_sel <- Supp_Table_9 |>
                     mutate(gene_status = case_when(

  # If the gene's logFC is equal to or higher than 2 AND its adjusted p-value is below the threshold, classify this
  # gene as upregulated.
  log2FC >= 1.83 & Padj <= 0.01 ~ "Upregulated",

  # If the gene's logFC is equal to or lower than -2 AND its adjusted p-value is below the threshold, classify this
  # gene as downregulated.
  log2FC <= -1.83 & Padj <= 0.01 ~ "Downregulated",

  # If the gene's absolute logFC value (module) is between -1.83 and 1.83, OR its adjusted p-value is above the threshold,
  # classify that gene as non-significant.
  abs(log2FC) < 1.83 | Padj > 0.01 ~ "Non-significant"
))

# Select only the 4 columns to plot
columns_to_graph <- Supp_Table_9_sel |>
                    select(Genes, log2FC, `log10(1/Padj)`, gene_status)

# Plot volcano
P <- columns_to_graph |>
  ggplot(aes(x = log2FC, y = `log10(1/Padj)`, color = gene_status)) +
  geom_point(shape = 19) +
  theme_classic() +
  theme(
    axis.title.x = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 13, face = "bold", color = "black"),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.position = "right",
    legend.justification = "top"
  ) +
  scale_x_reverse() +
  scale_color_manual(values = c("#2E8B57", "#808080", "#226E9C"), name = "Gene status") +
  geom_vline(xintercept = c(-1.83, 1.83), linetype = "dashed") +
  geom_text_repel(
    data = columns_to_graph[columns_to_graph$Genes %in% c("SOX10", "MITF", "CDS2", 
                                                          "CDIPT", "RASGRP3", "TFAP2A", 
                                                          "GNAQ", "RIC8A", "LCE2D", 
                                                          "PRKCE"), ],
    aes(label = Genes),
    color = "black",
    min.segment.length = 0,
    fontface = "italic", size = 4
  ) +
  geom_vline(xintercept = c(-1.83, 1.83), linetype = "dashed") +
  xlab(expression("log"[2] ~ "(Fold change)")) +
  ylab(expression("log"[10] ~ "(1 / Padj)")) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Rbase save PDF
pdf(here("results/plots/fig4_volcano.pdf"), width = 12, height = 7)
print(P)
dev.off()


# ggsave
ggsave(
  filename = here("results/plots/fig4_volcano.jpeg"),
  plot = P,
  device = "jpeg",
  width = 12,
  height = 7,
  units = "in",
  dpi = 600
)

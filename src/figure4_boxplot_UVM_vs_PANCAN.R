# 11th June 24. Here we take the file from sumeet.patiyal@nih.gov 
# and graph it for publication.
# Updated july 2024

# Load necessary libraries
library(ggplot2)
library(readr)
library(extrafont)
library(here)

# Load fonts (run only once)

loadfonts(device = "pdf")
# Use 'device = "pdf"' on non-Windows platforms

# Import the data
data <- read_csv(here("data/MAGeCK_UVM_Vs_Pancancer_2_input_file.csv"))

# Define the order of the genes
gene_order <- c("SOX10", "MITF", "CDS2", "CDIPT", "RASGRP3", "TFAP2A", 
                "GNAQ", "RIC8A", "LCE2D", "PRKCE")

# Plot boxplot
p <- ggplot(data, aes(x = factor(Gene, levels = gene_order), 
                      y = Rank, fill = Label)) +
  geom_boxplot(width = 1, outlier.alpha = 0.7) +
  scale_fill_manual(values = c("UVM" = "#226E9C", "Pan-Cancer" = "#AB1866")) +
  labs(title = "Gene Rank Comparison", x = "Gene", y = "Rank") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 14, face = "bold.italic", color = "black"),
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14),
    title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.ticks = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )
p

# Rbase save PDF
pdf(here("results/plots/MAGeCK_UVM_vs_Pancancer_generank10.pdf"), width = 10, height = 10)
plot(p)
dev.off()

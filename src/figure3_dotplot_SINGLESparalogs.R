# Graph the data. 2nd July 2024.

# Load libraries
library(ggplot2)
library(readr)
library(here)

# Load the data
whole_genome_screen_filtered <- read_csv(here("data/2nd_July_2024.csv"))

# Reverse the levels of Gene
whole_genome_screen_filtered$Gene <- factor(whole_genome_screen_filtered$Gene, levels = rev(sort(unique(whole_genome_screen_filtered$Gene))))

# saveRDS(whole_genome_screen_filtered, "whole_genome_screen_filtered_object_for_plotting.rds")

# Create the plot
DotPlot_SINGLESparalogs <- ggplot(whole_genome_screen_filtered) +
  geom_point(size = 2.2, aes(x = Beta, y = Gene), colour = "blue") +
  theme_classic() +
  theme(axis.line.x = element_line(linewidth = 0.2, color = "black")) +
  theme(axis.line.y = element_line(linewidth = 0.2, color = "black")) +
  labs(x = bquote("Beta"), y = "") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) +
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(nrow = 2, reverse = TRUE)) +
  theme(axis.text.y = element_text(face = "italic", color = "black", margin = margin(t = 0, b = 0))) +
  theme(axis.text.x = element_text(color = "black")) +
  theme(text = element_text(size = 20))
DotPlot_SINGLESparalogs

# Save PDF
pdf(here("results/plots/DotPlot_SINGLESparalogs10.pdf"), height = 7, width = 6.1)
print(DotPlot_SINGLESparalogs)
dev.off()

# Save tiff
ggsave(here("results/plots/DotPlot_SINGLESparalogs10.tiff"), units = "in", width = 6.1, height = 7, dpi = 300)

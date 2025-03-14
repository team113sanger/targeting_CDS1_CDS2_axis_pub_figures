# July 2024
# Figure 4 - Enrichment dot plot CDS2

library(tidyverse)
library(here)

# Load pathway results
CDS2 <- read_tsv(here("data/CDS2_Pathway.csv"))


# Order and plot
CDS2_ord <- CDS2 |> 
            mutate(pathway = fct_reorder(Pathway, Fold_Enrichment))

p <- ggplot(CDS2_ord, aes(x = Fold_Enrichment, y = pathway)) +
  labs(x = "Fold Enrichment", y = "", title = "Gene Set Enrichment Analysis") +
  geom_point(size = 5.9, color = "blue") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", size = 15, color = "black", vjust = 0.5),
    axis.title.y = element_text(lineheight = 2, color = "black", size = 14),
    plot.title = element_text(face = "bold", size = 18),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = 13, color = "black", face = "bold")
  ) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20))

# Rbase save PDF
pdf(here("results/plots/CDS2_enrichment17.pdf"), width = 17, height = 10)
plot(p)
dev.off()

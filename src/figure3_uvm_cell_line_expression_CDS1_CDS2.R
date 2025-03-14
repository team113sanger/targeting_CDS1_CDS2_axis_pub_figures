# Uveal cell lines' gene expression - boxplot
# Updated July 2024

# Load libraries
library(tidyverse)
library(tidylog)
library(here)


# Read expression data
raw_exp <- read_table(here("data/6460_htseq_tpm-ENSv103_NEW_RNA_Seq.txt"))

exp <- raw_exp %>%
  pivot_longer(3:ncol(.), names_to = "sample", values_to = "tpm") %>%
  rename(gene = external_gene_name, ensembl = ENSEMBL_GENE_ID) |>
  mutate(line = sub("_.*", "", sample)) |>
  mutate(rep = gsub("^[^_]*_([^_]*)_.*", "\\1", sample)) |>
  relocate(line, rep, sample, gene, ensembl, tpm) |>
  mutate(log2_tpm_plus_one = log2(tpm + 1))

# Mean
mean_exp <- exp |>
  group_by(line, gene, ensembl) |>
  summarise(mean_log2_tpm_plus_one = mean(log2_tpm_plus_one))

# Gene pairs
d <- mean_exp |>
  filter(gene %in% c("CDS1", "CDS2", "RIC8A", "RIC8B", "SPTSSA", "SPTSSB")) %>%
  mutate(gene = factor(gene, levels = c("CDS2", "CDS1", "RIC8A", "RIC8B", "SPTSSA", "SPTSSB"))) %>%
  mutate(gene_pair = case_when(
    gene == "CDS1" | gene == "CDS2" ~ "CDS1_CDS2",
    gene == "RIC8B" | gene == "RIC8A" ~ "RIC8A_RIC8B",
    gene == "SPTSSB" | gene == "SPTSSA" ~ "SPTSSA_SPTSSB",
  ), .before = gene) |>
  mutate(gene_n = if_else(gene %in% c("CDS1", "RIC8B", "SPTSSB"), "gene_2", "gene_1")) %>%
  mutate(gene_pair = factor(gene_pair, levels = c("CDS1_CDS2", "RIC8A_RIC8B", "SPTSSA_SPTSSB")))

# Plot
boxp <- d |>
  ggplot(aes(x = gene, y = mean_log2_tpm_plus_one, fill = gene_n)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "grey20") +
  geom_point(alpha = 0.9, color = "grey30") +
  theme_bw() +
  geom_line(aes(x = gene, y = mean_log2_tpm_plus_one, group = line), color = "grey20", alpha = 0.3) +
  facet_grid(~gene_pair, scales = "free") +
  scale_fill_discrete(type = c("#E6A667", "#558593")) + # option 2 for blue color #4682B4
  theme(legend.position = "none") +
  labs(x = "Gene", y = "Mean Expression Log2(TPM+1)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    # strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.75)
  ) +
  theme(axis.text.x = element_text(size = 20, face = "italic", color = "black")) +
  theme(axis.title.y = element_text(size = 18)) +
  scale_y_continuous(
    limits = c(0, 8), breaks = seq(0, 14, by = 2),
    minor_breaks = seq(0, 14, 0.5)
  )

boxp


# Save PDF
pdf(here("results/plots/uvm_cell_line_expression_boxplot.pdf"), width = 10, height = 7)
plot(boxp)
dev.off()

# Save tiff
ggsave(here("results/plots/uvm_cell_line_expression_boxplot.tiff"), units = "in", width = 10, height = 10, dpi = 300)

# TCGA uveal cell lines
# Gene expression boxplot CDS2-CDS1 RIC8A-RIC8B SPTSSA-SPTSSB
# Updated July 2024

# Load libraries
library(tidyverse)
library(ggplot2)
library(tidylog)
library(broom)
library(here)

# Load data
mRNA_Expression_RSEM_Batch_normalized_from_Illumina_HiSeq_RNASeqV2 <- read_delim(here("data/mRNA_expression_RSEM_batch_normalized_illumina_HiSeq_RNASeqV2.txt"),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)


# Load expression data and process
exp <- mRNA_Expression_RSEM_Batch_normalized_from_Illumina_HiSeq_RNASeqV2 |>
  as_tibble() |>
  mutate(STUDY_ID = toupper(sub("_.*", "", STUDY_ID))) %>%
  pivot_longer(3:ncol(.), names_to = "gene", values_to = "tpm") |>
  filter(gene != "MDM2") |>
  filter(gene != "MDM4") |>
  mutate(log_tpm_plus_one = log2(tpm + 1)) |>
  mutate(gene_pair = case_when(
    gene == "CDS1" | gene == "CDS2" ~ "CDS1_CDS2",
    gene == "RIC8A" | gene == "RIC8B" ~ "RIC8A_RIC8B",
    # gene == 'MDM2' | gene == 'MDM4' ~ 'MDM2_MDM4',
    gene == "SPTSSA" | gene == "SPTSSB" ~ "SPTSSA_SPTSSB",
  ), .before = gene) |>
  mutate(gene_n = if_else(gene %in% c("CDS1", "MDM4", "RIC8B", "SPTSSB"), "gene_2", "gene_1")) |>
  mutate(gene_pair = factor(gene_pair, levels = c("CDS1_CDS2", "RIC8A_RIC8B", "SPTSSA_SPTSSB"))) |>
  arrange(SAMPLE_ID, gene_pair, gene_n) %>%
  mutate(gene = factor(gene, levels = unique(x = .[["gene"]])))

# saveRDS(exp, "exp_object_for_plotting.rds")

# Plot
plot <- ggplot(exp, aes(x = gene, y = log_tpm_plus_one, fill = gene_n)) +
  geom_boxplot(width = 0.3, outlier.shape = NA, color = "grey20") +
  geom_point(alpha = 0.9, color = "grey30") +
  theme_bw() +
  geom_line(aes(x = gene, y = log_tpm_plus_one, group = SAMPLE_ID), color = "grey20", alpha = 0.1) +
  facet_grid(~gene_pair, scales = "free", space = "free") +
  scale_fill_discrete(type = c("#E6A667", "#558593")) +
  facet_grid(~gene_pair, scales = "free") +
  theme(legend.position = "none") +
  labs(x = "Gene", y = "Log2(RSEM+1)") +
  # labs(y = expression(log[2](TPM + 1))) +
  theme( # panel.grid.major = element_blank(),
    # panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  ) +
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2), minor_breaks = seq(0, 14, 0.5)) +
  theme(axis.text.x = element_text(face = "italic", color = "black")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    # strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text = element_blank(),
    axis.title.x = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  theme(axis.text.x = element_text(size = 20)) +
  theme(axis.text.y = element_text(size = 20, colour = "black")) +
  theme(axis.title.y = element_text(size = 20))
plot

# Save PDF
pdf(here("results/plots/tcga_uvm_exp_boxplot10.pdf"), height = 7, width = 10)
print(plot)
dev.off()

# Save tiff
ggsave(here("results/plots/tcga_uvm_exp_boxplot.tiff"), units = "in", width = 10, height = 10, dpi = 300)


# Perform t-test
exp |>
group_by(gene_pair, gene, gene_n) |>
summarise(mean_log_tpm_plus_one = mean(log_tpm_plus_one), 
          median_log_tpm_plus_one = median(log_tpm_plus_one))


exp |> 
  select(STUDY_ID, SAMPLE_ID, gene_pair, gene_n, log_tpm_plus_one) |>
  pivot_wider(names_from = gene_n, values_from = log_tpm_plus_one) |>
  group_by(gene_pair) %>%
  do(tidy(t.test(.$gene_1, .$gene_2, mu = 0, alt = "one.sided", alternative = "greater", paired = TRUE, conf.level = 0.95))) %>%
  mutate(p_bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
  mutate(fdr = p.adjust(p.value, method = "hochberg"))


exp |>
  group_by(gene) |>
  summarise(
    median = median(log_tpm_plus_one),
    min = min(log_tpm_plus_one),
    max = max(log_tpm_plus_one),
    IQR_low = quantile(log_tpm_plus_one, 0.25),
    IQR_high = quantile(log_tpm_plus_one, 0.75)
  )

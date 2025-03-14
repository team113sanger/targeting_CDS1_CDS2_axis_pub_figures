# Count the number of times a gene is in the list

library(tidyverse)
library(biomaRt)
library(here)


combined_singles_results_scaled <- read_delim(
    here("data/combined_singles_results.scaled.tsv"),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

CEGv2_copy <- read_delim(
    here("data/CEGv2.txt"),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

sig_single <- combined_singles_results_scaled |> 
              filter(is_depleted_mageck == 1)

gene_counts <- sig_single |>
  group_by(gene) |>
  summarize(count = n()) |>
  arrange(desc(count))

write_csv(gene_counts, here("results/tables/gene_counts.csv"))


# Use the Ensembl mart
ensembl <- useMart("ensembl", 
                   dataset = "hsapiens_gene_ensembl")

# List of genes to check
genes_to_check <- c(
  "ATP6V1B2", "BUB3", "CCND1", "CHORDC1", "CKAP5", "DYNLRB1", "MYC", "POLR2L", 
  "PRPF38A", "RPL14", "RPL8", "RRM1", "SF3B3", "ALDOA", "CDIPT", "CLNS1A", 
  "CLTC", "COPA", "HSCB","KIF11", "LONP1", "MED18", "MED20", "MRPL38", 
  "MTG2", "NAA10", "NAE1", "PPP4C", "PRPF8", "PSMB3", "RFC2", "RPL10A", 
  "RPL11", "RPL13", "RPL3", "RPL36", "RPP21", "RPS19", "RRM2", "SF1", "SF3B1", 
  "SF3B2", "SNRPD3", "TIMM10", "TOMM40", "TOP2A", "VCP", "ACTL6A", "AKIRIN2",
  "ATP5F1B", "CCT6A", "CTR9", "ELL", "ERH", "GAPDH", "HMGCS1", "INTS6", 
  "MRPL41", "PRKDC","RAN", "RIC8A", "RPS11", "RPS23", "RPS5", "SFPQ", 
  "SMC1A", "SNRPD2", "TUBB", "VMP1", "ABT1", "ATP5F1A", 
  "ATP6V1C1", "ATP6V1D", "CDK2", "CDK7", "COX6B1", 
  "DHX9", "EIF3B", "HUS1", "ILF3", "LSM4", "NDUFB4", "NKAP", "PPP2R1A",
   "RBBP4", "RNF113A", "RPS3", "SBNO1", "SF3A1",
  "SOD1", "TRAPPC3", "TUFM", "UQCRC1", "XAB2", "AHCYL1", "CDS2", "CPSF3", 
  "EIF4A3", "GTF3C1", "KIF4A", "MED12", "NXT1", "PMPCB", "POLR2B", "PSMA2", 
  "RNF40", "RPL19", "RPL31", "RPS8", "SF3B5", "SIN3A", "SMARCA4", 
  "TOMM22", "TRRAP", "WDR12", "ZNF407"
)

# Count the number of genes in the list
num_genes <- length(genes_to_check)

# 116. This number matches the gene-counts dataframe.

# Retrieve paralogs
paralog_data <- getBM(
  attributes = c("external_gene_name", "hsapiens_paralog_associated_gene_name"),
  filters = "external_gene_name",
  values = genes_to_check,
  mart = ensembl
)

# Filter genes that have paralogs
genes_with_paralogs <- paralog_data[paralog_data$hsapiens_paralog_associated_gene_name != "", ]

# Unique list of genes with paralogs
unique_genes_with_paralogs <- unique(genes_with_paralogs$external_gene_name)



sig_single_paralogs <- sig_single |> 
                       filter(sig_single$gene %in% unique_genes_with_paralogs)

sig_single_paralogs_6 <- sig_single_paralogs |>
  group_by(gene) |>
  filter(n() >= 6) |>
  ungroup()

# filter genes for those that are on the common essential gene list V2.

sig_single_paralogs_6_final <- sig_single_paralogs_6 |> 
                               filter(!sig_single_paralogs_6$gene %in% CEGv2_copy$GENE)

write_csv(sig_single_paralogs_6_final, here("results/tables/sig_single_paralogs_6_final.csv"))

SinglesResults_ParalogPairDerived <- sig_single_paralogs_6_final |>
  group_by(gene) |>
  mutate(avg_neg_lfc = mean(neg.lfc, na.rm = TRUE)) |>
  ungroup()

# take SinglesResults_ParalogPairDerived and manually add the screen or not screen column

SinglesResults_ParalogPairDerived_previous_screen <- read_csv(here("data/SinglesResults_ParalogPairDerived_previous_screen.csv"))
# the BINARY column is from Colms spreadsheet of hits.

write_csv(SinglesResults_ParalogPairDerived, here("results/tables/SinglesResults_ParalogPairDerived.csv"))

DotPlot_SINGLESparalogs <- ggplot(SinglesResults_ParalogPairDerived_previous_screen) +
  geom_point(size = 2.2, aes(x = neg.lfc, y = (reorder(gene, avg_neg_lfc, decreasing = TRUE)), color = factor(BINARY)), size = 2) + 
  theme_classic() + 
  theme(axis.line.x = element_line(linewidth = 0.2, color = "black")) + 
  theme(axis.line.y = element_line(linewidth = 0.2, color = "black")) + 
  scale_color_manual(values = c("#C52930", "grey"), 
  labels = c("Not screened or member of previous SL screen non-hit", 
            "Member of previous SL screen hit")) + 
            labs(x = bquote("Normalised log2 fold change"), y = "") + 
            theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin()) + 
            theme(legend.title = element_blank()) + 
            guides(colour = guide_legend(nrow = 2, reverse = TRUE)) + 
            theme(axis.text.y = element_text(face = "italic")) + 
            theme(text = element_text(family = "Helvetica", size = 14))
DotPlot_SINGLESparalogs
dev.off()

# the above figure contains all genes including some in the 
# from the non-paralog genes. These include:
# ALDOA, RRM2,LSM4,SOD1,VMP1, CDK2,KIF4A, NAE1,ZNF407,SFPQ, PRKDC.

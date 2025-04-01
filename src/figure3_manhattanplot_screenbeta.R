# Screen Beta distributions - Manhattan plot
# 8th July 2024
# Updated July 20 2024

library(tidyverse)
library(ggrepel)
library(here)
# library(biomaRt)
library(readxl)
# Import the data from upstream steps
# https://github.com/team113sanger/uveal_melanoma_CRISPR_downstream/blob/main/results/MAGeCK_gene_corrected_beta.tsv

# https://github.com/team113sanger/uveal_melanoma_CRISPR_downstream/blob/main/results/MAGeCK_gene_corrected_fdr.tsv

beta <- read_tsv("data/MAGeCK_gene_corrected_beta.tsv") |> 
        pivot_longer(cols = -c(genes), names_to = "ModelName", values_to = "Beta")

fdr <- read_tsv("data/MAGeCK_gene_corrected_fdr.tsv") |> 
        pivot_longer(cols = -c(Gene), names_to = "ModelName", values_to = "FDR")


# Fetch Kosuke Yusa's gene ID mapping 
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S2211124716313353-mmc2.xlsx", destfile = "data/1-s2.0-S2211124716313353-mmc2.xlsx")
hgnc_ids <- read_xlsx("data/1-s2.0-S2211124716313353-mmc2.xlsx", sheet = 2) |> 
dplyr::select(Gene, HGNC_ID) |> 
distinct()


hgnc_ids |> 
select(HGNC_ID) |> 
write_tsv("results/query_hgnc_ids.txt")

# Biomart down at time of wiritng - fetch gene names from the HGNC ID manually 
# Mart export 
gene_id_update_map <- hgnc_ids |> 
            left_join(read_tsv("data/mart_export.txt"), by = c("HGNC_ID" = "HGNC ID")) |> 
            dplyr::select(Gene, `Gene name`) |>
            distinct() |> 
            mutate(`Gene name`= case_when(is.na(`Gene name`) ~ Gene, TRUE ~ `Gene name`))

gene_id_update_map |> dplyr::rename("Updated_ID" = "Gene name")


mageck_output <- left_join(beta, fdr, by = c("genes" = "Gene", "ModelName")) |>
dplyr::rename("Gene" = "genes") |>
left_join(gene_id_update_map, by = "Gene") |>
select(-Gene) |>
rename("Gene" = `Gene name`) |> 
relocate(Gene, .before =  "ModelName")

write_tsv(mageck_output, "data/mageck_output_updated_ids.tsv")

mageck_output_sigs_ggplot <- mageck_output |>
  mutate(
    Significant =
      case_when(
        Beta >= 0.5 & FDR <= 0.05 ~ "1",
        Beta <= -0.5 & FDR <= 0.05 ~ "1", TRUE ~ "0"
      )
  )

CRISPRInferredCommonEssentials_24Q2 <- read_csv("https://figshare.com/ndownloader/files/46489597")

# Clean data
DepMAP_24Q2_essentials <- gsub(
  "\\s*\\(.*?\\)\\s*", "",
  CRISPRInferredCommonEssentials_24Q2$Essentials
)
# Write file
write_csv(DepMAP_24Q2_essentials |> 
          as_tibble(), col_names = FALSE,
          here("results/tables/DepMAP_24Q2_essentials.csv"))



hart_cegs <- read_tsv(here("data/CEGv2.txt")) |> 
left_join(gene_id_update_map, by = c("GENE" = "Gene")) |>
mutate(`Gene name`= case_when(is.na(`Gene name`) ~ GENE, TRUE ~ `Gene name`)) |> 
pull(`Gene name`)

# Filter mageck_output to keep genes that are not in the Essentials list
filtered_mageck_output <- mageck_output |>
                          filter(!(Gene %in% DepMAP_24Q2_essentials)) |>
                          filter(!(Gene %in% hart_cegs))

# Expand the ModelName column by separating rows, and create new rows with _NE appended
expanded_mageck_output <- filtered_mageck_output |>
  separate_rows(ModelName, sep = ";") |>
  mutate(ModelName = paste0(ModelName, "_NE"))

expanded_mageck_output <- expanded_mageck_output |> mutate(
  Significant =
    case_when(
      Beta >= 0.5 & FDR <= 0.05 ~ "1",
      Beta <= -0.5 & FDR <= 0.05 ~ "1", TRUE ~ "0"
    )
)

# Count how many significant genes
expanded_mageck_output |> 
count(Significant)

# hits
ne_hits <- expanded_mageck_output |>
  group_by(Gene, Significant) |>
  tally() |>
  filter(Significant == 1) |>
  filter(n >= 6) |>
  mutate(plot_shape = TRUE) |>
  select(Gene, plot_shape)

# top positives
top_positives <- expanded_mageck_output |>
  filter(Beta > 0) |>
  group_by(Gene, Significant) |>
  tally() |>
  filter(Significant == 1) |>
  filter(n >= 3) |>
  arrange(desc(n)) |>
  ungroup() |>
  slice_head(n = 10) |>
  mutate(top_recurrents = TRUE) |>
  select(Gene, top_recurrents)

# top negatives
top_negatives <- expanded_mageck_output |>
  filter(Beta < 0) |>
  group_by(Gene, Significant) |>
  tally() |>
  filter(Significant == 1) |>
  filter(n >= 6) |>
  arrange(desc(n)) |>
  ungroup() |>
  slice_head(n = 10) |>
  mutate(top_recurrents = TRUE) |>
  select(Gene, top_recurrents)

# bind into top recurrents
top_recurrents <- bind_rows(top_positives, top_negatives)

# Combine the original mageck_output with the expanded dataframe, keeping original entries
combined_output <- bind_rows(mageck_output_sigs_ggplot, expanded_mageck_output)

# Write the output to a new CSV file
write_csv(combined_output, here("results/tables/combined_mageck_output.csv"))

# Reading and preparing data for plotting


# Subsetting the significant genes data
sign_gene <- combined_output[combined_output$Significant == 1, ]

# Getting the frequency of a significant gene
freq_table <- sign_gene |>
  group_by(Gene) |>
  summarise(frequency = n()) |>
  arrange(desc(frequency))

# Adding a column to denote if the gene is significant in 3 or above modelnames
sign_gene$sign_in_3 <- ifelse(sign_gene$Gene %in% freq_table$Gene[freq_table$frequency >= 3], TRUE, FALSE)
# glimpse(sign_gene)

# Creating the manhattan plot

# Preprocess step
data <- combined_output |>
  left_join(ne_hits) |>
  left_join(top_recurrents) |>
  replace_na(list(plot_shape = FALSE, top_recurrents = FALSE)) |>
  mutate(plot_shape = case_when(!grepl(ModelName, pattern = "NE") ~ FALSE, TRUE ~ plot_shape)) |>
  mutate(gene_name = case_when((top_recurrents & plot_shape) ~ Gene, TRUE ~ NA))

# saveRDS(data, "data_object_for_plotting.rds")

data
# Plot
set.seed(42)
p <- ggplot(data, aes(x = ModelName, y = Beta)) +
  geom_point(aes(color = Beta, shape = factor(plot_shape)),
    position = position_jitter(0.2), stroke = NA
  ) +
  scale_color_gradient2(midpoint = 0, low = "#FF1A1A", mid = "white", high = "navyblue") +
 geom_text_repel(aes(label = gene_name),
                min.segment.length = 0,
                size = 2.0,
                force = 1.1,
                box.padding = 0.25,
                max.iter = 10000,         # Increase iterations for better placement
                max.overlaps = Inf
 ) +
  labs(x = "\nModel Name", y = "Beta", shape = "Hit in > 6 models") +
  geom_hline(yintercept = c(-0.5, 0.5), col = "black", linetype = "dashed", linewidth = 0.5) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -45, face = "bold", color = "black", vjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 11, face = "bold"),
    title = element_text(size = 11, face = "bold"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

p

# Save PDF
pdf(here("results/plots/figure3_manhattanplot10.pdf"), height = 10, width = 15)
print(p)
dev.off()

# Save PNG
png(here("results/plots/figure3_manhattanplot.png"), height = 35 * 350, width = 55 * 350, res = 1200)
print(p)
dev.off()

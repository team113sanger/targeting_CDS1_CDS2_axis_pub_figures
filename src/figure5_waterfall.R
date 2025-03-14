# 27th June 2024


library(ggplot2)
library(tidyverse)
library(here)
library(curl)

# Import the data.These are all from 23Q4 release of DepMap

# RNA_Seq data

download.file("https://plus.figshare.com/ndownloader/files/43347204", destfile = "data/OmicsExpressionProteinCodingGenesTPMLogp1(23Q4).csv")
OmicsExpressionProteinCodingGenesTPMLogp1_23Q4_ <- read_csv(here("data/OmicsExpressionProteinCodingGenesTPMLogp1(23Q4).csv"))

# CRISPR Dependency Data
download.file("https://plus.figshare.com/ndownloader/files/43347204", destfile = "data/CRISPRGeneEffect_23Q4.csv")

CRISPRGeneEffect_23Q4 <- read_csv(here("data/CRISPRGeneEffect_23Q4.csv"))

# Models file
Model_23Q4 <- read_csv("https://plus.figshare.com/ndownloader/files/43746708")

# Select the columns of interest

RNA_Seq_CDS1_23Q4_CDS1 <- OmicsExpressionProteinCodingGenesTPMLogp1_23Q4_ |> select(`...1`, `CDS1 (1040)`)
CRISPRGeneEffect_23Q4_CDS2 <- CRISPRGeneEffect_23Q4 |> 
                              select(`...1`, `CDS2 (8760)`)
Model_23Q4_2 <- Model_23Q4 |> 
                select(ModelID, DepmapModelType)

# add an extra column to the dataframe so that it can be merged

Model_23Q4_3 <- Model_23Q4_2 |> 
                mutate(...1 = ModelID)

# Merge the dataframes

combined_df_1 <- full_join(Model_23Q4_3, CRISPRGeneEffect_23Q4_CDS2, by = "...1")

final_combined <- full_join(combined_df_1, RNA_Seq_CDS1_23Q4_CDS1, by = "...1")

final_combined_minus_UM <- final_combined |> 
                            filter(!DepmapModelType == "UM")




write_csv(final_combined_minus_UM, here("results/tables/final_combined_minus_UM.csv"))

# Remove rows with NA values for plotting
filtered_data <- na.omit(final_combined_minus_UM[, c("CDS2 (8760)", "CDS1 (1040)")])

# filtered_data |> mutate("Colour" = case_when(`CDS1 (1040)` >0 & `CDS1 (1040)` <=1 ~ "A", .default = "X")) #`CDS1 (1040)` >1))


filtered_data_2 <- filtered_data |> 
                   mutate("Colour" = case_when(`CDS1 (1040)` > 0 & `CDS1 (1040)` <= 2 ~ "A", `CDS1 (1040)` > 2 & `CDS1 (1040)` <= 4 ~ "B", `CDS1 (1040)` > 4 & `CDS1 (1040)` <= 6 ~ "C", `CDS1 (1040)` > 6 & `CDS1 (1040)` <= 20 ~ "D", .default = "D")) # `CDS1 (1040)` >1))


# Read and preprocess data. This is the same dataframe as provided by Jen but I deleted the cell lines marked as "eye".

# Order the data by CDS2_chronos in descending order
CDS1_CDS2_df <- filtered_data_2[order(filtered_data_2$`CDS2 (8760)`, decreasing = TRUE), ]

# Map colors based on the 'Colour' column
color_map <- c(
  "A" = "#247198",
  "B" = "#A8C2B9",
  "C" = "#EACA80",
  "D" = "#f5af3d",
  "E" = "#CB6B79",
  "F" = "#CB6B79"
)
CDS1_CDS2_df$Color <- color_map[CDS1_CDS2_df$Colour]

# Create the waterfall plot using ggplot2
Waterfall <- ggplot(CDS1_CDS2_df, aes(x = reorder(row.names(CDS1_CDS2_df), -`CDS2 (8760)`), y = `CDS2 (8760)`, fill = Colour)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_map) +
  labs(y = "CDS2 Chronos score") +
  theme_classic()

# Display the plot
print(Waterfall)

# Save the plot with 600 DPI
ggsave(
  filename = here("results/plots/Waterfall.jpg"),
  plot = Waterfall,
  device = "jpeg",
  width = 15,
  height = 5,
  units = "in",
  dpi = 600
)

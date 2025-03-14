### NItesh SHriwash
### Date 4th March, 2024

library(ggplot2)
library(here)
source(here('src/helper.R'))

# reading the data for plotting
data <- read_csv(here("data/CDS_scores_by_cancer_type.csv"))
data <- na.omit(data) # removing NA values
data_long <- reshape2::melt(data) # melting data for plotting (long format)
colnames(data_long) <- c("cancer_type", "Gene", "value") # changing column names

# order of cancer types
cancer_order <- c(
  "KICH", "THCA", "COAD", "BRCA", "KIRP", "READ", "LUAD", "PRAD", "ESCA",
  "HNSC", "UCEC", "STAD", "CESC", "PAAD", "BLCA", "LUSC", "OV", "ACC", "KIRC",
  "CHOL", "THYM", "UCS", "TGCT", "MESO", "PCPG", "LGG", "GBM", "UVM", "SARC", "LIHC",
  "SKCM", "DLBC", "LAML"
)

# converting cancer type column to factor for the purpose of ordering of the x axis labels
data_long[["cancer_type"]] <- factor(data_long[["cancer_type"]], 
                                    levels = cancer_order) 

# creating the split violin plot
p <- ggplot(data_long, aes(x = cancer_type, y = value, fill = Gene)) +
  geom_split_violin(trim = FALSE, scale = "width", aes(color = Gene)) +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F) +
  labs(x = "Cancer Type", y = "log2(TPM+1)") +
  scale_fill_manual(values = c("darkviolet", "gold2"), name = "Gene") +
  scale_color_manual(values = c("darkviolet", "gold2"), name = "Gene") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    axis.title = element_text(size = 15, face = "bold"),
    axis.ticks = element_line(color = "black", linewidth = 1),
    legend.text = element_text(size = 13, face = "bold"),
    title = element_text(size = 15, face = "bold"),
    panel.border = element_rect(color = "black", fill = "transparent"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# saving the image
jpeg(here("results/plots/split_violin.jpg"), height = 25 * 350, width = 45 * 350, res = 1200)
print(p)
dev.off()

#  Targeting the CDS1/2 axis as a therapeutic strategy in uveal melanoma and pan-cancer Figures

## Overview

This repository contains a series of R scripts that were used in generating the figures for **Targeting the CDS1/2 axis as a therapeutic strategy in uveal melanoma and pan-cancer**. 
It includes ab comparison of gene essentiality scores from uveal melanoma cell lines (`data/MAGeCK_gene_corrected_beta.tsv`; generated using MAGECK-MLE) with pan-cancer samples derived from DepMap. 
 

## Analysis Steps
Users can reproduce the scripts bundled in this package by creating a [Github codespace](https://docs.github.com/en/codespaces/getting-started/quickstart) using the `.devcontainer.json` file that is bundled with these scripts.


## Orangisation
```
├── LICENSE
├── README.md
├── data
│   ├── 230820_PairedScreen_DotPlot_Rinput.csv
│   ├── 2nd_July_2024.csv
│   ├── 6460_htseq_tpm-ENSv103_NEW_RNA_Seq.txt
│   ├── AHCYCL1.csv
│   ├── CDS1.tsv
│   ├── CDS2.tsv
│   ├── CDS2_Pathway.csv
│   ├── CDS_scores_by_cancer_type.csv
│   ├── CEGv2.txt
│   ├── CRISPRInferredCommonEssentials_24Q2.csv
│   ├── MAGeCK_UVM_Vs_Pancancer_2_input_file.csv
│   ├── MAGeCK_gene_corrected_beta.tsv
│   ├── Model.csv
│   ├── Model_23Q4.csv
│   ├── SinglesResults_ParalogPairDerived_previous_screen.csv
│   ├── Supplementary Table 6.xlsx
│   ├── Supplementary Table 9.csv
│   ├── annotation.tsv
│   ├── combined_singles_results.scaled.tsv
│   ├── mRNA_expression_RSEM_batch_normalized_illumina_HiSeq_RNASeqV2.txt
│   ├── mageck_output.tsv
│   ├── pan_genes.csv
│   └── ranked_list.txt
├── renv
│   ├── activate.R
│   └── settings.json
├── renv.lock
├── results
│   ├── plots
│   │   ├── 2024_SinglesDotPlotParalogs25th_6_24(AHCYCL1_fix).jpeg
│   │   ├── DotPlot_SINGLESparalogs10.pdf
│   │   ├── DotPlot_SINGLESparalogs10.tiff
│   │   ├── MAGeCK_UVM_vs_Pancancer_generank10.pdf
│   │   ├── Waterfall.jpg
│   │   ├── fig4_volcano.jpeg
│   │   ├── fig4_volcano.pdf
│   │   ├── figure3_manhattanplot.png
│   │   ├── figure3_manhattanplot10.pdf
│   │   ├── split_violin.jpg
│   │   ├── tcga_uvm_exp_boxplot.tiff
│   │   ├── tcga_uvm_exp_boxplot10.pdf
│   │   ├── uvm_cell_line_expression_boxplot.pdf
│   │   └── uvm_cell_line_expression_boxplot.tiff
│   └── tables
│       ├── DepMAP_24Q2_essentials.csv
│       ├── SinglesResults_ParalogPairDerived.csv
│       ├── combined_mageck_output.csv
│       ├── final_combined_minus_UM.csv
│       ├── gene_counts.csv
│       └── sig_single_paralogs_6_final.csv
└── src
    ├── figure2_pairs.R
    ├── figure2_singles.R
    ├── figure3_dotplot_SINGLESparalogs.R
    ├── figure3_manhattanplot_screenbeta.R
    ├── figure3_tcga_uvm_expression.R
    ├── figure3_uvm_cell_line_expression_CDS1_CDS2.R
    ├── figure4_boxplot_UVM_vs_PANCAN.R
    ├── figure4_enrichmentplot_CDS2.R
    ├── figure4_volcanoplot.R
    ├── figure5_split_violin_plot.R
    ├── figure5_waterfall.R
    └── helper.R
```

## Dependencies
This analysis was conducted with R (4.2.3). All R dependencies for this project are detailed within the project `renv.lock` file.

## Contact
- David J Adams (da1@sanger.ac.uk)
- Jamie Billington (jb63@sanger.ac.uk)

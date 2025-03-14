#  Targeting the CDS1/2 axis as a therapeutic strategy in uveal melanoma and pan-cancer Figures

[![DOI](https://zenodo.org/badge/948453959.svg)](https://doi.org/10.5281/zenodo.15025727)

## Overview

This repository contains a series of R scripts that were used in generating the final figures for **Targeting the CDS1/2 axis as a therapeutic strategy in uveal melanoma and pan-cancer**. It takes data generated in the related analyses in https://github.com/team113sanger/Targeting-the-CDS1-2-axis-as-a-therapeutic-strategy-in-uveal-melanoma-and-pan-cancer as inputs for plotting.


## Analysis Steps
Users can reproduce the scripts bundled in this package by creating a [Github codespace](https://docs.github.com/en/codespaces/getting-started/quickstart) using the `.devcontainer.json` file that is bundled with these scripts.


## Orangisation
```
.
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
│   ├── MAGeCK_UVM_Vs_Pancancer_2_input_file.csv
│   ├── MAGeCK_gene_corrected_beta.tsv
│   ├── SinglesResults_ParalogPairDerived_previous_screen.csv
│   ├── Supplementary Table 6.xlsx
│   ├── Supplementary Table 9.csv
│   ├── combined_singles_results.scaled.tsv
│   ├── mRNA_expression_RSEM_batch_normalized_illumina_HiSeq_RNASeqV2.txt
│   ├── mageck_output.tsv
│   ├── pan_genes.csv
│   └── ranked_list.txt
├── renv
│   ├── activate.R
│   ├── library
│   ├── settings.json
│   └── staging
├── renv.lock
├── results
│   ├── plots
│   └── tables
└── src
    ├── figure2_pairs.R
    ├── figure2_singles.R
    ├── figure3_dotplot_SINGLESparalogs.R
    ├── figure3_manhattanplot_screenbeta.R
    ├── figure3_tcga_uvm_expression.R
    ├── figure3_uvm_cell_line_expression_CDS1_CDS2.R
    ├── figure4_enrichmentplot_CDS2.R
    ├── figure5_split_violin_plot.R
    ├── figure5_waterfall.R
    └── helper.R
```

## Dependencies
This analysis was conducted with R (4.2.3). All R dependencies for this project are detailed within the project `renv.lock` file.

## Contact
- David J Adams (da1@sanger.ac.uk)
- Jamie Billington (jb63@sanger.ac.uk)

# Bacterial Responses to Biocide Chemicals

Analysis code and data for the research paper:

> **A Broad Examination on Bacterial Responses to a Wide Range of Biocide Chemicals and Exploration of Potential 'Biosensor' Strategy**

## Background

This project investigates how 21 environmental bacterial isolates respond to approximately 190 biocide compounds (pesticides, herbicides, fungicides, and insecticides) using high-throughput 96-well microplate growth assays. The analysis pipeline quantifies growth via area under the curve (AUC) of OD600 readings, identifies chemicals with statistically significant effects, clusters chemicals by bacterial response profiles, and tests for phylogenetic signal in the phenotypic responses.

## Experimental Design

| Parameter | Value |
|-----------|-------|
| Bacterial strains | 21 isolates (IDs: 74, 85, 88, 100, 186, 322, 331, 333, 350, 353, 371, 374, 380, 390, 398, 436, 442, 448, 487, 527, 565) |
| Chemical libraries | PS1 and PS2, ~96 compounds each |
| Replicates | 3 technical replicates per strain-library combination |
| Plate format | 96-well microplate (rows A-H, columns 1-12) |
| Measurement | OD600, 30-min intervals over 72 hours |
| Control | DMSO solvent control wells |
| Position control | Plate layout rotation across replicates 2 and 3 |

### Strain Taxonomy

The 21 isolates span 4 phyla and 15 genera:

| Strain ID | Species | Phylum | Genus |
|-----------|---------|--------|-------|
| 74 | *Bacillus soli* | Firmicutes | *Neobacillus* |
| 85 | *Bosea massiliensis* | Proteobacteria | *Bosea* |
| 88 | *Prolinoborus fasciculus* | Proteobacteria | *Prolinoborus* |
| 100 | *Pseudomonas baetica* | Proteobacteria | *Pseudomonas* |
| 186 | *Paenibacillus castaneae* | Firmicutes | *Paenibacillus* |
| 322 | *Pseudomonas vancouverensis* | Proteobacteria | *Pseudomonas* |
| 331 | *Pseudomonas baetica* | Proteobacteria | *Pseudomonas* |
| 333 | *Pseudomonas mandelii* | Proteobacteria | *Pseudomonas* |
| 350 | *Serratia* sp. | Proteobacteria | *Serratia* |
| 353 | *Aeromonas* sp. | Proteobacteria | *Aeromonas* |
| 371 | *Rhizobium herbae* | Proteobacteria | *Rhizobium* |
| 374 | *Aeromonas sobria* | Proteobacteria | *Aeromonas* |
| 380 | *Pseudomonas peli* | Proteobacteria | *Pseudomonas* |
| 390 | *Microbacterium arborescens* | Actinobacteria | *Microbacterium* |
| 398 | *Sphingobium czechense* | Proteobacteria | *Sphingobium* |
| 436 | *Brevundimonas staleyi* | Proteobacteria | *Brevundimonas* |
| 442 | *Hymenobacter glaciei* | Bacteroidetes | *Siccationidurans* |
| 448 | *Carnobacterium gallinarum* | Firmicutes | *Carnobacterium* |
| 487 | *Aeromonas popoffii* | Proteobacteria | *Aeromonas* |
| 527 | *Arthrobacter humicola* | Actinobacteria | *Arthrobacter* |
| 565 | *Fodinibacter* sp. | Actinobacteria | *Oryzobacter* |

## Analysis Pipeline

The analysis is performed in `Project_Main_Code.R`, organized into the following stages:

1. **Data loading and AUC calculation** (Sections 1-2) -- Raw OD600 time series are cleaned, and growth is quantified by spline-based AUC integration over 1-73 hours.
2. **Chemical mapping and plate correction** (Sections 3-3.5) -- Chemical identities are mapped to wells, with plate position rotation corrected for replicates 2 and 3.
3. **AUC ratio normalization** (Section 4) -- AUC values are normalized to per-strain mean DMSO control values.
4. **Statistical testing** (Section 5) -- Dunnett's test identifies chemicals with significantly different growth compared to the DMSO control.
5. **Clustering analysis** (Sections 5.5-6.1) -- Hierarchical clustering (Euclidean distance, average linkage) and model-based clustering (Gaussian mixture models via mclust) group chemicals by bacterial response similarity.
6. **PCA ordination** (Section 7) -- Principal component analysis visualizes strain-level variation in response to chemicals.
7. **Phylogenetic-phenotypic correlation** (Section 8) -- Mantel test (Kendall correlation, 9999 permutations) assesses whether phylogenetically related strains show similar chemical response profiles.
8. **Phylogenetic signal** (Section 9) -- Pagel's lambda tests and ancestral state reconstruction (`contMap`) evaluate phylogenetic signal for individual chemical responses.

## Repository Structure

```
.
├── Project_Main_Code.R                  # Main analysis script (~880 lines)
├── REQUIREMENTS.txt                     # R package dependencies
├── CLAUDE.md                            # AI assistant project guidelines
│
├── data/                                # All data files
│   ├── ps1.csv                          # Chemical library 1 plate layout
│   ├── ps2.csv                          # Chemical library 2 plate layout
│   ├── chemicaldetails.csv              # Chemical metadata (family, target type)
│   ├── aiden-strain-taxonomy.csv        # Full strain taxonomy (Species through Genus)
│   ├── spline-fits.csv                  # Pre-computed spline AUC for all strains/wells
│   ├── chemical_table_final.csv         # Clustering output (Cluster, Target, Dir1, Dir2)
│   ├── phylo.io_n.nwk                   # Phylogenetic tree (Newick, 21 strains)
│   ├── phylo.io_new.nwk                 # Re-rooted phylogenetic tree
│   ├── phylo.io_all.nwk                 # Full phylogenetic tree (including excluded strains)
│   ├── phylo.io_wo306.nwk              # Phylogenetic tree without strain 306
│   ├── ROF_Isolates_DEC_2021_Aiden.xlsx # Bacterial strain collection information
│   ├── Taxo.xlsx                        # Taxonomy data
│   ├── replicate-plate-layouts.xlsx     # Plate layout designs
│   ├── small_pesticide_library_layout_visual.csv
│   ├── aiden-sequences.fa               # 16S rRNA sequences (Aiden's strains)
│   ├── tom-sequences.fa                 # 16S rRNA sequences (Tom's strains)
│   ├── community-growth-curves.csv      # Community-level growth data (32 MB, not in git)
│   ├── isolate-growth-curves.csv        # Isolate-level growth data (325 MB, not in git)
│   └── wrksp.RData                      # Saved R workspace with intermediate results (not in git)
│
├── growth-curve/                        # Raw 96-well OD600 time series
│   ├── {strain}_p_{library}_r_{replica}.txt
│   ├── _{strain}_p_{library}_r_{replica}_2.txt  # Non-standard files (strains 85, 398, 436)
│   └── README_DATA_FORMAT.md
│
├── legacy/                              # Earlier modular R scripts (development history)
│   ├── auc.R                            # AUC calculation module
│   ├── dmsotest.R                       # DMSO control testing module
│   ├── mantel-test.R                    # Mantel test module
│   ├── tree.R                           # Clustering & phylo tree visualization
│   ├── load.R                           # Data loading snippet
│   └── Copy.R                           # Earlier version of main script (~576 lines)
│
├── write-up/                            # Manuscript and presentation (not in git)
│   ├── Draft.docx, Draft_2.docx         # Manuscript drafts
│   ├── Draft-TS.pdf, Draft_2-TS.pdf     # Supervisor-annotated versions
│   ├── Draft_2.pdf, Draft_3.pdf         # PDF exports
│   ├── Presentation.pptx                # Defence/seminar presentation
│   └── Zhang_EEC_MSc_01740950.pdf       # Final submitted manuscript
│
└── references/                          # Reference literature (not in git)
    ├── Latent_functional_diversity_may_accelerate_microbi.pdf
    ├── RJ-2016-021.pdf                  # R Journal (mclust methodology)
    ├── Scrucca_2010_StatComp.pdf        # mclust statistical computing
    └── Smith_2023_biorxiv.pdf
```

## Requirements

### R (>= 4.0)

Install all dependencies:

```r
# CRAN packages
install.packages(c(
  "factoextra", "circlize", "dplyr", "tidyr", "plotly", "ggplot2",
  "ggdendro", "tibble", "permute", "mclust", "DescTools", "plot3D",
  "berryFunctions", "plotrix", "gridExtra", "MASS", "shape", "formattable",
  "lattice", "vegan", "ggrepel", "ape", "phytools", "glmm", "wrapr"
))

# Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# GitHub
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("houyunhuang/ggcor")
```

## Usage

1. Open `Project_Main_Code.R` in RStudio.
2. Update the `setwd()` paths at the top of the script to point to your local `growth-curve/` directory.
3. Run section by section (the script is designed for interactive execution, not batch mode).

**Note:** The Dunnett's test step (Section 5) is computationally intensive and may take several hours depending on hardware.

## Key Outputs

- Heatmaps of AUC ratios and Dunnett's test significance across all strain-chemical combinations
- Chemical clustering dendrograms and model-based clustering with dimension reduction plots
- PCA biplots showing strain-level response variation
- Mantel test correlation plots linking phylogenetic and phenotypic distances
- Phylogenetic trees with ancestral state reconstruction for selected chemicals

## Citation

If you use this code or data, please cite:

> A Broad Examination on Bacterial Responses to a Wide Range of Biocide Chemicals and Exploration of Potential 'Biosensor' Strategy

## License

This project is for academic research purposes.

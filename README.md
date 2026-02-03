# Bacterial Responses to Biocide Chemicals

A research project examining bacterial responses to a wide range of biocide chemicals and exploring potential 'biosensor' strategies.

## Project Overview

This project investigates how multiple bacterial strains respond to over 100 different chemical compounds (pesticides, herbicides, fungicides, etc.). The analysis includes:

- Growth curve measurements and AUC (Area Under Curve) calculations
- Chemical clustering and grouping analysis
- Principal Component Analysis (PCA)
- Phylogenetic analysis and Mantel tests
- Biosensor candidate identification

## Directory Structure

```
.
├── Project_Main_Code.R              # Main R analysis script (880 lines)
├── REQUIREMENTS.txt                 # R package dependencies
├── bacteria-chemical-lib/           # Bacterial and chemical library data
│   ├── ROF_Isolates_DEC_2021_Aiden.xlsx  # Bacterial strain information
│   ├── Taxo.xlsx                    # Taxonomy data
│   ├── ps1.csv                      # Pesticide library 1 layout
│   ├── ps2.csv                      # Pesticide library 2 layout
│   ├── chemicaldetails.csv          # Chemical compound details (target, family)
│   ├── replicate-plate-layouts.xlsx # Plate layout designs
│   └── small_pesticide_library_layout_visual.csv  # Chemical library layouts (visual)
├── growth-curve/                    # Growth curve experimental data
│   └── [growth curve data files]    # Format: {strain}_p_{library}_r_{replica}.txt
└── data/                            # Additional data files
    ├── aiden-strain-taxonomy.csv    # Strain taxonomy information (template)
    ├── spline-fits.csv              # Additional spline fit data (template)
    └── phylo.io_n.nwk               # Phylogenetic tree - Newick format (21 strains)
```

## Requirements

### R Packages

The following R packages are required:

**Data Processing:**
- dplyr, tidyr, tibble, wrapr

**Statistical Analysis & Clustering:**
- factoextra, ComplexHeatmap, mclust, DescTools, MASS, vegan, permute

**Phylogenetic Analysis:**
- ape, phytools, ggtree

**Visualization:**
- ggplot2, plotly, ggrepel, ggdendro, ggcor
- circlize, gridExtra, lattice
- plot3D, berryFunctions, plotrix, shape, formattable

**Other:**
- glmm

### Installation

Install all required packages in R:

```r
install.packages(c(
  "factoextra", "circlize", "dplyr", "tidyr", "plotly", "ggplot2",
  "ggdendro", "tibble", "permute", "mclust", "DescTools", "plot3D",
  "berryFunctions", "plotrix", "gridExtra", "MASS", "shape", "formattable",
  "lattice", "vegan", "ggrepel", "ape", "phytools", "glmm", "wrapr"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# GitHub packages
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("houyunhuang/ggcor")
```

## Data Files

### Input Data Format

**Growth curve files** (`{strain}_p_{library}_r_{replica}.txt`):
- Tab-separated text files
- Contains Time, Temperature, and OD600 readings for 96 wells (A1-H12)
- 72-hour growth measurements

**Chemical library files** (`ps1.csv`, `ps2.csv`):
- CSV format with row labels (A-H) and columns (1-12)
- Contains chemical compound names for each well position

**Strain taxonomy** (`aiden-strain-taxonomy.csv`):
- Contains Species and Strain.id columns
- Maps strain IDs to taxonomic information

**Phylogenetic tree** (`phylo.io_n.nwk`):
- Newick format phylogenetic tree
- Contains evolutionary relationships between 21 strains with branch lengths

## Usage

1. Set the working directory in `Project_Main_Code.R` to your data location
2. Ensure all required data files are in place
3. Run the R script section by section:
   - Section 1: Load packages and define functions
   - Section 2: Read data and calculate AUC values
   - Section 3: Integrate chemical data
   - Section 4: Calculate AUC ratios
   - Section 5: Clustering analysis
   - Section 6: PCA analysis
   - Section 7-9: Phylogenetic analysis

## Experimental Design

- **Strains**: 21 bacterial strains (IDs: 74, 85, 88, 100, 186, 322, 331, 333, 350, 353, 371, 374, 380, 390, 398, 436, 442, 448, 487, 527, 565)
- **Chemical Libraries**: 2 libraries (PS1, PS2) with ~96 compounds each
- **Replicates**: 3 technical replicates per condition
- **Plate Format**: 96-well microplates (rows A-H, columns 1-12)
- **Controls**: DMSO solvent controls and empty wells

## Output

The analysis generates:
- Growth curve visualizations
- Heatmaps of AUC values
- Chemical clustering dendrograms
- PCA biplots
- Phylogenetic trees with phenotypic mapping
- Mantel test results for phylogeny-phenotype correlation

## Citation

If you use this code or data, please cite the accompanying paper:
"A Broad Examination on Bacterial Responses to a Wide Range of Biocide Chemicals and Exploration of Potential 'Biosensor' Strategy"

## License

This project is for academic research purposes.

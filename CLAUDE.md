# CLAUDE.md - Bacterial Response to Biocides Project

## Project Overview

This is a **biological science research project** analyzing bacterial growth responses to biocide chemicals (pesticides, herbicides, fungicides, insecticides). The codebase is an R-based analysis pipeline (`Project_Main_Code.R`) processing 96-well microplate OD600 growth curve data for 21 bacterial strains (spanning 4 phyla: Proteobacteria, Firmicutes, Actinobacteria, Bacteroidetes; 15 genera) across ~190 chemical compounds, with downstream statistical and phylogenetic analyses.

## Scientific Rigor Requirements

### Absolute Rules

- **Never fabricate, infer, or speculate about experimental data or results.** Only state what the data explicitly shows.
- **Never modify raw data files** in `growth-curve/` or `data/` unless explicitly instructed. These are primary experimental records.
- **Never alter statistical test parameters** (significance thresholds, permutation counts, distance metrics) without explicit justification and user approval.
- **Preserve the exact sample sizes, strain IDs, and chemical names** as they appear in the original data. Do not rename, merge, or drop samples without documented scientific rationale.
- **Do not round or truncate numerical results** beyond what is standard for the measurement precision (OD600 readings to 2-3 decimal places; AUC values retain full precision).

### Scientific Communication

- Use precise biological and statistical terminology. Refer to organisms by their actual species names from `aiden-strain-taxonomy.csv` (e.g., *Pseudomonas baetica* strain 100, *Aeromonas sobria* strain 374, *Hymenobacter glaciei* strain 442). The old simplified labels ("Pseudomonas sp.", "Bacillus sp.", "Enterobacter sp.") are incorrect — the 21 strains span 15 genera across 4 phyla. Refer to chemicals by their standard names as listed in the data files.
- Distinguish clearly between **statistically significant** and **biologically meaningful** differences.
- When describing results, state the statistical test used, the test statistic, p-value, and effect size where applicable.
- Never claim causation from correlational analyses (e.g., Mantel test results show correlation, not causation).
- Acknowledge limitations: incomplete replicates, missing strains (302, 306 excluded from main analysis), position effects in microplates.

### Data Integrity

- The plate layout rotation scheme (Section 3.5) corrects for microplate position effects across replicates 2 and 3. Do not bypass or simplify this correction.
- DMSO is the solvent control. AUC ratios are normalized to per-strain mean DMSO AUC. Do not use alternative normalization without justification.
- Known data exclusions are documented in the code (e.g., strain 100, library 1, replica 3, column 8 — experimental error). Maintain these exclusions.
- Strains 85, 398, 436 have non-standard file naming (`_85_p_...`, `_398_p_...`, `_436_p_...`) and are loaded separately. Account for this when modifying data loading.
- Strains 331, 371, 74 use pre-computed spline fits from `spline-fits.csv` (contributed by a collaborator). Their raw growth curves are not in the repository.

## Code Architecture

### Main Script: `Project_Main_Code.R` (~880 lines)

The script is organized into sequential sections:

| Section | Purpose | Key Methods |
|---------|---------|-------------|
| 1 | Setup, package loading, utility functions | `remove.temp()`, `clean.time()`, `AreaUnderSpline()` |
| 2 | Data loading and AUC calculation | Spline integration over 1-73 hours |
| 2.5 | Growth curve visualization | ggplot2 faceted plots |
| 3 | Chemical library mapping | Well-to-chemical assignment for PS1/PS2 |
| 3.5 | Plate rotation correction | Position-effect correction for replicates 2 & 3 |
| 4 | AUC ratio calculation | Per-strain DMSO-normalized ratios |
| 5 | Statistical testing & heatmaps | Dunnett's test (vs DMSO control), ComplexHeatmap |
| 5.5 | Data filtering | Remove chemicals with no significant effect on any strain |
| 6 | Hierarchical clustering | Euclidean distance, average linkage, k=7 |
| 6.1 | Model-based clustering | mclust EM algorithm, Gaussian mixtures, MclustDR |
| 7 | PCA ordination | `prcomp()` on log10-transformed AUC ratios |
| 8 | Mantel test | Kendall correlation, 9999 permutations, phylo vs pheno distance |
| 9 | Phylogenetic signal | Pagel's lambda, `contMap`, phylomorphospace |

### Data Flow

```
Raw OD600 (.txt) --> remove temperature --> clean time format
    --> pivot to long format --> spline AUC (1-73h)
    --> map chemicals to wells --> correct plate rotation
    --> normalize to DMSO (AUC ratio)
    --> Dunnett's test --> filter significant chemicals
    --> clustering / PCA / Mantel test / phylogenetic signal
```

### Key Data Files

- `growth-curve/{strain}_p_{lib}_r_{rep}.txt` — Raw 96-well OD600 time series (tab-separated)
- `data/ps1.csv`, `data/ps2.csv` — Chemical-to-well mapping for two libraries
- `data/chemicaldetails.csv` — Chemical metadata (family, target organism type)
- `data/aiden-strain-taxonomy.csv` — Strain ID to species mapping
- `data/phylo.io_n.nwk` — Newick format phylogenetic tree (21 strains)
- `data/spline-fits.csv` — Pre-computed spline AUC for all strains and wells (6049 rows, includes compound metadata)
- `data/chemical_table_final.csv` — Clustering analysis output (Cluster assignment, Target, Dir1/Dir2 coordinates)
- `data/community-growth-curves.csv` — Community-level consolidated growth data (32 MB, git-ignored)
- `data/isolate-growth-curves.csv` — Isolate-level consolidated growth data (325 MB, git-ignored)
- `data/wrksp.RData` — Saved R workspace with intermediate computation results (git-ignored)
- `legacy/` — Earlier modular R scripts from development (auc.R, dmsotest.R, mantel-test.R, tree.R, etc.)

## Working Conventions

- The script uses hardcoded `setwd()` paths that must be updated per machine. The primary working directory should be set to the `growth-curve/` folder, with `data/` files referenced by relative paths from there.
- All strain IDs are character strings (e.g., `"88"`, `"100"`), not integers. Maintain this throughout.
- Chemical names use uppercase as they appear in the library CSV files (e.g., `"ABAMECTIN"`, `"DMSO"`). Do not change case.
- The analysis pipeline is designed to run section-by-section interactively in RStudio, not as a batch script.

## What Not To Do

- Do not add machine learning models or predictive analytics unless explicitly requested. This is a descriptive/exploratory study.
- Do not auto-generate figures or tables that might be mistaken for actual experimental results.
- Do not suggest alternative statistical methods without explaining the assumptions, trade-offs, and how they compare to the existing approach.
- Do not install or suggest packages beyond those listed in `REQUIREMENTS.txt` without user approval.
- Do not commit generated output files (plots, CSVs) to the repository — they are excluded by `.gitignore`.

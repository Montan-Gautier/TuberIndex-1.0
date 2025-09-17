# TuberIndex 1.0
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17077591.svg)](https://doi.org/10.5281/zenodo.17077591)

**Author:** Montan Gautier  
**Date:** September 16, 2025  

TuberIndex is the first open-access dataset compiling literature on truffles (*Tuberaceae*) and their ecological interactions.

The **TuberIndex 1.0 dataset** is available on Zenodo: [https://doi.org/10.5281/zenodo.17077591].
It consists of three UTF-8 encoded CSV files and two complementary files.  

---

## Dataset overview

- **493 documents** (French-language, 17th–21st century).  
- **3,508 interaction records**:  
  - 26 truffle taxa (*Tuber* genus, 24 species, 1 subspecies)  
  - 418 plant taxa (379 Angiosperms, 29 Gymnosperms, 2 Pteridophytes, 8 Plantae)  
  - 53 fungal taxa (19 Ascomycetes, 34 Basidiomycetes)  
- **1,283 taxon names** cited in the literature, each matched with a validated scientific name.

---

## R scripts

This repository provides two R scripts to ensure **transparency and reproducibility**:

- **`TuberIndex_Library.R`**  
  Processes raw bibliographic metadata into the standardized file `TuberIndex_Library.csv`.

- **`TuberIndex_Taxonomy.R`**  
  Standardizes common names, matches them with scientific names, and outputs the reference file `TuberIndex_Taxonomy.csv`.

Both scripts use open-source R packages and correspond to the raw data archived in Zenodo.

---

## Citation

If you use TuberIndex 1.0, please cite:  
Gautier M., Taschen E., Lescureux N., Richard F. (2025). *TuberIndex 1.0, a dataset of ecological interactions from five centuries of French literature on Tuberaceae*. Scientific Data (submitted). Zenodo. https://doi.org/10.5281/zenodo.17077591

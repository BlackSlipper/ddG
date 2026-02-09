# ddG
Code Repository for ddG paper

# Analysis of Protein Stability and Amino Acid Transitions in Cancer Genomics
This repository contains Python scripts used for identifying significant protein stability changes ($\Delta\Delta G$) and analyzing amino acid substitution patterns in cancer-related datasets (COSMIC) versus population-based datasets (gnomAD).

### 1. Project Overview
The analysis is focused on:
- Gene-Specific Thresholding (GST): Defining stability thresholds for each protein based on the distribution of variants in the general population (gnomAD) to identify extreme outliers in cancer (COSMIC).
- Amino Acid Transition Matrix: Categorizing amino acids into four chemical classes (Nonpolar, Polar, Positive, Negative) and quantifying the frequency of class-to-class transitions.

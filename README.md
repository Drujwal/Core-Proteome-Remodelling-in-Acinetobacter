# Core Proteome Remodelling in Acinetobacter

This repository contains the complete bioinformatics pipeline and statistical suite for the phylogenomic analysis of 67 *Acinetobacter* species. 

## 📌 Project Overview
Our research identifies a "Biophysical Bottleneck" and "Mass Compression" within resistant lineages of *Acinetobacter*. This repository provides the tools to replicate the transition from raw DNA sequences to high-confidence statistical validations (p < 0.0001).

## 📂 Repository Structure
- `core.py`: Extracts conserved core genes from 67 *Acinetobacter* assemblies.
- `stitch.py`: Concatenates core gene alignments into a single phylogenomic super-matrix.
- `translate.py`: Translates DNA sequences to proteins using **Bacterial Genetic Code (Table 11)**.
- `biophysical_stats.py`: Analyzes pI, MW, and CDS loading across Resistant (R) vs Sensitive (S) cohorts.
- `categorical_analysis.py`: Performs Fisher’s Exact Test for genomic resistome partitioning.
- `evolutionary_dynamics.py`: Automated suite for confirming Biophysical Bottlenecks and Mass Compression.
- `point_results.csv`: Complete biophysical and genomic dataset for the 67 isolates.
- `requirements.txt`: Python library dependencies (Pandas, Scipy, Numpy).

## 🚀 How to Run the Analysis

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/Drujwal/Core-Proteome-Remodelling-in-Acinetobacter.git](https://github.com/Drujwal/Core-Proteome-Remodelling-in-Acinetobacter.git)
   2. Install Dependencies:
   pip install -r requirements.txt
   3. Step 1: Translate Sequences (DNA to Protein):
   python translate.py
   4.Step 2: Run Statistical Suite:
   python evolutionary_dynamics.py
   Key Validations IncludedIsoelectric Point (pI) Bottleneck: Confirmed via Levene’s Test for variance compression in resistant clades.Molecular Weight (MW) Compression: Validated via Mann-Whitney U test (reduction in average mass).Resistome Partitioning: Validated via Fisher’s Exact Test (p < 0.0001).
   ## 📂 Data Availability
The cds for the 67 *Acinetobacter* isolates (used as input for `core.py` and `stitch.py`) are retrieved from the **EnsemblBacteria (http://bacteria.ensembl.org/index.html)**. 

To maintain repository performance, the large-scale sequence alignments and the final `point_results.csv` are available from the **corresponding author upon reasonable request** or will be released upon formal publication.
   🎓 AffiliationDepartment of Biochemistry, School of Bioengineering and Biosciences,Lovely Professional University (LPU), Punjab, India.

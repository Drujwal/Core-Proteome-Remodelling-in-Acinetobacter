from pathlib import Path
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd

# ============================================================
# USER SETTINGS
# ============================================================

FOLDER = Path(
    r"C:\Users\xenou\OneDrive\Desktop\Major REVISION AcinetobActer\CALCULATE"
)

# Input files generated from the proteome reconstruction script
GENOME_SUMMARY_CSV = FOLDER / "RECALCULATED_Genome_Summary.csv"

# File containing ARG phenotype mapping (must contain genome identifier and ARG status)
# If your phenotype metadata is in a separate file, specify its path here:
METADATA_CSV = FOLDER / "ARG_Phenotypes.csv"

# Outgroup species name to exclude
OUTGROUP_NAME = "Alkanindiges hydrocarboniclasticus"

# Output summary file
STAT_RESULTS_TXT = FOLDER / "MULTIVARIABLE_REGRESSION_RESULTS.txt"


# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

if not GENOME_SUMMARY_CSV.exists():
    raise FileNotFoundError(f"Could not find: {GENOME_SUMMARY_CSV}")

genome_df = pd.read_csv(GENOME_SUMMARY_CSV)

# Load metadata if separate, or ensure 'ARG_Status' exists
if METADATA_CSV.exists():
    metadata_df = pd.read_csv(METADATA_CSV)
    # Merge on Genome_File or Species column
    df = pd.merge(genome_df, metadata_df, on="Genome_File", how="inner")
else:
    df = genome_df.copy()

# Ensure 'ARG_Status' column exists (0 for ARG-negative, 1 for ARG-positive)
if "ARG_Status" not in df.columns:
    raise KeyError(
        "Column 'ARG_Status' not found. Please ensure your dataset includes an 'ARG_Status' "
        "column with binary values (1 = ARG-positive, 0 = ARG-negative) or group labels."
    )

# Filter out the outgroup
df_filtered = df[
    ~df["Species"].str.contains(OUTGROUP_NAME, case=False, na=False)
].copy()

print(f"Total genomes after excluding outgroup: {len(df_filtered)}")


# ============================================================
# MULTIVARIABLE LINEAR REGRESSION
# ============================================================

# Dependent variable: Avg_MolWt_Da
# Predictor: C(ARG_Status) [Categorical: 0 vs 1]
# Covariate: Total_Proteins
model_formula = "Avg_MolWt_Da ~ C(ARG_Status, Treatment(reference=0)) + Total_Proteins"

model = smf.ols(formula=model_formula, data=df_filtered).fit()


# ============================================================
# EXTRACT KEY STATISTICS
# ============================================================

# Get parameters for ARG_Status[T.1]
arg_coef = model.params["C(ARG_Status, Treatment(reference=0))[T.1]"]
arg_pvalue = model.pvalues["C(ARG_Status, Treatment(reference=0))[T.1]"]
conf_int = model.conf_int().loc["C(ARG_Status, Treatment(reference=0))[T.1]"]

summary_text = f"""
============================================================
MULTIVARIABLE LINEAR REGRESSION ANALYSIS
============================================================
Formula: {model_formula}
Number of observations: {int(model.nobs)}

KEY FINDINGS (ARG-Positive vs. ARG-Negative):
------------------------------------------------------------
Adjusted Effect Size (Mean MW Diff) : {arg_coef:.2f} Da
95% Confidence Interval             : [{conf_int[0]:.2f}, {conf_int[1]:.2f}]
p-value                            : {arg_pvalue:.4f}
============================================================

FULL REGRESSION MODEL SUMMARY:
{model.summary()}
"""

print(summary_text)

# Save text summary to disk
with open(STAT_RESULTS_TXT, "w") as f:
    f.write(summary_text)

print(f"Results saved to: {STAT_RESULTS_TXT}")
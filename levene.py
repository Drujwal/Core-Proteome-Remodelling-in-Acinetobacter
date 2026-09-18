import pandas as pd
from scipy.stats import levene
from pathlib import Path

# ============================================================
# LEVENE'S TEST FOR EQUALITY OF VARIANCES
# Dataset: Acinetobacter ARG-positive (+) vs ARG-negative (-)
# Outgroup has already been removed.
#
# Variables tested:
# 1. Avg_pI
# 2. Avg_MolWt_Da
# 3. Total_Proteins
# ============================================================

# Folder containing this script and Excel file
folder = Path(__file__).resolve().parent

# Input Excel file
input_file = folder / "Tests.xlsx"

# Output result file
output_file = folder / "Levene_results.txt"

# ------------------------------------------------------------
# Read Excel file
# ------------------------------------------------------------

df = pd.read_excel(input_file)

# Remove accidental spaces from column names
df.columns = df.columns.str.strip()

# ------------------------------------------------------------
# Check required columns
# ------------------------------------------------------------

required_columns = [
    "Species",
    "ARG Phenotype",
    "Avg_pI",
    "Avg_MolWt_Da",
    "Total_Proteins"
]

missing = [col for col in required_columns if col not in df.columns]

if missing:
    print("ERROR: The following required columns are missing:")
    for col in missing:
        print(f"  - {col}")

    print("\nColumns actually found:")
    print(list(df.columns))

    input("\nPress Enter to exit...")
    raise SystemExit

# ------------------------------------------------------------
# Clean ARG phenotype
# ------------------------------------------------------------

df["ARG Phenotype"] = (
    df["ARG Phenotype"]
    .astype(str)
    .str.strip()
)

# Keep only ARG-positive (+) and ARG-negative (-)
df = df[df["ARG Phenotype"].isin(["+", "-"])].copy()

# ------------------------------------------------------------
# Variables to test
# ------------------------------------------------------------

variables = [
    "Avg_pI",
    "Avg_MolWt_Da",
    "Total_Proteins"
]

# ------------------------------------------------------------
# Convert variables to numeric
# ------------------------------------------------------------

for variable in variables:
    df[variable] = pd.to_numeric(
        df[variable],
        errors="coerce"
    )

# ------------------------------------------------------------
# Function for descriptive statistics
# ------------------------------------------------------------

def descriptive_stats(series):
    return {
        "n": len(series),
        "mean": series.mean(),
        "median": series.median(),
        "Q1": series.quantile(0.25),
        "Q3": series.quantile(0.75),
        "IQR": series.quantile(0.75) - series.quantile(0.25),
        "SD": series.std(ddof=1),
        "variance": series.var(ddof=1),
        "minimum": series.min(),
        "maximum": series.max()
    }

# ------------------------------------------------------------
# Run Levene's test for all three variables
# ------------------------------------------------------------

results = {}

for variable in variables:

    # Remove missing values for this variable
    data = df.dropna(subset=[variable]).copy()

    # Separate ARG-positive and ARG-negative groups
    arg_positive = data.loc[
        data["ARG Phenotype"] == "+",
        variable
    ]

    arg_negative = data.loc[
        data["ARG Phenotype"] == "-",
        variable
    ]

    # Levene's test using median-centered version
    statistic, p_value = levene(
        arg_positive,
        arg_negative,
        center="median"
    )

    # Store results
    results[variable] = {
        "positive": arg_positive,
        "negative": arg_negative,
        "statistic": statistic,
        "p_value": p_value,
        "positive_stats": descriptive_stats(arg_positive),
        "negative_stats": descriptive_stats(arg_negative)
    }

# ------------------------------------------------------------
# Write results to text file
# ------------------------------------------------------------

with open(output_file, "w", encoding="utf-8") as f:

    f.write("LEVENE'S TEST RESULTS\n")
    f.write("=" * 75 + "\n\n")

    f.write("DATASET\n")
    f.write("-" * 75 + "\n")
    f.write(f"Input file: {input_file.name}\n")
    f.write(f"Total isolates analyzed: {len(df)}\n")
    f.write("Outgroup: Removed before analysis\n\n")

    # --------------------------------------------------------
    # Group sizes
    # --------------------------------------------------------

    f.write("GROUP SIZES\n")
    f.write("-" * 75 + "\n")
    f.write(
        f"ARG-positive (+): "
        f"{sum(df['ARG Phenotype'] == '+')}\n"
    )
    f.write(
        f"ARG-negative (-): "
        f"{sum(df['ARG Phenotype'] == '-')}\n"
    )

    # --------------------------------------------------------
    # Results for each variable
    # --------------------------------------------------------

    for variable in variables:

        result = results[variable]

        f.write("\n\n")
        f.write("=" * 75 + "\n")
        f.write(f"VARIABLE: {variable}\n")
        f.write("=" * 75 + "\n\n")

        # ARG-positive statistics
        f.write("ARG-POSITIVE (+)\n")
        f.write("-" * 40 + "\n")

        for key, value in result["positive_stats"].items():
            f.write(f"{key}: {value:.6f}\n")

        # ARG-negative statistics
        f.write("\nARG-NEGATIVE (-)\n")
        f.write("-" * 40 + "\n")

        for key, value in result["negative_stats"].items():
            f.write(f"{key}: {value:.6f}\n")

        # Levene test
        f.write("\nLEVENE'S TEST\n")
        f.write("-" * 40 + "\n")
        f.write("Test: Levene's test for equality of variances\n")
        f.write("Center: median\n")
        f.write(
            f"Levene statistic: "
            f"{result['statistic']:.6f}\n"
        )
        f.write(
            f"p-value: "
            f"{result['p_value']:.10g}\n"
        )

        # Interpretation
        f.write("\nINTERPRETATION\n")
        f.write("-" * 40 + "\n")

        if result["p_value"] < 0.05:
            f.write(
                "The variances differ significantly between the "
                "ARG-positive and ARG-negative groups (p < 0.05).\n"
            )
        else:
            f.write(
                "No statistically significant difference in variance "
                "was detected between the ARG-positive and ARG-negative "
                "groups (p >= 0.05).\n"
            )

    f.write("\n\n")
    f.write("=" * 75 + "\n")
    f.write("END OF LEVENE'S TEST ANALYSIS\n")
    f.write("=" * 75 + "\n")

# ------------------------------------------------------------
# Print summary to Command Prompt
# ------------------------------------------------------------

print("=" * 75)
print("LEVENE'S TEST ANALYSIS COMPLETED")
print("=" * 75)

print(
    f"ARG-positive (+): "
    f"{sum(df['ARG Phenotype'] == '+')}"
)

print(
    f"ARG-negative (-): "
    f"{sum(df['ARG Phenotype'] == '-')}"
)

print("\nRESULTS:")

for variable in variables:

    result = results[variable]

    print("\n" + "-" * 60)
    print(variable)
    print("-" * 60)

    print(
        f"Levene statistic: "
        f"{result['statistic']:.6f}"
    )

    print(
        f"p-value: "
        f"{result['p_value']:.10g}"
    )

    if result["p_value"] < 0.05:
        print("Result: SIGNIFICANT (p < 0.05)")
    else:
        print("Result: NOT SIGNIFICANT (p >= 0.05)")

print("\n" + "=" * 75)
print(f"Complete results saved to:")
print(output_file)
print("=" * 75)

input("\nPress Enter to exit...")
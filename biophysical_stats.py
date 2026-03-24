import pandas as pd
from scipy import stats

# Load your dataset (ensure 'acinetobacter_data.csv' is in the same folder)
df = pd.read_csv('acinetobacter_data.csv')

# Split into cohorts based on your Red (Resistant) and Green (Sensitive) markers
res = df[df['Cohort'] == 'Resistant']
sens = df[df['Cohort'] == 'Sensitive']

print("--- Biophysical Statistical Analysis ---")

# --- Test 1: Levene’s Test (pI Variance / Bottleneck) ---
stat_pI, p_levene = stats.levene(res['pI'], sens['pI'])
print(f"Levene’s Test (pI Variance) p-value: {p_levene}")

# --- Test 2: Mann-Whitney U (Molecular Weight) ---
stat_mw, p_mw = stats.mannwhitneyu(res['MW'], sens['MW'], alternative='two-sided')
print(f"Mann-Whitney U (Molecular Weight) p-value: {p_mw}")

# --- Test 3: Mann-Whitney U (Total Protein Count) ---
stat_prot, p_prot = stats.mannwhitneyu(res['Protein_Count'], sens['Protein_Count'], alternative='two-sided')
print(f"Mann-Whitney U (Protein Count) p-value: {p_prot}")
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from scipy import stats
import os

# --- SETUP ---
folder_path = r'D:\PROTEIN'
target_file = 'point_results.csv'
full_path = os.path.join(folder_path, target_file)

def run_q1_stats():
    print("--- Acinetobacter Evolutionary Dynamics: Statistical Suite ---")
    
    if not os.path.exists(full_path):
        print(f"ERROR: Could not find '{target_file}' in '{folder_path}'")
        return

    try:
        # Load CSV (handles various delimiters)
        df = pd.read_csv(full_path, sep=None, engine='python')
        df.columns = df.columns.str.strip()

        # Filtering based on 'R' (Resistant) and 'S' (Sensitive)
        res = df[df['Phenotype'].str.strip().str.upper() == 'R']
        sen = df[df['Phenotype'].str.strip().str.upper() == 'S']
        
        if len(res) == 0 or len(sen) == 0:
            print(f"ERROR: No groups found. Phenotypes in file are: {df['Phenotype'].unique()}")
            return
            
    except Exception as e:
        print(f"ERROR: {e}")
        return

    # --- METRICS TO TEST ---
    metrics = [
        ('Avg_pI', 'Global Isoelectric Point (pI)'),
        ('Avg_MolWt_Da', 'Molecular Weight (MW)'),
        ('Total_Proteins', 'Total Protein Count (CDS)')
    ]

    print(f"Analyzing: {len(res)} Resistant (R) vs {len(sen)} Sensitive (S) isolates.\n")

    for col_name, label in metrics:
        r_data = res[col_name].dropna().values
        s_data = sen[col_name].dropna().values

        # Mann-Whitney U (Tests shift in the Median)
        u_stat, p_u = stats.mannwhitneyu(r_data, s_data, alternative='two-sided')
        
        # Levene's/Brown-Forsythe (Tests difference in Variance/Spread)
        lev_stat, p_lev = stats.levene(r_data, s_data, center='median')

        print(f"== {label} ==")
        print(f"   R Median: {np.median(r_data):.2f} | S Median: {np.median(s_data):.2f}")
        print(f"   Mann-Whitney p-value: {p_u:.4f}")
        print(f"   Levene Variance p-value: {p_lev:.4f}")
        
        if label == 'Global Isoelectric Point (pI)' and p_lev < 0.05:
            print("   >>> BIOPHYSICAL BOTTLENECK CONFIRMED (Significant pI Compression)")
        
        if label == 'Molecular Weight (MW)' and p_u < 0.05:
            print("   >>> MASS COMPRESSION CONFIRMED (Significant MW Reduction)")

        if label == 'Total_Proteins' and p_u < 0.05:
            print("   >>> GENOMIC EXPANSION CONFIRMED (Significant CDS Increase)")
        print("-" * 55)

if __name__ == "__main__":
    run_q1_stats()
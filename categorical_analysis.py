import numpy as np
from scipy import stats

# --- Fisher’s Exact Test (Resistome Gap) ---
# Table: [Clinical_Resistant, Clinical_Sensitive], [Env_Resistant, Env_Sensitive]
# Replace A, B, C, D with your real species counts from the 67 isolates
A, B, C, D = 38, 0, 0, 29 

contingency_table = [[A, B], [C, D]]

odds_ratio, p_value = stats.fisher_exact(contingency_table)

print("--- Fisher’s Exact Test Results ---")
print(f"Table: {contingency_table}")
print(f"P-value: {p_value}")

if p_value < 0.05:
    print("Result: Statistically significant genomic divide between niches.")
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.spatial.distance import braycurtis
import os

# ── Sample metadata ──────────────────────────────────────
samples = {
    'SRR19064316': 'PD',
    'SRR19064317': 'PD',
    'SRR19064453': 'HC',
    'SRR19064454': 'HC',
    'SRR19064809': 'HC',
}

bracken_dir = os.path.expanduser('~/metagenomics_project/bracken_output')
results_dir = os.path.expanduser('~/metagenomics_project/results')
os.makedirs(results_dir, exist_ok=True)

# ── Load bracken files ────────────────────────────────────
def load_bracken(srr):
    path = os.path.join(bracken_dir, f'{srr}_bracken.txt')
    df = pd.read_csv(path, sep='\t')
    df = df[['name', 'fraction_total_reads']].set_index('name')
    df.columns = [srr]
    return df

print("Loading bracken files...")
dfs = [load_bracken(srr) for srr in samples.keys()]
abundance = pd.concat(dfs, axis=1).fillna(0)
print(f"Loaded {len(abundance)} species across {len(samples)} samples")

# ── Alpha diversity (Shannon) ─────────────────────────────
def shannon(row):
    row = row[row > 0]
    return -np.sum(row * np.log(row))

print("\nCalculating alpha diversity...")
alpha = pd.DataFrame({
    'SRR': list(samples.keys()),
    'Group': list(samples.values()),
    'Shannon': [shannon(abundance[srr]) for srr in samples.keys()]
})
print(alpha.to_string(index=False))

# ── Save alpha diversity stats ────────────────────────────
alpha.to_csv(f'{results_dir}/alpha_diversity.csv', index=False)

pd_shannon = alpha[alpha['Group']=='PD']['Shannon'].values
hc_shannon = alpha[alpha['Group']=='HC']['Shannon'].values
stat, p = mannwhitneyu(pd_shannon, hc_shannon, alternative='two-sided')
print(f"\nMann-Whitney U test PD vs HC Shannon: p={p:.4f}")

with open(f'{results_dir}/diversity_stats.txt', 'w') as f:
    f.write("ALPHA DIVERSITY — SHANNON INDEX\n")
    f.write("="*40 + "\n")
    f.write(alpha.to_string(index=False))
    f.write(f"\n\nPD mean Shannon: {pd_shannon.mean():.4f}")
    f.write(f"\nHC mean Shannon: {hc_shannon.mean():.4f}")
    f.write(f"\nMann-Whitney U p-value: {p:.4f}\n")

# ── Plot alpha diversity ──────────────────────────────────
plt.figure(figsize=(6,5))
sns.boxplot(data=alpha, x='Group', y='Shannon', palette={'PD':'#E74C3C','HC':'#2ECC71'})
sns.stripplot(data=alpha, x='Group', y='Shannon', color='black', size=8, jitter=False)
plt.title('Alpha Diversity (Shannon Index)\nPD vs Healthy Controls')
plt.ylabel('Shannon Index')
plt.xlabel('')
plt.tight_layout()
plt.savefig(f'{results_dir}/alpha_diversity.png', dpi=150)
plt.close()
print("Alpha diversity plot saved.")

# ── Top 15 species comparison ─────────────────────────────
pd_samples = [s for s,g in samples.items() if g=='PD']
hc_samples = [s for s,g in samples.items() if g=='HC']

pd_mean = abundance[pd_samples].mean(axis=1)
hc_mean = abundance[hc_samples].mean(axis=1)

comparison = pd.DataFrame({
    'PD_mean': pd_mean,
    'HC_mean': hc_mean
})
comparison['log2FC'] = np.log2((comparison['PD_mean']+1e-6) / (comparison['HC_mean']+1e-6))
comparison = comparison.sort_values('log2FC')

top_depleted = comparison.head(10)
top_enriched = comparison.tail(10)
top15 = pd.concat([top_depleted, top_enriched])

plt.figure(figsize=(10,8))
colors = ['#2ECC71' if x < 0 else '#E74C3C' for x in top15['log2FC']]
plt.barh(top15.index, top15['log2FC'], color=colors)
plt.axvline(0, color='black', linewidth=0.8)
plt.xlabel('Log2 Fold Change (PD / HC)')
plt.title('Top Depleted & Enriched Species in PD vs HC')
plt.tight_layout()
plt.savefig(f'{results_dir}/taxa_comparison.png', dpi=150)
plt.close()
print("Taxa comparison plot saved.")

print(f"\n✅ All results saved to {results_dir}/")
print("Files: alpha_diversity.csv, diversity_stats.txt, alpha_diversity.png, taxa_comparison.png")

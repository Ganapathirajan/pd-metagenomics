import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import braycurtis
from scipy.spatial.distance import squareform
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

# ── Bray-Curtis distance matrix ───────────────────────────
print("Calculating Bray-Curtis distances...")
srr_list = list(samples.keys())
n = len(srr_list)
dist_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        dist_matrix[i][j] = braycurtis(
            abundance[srr_list[i]],
            abundance[srr_list[j]]
        )

dist_df = pd.DataFrame(dist_matrix,
    index=srr_list, columns=srr_list)
dist_df.to_csv(f'{results_dir}/bray_curtis_matrix.csv')
print("\nBray-Curtis Distance Matrix:")
print(dist_df.round(3).to_string())

# ── PCoA ─────────────────────────────────────────────────
print("\nRunning PCoA...")
# Double centre the distance matrix
D = dist_matrix ** 2
n = D.shape[0]
H = np.eye(n) - np.ones((n,n))/n
B = -0.5 * H @ D @ H

eigenvalues, eigenvectors = np.linalg.eigh(B)
idx = np.argsort(eigenvalues)[::-1]
eigenvalues = eigenvalues[idx]
eigenvectors = eigenvectors[:, idx]

# Take positive eigenvalues only
pos = eigenvalues > 0
coords = eigenvectors[:, pos] * np.sqrt(eigenvalues[pos])

explained = eigenvalues[pos] / eigenvalues[pos].sum() * 100

# ── Plot PCoA ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8,6))

colors = {'PD': '#E74C3C', 'HC': '#2ECC71'}
markers = {'PD': 'o', 'HC': 's'}

for i, srr in enumerate(srr_list):
    group = samples[srr]
    ax.scatter(coords[i,0], coords[i,1],
               c=colors[group],
               marker=markers[group],
               s=150, zorder=5,
               label=group)
    ax.annotate(srr[-3:],
                (coords[i,0], coords[i,1]),
                textcoords="offset points",
                xytext=(8,5), fontsize=8)

# Remove duplicate legend entries
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(),
          title='Group', fontsize=11)

ax.set_xlabel(f'PC1 ({explained[0]:.1f}%)', fontsize=12)
ax.set_ylabel(f'PC2 ({explained[1]:.1f}%)', fontsize=12)
ax.set_title('Beta Diversity — PCoA (Bray-Curtis)\nPD vs Healthy Controls', fontsize=13)
ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
ax.axvline(0, color='gray', linewidth=0.5, linestyle='--')
plt.tight_layout()
plt.savefig(f'{results_dir}/beta_diversity_pcoa.png', dpi=150)
plt.close()
print("PCoA plot saved.")
print(f"\nPC1 explains: {explained[0]:.1f}%")
print(f"PC2 explains: {explained[1]:.1f}%")
print("\n✅ Beta diversity complete.")

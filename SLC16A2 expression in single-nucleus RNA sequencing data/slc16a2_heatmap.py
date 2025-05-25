#MIT License

#Copyright (c) 2025 Jonas Rotter

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

adata1 = sc.read_h5ad("data/9082ad42-b5ba-449f-ad10-1b988ac79eaa.h5ad")
adata2 = sc.read_h5ad("data/f9ecb4ba-b033-4a93-b794-05e262dc1f59.h5ad")

def compute_fraction_expressing(adata, gene='SLC16A2'):
    # Confirm gene exists
    if gene not in adata.var['Gene'].values:
        raise ValueError(f"Gene {gene} not found in adata.var['Gene']")

    # Boolean mask and extract expression values
    gene_mask = adata.var['Gene'] == gene
    expr = adata[:, gene_mask].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = expr.flatten()

    # Binary expression status
    adata.obs['expressing_SLC16A2'] = expr > 1

    # Group by cell_type and tissue
    grouped = adata.obs.groupby(['cell_type', 'tissue']).agg(
        total_cells=('expressing_SLC16A2', 'size'),
        expressing_cells=('expressing_SLC16A2', 'sum')
    )

    # Compute fraction expressing
    grouped['fraction_expressing'] = grouped['expressing_cells'] / grouped['total_cells']

    # Reshape
    return grouped['fraction_expressing'].unstack(fill_value=0)

# Compute fraction matrices for each adata
fraction_matrix1 = compute_fraction_expressing(adata1)
fraction_matrix2 = compute_fraction_expressing(adata2)

# Combine the two matrices by summing numerator and denominator before recomputing the fraction
def merge_fraction_matrices(adata1, adata2):
    def get_counts(adata):
        gene_mask = adata.var['Gene'] == 'SLC16A2'
        expr = adata[:, gene_mask].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        expr = expr.flatten()
        adata.obs['expressing_SLC16A2'] = expr > 1
        return adata.obs.groupby(['cell_type', 'tissue']).agg(
            total_cells=('expressing_SLC16A2', 'size'),
            expressing_cells=('expressing_SLC16A2', 'sum')
        )

    counts1 = get_counts(adata1)
    counts2 = get_counts(adata2)

    # Combine counts
    merged = counts1.add(counts2, fill_value=0)

    # Recalculate fraction
    merged['fraction_expressing'] = merged['expressing_cells'] / merged['total_cells']
    return merged['fraction_expressing'].unstack(fill_value=0)

# Final merged fraction matrix
fraction_matrix = merge_fraction_matrices(adata1, adata2)

# Convert to percentage
percent_matrix = fraction_matrix * 100

# Capitalize and sort cell types and tissue names
percent_matrix.index = percent_matrix.index.str.capitalize()
percent_matrix.columns = percent_matrix.columns.str.capitalize()
percent_matrix = percent_matrix.sort_index().sort_index(axis=1)

# Step 2: Custom colormap
colors = [
    (0.0, '#f7fbff'),    # Lightest blue
    (0.3, '#6baed6'),    # Medium blue (at 30%)
    (1.0, '#cb181d')     # Red for high outliers (>30%)
]
custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', colors)

# Step 3: Normalization clipped at 30%
class ClippedNormalize(Normalize):
    def __call__(self, value, clip=None):
        value = np.clip(value, 0, 30) + (value > 30) * (value - 30)
        return super().__call__(value, clip)

# Step 4: Plot
plt.figure(figsize=(14, 10))
ax = sns.heatmap(
    percent_matrix,
    cmap=custom_cmap,
    norm=ClippedNormalize(vmin=0, vmax=100),
    annot=True,
    fmt=".1f",
    linewidths=0.5,
    cbar_kws={'label': 'Percent Expressing SLC16A2'}
)

ax.set_title("SLC16A2 expression (% expressing cells) in the human brain", pad=20, fontsize=22, weight='bold')
ax.set_xlabel("Brain region", fontsize=20, weight='bold')
ax.set_ylabel("Cell type", fontsize=20, weight='bold')
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)

plt.tight_layout()
plt.savefig("heatmap_scl16a2_expression_gt1.png")
plt.savefig("heatmap_scl16a2_expression_gt1.svg")
plt.show()
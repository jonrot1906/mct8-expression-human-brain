#!/usr/bin/env python3
"""SLC16A2 (MCT8) expression in the human brain: main-text figure 7 and robustness checks.

Outputs
    main_figure/            panels a-c as PDF/SVG/PNG plus source_data_*.csv
    robustness_threshold/   raw>1 vs log2(CP10k+1)>1 vs depth-explicit Pearson residuals
    robustness_pseudobulk/  fine-cluster pseudobulk CPM vs per-cell detection rate

MIT License. Copyright (c) 2026 Jonas Rotter.
"""

import gc
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

try:
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
except ImportError:
    sys.exit("Missing dependency: pip install abc_atlas_access\n"
             "(see https://github.com/AllenInstitute/abc_atlas_access )")


# ==================================================================================================
# Configuration
# ==================================================================================================
DOWNLOAD_BASE = Path('XXX')                          # local ABC Atlas cache directory
OUT_DIR = Path('revision')
PINNED_MANIFEST = 'releases/20260415/manifest.json'   # pin for a reproducible published run

EXPR_DIRECTORY = 'WHB-10Xv3'
TAXONOMY_DIRECTORY = 'WHB-taxonomy'
CELL_METADATA_FILE = 'cell_metadata'
MEMBERSHIP_FILE = 'cluster_to_cluster_annotation_membership'
EXPR_MATRICES = {'Neuronal': 'WHB-10Xv3-Neurons/raw',
                 'Non-neuronal': 'WHB-10Xv3-Nonneurons/raw'}

# Column names in the pinned manifest. Verified at start-up by require_columns(), which prints the
# available columns and exits if the manifest ever renames one.
CELL_LABEL_COL = 'cell_label'
REGION_COL = 'anatomical_division_label'
EMBED_COLS = ('x', 'y')
CLUSTER_ALIAS_COL = 'cluster_alias'
GENE_SYMBOL_COL = 'gene_symbol'                  # in .var of the expression matrices
TERM_SET_COL = 'cluster_annotation_term_set_name'
TERM_NAME_COL = 'cluster_annotation_term_name'

GENE = 'SLC16A2'
EXPR_COL = 'SLC16A2'                             # raw counts per nucleus
POS_COL = 'SLC16A2_pos'                          # raw count > EXPR_THRESHOLD
EXPR_THRESHOLD = 1
MIN_CELLS = 20
DROP_UNASSIGNED_FROM_STATS = True                # drop the atlas "no NT assigned" bucket from panel c

PREFERRED_REGION_ORDER = [
    'cerebral cortex', 'hippocampal formation', 'cerebral nuclei', 'thalamic complex',
    'hypothalamus', 'midbrain', 'pons', 'cerebellum', 'myelencephalon', 'spinal cord',
]
TELENCEPHALIC_REGIONS = ['cerebral cortex', 'hippocampal formation', 'cerebral nuclei']

# PV/SST/VIP is not part of the atlas neurotransmitter term set, so it is called from canonical
# markers within atlas-GABAergic neurons only for comparison.
INTERNEURON_MARKER_GENES = {'PV': ['PVALB'], 'SST': ['SST'], 'VIP': ['VIP']}
GABA_SUBTYPE_ORDER = ['PV', 'SST', 'VIP', 'Other GABAergic']

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 12
mpl.rcParams['pdf.fonttype'] = 42               
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['axes.linewidth'] = 1.0
sns.set_style('white')

MAIN_DIR = OUT_DIR / 'main_figure'
THRESHOLD_DIR = OUT_DIR / 'robustness_threshold'
PSEUDOBULK_DIR = OUT_DIR / 'robustness_pseudobulk'
for d in (MAIN_DIR, THRESHOLD_DIR, PSEUDOBULK_DIR):
    d.mkdir(parents=True, exist_ok=True)


def require_columns(available, needed, what):
    missing = [c for c in needed if c not in available]
    if missing:
        sys.exit(f"{what}: expected column(s) {missing} not found. Available: {sorted(available)}")


def save_csv(df, out_dir, name, index=False):
    path = out_dir / name
    df.to_csv(path, index=index)
    print(f'Saved {path}')


def save_figure(fig, out_dir, name, dpi=300):
    base = out_dir / name
    for ext in ('pdf', 'svg', 'png'):
        fig.savefig(f'{base}.{ext}', dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {base}.pdf / .svg / .png')


# ==================================================================================================
# 1. Atlas download, metadata and taxonomy
# ==================================================================================================
DOWNLOAD_BASE.mkdir(parents=True, exist_ok=True)
abc = AbcProjectCache.from_cache_dir(DOWNLOAD_BASE)
abc.load_manifest(PINNED_MANIFEST)
print(f'[abc] manifest: {abc.current_manifest}')

available_matrices = set(abc.list_expression_matrix_files(EXPR_DIRECTORY))
missing_matrices = [f for f in EXPR_MATRICES.values() if f not in available_matrices]
if missing_matrices:
    sys.exit(f'Expression matrices {missing_matrices} not in {EXPR_DIRECTORY}. '
             f'Available: {sorted(available_matrices)}')

cell_meta = abc.get_metadata_dataframe(directory=EXPR_DIRECTORY, file_name=CELL_METADATA_FILE)
if cell_meta.index.name != CELL_LABEL_COL:
    require_columns(set(cell_meta.columns), [CELL_LABEL_COL], 'cell metadata')
    cell_meta = cell_meta.set_index(CELL_LABEL_COL)
require_columns(set(cell_meta.columns), [REGION_COL, CLUSTER_ALIAS_COL, *EMBED_COLS], 'cell metadata')
print(f'[abc] cell metadata: {len(cell_meta):,} nuclei')

membership = abc.get_metadata_dataframe(directory=TAXONOMY_DIRECTORY, file_name=MEMBERSHIP_FILE)
require_columns(set(membership.columns), [CLUSTER_ALIAS_COL, TERM_SET_COL, TERM_NAME_COL], 'taxonomy')


def term_set_series(term_set):
    """cluster_alias -> term name for one taxonomy term set."""
    hit = membership[TERM_SET_COL].astype(str).str.lower() == term_set.lower()
    if not hit.any():
        avail = sorted(membership[TERM_SET_COL].astype(str).unique())
        sys.exit(f"Term set '{term_set}' not found. Available: {avail}")
    sub = membership.loc[hit, [CLUSTER_ALIAS_COL, TERM_NAME_COL]].drop_duplicates(CLUSTER_ALIAS_COL)
    return sub.set_index(CLUSTER_ALIAS_COL)[TERM_NAME_COL]


nt_by_alias = term_set_series('neurotransmitter')
supercluster_by_alias = term_set_series('supercluster')
cluster_name_by_alias = term_set_series('cluster')   # finest level, the pseudobulk unit
print(f'[abc] taxonomy: {nt_by_alias.nunique()} neurotransmitter types, '
      f'{supercluster_by_alias.nunique()} superclusters, {cluster_name_by_alias.nunique():,} clusters')


# ==================================================================================================
# 2. Expression matrices and per-nucleus annotation
# ==================================================================================================
def get_expression(adata, gene):
    """One gene's column of .X as a dense 1-D vector: the only place the matrix is densified."""
    mask = (adata.var[GENE_SYMBOL_COL] == gene).to_numpy()
    if mask.sum() == 0:
        raise ValueError(f"Gene {gene} not found in .var['{GENE_SYMBOL_COL}']")
    expr = adata[:, mask].X
    if hasattr(expr, 'toarray'):
        expr = expr.toarray()
    return np.asarray(expr).ravel()


def attach_annotations(adata, compartment):
    meta = cell_meta.reindex(adata.obs_names)
    n_missing = int(meta[CLUSTER_ALIAS_COL].isna().sum())
    if n_missing:
        print(f'[{compartment}] WARNING: {n_missing:,}/{adata.n_obs:,} nuclei have no metadata match '
              '(matrix barcodes and cell_label must share one convention).')
    alias = meta[CLUSTER_ALIAS_COL]
    adata.obs['tissue'] = meta[REGION_COL].astype(str).str.lower().values
    adata.obs['cluster_alias'] = alias.values
    adata.obs['supercluster'] = alias.map(supercluster_by_alias).astype('object').values
    adata.obs['neurotransmitter_type'] = alias.map(nt_by_alias).astype('object').values
    adata.obsm['X_umap'] = meta[list(EMBED_COLS)].to_numpy(dtype=float)


adatas = {}
for compartment, fname in EXPR_MATRICES.items():
    path = abc.get_file_path(directory=EXPR_DIRECTORY, file_name=fname)
    adata = sc.read_h5ad(path)
    adata.raw = None
    adata.layers.clear()
    gc.collect()
    require_columns(set(adata.var.columns), [GENE_SYMBOL_COL], f'{compartment} .var')
    slc = get_expression(adata, GENE)
    probe = slc[:5000]
    if not np.allclose(probe, np.round(probe)):
        sys.exit(f'{compartment}: {GENE} values are not integers, so this is not the /raw matrix. '
                 'The positivity threshold is defined on raw counts.')
    adata.obs[EXPR_COL] = slc
    adata.obs[POS_COL] = slc > EXPR_THRESHOLD
    attach_annotations(adata, compartment)
    adatas[compartment] = adata
    print(f'[load] {compartment}: {adata.n_obs:,} nuclei x {adata.n_vars:,} genes')

adata_neuronal, adata_nonneuronal = adatas['Neuronal'], adatas['Non-neuronal']
neu_obs, nonneu_obs = adata_neuronal.obs, adata_nonneuronal.obs


# ==================================================================================================
# 3. Neurotransmitter labels and GABAergic PV/SST/VIP subtypes
# ==================================================================================================
neu_obs['neurotransmitter_type'] = neu_obs['neurotransmitter_type'].fillna('Unassigned')
is_gaba = neu_obs['neurotransmitter_type'].astype(str).str.contains('GABA').to_numpy()


def marker_mean(adata, genes):
    present = [g for g in genes if (adata.var[GENE_SYMBOL_COL] == g).any()]
    missing = sorted(set(genes) - set(present))
    if missing:
        print(f'[subtype] marker genes not found, excluded: {missing}')
    if not present:
        return None
    return np.vstack([get_expression(adata, g) for g in present]).mean(axis=0)


def argmax_assign(adata, marker_dict, restrict_mask):
    cols = {cls: m for cls, genes in marker_dict.items()
            if (m := marker_mean(adata, genes)) is not None}
    scores = pd.DataFrame(cols, index=adata.obs.index).loc[restrict_mask]
    return scores.idxmax(axis=1).where(scores.max(axis=1) > 0)   # NaN where no marker is expressed


best_sub = argmax_assign(adata_neuronal, INTERNEURON_MARKER_GENES, is_gaba)
gaba_subtype = pd.Series('Not applicable', index=neu_obs.index, dtype=object)
gaba_subtype.loc[is_gaba] = best_sub.reindex(neu_obs.index[is_gaba]).fillna('Other GABAergic').values
neu_obs['gabaergic_subtype'] = pd.Categorical(
    gaba_subtype.values, categories=GABA_SUBTYPE_ORDER + ['Not applicable'], ordered=True)
print(neu_obs['neurotransmitter_type'].value_counts())


# ==================================================================================================
# 4. Shared aggregation helpers
# ==================================================================================================
def region_order(columns):
    lower = {str(c).lower(): c for c in columns}
    ordered = [lower[r] for r in PREFERRED_REGION_ORDER if r in lower]
    return ordered + sorted(c for c in columns if c not in ordered)


def group_region_matrices(obs, group_col, region_col='tissue'):
    """(% positive, n cells) matrices of group x region; percentages masked below MIN_CELLS."""
    g = obs.groupby([group_col, region_col], observed=True).agg(
        n_cells=(POS_COL, 'size'), n_pos=(POS_COL, 'sum'))
    n_mat = g['n_cells'].unstack()
    pct_mat = (100 * g['n_pos'] / g['n_cells']).unstack().where(n_mat >= MIN_CELLS)
    cols = region_order(n_mat.columns)
    return pct_mat.reindex(columns=cols), n_mat.reindex(columns=cols)


def long_counts(obs, group_col, region_col='tissue'):
    """Source numbers behind a panel: one row per (group x region)."""
    g = (obs.groupby([group_col, region_col], observed=True)
            .agg(n_cells=(POS_COL, 'size'), n_pos=(POS_COL, 'sum'))
            .reset_index())
    g['pct_positive'] = 100.0 * g['n_pos'] / g['n_cells']
    return g.rename(columns={group_col: 'cell_group', region_col: 'region'})


def collapse_primary_nt(label):
    """Collapse a combinatorial (co-release) atlas label such as 'DA VGLUT2' or 'GABA HDC' to one
    primary class, and flag glutamate co-release. Priority: HDC > DA > CHOL > SER > GLY > GABA >
    VGLUT*, so a modulatory transmitter wins over the fast one it is co-released with."""
    s = str(label)
    if s.strip().lower() in ('nan', 'none', '', 'unassigned'):
        return 'Unassigned', False
    toks = [t for t in re.split(r'[\s/_,;-]+', s.upper()) if t]

    def has(*keys):
        return any(any(t == k or t.startswith(k) for k in keys) for t in toks)

    has_vglut = any(t.startswith('VGLUT') for t in toks) or has('GLUT', 'SLC17')
    if has('HDC', 'HIST'):                            prim = 'Histaminergic'
    elif any(t == 'DA' for t in toks) or has('DOPA'): prim = 'Dopaminergic'
    elif has('CHOL', 'ACH'):                          prim = 'Cholinergic'
    elif has('SER', '5HT'):                           prim = 'Serotonergic'
    elif has('GLY'):                                  prim = 'Glycinergic'
    elif has('GABA'):                                 prim = 'GABAergic'
    elif has_vglut:                                   prim = 'Glutamatergic'
    else:                                             prim = 'Unassigned'
    return prim, bool(has_vglut and prim != 'Glutamatergic')


# ==================================================================================================
# 5. Main-text figure: panels a, b, c
# ==================================================================================================
MIN_LABEL_CELLS = 200 
LOW_N_CLASS = 200
NODATA_COLOR = '#f5f5f5'

EXPR_CMAP = sns.color_palette('flare', as_cmap=True)
HEAT_CMAP = LinearSegmentedColormap.from_list(
    'blues', [(0.00, '#f7fbff'), (0.50, '#4292c6'), (1.00, '#08306b')])
HEAT_CMAP.set_bad(NODATA_COLOR)

REGION_COLS = region_order(pd.unique(
    pd.concat([neu_obs['tissue'], nonneu_obs['tissue']]).astype(str)))

HM_FS_ANNOT, HM_FS_TICK, HM_FS_AXIS, HM_FS_NCELL, HM_FS_CBAR = 12, 13, 14, 12, 12
HM_CELL_W, HM_CELL_H = 0.60, 0.46
HM_MARGIN_W, HM_MARGIN_H = 3.6, 1.5
HM_REF_WIDTH = len(REGION_COLS) * HM_CELL_W + HM_MARGIN_W


def heatmap_figure(nrows, ncols):
    fig_w = ncols * HM_CELL_W + HM_MARGIN_W
    fig = plt.figure(figsize=(fig_w, nrows * HM_CELL_H + HM_MARGIN_H))
    gs = fig.add_gridspec(1, 3, width_ratios=[ncols * HM_CELL_W, 0.75, 0.3], wspace=0.03)
    return (fig, fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]),
            fig_w / HM_REF_WIDTH)


def draw_heatmap(ax, ax_side, cax, pct_mat, n_mat, row_order, row_labels=None, ylabel='',
                 xlabel='Brain region', fs_scale=1.0):
    mat = pct_mat.reindex(index=row_order)
    nmat = n_mat.reindex(index=row_order)
    vmax = max(np.ceil(np.nanmax(mat.values)), 1.0)
    sns.heatmap(mat, cmap=HEAT_CMAP, vmin=0, vmax=vmax, annot=True, fmt='.1f',
                annot_kws={'fontsize': HM_FS_ANNOT * fs_scale}, linewidths=0.5, linecolor='white',
                ax=ax, cbar_ax=cax, cbar_kws={'label': '% expressing SLC16A2'})
    cax.yaxis.label.set_size(HM_FS_TICK * fs_scale)
    cax.tick_params(labelsize=HM_FS_CBAR * fs_scale)
    ax.set_facecolor(NODATA_COLOR)
    ax.set_xticklabels([str(c).capitalize() for c in mat.columns], rotation=45, ha='right',
                       rotation_mode='anchor', fontsize=HM_FS_TICK * fs_scale)
    ax.set_yticklabels(row_labels if row_labels is not None else list(row_order), rotation=0,
                       fontsize=HM_FS_TICK * fs_scale)
    ax.set_xlabel(xlabel, fontsize=HM_FS_AXIS * fs_scale, weight='bold')
    ax.set_ylabel(ylabel, fontsize=HM_FS_AXIS * fs_scale, weight='bold')

    ax_side.set_ylim(ax.get_ylim())
    ax_side.set_xlim(0, 1)
    ax_side.axis('off')
    ax_side.text(0.0, -0.4, 'n cells', fontsize=HM_FS_NCELL * fs_scale, weight='bold',
                 va='bottom', ha='left')
    totals = nmat.sum(axis=1)
    for i, r in enumerate(row_order):
        v = totals.get(r, np.nan)
        ax_side.text(0.0, i + 0.5, f'{int(v):,}' if pd.notna(v) else '-',
                     fontsize=HM_FS_NCELL * fs_scale, va='center', ha='left')


# ---- Panel a: UMAPs of both compartments, SLC16A2+ nuclei over grey negatives ---------------------
def label_xy(xy_sub):
    """Anchor a callout at the cloud's density peak. The median lands in the empty gap whenever a
    class is split across two UMAP blobs."""
    if len(xy_sub) < 50:
        return float(np.median(xy_sub[:, 0])), float(np.median(xy_sub[:, 1]))
    hist, xe, ye = np.histogram2d(xy_sub[:, 0], xy_sub[:, 1], bins=40)
    ix, iy = np.unravel_index(np.argmax(hist), hist.shape)
    inbin = ((xy_sub[:, 0] >= xe[ix]) & (xy_sub[:, 0] < xe[ix + 1]) &
             (xy_sub[:, 1] >= ye[iy]) & (xy_sub[:, 1] < ye[iy + 1]))
    if inbin.sum() == 0:
        return float(np.median(xy_sub[:, 0])), float(np.median(xy_sub[:, 1]))
    return float(np.median(xy_sub[inbin, 0])), float(np.median(xy_sub[inbin, 1]))


def place_cluster_labels(ax, labels_xy, fontsize=8):
    """Leader-line callouts, pushed radially outward from the centroid of all anchors."""
    if not labels_xy:
        return
    xs = np.array([p[1][0] for p in labels_xy])
    ys = np.array([p[1][1] for p in labels_xy])
    cx0, cy0 = float(np.median(xs)), float(np.median(ys))
    reach = 0.16 * max(float(np.ptp(xs)), float(np.ptp(ys)), 1e-6)
    for text, (ax_, ay_) in labels_xy:
        dx, dy = ax_ - cx0, ay_ - cy0
        d = np.hypot(dx, dy) or 1.0
        tx, ty = ax_ + reach * dx / d, ay_ + reach * dy / d
        ax.annotate(text, xy=(ax_, ay_), xytext=(tx, ty), fontsize=fontsize, weight='bold',
                    ha=('left' if tx >= ax_ else 'right'), va='center', zorder=6,
                    arrowprops=dict(arrowstyle='-', lw=0.5, color='0.4', shrinkA=0, shrinkB=2),
                    path_effects=[pe.withStroke(linewidth=2.0, foreground='white')])


def panel_a():
    """Two UMAPs (~180 mm full width): every nucleus plotted, positives drawn on top and coloured by
    raw count. Point clouds are rasterised, all text stays vector. Non-neuronal callouts carry the
    atlas supercluster; the WHB cell metadata has no finer per-nucleus cell-type column."""
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.6))

    pos_counts = np.concatenate([ad.obs[EXPR_COL].to_numpy()[ad.obs[POS_COL].to_numpy()]
                                 for ad in (adata_neuronal, adata_nonneuronal)])
    vmax = max(float(np.percentile(pos_counts, 99)), EXPR_THRESHOLD + 1)
    norm = Normalize(vmin=EXPR_THRESHOLD, vmax=vmax)

    rows, mappable = [], None
    for ax, ad, comp in [(axes[0], adata_neuronal, 'Neuronal'),
                         (axes[1], adata_nonneuronal, 'Non-neuronal')]:
        xy = ad.obsm['X_umap']
        expr = ad.obs[EXPR_COL].to_numpy()
        pos = ad.obs[POS_COL].to_numpy()
        ax.scatter(xy[~pos, 0], xy[~pos, 1], s=0.4, c='#d9d9d9', linewidths=0,
                   rasterized=True, zorder=1)
        mappable = ax.scatter(xy[pos, 0], xy[pos, 1], s=1.0, c=expr[pos], cmap=EXPR_CMAP, norm=norm,
                              linewidths=0, rasterized=True, zorder=2)
        ax.set_axis_off()
        ax.set_title(comp, fontsize=11, weight='bold')

        if comp == 'Neuronal':
            ax.text(np.median(xy[:, 0]), np.median(xy[:, 1]), 'Neuron', fontsize=10, weight='bold',
                    ha='center', va='center', zorder=5,
                    path_effects=[pe.withStroke(linewidth=2.5, foreground='white')])
            rows.append(('Neuronal', 'Neuron', int(pos.size), int(pos.sum())))
        else:
            ct = ad.obs['supercluster'].astype(str)
            labels_xy = []
            for name, count in ct.value_counts().items():
                if name in ('nan', 'None', 'Unassigned', 'unknown') or count < MIN_LABEL_CELLS:
                    continue
                m = (ct == name).to_numpy()
                labels_xy.append((name, label_xy(xy[m])))
                rows.append(('Non-neuronal', name, int(m.sum()), int(pos[m].sum())))
            place_cluster_labels(ax, labels_xy)
            print(f'[panel a] {len(labels_xy)} non-neuronal cell types labelled '
                  f'(>= {MIN_LABEL_CELLS} nuclei).')

    cb = fig.colorbar(mappable, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label('SLC16A2 (raw counts)', fontsize=8)
    cb.ax.tick_params(labelsize=7)

    df = pd.DataFrame(rows, columns=['compartment', 'cell_type', 'n_cells', 'n_pos'])
    df['pct_positive'] = 100.0 * df['n_pos'] / df['n_cells']
    save_csv(df, MAIN_DIR, 'source_data_panelA_umap_celltype.csv')
    print(f'[panel a] colour scale: raw counts {EXPR_THRESHOLD}..{vmax:.1f} (99th percentile).')
    save_figure(fig, MAIN_DIR, 'panelA_umap_expression', dpi=600)


# ---- Panel b: % SLC16A2-positive, cell type x brain region ---------------------------------------
def panel_b():
    """One pooled 'Neuron' row (neuronal detail is panel c) above every non-neuronal supercluster."""
    nn = nonneu_obs[['tissue', POS_COL]].copy()
    nn['group'] = nonneu_obs['supercluster'].astype(str)
    ne = neu_obs[['tissue', POS_COL]].copy()
    ne['group'] = 'Neuron'
    obs = pd.concat([nn[['group', 'tissue', POS_COL]], ne[['group', 'tissue', POS_COL]]],
                    ignore_index=True).astype({'group': str, 'tissue': str})

    pct, n = group_region_matrices(obs, 'group')
    order = (['Neuron'] if 'Neuron' in pct.index else []) + \
            sorted([r for r in pct.index if r != 'Neuron'], key=lambda r: str(r).lower())
    save_csv(long_counts(obs, 'group'), MAIN_DIR, 'source_data_panelB_celltype_region.csv')

    pct = pct.reindex(columns=REGION_COLS)
    n = n.reindex(columns=REGION_COLS)
    fig, ax, ax_side, cax, fs_scale = heatmap_figure(len(order), len(REGION_COLS))
    draw_heatmap(ax, ax_side, cax, pct, n, order, ylabel='Cell type', fs_scale=fs_scale)
    save_figure(fig, MAIN_DIR, 'panelB_celltype_region')


# ---- Panel c: % SLC16A2-positive, primary neurotransmitter class x brain region ------------------
def panel_c():
    """The GABAergic PV/SST/VIP/Other breakdown is exported to CSV rather than drawn, restricted to
    the telencephalic regions where that interneuron taxonomy applies."""
    nt = neu_obs[['tissue', POS_COL]].copy()
    nt['nt_primary'] = [collapse_primary_nt(l)[0] for l in neu_obs['neurotransmitter_type'].astype(str)]
    nt = nt.astype({'nt_primary': str, 'tissue': str})

    map_rows = [(lab, *collapse_primary_nt(lab), int(cnt))
                for lab, cnt in neu_obs['neurotransmitter_type'].astype(str).value_counts().items()]
    map_df = pd.DataFrame(map_rows,
                          columns=['atlas_label', 'primary_class', 'glutamate_corelease', 'n_cells'])
    save_csv(map_df.sort_values(['primary_class', 'n_cells'], ascending=[True, False]),
             MAIN_DIR, 'source_data_panelC_nt_label_mapping.csv')
    unmatched = map_df.loc[(map_df['primary_class'] == 'Unassigned') &
                           (~map_df['atlas_label'].str.lower().isin(['unassigned', 'nan', 'none', '']))]
    if len(unmatched):
        print('[panel c] WARNING: atlas NT labels that fell through to Unassigned (audit these):')
        for _, r in unmatched.iterrows():
            print(f"    {r['atlas_label']!r:<28s} n={r['n_cells']:,}")
    save_csv(long_counts(nt, 'nt_primary'), MAIN_DIR, 'source_data_panelC_primary_nt_region.csv')

    pooled = (nt.groupby('nt_primary', observed=True)
                .agg(n_cells=(POS_COL, 'size'), n_pos=(POS_COL, 'sum')).reset_index())
    pooled['pct_positive'] = 100.0 * pooled['n_pos'] / pooled['n_cells']
    pooled['low_confidence'] = pooled['n_cells'] < LOW_N_CLASS
    save_csv(pooled, MAIN_DIR, 'source_data_panelC_primary_nt_pooled.csv')
    low_n = set(pooled.loc[pooled['low_confidence'], 'nt_primary'])
    if low_n:
        print(f'[panel c] pooled n < {LOW_N_CLASS}, marked "*": {sorted(low_n)}')

    pct, n = group_region_matrices(nt, 'nt_primary')
    order = [c for c in pooled.sort_values('n_cells', ascending=False)['nt_primary'] if c in pct.index]
    if DROP_UNASSIGNED_FROM_STATS:
        order = [c for c in order if c != 'Unassigned']

    gsub = neu_obs.loc[neu_obs['gabaergic_subtype'].astype(str) != 'Not applicable',
                       ['gabaergic_subtype', 'tissue', POS_COL]].astype(
                           {'gabaergic_subtype': str, 'tissue': str})
    sub_pct, sub_n = group_region_matrices(gsub, 'gabaergic_subtype')
    tel_cols = [c for c in pct.columns if c in TELENCEPHALIC_REGIONS]
    sub_order = [s for s in GABA_SUBTYPE_ORDER if s in sub_pct.index]
    export = sub_pct.reindex(index=sub_order, columns=tel_cols)
    export.index.name = 'gaba_subtype'
    export = export.reset_index()
    for c in tel_cols:
        export[f'{c} (n cells)'] = sub_n[c].reindex(sub_order).values
    save_csv(export, MAIN_DIR, 'source_data_panelC_gaba_subrows_region_pct.csv')

    pct = pct.reindex(index=order, columns=REGION_COLS)
    n = n.reindex(index=order, columns=REGION_COLS)
    fig, ax, ax_side, cax, fs_scale = heatmap_figure(len(order), len(REGION_COLS))
    draw_heatmap(ax, ax_side, cax, pct, n, order,
                 row_labels=[f'{c} *' if c in low_n else c for c in order],
                 ylabel='Neurotransmitter class', fs_scale=fs_scale)
    save_figure(fig, MAIN_DIR, 'panelC_neurotransmitter_region')


panel_a()
panel_b()
panel_c()
print(f"\nMain-text figure panels and source data in '{MAIN_DIR}/'. Positivity: raw {GENE} count > "
      f'{EXPR_THRESHOLD} (>= {EXPR_THRESHOLD + 1} counts), masked below {MIN_CELLS} nuclei per '
      '(group x region).')


# ==================================================================================================
# 6. Robustness check 1: threshold sensitivity and sequencing depth
# ==================================================================================================
# The raw>1 floor could in principle track sequencing depth, since a deeper nucleus reaches >= 2
# counts more easily. Per group we therefore report library size (the confound itself), % positive
# under the depth-sliding threshold log2(CP10k+1) > 1, and the mean Pearson residual (Lause et al.
# 2021), which puts depth into the null model instead. Concordant rankings across all three mean the
# raw>1 ranking is neither threshold-dependent nor a depth artefact.
LOGNORM_THRESHOLD = 1.0    # BrainPalmSeq normalisation: CP10k, log2(x+1), expressing := > 1
PEARSON_THETA = 100.0      # scanpy default overdispersion
PEARSON_POS_RESID = 0.0    # depth-relative positivity call: more SLC16A2 than depth predicts


def add_lognorm(adata):
    total = np.asarray(adata.X.sum(axis=1)).ravel().astype(float)   # sparse row sums, never densified
    slc = adata.obs[EXPR_COL].to_numpy().astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        cp10k = np.where(total > 0, slc / total * 1e4, 0.0)
    adata.obs['total_counts'] = total
    adata.obs['SLC16A2_lognorm'] = np.log2(1.0 + cp10k)
    adata.obs['SLC16A2_pos_lognorm'] = adata.obs['SLC16A2_lognorm'].to_numpy() > LOGNORM_THRESHOLD


def add_pearson_residuals(adatas):
    """Normalising the full matrix over millions of nuclei is impractical, so scanpy is handed a
    two-column [SLC16A2, rest-of-library] surrogate: a gene's Pearson residual depends only on the
    cell total and the gene's global rate, so SLC16A2's residual is identical to float precision.
    One null fit over both compartments keeps the residuals comparable."""
    slc = np.concatenate([a.obs[EXPR_COL].to_numpy().astype(np.float32) for a in adatas])
    tot = np.concatenate([a.obs['total_counts'].to_numpy().astype(np.float32) for a in adatas])
    surrogate = sc.AnnData(np.column_stack([slc, np.clip(tot - slc, 0.0, None)]))
    z = sc.experimental.pp.normalize_pearson_residuals(
        surrogate, theta=PEARSON_THETA, check_values=False, inplace=False)['X'][:, 0]
    print(f'[robustness] Pearson residuals (theta={PEARSON_THETA:g}) over {surrogate.n_obs:,} nuclei; '
          f'{GENE} = {slc.sum() / tot.sum():.3e} of all pooled UMIs.')
    i = 0
    for a in adatas:
        a.obs['SLC16A2_pearson'] = np.asarray(z[i:i + a.n_obs], dtype=float)
        i += a.n_obs


add_lognorm(adata_neuronal)
add_lognorm(adata_nonneuronal)
add_pearson_residuals([adata_neuronal, adata_nonneuronal])


def compare_metrics(df):
    """df columns: group, pos_raw, pos_ln, resid, total."""
    g = (df.groupby('group', observed=True)
           .agg(n_cells=('pos_raw', 'size'), median_total_counts=('total', 'median'),
                mean_total_counts=('total', 'mean'), n_pos_raw=('pos_raw', 'sum'),
                n_pos_lognorm=('pos_ln', 'sum'), mean_pearson=('resid', 'mean'),
                median_pearson=('resid', 'median'),
                n_pos_pearson=('resid', lambda s: int((s > PEARSON_POS_RESID).sum())))
           .reset_index())
    for metric in ('raw', 'lognorm', 'pearson'):
        g[f'pct_pos_{metric}'] = 100.0 * g[f'n_pos_{metric}'] / g['n_cells']
    g['pct_diff_lognorm_minus_raw'] = g['pct_pos_lognorm'] - g['pct_pos_raw']
    g['pct_diff_pearson_minus_raw'] = g['pct_pos_pearson'] - g['pct_pos_raw']
    return g.sort_values('pct_pos_raw', ascending=False).reset_index(drop=True)


def metric_frame(obs, group):
    return pd.DataFrame({'group': group,
                         'pos_raw': obs[POS_COL].to_numpy(),
                         'pos_ln': obs['SLC16A2_pos_lognorm'].to_numpy(),
                         'resid': obs['SLC16A2_pearson'].to_numpy(),
                         'total': obs['total_counts'].to_numpy()})


ct_cmp = compare_metrics(pd.concat([
    metric_frame(neu_obs, 'Neuron'),
    metric_frame(nonneu_obs, nonneu_obs['supercluster'].astype(str).to_numpy()),
], ignore_index=True))
save_csv(ct_cmp, THRESHOLD_DIR, 'robustness_celltype_raw_vs_lognorm.csv')

nt_df = metric_frame(neu_obs, [collapse_primary_nt(l)[0]
                               for l in neu_obs['neurotransmitter_type'].astype(str)])
if DROP_UNASSIGNED_FROM_STATS:
    nt_df = nt_df[nt_df['group'] != 'Unassigned']
nt_cmp = compare_metrics(nt_df)
save_csv(nt_cmp, THRESHOLD_DIR, 'robustness_neurotransmitter_raw_vs_lognorm.csv')


def spearman_raw_vs_lognorm(cmp):
    sub = cmp[cmp['n_cells'] >= MIN_CELLS]
    if len(sub) < 3:
        return np.nan, len(sub)
    return float(sub['pct_pos_raw'].corr(sub['pct_pos_lognorm'], method='spearman')), len(sub)


print(f'\n[robustness] threshold sensitivity (raw > {EXPR_THRESHOLD} vs '
      f'log2(CP10k+1) > {LOGNORM_THRESHOLD}) and depth-explicit Pearson residuals:')
for name, cmp in [('cell type', ct_cmp), ('neurotransmitter class', nt_cmp)]:
    rho, n_groups = spearman_raw_vs_lognorm(cmp)
    sub = cmp[cmp['n_cells'] >= MIN_CELLS]
    enough = len(sub) >= 3
    rho_resid = float(sub['pct_pos_raw'].corr(sub['mean_pearson'], method='spearman')) if enough else np.nan
    rho_pos = float(sub['pct_pos_raw'].corr(sub['pct_pos_pearson'], method='spearman')) if enough else np.nan
    dmin, dmax = sub['median_total_counts'].min(), sub['median_total_counts'].max()
    dfold = (dmax / dmin) if dmin and np.isfinite(dmin) and dmin > 0 else np.nan
    print(f'  by {name:<22s} rho(raw%, lognorm%) = {rho:.3f} | rho(raw%, mean resid) = {rho_resid:.3f} '
          f'| rho(raw%, Pearson% [>{PEARSON_POS_RESID:g}]) = {rho_pos:.3f} over {n_groups} groups | '
          f'median depth {dmin:,.0f}-{dmax:,.0f} UMIs ({dfold:.1f}x) | max |delta%| '
          f"lognorm={sub['pct_diff_lognorm_minus_raw'].abs().max():.1f} / "
          f"Pearson={sub['pct_diff_pearson_minus_raw'].abs().max():.1f} pts")
print('  The two checks probe the depth worry from opposite sides: lognorm% slides the threshold with\n'
      '  depth, the Pearson residual puts depth into the null model. Rho near 1 in both means the\n'
      '  raw>1 ranking is not a depth artefact.')


def threshold_figure():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 4.6), gridspec_kw={'width_ratios': [1.1, 1.0]})
    lim = 1.0
    for cmp, marker, color in [(ct_cmp, 'o', '#3182bd'), (nt_cmp, '^', '#e6550d')]:
        sub = cmp[cmp['n_cells'] >= MIN_CELLS]
        ax1.scatter(sub['pct_pos_raw'], sub['pct_pos_lognorm'], marker=marker, s=34, color=color,
                    edgecolors='white', linewidths=0.5, zorder=3)
        lim = max(lim, sub[['pct_pos_raw', 'pct_pos_lognorm']].to_numpy().max() * 1.08)
        for _, r in sub.iterrows():
            ax1.annotate(str(r['group'])[:18], (r['pct_pos_raw'], r['pct_pos_lognorm']), fontsize=6,
                         xytext=(2, 2), textcoords='offset points', zorder=4)
    ax1.plot([0, lim], [0, lim], ls='--', color='#888888', lw=1, zorder=1)
    rho_ct, _ = spearman_raw_vs_lognorm(ct_cmp)
    rho_nt, _ = spearman_raw_vs_lognorm(nt_cmp)
    ax1.set_xlim(0, lim)
    ax1.set_ylim(0, lim)
    ax1.set_xlabel(f'% positive, raw count > {EXPR_THRESHOLD}', fontsize=9, weight='bold')
    ax1.set_ylabel(f'% positive, log2(CP10k+1) > {LOGNORM_THRESHOLD}', fontsize=9, weight='bold')
    ax1.set_title(f'Threshold agreement\nSpearman rho: cell type {rho_ct:.3f}, NT class {rho_nt:.3f}',
                  fontsize=9, weight='bold', loc='left')
    ax1.legend(handles=[Line2D([0], [0], marker='o', color='#3182bd', ls='', label='cell type'),
                        Line2D([0], [0], marker='^', color='#e6550d', ls='', label='NT class'),
                        Line2D([0], [0], ls='--', color='#888888', label='y = x')],
               fontsize=7, frameon=False, loc='lower right')
    ax1.tick_params(labelsize=8)
    sns.despine(ax=ax1)

    dep = ct_cmp[ct_cmp['n_cells'] >= MIN_CELLS].sort_values('median_total_counts')
    ypos = np.arange(len(dep))
    ax2.barh(ypos, dep['median_total_counts'].values, color='#756bb1', edgecolor='white', linewidth=0.4)
    ax2.set_yticks(ypos)
    ax2.set_yticklabels([str(g)[:22] for g in dep['group']], fontsize=6.5)
    ax2.set_xlabel('Median total UMIs per nucleus', fontsize=9, weight='bold')
    ax2.set_title('Sequencing depth by cell type', fontsize=9, weight='bold', loc='left')
    ax2.tick_params(axis='x', labelsize=8)
    sns.despine(ax=ax2)

    fig.tight_layout()
    save_figure(fig, THRESHOLD_DIR, 'robustness_threshold_depth')


threshold_figure()


# ==================================================================================================
# 7. Robustness check 2: pseudobulk per fine cluster
# ==================================================================================================
# As single-nucleus detection is a sampling process, so the per-cell rate the panels report under-counts
# a low-abundance transcript. The standard remedy (Crowell et al. 2020, Squair et al. 2021) is to pool 
# alike cells, where effective depth reaches millions of UMIs. Pooling at the finest atlas level ('cluster'),
# pseudobulk CPM = 1e6 * sum(SLC16A2 counts) / sum(total UMIs), expressing := > threshold
# is compared against the per-cell raw>1 rate. Support for raw>1 means: concordant ranking, no
# cluster silent by pseudobulk yet positive by raw>1, and a visible dropout gap.
PSEUDOBULK_CPM_THRESHOLD = 1.0

alias_name = {str(k): str(v) for k, v in cluster_name_by_alias.items()}


def pseudobulk_frame(obs, celltype_group, nt_group):
    return pd.DataFrame({
        'cluster_alias': obs['cluster_alias'].astype(str).to_numpy(),
        'celltype_group': celltype_group,
        'nt_group': nt_group,
        'slc': obs[EXPR_COL].to_numpy().astype(float),
        'total': obs['total_counts'].to_numpy().astype(float),
        'pos_raw': obs[POS_COL].to_numpy(),
    })


pb_cells = pd.concat([
    pseudobulk_frame(neu_obs, 'Neuron',
                     np.array([collapse_primary_nt(l)[0]
                               for l in neu_obs['neurotransmitter_type'].astype(str)])),
    pseudobulk_frame(nonneu_obs, nonneu_obs['supercluster'].astype(str).to_numpy(), 'Non-neuronal'),
], ignore_index=True)


def pseudobulk_aggregate(cells, key):
    """Pooled CPM (expression level) beside the per-cell raw>1 rate (detection prevalence)."""
    g = (cells.groupby(key, observed=True)
              .agg(n_cells=('pos_raw', 'size'), slc_sum=('slc', 'sum'), total_sum=('total', 'sum'),
                   n_pos_raw=('pos_raw', 'sum'), median_total=('total', 'median')))
    g['pseudobulk_cpm'] = 1e6 * g['slc_sum'] / g['total_sum'].where(g['total_sum'] > 0)
    g['pseudobulk_expressed'] = g['pseudobulk_cpm'] > PSEUDOBULK_CPM_THRESHOLD
    g['pct_pos_raw'] = 100.0 * g['n_pos_raw'] / g['n_cells']
    return g.reset_index()


cl = pseudobulk_aggregate(pb_cells, 'cluster_alias')
cl = cl.merge(pb_cells.drop_duplicates('cluster_alias')[['cluster_alias', 'celltype_group', 'nt_group']],
              on='cluster_alias', how='left')
cl['cluster_name'] = cl['cluster_alias'].map(alias_name).fillna('cluster ' + cl['cluster_alias'].astype(str))
cl = cl[cl['n_cells'] >= MIN_CELLS].sort_values('pseudobulk_cpm', ascending=False).reset_index(drop=True)
save_csv(cl[['cluster_alias', 'cluster_name', 'celltype_group', 'nt_group', 'n_cells', 'median_total',
             'pseudobulk_cpm', 'pseudobulk_expressed', 'pct_pos_raw']],
         PSEUDOBULK_DIR, 'pseudobulk_by_cluster.csv')


def group_pseudobulk(cells, key, drop=None):
    g = pseudobulk_aggregate(cells, key)
    if drop:
        g = g[~g[key].isin(drop)]
    return g[g['n_cells'] >= MIN_CELLS].sort_values('pseudobulk_cpm', ascending=False).reset_index(drop=True)


ct_pb = group_pseudobulk(pb_cells, 'celltype_group')
nt_pb = group_pseudobulk(pb_cells, 'nt_group',
                         drop=['Non-neuronal', 'Unassigned'] if DROP_UNASSIGNED_FROM_STATS
                         else ['Non-neuronal'])
save_csv(ct_pb, PSEUDOBULK_DIR, 'pseudobulk_by_celltype.csv')
save_csv(nt_pb, PSEUDOBULK_DIR, 'pseudobulk_by_neurotransmitter.csv')

rho_cluster = float(cl['pseudobulk_cpm'].corr(cl['pct_pos_raw'], method='spearman'))
rho_ct = float(ct_pb['pseudobulk_cpm'].corr(ct_pb['pct_pos_raw'], method='spearman')) if len(ct_pb) >= 3 else np.nan
rho_nt = float(nt_pb['pseudobulk_cpm'].corr(nt_pb['pct_pos_raw'], method='spearman')) if len(nt_pb) >= 3 else np.nan
silent_but_detected = int((cl.loc[~cl['pseudobulk_expressed'], 'pct_pos_raw'] > 5).sum())
expressing = cl[cl['pseudobulk_expressed']]
flips = {t: int((cl['pseudobulk_cpm'] > t).sum()) for t in (0.5, 1.0, 2.0)}

print(f"\n[robustness] fine-cluster pseudobulk (n >= {MIN_CELLS} nuclei): {len(cl):,} clusters, "
      f"{int(cl['pseudobulk_expressed'].sum()):,} expressing (CPM > {PSEUDOBULK_CPM_THRESHOLD:g}).")
print(f'  (1) rank agreement pseudobulk CPM vs per-cell raw>1: Spearman rho = {rho_cluster:.3f} '
      f'over {len(cl):,} clusters.')
print(f"  (2) clusters silent by pseudobulk yet >5% positive by raw>1: {silent_but_detected} "
      f"(raw>1 is {'not ' if silent_but_detected == 0 else ''}inventing signal).")
if len(expressing):
    print(f"  (3) dropout gap across the {len(expressing):,} expressing clusters: per-cell detection "
          f"median {expressing['pct_pos_raw'].median():.1f}% "
          f"(range {expressing['pct_pos_raw'].min():.1f}-{expressing['pct_pos_raw'].max():.1f}%).")
print(f'  expressing clusters at CPM > {{0.5, 1, 2}}: {flips[0.5]:,}, {flips[1.0]:,}, {flips[2.0]:,}.')
print(f'  roll-up rank agreement: cell type rho = {rho_ct:.3f}, NT class rho = {rho_nt:.3f}.')


def pseudobulk_figure():
    floor = 1e-2       # log axis needs a positive floor for the zero-CPM clusters
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.7), gridspec_kw={'width_ratios': [1.1, 1.0]})
    cpm = cl['pseudobulk_cpm'].to_numpy(dtype=float).copy()
    cpm[~np.isfinite(cpm)] = 0.0
    ax1.scatter(cl['pct_pos_raw'], np.clip(cpm, floor, None), s=13,
                c=np.where(cl['pseudobulk_expressed'].to_numpy(), '#cb181d', '#bdbdbd'),
                edgecolors='white', linewidths=0.25, zorder=3)
    ax1.axhline(PSEUDOBULK_CPM_THRESHOLD, ls='--', color='#333333', lw=1, zorder=2)
    ax1.set_yscale('log')
    ax1.set_xlabel(f'% positive, raw count > {EXPR_THRESHOLD} (per-cell detection)', fontsize=9,
                   weight='bold')
    ax1.set_ylabel('Pseudobulk SLC16A2 (CPM, log)', fontsize=9, weight='bold')
    ax1.set_title(f'Fine clusters: expression vs detection\nSpearman rho = {rho_cluster:.3f} '
                  f'({len(cl):,} clusters)', fontsize=9, weight='bold', loc='left')
    ax1.text(0.98, 0.03, f'red = expressed (CPM > {PSEUDOBULK_CPM_THRESHOLD:g})',
             transform=ax1.transAxes, fontsize=6.5, ha='right', va='bottom', color='#cb181d')
    ax1.tick_params(labelsize=8)
    sns.despine(ax=ax1)

    for cmp, keycol, marker, color, label in [(ct_pb, 'celltype_group', 'o', '#3182bd', 'cell type'),
                                              (nt_pb, 'nt_group', '^', '#e6550d', 'NT class')]:
        ax2.scatter(cmp['pct_pos_raw'], np.clip(cmp['pseudobulk_cpm'], floor, None), marker=marker,
                    s=34, color=color, edgecolors='white', linewidths=0.5, label=label, zorder=3)
        for _, r in cmp.iterrows():
            ax2.annotate(str(r[keycol])[:16], (r['pct_pos_raw'], max(r['pseudobulk_cpm'], floor)),
                         fontsize=6, xytext=(2, 2), textcoords='offset points', zorder=4)
    ax2.set_yscale('log')
    ax2.set_xlabel(f'% positive, raw count > {EXPR_THRESHOLD}', fontsize=9, weight='bold')
    ax2.set_ylabel('Group pseudobulk SLC16A2 (CPM, log)', fontsize=9, weight='bold')
    ax2.set_title(f'Panel groupings\nrho: cell type {rho_ct:.3f}, NT class {rho_nt:.3f}',
                  fontsize=9, weight='bold', loc='left')
    ax2.legend(fontsize=7, frameon=False, loc='lower right')
    ax2.tick_params(labelsize=8)
    sns.despine(ax=ax2)

    fig.tight_layout()
    save_figure(fig, PSEUDOBULK_DIR, 'pseudobulk_vs_raw')


pseudobulk_figure()

print(f'\nDone. Data: ABC WHB atlas, manifest {abc.current_manifest}.')

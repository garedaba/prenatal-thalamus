"""
04 | Expression of age-related thalamic genes in single cells

Total expression of the genes in the transcriptomic age model is calculated
for each cell in a single-cell atlas of the second-trimester human thalamus
(Kim et al. 2023; n = 137,007 cells, 20 clusters) and displayed on the
precomputed UMAP embedding and by cell cluster.

Inputs : outputs from 01; config.KIM_META, KIM_UMAP, KIM_FEATURES, KIM_MATRIX
Outputs: outputs/age_model/kim2023_cell_scores_{tag}.csv
Figures: Fig. 2G

Results: "Gene expression indexes relative maturation of thalamic nuclei at mid-gestation"
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.io import mmread

from prenatal_thalamus import config
from prenatal_thalamus.utils import report, require, savefig

TAG = config.tag()

CLUSTER_COLOURS = {'EN1': '#1d84c2', 'EN2': '#d87129', 'IN1': '#3196d3', 'IN2': '#1b8c6a',
                   'IN3': '#dac729', 'IN4': '#3b779f', 'IN5': '#a34f20', 'IN6': '#b54b85',
                   'PC': '#f4ee73', 'EC': '#f6f08e', 'MG': '#ffd06a', 'FB': '#7f7b7c',
                   'DIV': '#14a479', 'EP': '#d58ab2', 'AC1': '#eaa042', 'AC2': '#6bb9e8',
                   'GL': '#b58531', 'IPC': '#5a5657', 'OL': '#95cbea', 'OPC': '#75c7ab'}


def load_kim():
    """Cell metadata with UMAP coordinates, gene symbols and the genes x cells count matrix."""
    require(config.KIM_META, config.KIM_UMAP, config.KIM_FEATURES, config.KIM_MATRIX)
    meta = pd.read_csv(config.KIM_META, delimiter='\t')
    umap = pd.read_csv(config.KIM_UMAP, delimiter='\t', header=None)
    meta['umap_x'] = meta['Cell'].map(dict(zip(umap[0], umap[1])))
    meta['umap_y'] = meta['Cell'].map(dict(zip(umap[0], umap[2])))
    features = pd.read_csv(config.KIM_FEATURES, delimiter='\t', header=None)[1]
    matrix = mmread(config.KIM_MATRIX).tocsr()
    return meta, features, matrix


def total_expression(matrix, features, genes):
    """Sum of expression over a gene set for each cell."""
    rows = np.where(features.isin(genes))[0]
    return np.asarray(matrix[rows].sum(axis=0)).ravel(), len(rows)


def main():
    print('\n=== 04 | Age-model genes in single cells (Kim et al. 2023) ===')
    genes = np.loadtxt(config.AGE_MODEL_DIR / f'thalamic_genes_{TAG}.txt', dtype=str)
    meta, features, matrix = load_kim()
    meta['score'], n_found = total_expression(matrix, features, genes)
    report('cells', f'{len(meta):,}', paper='137,007')
    report('age-model genes found', f'{n_found} / {len(genes)}')
    meta[['Cell', 'clusters', 'score']].to_csv(config.AGE_MODEL_DIR / f'kim2023_cell_scores_{TAG}.csv', index=False)


    print('\nExpression of age-model genes by cell cluster (Fig. 2G)')

    # Fig. 2G: UMAP coloured by total expression, inset by cluster
    fig, ax = plt.subplots(figsize=(8, 10))
    sns.scatterplot(x='umap_x', y='umap_y', hue='score', palette='Blues', data=meta, s=10, alpha=0.8,
                    linewidth=0, legend=False, ax=ax)
    inset = fig.add_axes([0.73, 0.56, 0.3, 0.3])
    sns.scatterplot(x='umap_x', y='umap_y', hue='clusters', hue_order=list(CLUSTER_COLOURS), palette=CLUSTER_COLOURS,
                    data=meta, s=5, alpha=0.5, linewidth=0, legend=False, ax=inset)
    for a in (ax, inset):
        a.set_axis_off()
    savefig(fig, 'kim2023-umap', panel='fig2g')

    # Fig. 2G (right): distribution of total expression within each cluster
    g = sns.FacetGrid(meta, row='clusters', hue='clusters', row_order=list(CLUSTER_COLOURS),
                      hue_order=list(CLUSTER_COLOURS), palette=CLUSTER_COLOURS, aspect=5, height=0.5)
    g.map(sns.kdeplot, 'score', bw_adjust=0.5, fill=True, alpha=1, linewidth=1.5, clip_on=False)
    for ax_, cluster in zip(g.axes.flat, CLUSTER_COLOURS):
        ax_.text(1, 0, cluster, color=CLUSTER_COLOURS[cluster], ha='left', va='center', transform=ax_.transAxes)
        ax_.set_facecolor('none')
    g.figure.subplots_adjust(hspace=-0.55)
    g.set_titles('')
    g.set(yticks=[], ylabel='')
    g.despine(bottom=True, left=True)
    g.set_xlabels('expression')
    savefig(g.figure, 'kim2023-cluster-distributions', panel='fig2g')


if __name__ == '__main__':
    main()

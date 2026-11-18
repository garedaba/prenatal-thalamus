"""
06 | Model fit to empirical node degree

1. Expected node degree under the best performing model (step 05) is compared
   with empirical degree, and tested against n = 1000 surrogate timing maps
   with matched spatial autocorrelation (BrainSMASH).
2. The model is thresholded to the empirical density and compared at the edge
   level: edge length distributions and overlap (Jaccard index) of short,
   mid-range and long edges, relative to n = 500 degree-preserving random
   networks.

Inputs : as step 05; config.GEODESIC_DISTANCE; outputs/generative_model/best_parameters.csv
Outputs: outputs/generative_model/
           best_model_degree.csv           - empirical and simulated degree per node
           spatial_nulls_timing.csv        - model fit for each surrogate map
           edge_overlap.csv                - Jaccard index by edge length (model and random nulls)
Figures: Fig. 4C, 4D; fig. S7

Results: "A generative model of thalamocortical constraints on the developing human connectome"
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from bct import randmio_und_connected, threshold_proportional
from matplotlib_venn import venn2
from scipy.stats import gaussian_kde, pearsonr
from tqdm import tqdm

from prenatal_thalamus import config
from prenatal_thalamus.network import (edges_by_length, generate_network, jaccard_index, load_network_data,
                                       spatial_null_test, spatial_surrogates)
from prenatal_thalamus.surface import plot_parcels
from prenatal_thalamus.utils import outdir, permutation_p, report, savefig

OUT = outdir(config.MODEL_DIR)


def degree_fit(data, eta, beta):
    """Fig. 4C, 4D: simulated vs empirical degree, tested against spatial nulls."""
    print(f'\nSimulated vs empirical degree, tested against {config.N_SPATIAL_NULLS} spatial nulls (Fig. 4D)')
    D, T, degree = data['D'], data['T_timing'], data['degree']
    k, _ = generate_network(D, T, eta, beta)
    pd.DataFrame({'node': np.arange(len(k)), 'empirical': degree, 'simulated': k}).to_csv(
        OUT / 'best_model_degree.csv', index=False)

    surrogates = spatial_surrogates(T)
    r, r_null, p = spatial_null_test(D, T, degree, eta, beta, surrogates)
    pd.DataFrame({'r': r_null}).to_csv(OUT / 'spatial_nulls_timing.csv', index=False)
    ci = pearsonr(k, degree).confidence_interval()   # 95% CI, Fisher z (assumes independent nodes)
    report('simulated vs empirical degree', f'r = {r:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}], p_spatial = {p:.3f}',
           paper='r = 0.41, p = 0.003')

    # Fig. 4D
    fig, ax = plt.subplots(1, 2, figsize=(9, 4), width_ratios=(0.35, 0.65))
    sns.regplot(x=k, y=degree, color='grey', scatter_kws={'facecolor': '#c24245', 'edgecolor': 'grey', 's': 60}, ax=ax[0])
    ax[0].set_xlabel('simulated degree')
    ax[0].set_ylabel('empirical degree')
    ax[0].set_title('model fit', loc='left')
    sns.histplot(r_null, color='#c24245', ax=ax[1])
    ax[1].axvline(r, color='black')
    ax[1].text(r, ax[1].get_ylim()[1] * 0.9, f'p = {p:.3f} ', ha='right')
    ax[1].set_xlabel('r')
    ax[1].set_title('null distribution', loc='left')
    fig.tight_layout()
    savefig(fig, 'model-fit', panel='fig4d')
    return k


def edge_overlap(data, eta, beta):
    """Fig. S7: edge length distributions and edge overlap with the empirical network."""
    print(f'\nEdge overlap with empirical network vs {config.N_RANDOM_NETWORKS} degree-preserving random networks (fig. S7)')
    D = data['D']
    _, P = generate_network(D, data['T_timing'], eta, beta)
    empirical = (data['connectivity'] > 0).astype(int)
    synthetic = (threshold_proportional(P, config.NETWORK_DENSITY) > 0).astype(int)

    # edge length bins from the empirical network (quartiles)
    bounds = np.percentile(D[empirical == 1], [25, 75])
    emp_bins = edges_by_length(empirical, D, bounds)
    syn_bins = edges_by_length(synthetic, D, bounds)

    # degree-preserving random networks
    rng_state = np.random.get_state()
    np.random.seed(config.RANDOM_SEED)
    random_nets = [randmio_und_connected(empirical, 5)[0] for _ in tqdm(range(config.N_RANDOM_NETWORKS), desc='  random networks')]
    np.random.set_state(rng_state)

    rows, null = [], {}
    for length in emp_bins:
        observed = jaccard_index(syn_bins[length], emp_bins[length])
        null[length] = [jaccard_index(edges_by_length(R, D, bounds)[length], emp_bins[length]) for R in random_nets]
        p = permutation_p(observed, null[length])
        rows.append([length, observed, np.mean(null[length]), p])
        report(f'Jaccard, {length} edges', f'{observed:.2f} (random: {np.mean(null[length]):.2f}), p = {p:.3f}',
               paper='p = 0.002')
    pd.DataFrame(rows, columns=['edge_length', 'jaccard', 'jaccard_random_mean', 'p']).to_csv(OUT / 'edge_overlap.csv', index=False)
    report('Jaccard, all edges', f'{jaccard_index(synthetic, empirical):.2f}')

    # edge length CDFs
    x = np.linspace(0, 1, 100)

    def cdf(A):
        density = gaussian_kde(D[A > 0])(x)
        return np.cumsum(density / density.sum())

    fig, ax = plt.subplot_mosaic('AABB\nAACC\nAADD', figsize=(10, 5))
    for R in random_nets:
        ax['A'].plot(x, cdf(R), color='grey', alpha=0.1, zorder=0)
    ax['A'].plot(x, cdf(empirical), color='black', lw=2, label='empirical')
    ax['A'].plot(x, cdf(synthetic), color='#851c29', lw=2, label='model')
    ax['A'].set_xlabel(r'$D_{ij}$')
    ax['A'].set_ylabel('CDF')
    ax['A'].legend()
    ax['A'].set_title('edge length distribution', loc='left')

    for key, length, colour in zip('BCD', ['long', 'mid-range', 'short'], ['#c24245', '#851c29', '#441016']):
        ax[key].plot(x, gaussian_kde(null[length])(x), color=colour)
        ax[key].axvline(jaccard_index(syn_bins[length], emp_bins[length]), color=colour, lw=2)
        ax[key].set_xlim(0, 0.8)
        ax[key].set_yticks([])
        ax[key].set_title(f'{length} edges', loc='left')
    ax['D'].set_xlabel('Jaccard index')
    fig.tight_layout()
    savefig(fig, 'edge-length-overlap', panel='figS7')

    # overlap of edges in each length bin
    fig, ax = plt.subplots(3, 1, figsize=(3, 9))
    for a, length in zip(ax, ['short', 'mid-range', 'long']):
        e, s = emp_bins[length], syn_bins[length]
        venn2(subsets=((e * (1 - s)).sum(), (s * (1 - e)).sum(), (e * s).sum()),   # counts over full symmetric matrix
              set_labels=['empirical', 'model'], set_colors=['black', '#c24245'], alpha=0.8, ax=a)
        a.set_title(f'{length} edges', loc='left')
    savefig(fig, 'edge-overlap-venn', panel='figS7')


def main():
    print('\n=== 06 | Model fit to the empirical network ===')
    data = load_network_data()
    best = pd.read_csv(config.MODEL_DIR / 'best_parameters.csv').set_index('density').loc[config.NETWORK_DENSITY]
    report('best parameters (step 05)', f'η = {best["eta"]:.2f}, β = {best["beta"]:.2f}')

    k = degree_fit(data, best['eta'], best['beta'])
    edge_overlap(data, best['eta'], best['beta'])

    # Fig. 4C: surface maps (VTK rendering)
    print()
    plot_parcels(data['degree'], data['parc'], 'empirical-degree', panel='fig4c')
    plot_parcels(k, data['parc'], 'simulated-degree', panel='fig4c')


if __name__ == '__main__':
    main()

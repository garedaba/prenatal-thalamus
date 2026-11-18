"""
05 | Probabilistic generative model of cortical node degree

Connection probability between cortical nodes i and j depends on their
distance (D_ij) and the difference in estimated timing of thalamic
innervation (T_i):

    P(i,j) ∝ exp(-η·D_ij - β·|T_i - T_j|)

Expected node degree is scaled to the empirical network density (25%).
Model fit is the Pearson correlation between expected and empirical degree,
evaluated over a 50 x 50 grid of η and β.

Inputs : config.CORTICAL_CONNECTIVITY, EUCLIDEAN_DISTANCE, CORTICAL_PARCELLATION,
         CORTICAL_INNERVATION_TIMING, DHCP_SURFACE_LEFT
Outputs: outputs/generative_model/
           sweep_timing_density{d}.csv     - model fit over the parameter grid
           best_parameters.csv             - best η, β at each density
Figures: Fig. 4A, 4B; fig. S8

Results: "A generative model of thalamocortical constraints on the developing human connectome"
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

from prenatal_thalamus import config
from prenatal_thalamus.network import best_parameters, load_network_data, parameter_sweep, plot_sweep
from prenatal_thalamus.surface import plot_parcels
from prenatal_thalamus.utils import outdir, report, savefig

OUT = outdir(config.MODEL_DIR)


def innervation_timing(data):
    """Fig. 4A: timing of thalamic innervation vs distance to primary sensory cortex."""
    print('\nTiming of thalamic innervation (Fig. 4A)')
    timing = 1 - data['T_timing']   # high = earlier innervation (more mature nuclei)

    res = pearsonr(data['T_sa'], timing)
    ci = res.confidence_interval()   # 95% CI, Fisher z-transform
    report('timing vs distance to primary cortex',
           f'r = {res.statistic:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}], p = {res.pvalue:.3g}',
           paper='r = -0.70, p < 0.001')
    res = pearsonr(data['degree'], timing)
    ci = res.confidence_interval()
    report('timing vs node degree',
           f'r = {res.statistic:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}], p = {res.pvalue:.2f}',
           paper='r = 0.11, p = 0.28')

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.regplot(x=data['T_sa'], y=timing, color='grey',
                scatter_kws={'s': 80, 'edgecolor': 'grey', 'facecolor': '#ffd788'}, ax=ax)
    ax.set_xlabel('distance to primary')
    ax.set_ylabel('thalamocortical timing')
    savefig(fig, 'timing-vs-distance', panel='fig4a')

    # surface maps (VTK rendering)
    plot_parcels(timing, data['parc'], 'thalamocortical-timing', panel='fig4a', cmap='RdBu')
    plot_parcels(data['T_sa'], data['parc'], 'distance-to-primary', panel='fig4a', cmap='Spectral_r')


def main():
    print('\n=== 05 | Generative model: P(i,j) ∝ exp(-η·D_ij - β·|T_i - T_j|) ===')

    data = load_network_data()
    print('\nParameter sweep: fit to empirical node degree (Fig. 4B, fig. S8)')
    best = []
    for density in [config.NETWORK_DENSITY] + config.DENSITIES_SUPP:
        d = data if density == config.NETWORK_DENSITY else load_network_data(density=density)
        sweep = parameter_sweep(d['D'], d['T_timing'], d['degree'], density=density)
        sweep.to_csv(OUT / f'sweep_timing_density{density}.csv', index=False)
        b = best_parameters(sweep)
        best.append([density, b['eta'], b['beta'], b['r']])

        main_analysis = density == config.NETWORK_DENSITY
        report(f'best model, density {density:.0%}', f'η = {b["eta"]:.2f}, β = {b["beta"]:.2f}, r = {b["r"]:.2f}',
               paper='η = 1.43, β = 8.17, r = 0.41' if main_analysis else None)

        fig, ax = plt.subplots(figsize=(6, 5))
        plot_sweep(sweep, ax)
        ax.set_title(f'network density = {density:.0%}', loc='left')
        savefig(fig, f'parameter-sweep-density{density}', panel='fig4b' if main_analysis else 'figS8')

    pd.DataFrame(best, columns=['density', 'eta', 'beta', 'r']).to_csv(OUT / 'best_parameters.csv', index=False)

    innervation_timing(data)


if __name__ == '__main__':
    main()

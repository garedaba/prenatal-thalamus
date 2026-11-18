"""
07 | Alternative models: distance only and sensory-association axis

1. Distance only (exponential distance rule): β = 0, sweep over η.
2. Sensory-association (S-A) axis: node attribute T replaced with the minimum
   distance of each node to primary sensory cortex (S1, A1, V1). Full
   parameter sweep and spatial null test, as for the timing model.

Inputs : as step 05; config.GEODESIC_DISTANCE
Outputs: outputs/generative_model/
           sweep_distance_only.csv
           sweep_sa_axis.csv
           spatial_nulls_sa_axis.csv
           model_comparison.csv
Figures: fig. S9

Results: "A generative model of thalamocortical constraints on the developing human connectome"
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr

from prenatal_thalamus import config
from prenatal_thalamus.network import (best_parameters, generate_network, load_network_data, parameter_sweep,
                                       plot_sweep, spatial_null_test, spatial_surrogates)
from prenatal_thalamus.surface import plot_parcels
from prenatal_thalamus.utils import outdir, report, savefig

OUT = outdir(config.MODEL_DIR)


def fit_ci(D, T, degree, eta, beta):
    """95% CI (Fisher z) for the correlation between expected and empirical degree."""
    k, _ = generate_network(D, T, eta, beta)
    return pearsonr(k, degree).confidence_interval()


def main():
    print('\n=== 07 | Alternative models: distance only and S-A axis ===')
    data = load_network_data()
    D, degree = data['D'], data['degree']
    comparison = []

    # thalamocortical timing model (step 05)
    timing = pd.read_csv(config.MODEL_DIR / 'best_parameters.csv').set_index('density').loc[config.NETWORK_DENSITY]
    ci = fit_ci(D, data['T_timing'], degree, timing['eta'], timing['beta'])
    comparison.append(['thalamocortical timing', timing['eta'], timing['beta'], timing['r'], ci.low, ci.high])

    # 1. distance only ---------------------------------------------------------
    print('\nDistance only (β = 0)')
    sweep = parameter_sweep(D, data['T_sa'], degree, betas=[0.0])   # T has no effect when β = 0
    sweep.to_csv(OUT / 'sweep_distance_only.csv', index=False)
    best = best_parameters(sweep)
    ci = fit_ci(D, data['T_sa'], degree, best['eta'], 0.0)
    report('best distance-only model', f'η = {best["eta"]:.2f}, r = {best["r"]:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}]',
           paper='r = 0.21')
    comparison.append(['distance only', best['eta'], 0.0, best['r'], ci.low, ci.high])

    # 2. sensory-association axis ---------------------------------------------
    print('\nS-A axis: T = distance to primary sensory cortex (fig. S9)')
    sweep = parameter_sweep(D, data['T_sa'], degree)
    sweep.to_csv(OUT / 'sweep_sa_axis.csv', index=False)
    best = best_parameters(sweep)
    r, r_null, p = spatial_null_test(D, data['T_sa'], degree, best['eta'], best['beta'],
                                     spatial_surrogates(data['T_sa']))
    pd.DataFrame({'r': r_null}).to_csv(OUT / 'spatial_nulls_sa_axis.csv', index=False)
    ci = fit_ci(D, data['T_sa'], degree, best['eta'], best['beta'])
    report('best S-A axis model',
           f'η = {best["eta"]:.2f}, β = {best["beta"]:.2f}, r = {r:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}], p_spatial = {p:.3f}',
           paper='r = 0.27, p = 0.056')
    comparison.append(['S-A axis', best['eta'], best['beta'], r, ci.low, ci.high])

    # comparison ----------------------------------------------------------------
    print('\nModel comparison')
    comparison = pd.DataFrame(comparison, columns=['model', 'eta', 'beta', 'r', 'ci95_low', 'ci95_high'])
    comparison.to_csv(OUT / 'model_comparison.csv', index=False)
    print(comparison.round(2).to_string(index=False))
    print()

    # fig. S9
    k, _ = generate_network(D, data['T_sa'], best['eta'], best['beta'])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), width_ratios=(0.55, 0.45))
    plot_sweep(sweep, ax[0])
    sns.regplot(x=k, y=degree, color='grey', scatter_kws={'facecolor': '#c24245', 'edgecolor': 'grey', 's': 60}, ax=ax[1])
    ax[1].set_xlabel('simulated degree')
    ax[1].set_ylabel('empirical degree')
    ax[1].set_title(f'r = {r:.2f}, p = {p:.3f}', loc='left')
    fig.tight_layout()
    savefig(fig, 'sa-axis-model', panel='figS9')

    # surface maps (VTK rendering)
    plot_parcels(degree, data['parc'], 'sa-axis-empirical-degree', panel='figS9')
    plot_parcels(k, data['parc'], 'sa-axis-simulated-degree', panel='figS9')


if __name__ == '__main__':
    main()

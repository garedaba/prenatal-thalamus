"""
08 | Synthetic deafferentation of primary visual cortex

The thalamic timing constraint is removed for nodes in primary visual cortex
(β = 0 for all connections to/from V1 nodes), returning these nodes to an
exponential distance rule. Expected connectivity of V1 nodes to the rest of
the cortex is compared before and after deafferentation.

Inputs : as step 05; outputs/generative_model/best_parameters.csv
Outputs: outputs/generative_model/ablation_v1.csv
Figures: fig. S10b

Results: "Removing thalamocortical constraints mirrors findings from experimental ablation studies"
"""

import numpy as np
import pandas as pd

from prenatal_thalamus import config
from prenatal_thalamus.network import generate_network, load_network_data
from prenatal_thalamus.surface import plot_parcels
from prenatal_thalamus.utils import outdir, report

OUT = outdir(config.MODEL_DIR)


def main():
    print('\n=== 08 | Synthetic deafferentation of primary visual cortex (fig. S10) ===')
    data = load_network_data()
    best = pd.read_csv(config.MODEL_DIR / 'best_parameters.csv').set_index('density').loc[config.NETWORK_DENSITY]
    v1 = config.V1_NODES

    report('best parameters (step 05)', f'η = {best["eta"]:.2f}, β = {best["beta"]:.2f}')
    report('V1 nodes (0-indexed)', ', '.join(map(str, v1)))

    k_before, P_before = generate_network(data['D'], data['T_timing'], best['eta'], best['beta'])
    k_after, P_after = generate_network(data['D'], data['T_timing'], best['eta'], best['beta'], ablate=v1)

    # mean expected connectivity from V1 nodes to every cortical node
    before = P_before[v1].mean(0)
    after = P_after[v1].mean(0)
    diff = after - before
    pd.DataFrame({'node': np.arange(len(before)), 'v1_connectivity_before': before,
                  'v1_connectivity_after': after, 'degree_before': k_before, 'degree_after': k_after}
                 ).to_csv(OUT / 'ablation_v1.csv', index=False)

    print()

    mask = np.isin(np.arange(config.N_NODES), v1).astype(float)
    plot_parcels(mask, data['parc'], 'deafferented-cortex', panel='figS10b', cmap='Reds')
    plot_parcels(before, data['parc'], 'v1-connectivity-before', panel='figS10b')
    plot_parcels(after, data['parc'], 'v1-connectivity-after', panel='figS10b')
    plot_parcels(diff, data['parc'], 'v1-connectivity-difference', panel='figS10b')


if __name__ == '__main__':
    main()

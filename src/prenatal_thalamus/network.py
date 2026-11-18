"""
Cortical network data and the probabilistic generative model (steps 05-08).

Model:  P(i,j) ∝ exp(-η·D_ij - β·|T_i - T_j|)

D_ij : Euclidean distance between cortical nodes (scaled to [0, 1])
T_i  : node attribute (estimated timing of thalamic innervation, or distance
       to primary sensory cortex)
η    : distance penalty;  β : attribute (dis)similarity penalty
η > 0, β = 0  -> exponential distance rule (EDR)
η = 0, β > 0  -> purely homophilic attachment
"""

import shutil

import numpy as np
import pandas as pd
from bct import degrees_und, threshold_proportional
from brainsmash.mapgen.base import Base
from scipy.spatial.distance import jaccard
from scipy.stats import pearsonr

from . import config
from .utils import rank_normalise, require


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
def copy_matlab_export(path, name=None):
    """Copy an input from the MATLAB export folder (matlab_analysis/Data4GenMdl) if it is not in DATA_DIR."""
    src = config.MATLAB_EXPORT_DIR / (name or path.name)
    if not path.exists() and src.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src, path)
        print(f'  copied {src} -> {path}')


def load_network_data(density=config.NETWORK_DENSITY, load_timing=True):
    """
    Load the empirical cortical network and node attributes.

    Returns a dict with:
      connectivity : thresholded within-hemisphere network (N x N), averaged
                     across hemispheres
      degree       : empirical binary node degree (N,)
      D            : Euclidean distance, scaled to [0, 1] (N x N)
      T_sa         : distance to primary sensory cortex, rank-normalised (N,)
      T_timing     : estimated timing of thalamic innervation, rank-normalised (N,)
      parc         : vertex-wise cortical parcellation (left hemisphere)
    """
    N = config.N_NODES
    if load_timing:
        copy_matlab_export(config.CORTICAL_INNERVATION_TIMING)
    paths = [config.CORTICAL_CONNECTIVITY, config.EUCLIDEAN_DISTANCE, config.CORTICAL_PARCELLATION]
    if load_timing:
        paths.append(config.CORTICAL_INNERVATION_TIMING)
    require(*paths)

    # within-hemisphere connections, averaged across hemispheres
    connectivity = np.loadtxt(config.CORTICAL_CONNECTIVITY, delimiter=',')
    connectivity = (connectivity[:N, :N] + connectivity[N:, N:]) / 2
    connectivity = threshold_proportional(connectivity, density)
    degree = degrees_und(connectivity)

    # Euclidean distance (left hemisphere), scaled to [0, 1]
    D = np.loadtxt(config.EUCLIDEAN_DISTANCE, delimiter=',')[:N, :N]
    D = D / D.max()

    # sensory-association axis: minimum distance to any primary sensory node
    T_sa = rank_normalise(D[config.PRIMARY_SENSORY_NODES, :].min(0))

    data = dict(connectivity=connectivity, degree=degree, D=D, T_sa=T_sa,
                parc=np.loadtxt(config.CORTICAL_PARCELLATION, delimiter=','))

    if load_timing:
        timing = np.loadtxt(config.CORTICAL_INNERVATION_TIMING, delimiter=',')
        data['T_timing'] = rank_normalise(-timing)

    return data


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
def generate_network(D, T, eta, beta, density=config.NETWORK_DENSITY, ablate=None):
    """
    Expected connection weights under the generative model.

    Parameters
    ----------
    D : (N, N) distance matrix
    T : (N,) node attribute
    eta, beta : distance and attribute penalties
    density : network density used to scale weights to the expected number of edges
    ablate : list of node indices for which the attribute term is removed
             (β set to 0 for all connections to/from these nodes)

    Returns
    -------
    k : (N,) expected node degree
    P : (N, N) expected connection weights, scaled so that
        sum_{i<j} P_ij = number of edges at the given density
    """
    N = len(T)
    T_diff = np.abs(T[:, None] - T[None, :])
    logP = -eta * D - beta * T_diff

    if ablate is not None:
        for i in ablate:
            logP[i, :] = -eta * D[i, :]
            logP[:, i] = -eta * D[:, i]

    np.fill_diagonal(logP, -np.inf)
    P = np.exp(logP)
    P = (P + P.T) / 2

    # scale to expected number of edges: E[k_i] = m * s_i / sum_{a<b} w_ab
    m = round((N * (N - 1) // 2) * density)
    P = P * m / np.triu(P, 1).sum()

    return P.sum(0), P


def model_fit(D, T, degree, eta, beta, density=config.NETWORK_DENSITY):
    """Pearson correlation between expected and empirical node degree."""
    k, _ = generate_network(D, T, eta, beta, density)
    if np.allclose(k, k[0]):
        return np.nan   # all edges equally likely
    return pearsonr(k, degree)[0]


def parameter_sweep(D, T, degree, etas=config.ETA_GRID, betas=config.BETA_GRID,
                    density=config.NETWORK_DENSITY):
    """Evaluate model fit over a grid of η and β values."""
    results = [[eta, beta, model_fit(D, T, degree, eta, beta, density)]
               for eta in etas for beta in betas]
    return pd.DataFrame(results, columns=['eta', 'beta', 'r'])


def best_parameters(sweep):
    """Row of a parameter sweep with the highest model fit."""
    return sweep.loc[sweep['r'].idxmax()]


def plot_sweep(sweep, ax):
    """Heatmap of model fit (r) over η and β, with the best parameters marked."""
    import seaborn as sns
    grid = sweep.pivot(index='beta', columns='eta', values='r')
    sns.heatmap(grid, cmap='OrRd_r', cbar_kws={'label': 'r'}, ax=ax)
    best = best_parameters(sweep)
    ax.scatter(grid.columns.get_loc(best['eta']) + 0.5, grid.index.get_loc(best['beta']) + 0.5,
               s=100, edgecolor='black', facecolor='#f54927', zorder=10)
    ticks = np.linspace(0, len(grid) - 1, 5, dtype=int)
    ax.set_xticks(ticks + 0.5, np.round(grid.columns[ticks], 2), rotation=0)
    ax.set_yticks(ticks + 0.5, np.round(grid.index[ticks], 2))
    ax.invert_yaxis()
    ax.set_xlabel(r'$\eta$')
    ax.set_ylabel(r'$\beta$')


# ---------------------------------------------------------------------------
# Null models
# ---------------------------------------------------------------------------
def spatial_surrogates(T, n=config.N_SPATIAL_NULLS, seed=config.RANDOM_SEED):
    """Surrogate maps with matched spatial autocorrelation (BrainSMASH)."""
    copy_matlab_export(config.GEODESIC_DISTANCE, name='Left_geo_distance.csv')
    require(config.GEODESIC_DISTANCE)
    with open(config.GEODESIC_DISTANCE) as f:      # whitespace- or comma-delimited
        delimiter = ',' if ',' in f.readline() else None
    geodesic = np.loadtxt(config.GEODESIC_DISTANCE, delimiter=delimiter)
    base = Base(x=T, D=geodesic, resample=True, seed=seed)
    return base(n=n)


def spatial_null_test(D, T, degree, eta, beta, surrogates, density=config.NETWORK_DENSITY):
    """
    Compare model fit with that obtained from surrogate attribute maps at the
    same parameters. Returns observed r, null r values and p-value.
    """
    r_obs = model_fit(D, T, degree, eta, beta, density)
    r_null = np.array([model_fit(D, s, degree, eta, beta, density) for s in surrogates])
    p = (np.sum(r_null > r_obs) + 1) / (len(r_null) + 1)
    return r_obs, r_null, p


# ---------------------------------------------------------------------------
# Edge-level comparison (fig. S7)
# ---------------------------------------------------------------------------
def jaccard_index(a, b):
    """Jaccard similarity between two binary adjacency matrices."""
    return 1 - jaccard((a != 0).ravel(), (b != 0).ravel())


def edges_by_length(A, D, bounds):
    """Split a binary network into short, mid-range and long edges."""
    lo, hi = bounds
    return {'short': A * (D < lo),
            'mid-range': A * ((D > lo) & (D < hi)),
            'long': A * (D > hi)}

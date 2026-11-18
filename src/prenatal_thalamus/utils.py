"""General helpers: input checks, reporting and figure output."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from . import config


def require(*paths):
    """Raise an informative error if any required input is missing."""
    missing = [str(p) for p in paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            'Required input(s) not found:\n  ' + '\n  '.join(missing) +
            '\nSee README.md (Data) and src/prenatal_thalamus/config.py for sources.')


def outdir(path):
    """Create (if needed) and return an output directory."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def report(label, value, paper=None):
    """Print a result alongside the value reported in the paper."""
    line = f'  {label:<42} {value}'
    if paper is not None:
        line += f'   [paper: {paper}]'
    print(line)


def savefig(fig, name, panel=None):
    """Save a figure panel to FIGURE_DIR as PNG (300 dpi)."""
    outdir(config.FIGURE_DIR)
    fname = f'{panel}_{name}.png' if panel else f'{name}.png'
    fig.savefig(config.FIGURE_DIR / fname, dpi=300, facecolor='white', bbox_inches='tight')
    plt.close(fig)
    print(f'  -> figures/{fname}')


def rank_normalise(x):
    """Rank-transform to (0, 1]."""
    x = np.asarray(x)
    return np.argsort(np.argsort(x)) / len(x) + 1e-8


def permutation_p(observed, null):
    """One-sided permutation p-value: P(null >= observed)."""
    null = np.asarray(null)
    return (np.sum(null >= observed) + 1) / (len(null) + 1)

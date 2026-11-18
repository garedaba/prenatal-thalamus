"""
Display parcel-wise values on the dHCP 40-week cortical surface.

Set environment variable PRENATAL_THALAMUS_SURFACES=0 to skip surface figures  
"""

import numpy as np

from . import config
from .utils import require, savefig


def parcels_to_vertices(values, parc):
    """Map one value per parcel (1..N) to surface vertices (0 = medial wall)."""
    lookup = np.r_[np.nan, np.asarray(values, dtype=float)]
    return lookup[parc.astype(int)]


def plot_parcels(values, parc, name, panel=None, cmap='RdGy_r', color_range=None, title=None):
    """Plot parcel values on the lateral and medial left hemisphere and save."""
    if not config.PLOT_SURFACES:
        print(f'  -> skipped {panel}_{name} (PRENATAL_THALAMUS_SURFACES=0)')
        return
    from surfplot import Plot   # imported here so analyses run without VTK

    require(config.DHCP_SURFACE_LEFT)
    p = Plot(str(config.DHCP_SURFACE_LEFT), zoom=1.25)
    p.add_layer({'left': parcels_to_vertices(values, parc)}, cbar=True, cmap=cmap,
                color_range=color_range)
    p.add_layer({'left': parc}, cbar=False, cmap='gray', as_outline=True)
    fig = p.build()
    if title:
        fig.axes[0].set_title(title)
    savefig(fig, name, panel)

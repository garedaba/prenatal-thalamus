"""
Paths and analysis parameters shared by all analysis steps.

All input data are read from DATA_DIR. No data are distributed with this
repository: download each dataset from its original source (see README.md) and
place it at the path given below. Entries marked PLACEHOLDER must be supplied
by the user before the corresponding step will run.

DATA_DIR defaults to ./data and can be overridden with the environment variable
PRENATAL_THALAMUS_DATA.
"""

import os
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get('PRENATAL_THALAMUS_DATA', ROOT / 'data'))
OUTPUT_DIR = ROOT / 'outputs'      # tables, models, NIfTI images
FIGURE_DIR = ROOT / 'figures'      # figure panels
MATLAB_DIR = ROOT / 'matlab_analysis'
MATLAB_EXPORT_DIR = MATLAB_DIR / 'Data4GenMdl'   # network inputs written by GetData4GenMdl.m
MATLAB_DATA_DIR = MATLAB_DIR / 'data'            # inputs read by CompileData.m

# Cortical surface figures require VTK rendering; set PRENATAL_THALAMUS_SURFACES=0 to skip them
PLOT_SURFACES = os.environ.get('PRENATAL_THALAMUS_SURFACES', '1') != '0'

# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------

# --- Preprocessed prenatal microarray (BrainSpan LMD, 4 specimens, 15-21 PCW)
# Microarray samples matched to μBrain thalamic nuclei, absent calls removed,
# one probe per gene selected by differential stability (DS) calculated over
# thalamic nuclei, and genes with DS <= 0.3 removed (Methods: Microarray data).
# `thalamus_ratio` is the per-gene ratio of mean expression in thalamus to mean
# expression in all other forebrain samples. Stricter DS cutoffs (fig. S1) are
# applied to the DS column.
# Source: https://zenodo.org/records/24297306 (raw data: https://www.brainspan.org/).
# Not interchangeable with the cortical files in the same record (.v1, .v2), which define DS differently.
#
# PLACEHOLDER: thalamic nuclei, mean expression per nucleus x specimen x gene
#   (columns used: thalamic_nucleus, specimen_id, gene, expression, DS, thalamus_ratio)
MICROARRAY_THALAMUS = DATA_DIR / 'microarray' / 'uBrain-processed-lmd-data.thalamus.csv'

# --- Prenatal thalamic atlas -----------------------------------------------
# μBrain thalamic atlas. Source: https://zenodo.org/records/24297306
# PLACEHOLDER: label image of thalamic nuclei on the dHCP 40-week template (0.15 mm)
# and its lookup table (no header; one '<label id> <nucleus name>' per line).
THALAMUS_LABELS_40W = DATA_DIR / 'thalamus_atlas' / 'uBrain-thalamus-labels-to-extdhcp40w-hires.nii.gz'
THALAMUS_LABELS_LUT = DATA_DIR / 'thalamus_atlas' / 'thalamus-nucleus-labels.txt'

# --- Preprocessed bulk RNA-seq (PsychENCODE / BrainSpan) -------------------
# RPKM (CQN + ComBat normalised), protein-coding genes, specimens up to 1 year
# postnatal, long format; thalamic samples labelled region = 'THALAMUS'
# (Methods: RNA-seq data). Source data: http://development.psychencode.org/
# PLACEHOLDER: (columns used: symbol, sample, region, rpkm, age [post-conceptional days])
RNASEQ_LONG = DATA_DIR / 'rnaseq' / 'PsychENCODE-bulk-RPKM-data.csv'

# --- Single-cell RNA-seq: Kim et al. 2023 ----------------------------------
# PLACEHOLDER: source: https://cells.ucsc.edu/?ds=dev-thal
KIM_DIR = DATA_DIR / 'scrna' / 'kim2023'
KIM_META = KIM_DIR / 'meta.tsv'
KIM_UMAP = KIM_DIR / 'UMAP.coords.tsv'
KIM_FEATURES = KIM_DIR / 'features.tsv'
KIM_MATRIX = KIM_DIR / 'matrix.mtx.gz'

# --- Mouse FO / HO marker genes: Lo Giudice et al. 2024 --------------------
# PLACEHOLDER: FO and HO marker gene lists (single column 'gene', mouse symbols)
# from Table S8 of Lo Giudice et al. Development 2024.
LOGIUDICE_FO = DATA_DIR / 'mouse' / 'logiudice2024_FO.csv'
LOGIUDICE_HO = DATA_DIR / 'mouse' / 'logiudice2024_HO.csv'

# PLACEHOLDER: Ensembl BioMart human-mouse ortholog table
# (columns: 'Gene name', 'Mouse gene name').
HUMAN_MOUSE_ORTHOLOGS = DATA_DIR / 'mouse' / 'human_mouse_orthologs.txt'

# --- Mouse thalamocortical hierarchy: Harris et al. 2019 -------------------
# PLACEHOLDER: Supplementary Table (41586_2019_1716_MOESM10_ESM.xlsx), sheet 2.
HARRIS_HIERARCHY = DATA_DIR / 'mouse' / '41586_2019_1716_MOESM10_ESM.xlsx'

# --- Neonatal connectivity (dHCP) ------------------------------------------
# PLACEHOLDER: group-average cortico-cortical connectivity (200 x 200; nodes
# 0-99 left hemisphere, 100-199 right hemisphere), generated from dHCP
# tractography (see Methods).
CONNECTIVITY_DIR = DATA_DIR / 'connectivity'
CORTICAL_CONNECTIVITY = CONNECTIVITY_DIR / 'Connectivity.csv'

# PLACEHOLDER: Euclidean distance between parcel centroids (200 x 200).
EUCLIDEAN_DISTANCE = CONNECTIVITY_DIR / 'Distance_euclidean.csv'

# PLACEHOLDER: geodesic distance between left hemisphere parcels (100 x 100),
# used to generate spatially autocorrelated surrogate maps. Written by
# matlab_analysis/GetData4GenMdl.m (as Left_geo_distance.csv); copied from
# MATLAB_EXPORT_DIR if missing.
GEODESIC_DISTANCE = CONNECTIVITY_DIR / 'Left_geo_distance.txt'

# PLACEHOLDER: 100-parcel cortical parcellation on the dHCP 40-week left
# hemisphere 32k surface (one value per vertex; 0 = medial wall).
CORTICAL_PARCELLATION = CONNECTIVITY_DIR / 'parc.csv'

# PLACEHOLDER: estimated timing of thalamic innervation of each cortical node
# (100 values): maturation-weighted thalamocortical connectivity, derived from
# nucleus maturation (step 02) and neonatal thalamocortical tractography
# (Methods). Calculated for right hemisphere parcels; the template is symmetric,
# so parcel order matches the left hemisphere. Written by
# matlab_analysis/GetData4GenMdl.m; copied from MATLAB_EXPORT_DIR if missing.
CORTICAL_INNERVATION_TIMING = CONNECTIVITY_DIR / 'cortical_innervation_timing.csv'

# dHCP 40-week symmetric surface template (very inflated, 32k); included in this
# repository. Source: https://brain-development.org/brain-atlases/
DHCP_SURFACE_LEFT = DATA_DIR / 'dhcp_surfaces' / 'week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii'

# ---------------------------------------------------------------------------
# Analysis parameters
# ---------------------------------------------------------------------------
RANDOM_SEED = 101

# Microarray specimens and age (PCW)
SPECIMEN_AGE = {'H376.IIIA.02': 15.0,
                'H376.IIIB.02': 16.0,
                'H376.IV.03': 21.0,
                'H376.IV.02': 21.57}
SPECIMENS = list(SPECIMEN_AGE)

# Gene selection (step 01)
DS_CUTOFF = 0.3               # differential stability threshold
DE_CUTOFF = 1.3               # ratio of mean expression in thalamus vs rest of brain
MAX_MISSING = 0.75            # drop genes missing in > 75% of nuclei in any specimen
EXCLUDED_NUCLEI = ['zona-incerta', 'reticular', 'lateral-habenular',
                   'medial-habenular', 'paraventricular']   # ventral thalamus / epithalamus

# Transcriptomic age model (step 01)
RIDGE_ALPHA = 1.0
CV_REPEATS = 1000
CV_TEST_SIZE = 0.5

# FO / HO designation (Fig. 2A)
NUCLEUS_ORDER = {'medial-geniculate': 'FO',
                 'lateral-geniculate': 'FO',
                 'ventral-posterior-lateral': 'FO',
                 'ventral-posterior-medial': 'FO',
                 'anterodorsal': 'HO',
                 'anteromedial': 'HO',
                 'anteroventral': 'HO',
                 'lateral-posterior': 'HO',
                 'mediodorsal': 'HO',
                 'mediodorsal (magnocellular)': 'HO',
                 'pulvinar (anterior)': 'HO',
                 'pulvinar (inferior)': 'HO',
                 'pulvinar (lateral)': 'HO',
                 'pulvinar (medial)': 'HO'}

# Cortical network and generative model (steps 05-08)
N_NODES = 100                 # cortical nodes per hemisphere
NETWORK_DENSITY = 0.25        # proportion of strongest edges retained
DENSITIES_SUPP = [0.15, 0.20, 0.30]   # fig. S8
ETA_GRID = np.linspace(0.0, 10, 50)    # distance penalty
BETA_GRID = np.linspace(0.01, 10, 50)  # attribute (dis)similarity penalty
N_SPATIAL_NULLS = 1000
N_RANDOM_NETWORKS = 500       # degree-preserving nulls for edge overlap (fig. S7)

# Nodes in primary sensory cortex (0-indexed, left hemisphere parcels).
# CHECK: originally given 1-indexed as [15, 42, 90, 2, 62, 86, 8, 98]; only
# the V1 nodes (7, 97) are identified explicitly in the original code.
PRIMARY_SENSORY_NODES = [14, 41, 89, 1, 61, 85, 7, 97]   # S1, A1, V1
V1_NODES = [7, 97]

# ---------------------------------------------------------------------------
# Derived outputs (written by earlier steps, read by later ones)
# ---------------------------------------------------------------------------
AGE_MODEL_DIR = OUTPUT_DIR / 'age_model'
MATURATION_DIR = OUTPUT_DIR / 'maturation'
MODEL_DIR = OUTPUT_DIR / 'generative_model'


def tag(ds=None, de=None):
    """Suffix identifying the gene selection parameters used for an output."""
    ds = DS_CUTOFF if ds is None else ds
    de = DE_CUTOFF if de is None else de
    return f'ds{ds}_de{de}'

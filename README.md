# prenatal-thalamus

Analysis code for:

> Oldham S, Yang JY, Lautarescu A, Bonthrone AF, Cruddas J, Tournier J-D, Batalle D & Ball G.
> **Thalamocortical constraints on areal connectivity in the developing human brain.**
> *Science Advances*, 2026. DOI to be added on publication.

## Summary

Thalamic nuclei acquire distinct molecular identities during gestation, with first-order 
(FO) sensory relays maturing before higher-order (HO) nuclei. Correspondingly thalamic 
afferents innervate the cortex in a precise order. Recent models propose that thalamic 
input to primary cortex constrains the development of wider cortical networks, directing 
highly connected hubs to form in association cortex.

This study combines post-mortem measure of gene expression, neonatal neuroimaging and network
modelling to test how thalamic maturation relates to the emergence of cortical 
network organisation:

1. **Thalamic maturation.** A transcriptomic signature of age, trained on bulk RNA-seq
   from the developing thalamus, is applied to microarray data from 30 thalamic nuclei at 
   15–21 post-conception weeks. Predicted maturation follows a lateral-to-medial gradient. 
   FO nuclei are more mature than HO nuclei at mid-gestationwith age-related genes expressed 
   by excitatory thalamic neurons.
2. **Thalamocortical connectivity.** In neonatal diffusion MRI (dHCP, n = 182),
   later-maturing nuclei connect more widely across the cortex but do not preferentially
   target cortical hubs.
3. **Generative network model.** We build a model of cortico-cortical connectivity that balances
   wiring distance and differences in the estimated timing of thalamic innervation to capture 
   variation in the spatial distribution of neonatal node degree. Removing the timing constraint 
   in primary visual cortex mirrors the macroscale effects of experimental deafferentation.


## Repository structure

```
analysis/                  numbered analysis steps, run in order
src/prenatal_thalamus/     shared code
    config.py              all data paths and analysis parameters
    expression.py          gene selection, age model, maturation index
    network.py             network data, generative model, null models
    surface.py             plotting on the dHCP cortical surface
    utils.py               input checks, reporting, figure output
matlab_analysis/           MATLAB: thalamocortical connectivity (Fig. 3) and network model inputs
data/                      shared data resources
environment.yml            conda environment
```

All paths and parameters are defined in `src/prenatal_thalamus/config.py`. Outputs are
written to `outputs/` and figure panels to `figures/`. Panel file names are prefixed with
the figure they belong to (for example `fig2c_...`, `figS7_...`).

## Analysis steps

| Step | Script | Analysis | Figures |
|---|---|---|---|
| 01 | `01_age_model.py` | Thalamus-enriched genes; ridge model of tissue age; split-half cross-validation | 2B |
| 02 | `02_nucleus_maturation.py` | Age model applied to thalamic nuclei; maturation index; robustness checks | 2C, 2D, 2E, S1–S3 |
| 03 | `03_fo_ho_hierarchy.py` | FO vs HO nuclei; mouse feedforward/feedback hierarchy; FO/HO marker genes | 2E, 2F |
| 04 | `04_single_cell_expression.py` | Expression of age-model genes in a second-trimester thalamus scRNA-seq atlas | 2G |
| 05 | `05_generative_model.py` | Generative model: timing of innervation vs distance; η–β parameter sweep | 4A, 4B, S8 |
| 06 | `06_model_fit.py` | Simulated vs empirical degree, spatial nulls; edge length and edge overlap | 4C, 4D, S7 |
| 07 | `07_null_models.py` | Distance-only and sensory-association axis models | S9 |
| 08 | `08_ablation.py` | Synthetic deafferentation of primary visual cortex | S10 |

Each script reports its key results next to the values reported in the paper.

The thalamocortical connectivity analyses (Fig. 3) are in MATLAB, in `matlab_analysis/` (see
[`matlab_analysis/README.md`](matlab_analysis/README.md)). They also export the network inputs to
steps 05–08, including the estimated timing of cortical innervation (Fig. 4A).

The following are not included in this repository:

- Construction of the 3D μBrain thalamic atlas and its alignment to the dHCP fetal and
  neonatal templates (Fig. 1). See Ball et al. 2024, *Nat Commun* 15:9685. The resulting
  thalamic atlas files are available at https://zenodo.org/records/24297306.


## Installation

```bash
conda env create -f environment.yml
conda activate prenatal-thalamus
```

## Data

**For open data not included in this repository.** Download each dataset from its original
source and place it under `data/`, or set `PRENATAL_THALAMUS_DATA` to another location.
The expected path of every file is listed in `src/prenatal_thalamus/config.py`. 

| Dataset | Used in | Source |
|---|---|---|
| Preprocessed prenatal microarray data (BrainSpan LMD, 4 specimens, 15–21 PCW), matched to μBrain thalamic nuclei(`uBrain-processed-lmd-data.thalamus.csv`)¹ | 01, 02 | https://zenodo.org/records/24297306 (raw data: https://www.brainspan.org/) |
| Preprocessed bulk tissue RNA-seq (PsychENCODE/BrainSpan RPKM) | 01, 02 | http://development.psychencode.org/ |
| μBrain thalamic atlas: nucleus labels on the dHCP 40-week template, and lookup table | 02 | https://zenodo.org/records/24297306 |
| Ensembl BioMart human–mouse orthologs | 03 | https://www.ensembl.org/biomart |
| Mouse thalamocortical hierarchy (Harris et al. 2019, *Nature*, Supplementary Table) | 03 | https://doi.org/10.1038/s41586-019-1716-z |
| Mouse FO/HO marker genes (Lo Giudice et al. 2024, *Development*, Table S8) | 03 | https://doi.org/10.1242/dev.202764 |
| Developing human thalamus scRNA-seq (Kim et al. 2023, *Science*) | 04 | https://cells.ucsc.edu/?ds=dev-thal |
| dHCP neonatal diffusion MRI (third release) | `matlab_analysis/` | NIMH Data Archive, collection 3955 (requires a data use certification; https://biomedia.github.io/dHCP-release-notes/) |
| dHCP 40-week neonatal surface template (included in `data/dhcp_surfaces/`) | 05–08 | https://brain-development.org/brain-atlases/ |

¹ Differential stability (DS) for this file was calculated over thalamic nuclei only, across
all specimen pairs. The cortical microarray files in the same record
(`uBrain-processed-lmd-data.v1.csv`, Ball et al. 2024; `uBrain-processed-lmd-data.v2.csv`,
doi:10.1038/s42003-025-08962-4) define DS differently and select different probes for some
genes, so they cannot be used here. `processing-pipelines.pdf` in the record explains the
differences.

The cortical network inputs to steps 05–08 are the group connectivity matrix, the
inter-node distances, the 100-node cortical parcellation and the estimated timing of
thalamic innervation. These group-level files are included in `data/connectivity/`, so steps
05–08 run without the dHCP imaging data. They were derived from dHCP data as described in the
Methods and exported by `matlab_analysis/GetData4GenMdl.m`.

### Expected data layout

```
data/
├── microarray/                                  # Zenodo 24297306
│   └── uBrain-processed-lmd-data.thalamus.csv
├── rnaseq/
│   └── PsychENCODE-bulk-RPKM-data.csv
├── thalamus_atlas/                              # Zenodo 24297306
│   ├── uBrain-thalamus-labels-to-extdhcp40w-hires.nii.gz
│   └── thalamus-nucleus-labels.txt
├── scrna/kim2023/                               # download from UCSC Cell Browser
│   ├── meta.tsv
│   ├── UMAP.coords.tsv
│   ├── features.tsv
│   └── matrix.mtx.gz
├── mouse/
│   ├── logiudice2024_FO.csv                     # from Lo Giudice et al. 2024, Table S8
│   ├── logiudice2024_HO.csv
│   ├── human_mouse_orthologs.txt
│   └── 41586_2019_1716_MOESM10_ESM.xlsx         # Download from Nature publication
├── connectivity/                                # included in this repository
│   ├── Connectivity.csv
│   ├── Distance_euclidean.csv
│   ├── Left_geo_distance.txt
│   ├── parc.csv
│   └── cortical_innervation_timing.csv
└── dhcp_surfaces/                               # included in this repository
    └── week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii
```

| File | Used by | Contents |
|---|---|---|
| `microarray/uBrain-processed-lmd-data.thalamus.csv` | 01, 02 | Mean expression per thalamic nucleus, specimen and gene (DS > 0.3), with each gene's DS and `thalamus_ratio` (mean expression in thalamus / rest of forebrain) |
| `rnaseq/PsychENCODE-bulk-RPKM-data.csv` | 01, 02 | Bulk RNA-seq RPKM, long format, with sample age |
| `thalamus_atlas/uBrain-thalamus-labels-to-extdhcp40w-hires.nii.gz` | 02 | Thalamic nucleus labels on the dHCP 40-week template (0.15 mm) |
| `thalamus_atlas/thalamus-nucleus-labels.txt` | 02 | Label lookup table (no header; label ID and nucleus name per line) |
| `scrna/kim2023/*` | 04 | Kim et al. 2023 cell metadata, UMAP coordinates, gene names and count matrix |
| `mouse/logiudice2024_{FO,HO}.csv` | 03 | Mouse FO and HO marker genes from Lo Giudice et al. 2024, Table S8 (column `gene`) |
| `mouse/human_mouse_orthologs.txt` | 03 | Human–mouse orthologs (`Gene name`, `Mouse gene name`) |
| `mouse/41586_2019_1716_MOESM10_ESM.xlsx` | 03 | Harris et al. 2019 hierarchy scores |
| `connectivity/Connectivity.csv` | 05–08 | Group cortico-cortical connectivity (200 × 200) |
| `connectivity/Distance_euclidean.csv` | 05–08 | Euclidean distance between nodes (200 × 200) |
| `connectivity/Left_geo_distance.txt` | 06, 07 | Geodesic distance between left hemisphere nodes (100 × 100) |
| `connectivity/parc.csv` | 05–08 | 100-node parcellation, one value per vertex of the left 32k surface |
| `connectivity/cortical_innervation_timing.csv` | 05–08 | Estimated timing of thalamic innervation, one value per node |
| `dhcp_surfaces/week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii` | 05–08 | Surface used to display cortical maps |

The full list of required columns is given in `src/prenatal_thalamus/config.py`.

## Running the analysis

From the repository root, with the environment active:

```bash
python analysis/01_age_model.py
python analysis/02_nucleus_maturation.py
python analysis/03_fo_ho_hierarchy.py
python analysis/04_single_cell_expression.py
python analysis/05_generative_model.py
python analysis/06_model_fit.py
python analysis/07_null_models.py
python analysis/08_ablation.py
```

Main parameters (`config.py`):

| Parameter | Value | Description |
|---|---|---|
| `DS_CUTOFF` | 0.3 | Differential stability threshold for microarray probes |
| `DE_CUTOFF` | 1.3 | Ratio of mean expression in thalamus vs other brain regions |
| `RIDGE_ALPHA` | 1.0 | Ridge regularisation |
| `CV_REPEATS` | 1000 | Split-half cross-validation repeats |
| `NETWORK_DENSITY` | 0.25 | Proportion of strongest cortico-cortical edges retained |
| `ETA_GRID`, `BETA_GRID` | 50 values in [0, 10] and [0.01, 10] | Model parameter sweep |
| `N_SPATIAL_NULLS` | 1000 | BrainSMASH surrogate maps |
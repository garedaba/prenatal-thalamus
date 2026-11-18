# Neonatal thalamocortical connectivity: MATLAB analyses

MATLAB code accompanying the **prenatal-thalamus** project:

> Oldham S, Yang JY, Lautarescu A, Bonthrone AF, Cruddas J, Tournier J-D, Batalle D & Ball G.  
> *Thalamocortical constraints on areal connectivity in the developing human brain.*

This component combines neonatal dHCP structural connectivity with nucleus-level estimates of prenatal thalamic maturation. It analyses how maturation relates to the distribution of thalamocortical connectivity, produces figure panels, and exports connectivity, distance and innervation estimates for the generative network model.

It covers the thalamocortical analyses associated with **Fig. 3** and the preparation of the **Fig. 4** model inputs. 

## Before running

The scripts use relative paths and shared MATLAB workspace variables. Run them from this directory, in the order below, and add the `./data` folder to the MATLAB path. Supply the required datasets and ensure the helper functions listed below are available.

## Requirements

| Requirement | Used for |
|---|---|
| MATLAB | Scripts, tables, graphs, matrix operations and figure generation.  |
| Statistics and Machine Learning Toolbox | Functions including `corr`, `pdist`, `squareform`, `zscore`, and the `nan*` functions used in the scripts. |
| [gifti toolbox](https://www.gllmflndn.com/software/matlab/gifti/) | Reading cortical surfaces and the medial-wall mask.  |
| [plotSurfaceROIBoundary](https://github.com/StuartJO/plotSurfaceROIBoundary) | Optional cortical surface maps in `MakeThalConnCortFigs.m`. |
| ImageMagick `magick` executable | Optional montage assembly in `MakeAllTractsFigs.m` and `MakeThalConnCortFigs.m`. Make sure to edit hard-coded Windows executable paths for your system. |

The included `.tck` reader is sufficient for loading the example tractogram; MATLAB plotting scripts do not call MRtrix commands.

Add this directory, `./data` and any external MATLAB helpers to your MATLAB path (including the gifti and plotSurfaceROIBoundary toolboxes). Including subfolders of `./data` makes shared resources in that folder tree available:

```matlab
cd('/path/to/matlab_code')
addpath(pwd)
addpath(genpath(fullfile(pwd,'data')))

```

## Required data

Keep the `./data` folder on the MATLAB path throughout the analysis. Inputs referenced by explicit relative paths, such as `./rand100_conns/`, must still be located at those paths relative to the working directory, or have their paths updated in the relevant scripts.

| Filename or path relative to this directory | Contents | Required by |
|---|---|---|
| `SubTCedits_mats/*_rand100.csv` | Headerless numeric connectome for each subject/session, with 200 cortical nodes followed by up to 40 thalamic labels. | `CompileData.m`, through to `Get_dHCP_Indv_TC.m` |
| `rand100_conns/*_rand100.csv` | Headerless 200 × 200 cortico-cortical connectome for each subject/session. | `CompileData.m`, through to `GetIndvRand100.m` |
| `rand100_conns/*_mu.txt` | One scalar SIFT2 μ value for each corresponding cortico-cortical connectome. | `GetIndvRand100.m` |
| `thalamic_all-regions_lut.csv` | Atlas lookup table containing `label_id` and `labels`; streamline figures also require `R`, `G`, `B` in the 0–255 range. | Compilation and streamline figures |
| `thalamic_nucleus_maturation_1.3thalamic_genes_0.3DS.csv` | Nucleus/specimen maturation estimates, with the first column identifying nucleus labels, and columns `predicted_age` and `maturation_index`. The first column is renamed `labels` before joining to the lookup table. | `CompileData.m` |
| `dHCP_wk40_100random_parc.txt` | **Included:** 32,492 vertex labels; 0 = unassigned, 1–100 = cortical parcels. | Figure analysis, cortical maps and model export |
| `NucRandPerm.mat` | **Included:** `NucPerm`, a 30 × 1000 matrix of 1-based nucleus permutation indices. | Main figure analysis |
| `dHCP40wk_surface/week-40_hemi-right_space-dhcpSym_dens-32k_vinflated.surf.gii` | Right very-inflated surface with vertex ordering matching the parcellation. | Main figure analysis and cortical maps |
| `dHCP40wk_surface/week-40_hemi-{left,right}_space-dhcpSym_dens-32k_wm.surf.gii` | Left and right white-matter surfaces. | Model export |
| `DemographicDataQCpass.mat` | MATLAB table `DemographicDataQCpass`, with `ScanAge`, `BirthAge` in decimal weeks and `sex` as strings/cells/categorical values including `female`. | Optional demographic summary |

The parent README identifies dHCP neonatal diffusion MRI and template surfaces as the imaging sources, and the μBrain thalamic atlas and upstream transcriptomic analysis as the sources of nucleus labels and maturation estimates. Use the acquisition, processing and cohort criteria described in the parent project.

### Subject filenames and node ordering

For a subject/session prefix such as `CC01236XX16_155830`, the cortical input pair is:

```text
rand100_conns/CC01236XX16_155830_rand100.csv
rand100_conns/CC01236XX16_155830_mu.txt
```

Keep only intended input CSVs in the two connectome directories. The loaders discover all `*.csv` files and expect names ending in `_rand<number>.csv`; the cortical loader then explicitly opens `<prefix>_rand100.csv`.

The thalamocortical extraction assumes this matrix ordering:

| Matrix indices | Meaning in these scripts |
|---|---|
| `1:100` | First cortical hemisphere, treated as left in the model export |
| `101:200` | Right cortical hemisphere; used for the thalamocortical analyses |
| `201:240` | Thalamic atlas IDs `1:40` |

Both connectome loaders calculate `A = (a + a') .* ~eye(length(a))`. Inputs must use a convention appropriate for this operation: a matrix that is already symmetric will have its off-diagonal values doubled.

The thalamic loader pads absent **trailing** nucleus rows with zeros. It does not remap arbitrary missing labels, so retain atlas index positions when constructing inputs. The lookup-table name extraction also assumes one leading background row followed by labels 1–40 in order.

The scripts do not perform cohort selection or imaging QC. The analysed sample is determined by the files present. The parent project describes a cohort of 182 neonates, but this number is not enforced here; ensure both connectome directories contain the intended cohort.

## Running the main workflow

From this directory, with the inputs and dependencies available:

```matlab
% Set up the MATLAB path.
addpath(pwd)
addpath(genpath(fullfile(pwd,'data')))

% 1. Compile connectivity and nucleus maturation data.
CompileData

% 2. Export inputs for the Python generative model.
% Uses variables left in the workspace by CompileData.
GetData4GenMdl

% 3. Analyse maturation-connectivity relationships and draw Fig. 3 panels.
PerformThalGramAnalysis_Fig

% 3. Draw Fig. 3 panel c.
GetThalMatuExamples
```

To restart from compiled data without reloading the individual connectomes:

```matlab
addpath(pwd)
addpath(genpath(fullfile(pwd,'data')))
load('Data4Analysis.mat')
GetData4GenMdl
PerformThalGramAnalysis_Fig
GetThalMatuExamples
```

The main figure analysis handles nucleus permutations.

## What the analysis does

### 1. Compile connectivity and maturation

`CompileData.m` calls the two individual-connectome loaders and constructs `Data4Analysis.mat`.

- **Thalamocortical connectivity:** extract thalamic rows and right cortical columns, then divide each nucleus's connectivity by its total cortical connectivity within each subject. Undefined proportions from zero-total rows become zero before group averaging.
- **Nucleus selection:** remove atlas IDs **10, 11, 12, 29 and 30**, leaving 35 nucleus rows. Apply the same exclusions to connectivity, nucleus names and maturation estimates.
- **Group thalamocortical matrix:** average the subject-level proportions, then apply `BF_NormalizeMatrix(...,'scaledSigmoid')`. This normalises **each cortical column across nuclei**, using a z-score, sigmoid and rescaling to [0,1]. The resulting `norm_TCprop` is a transformed relative-connectivity matrix; its rows no longer sum to one.
- **Cortico-cortical connectivity:** multiply each subject's symmetrised matrix by its μ scalar. Keep edges present in **strictly more than 30%** of subjects, average each edge over subjects in which it is present, and retain approximately the strongest **20%** of all possible undirected edges.
- **Cortical degree:** count positive edges in the thresholded group matrix. `RDeg` contains degree for nodes 101–200.
- **Nucleus maturation:** join maturation records to atlas labels and average `predicted_age` and `maturation_index` across matching records for each atlas ID.

`strength_threshold.m` retains all edges tied at the cutoff. The realised density can therefore differ from the requested value, and can be lower if too few positive edges survive the consistency filter.

### 2. Export model inputs

`GetData4GenMdl.m` selects a representative vertex for each parcel by minimising its mean within-parcel shortest-path distance on the white-matter mesh. It calculates inter-parcel graph distances and Euclidean distances between representative coordinates.

It also estimates a cortical value for predicted age and maturation index by weighting nucleus estimates with `norm_TCprop`:

These estimates are indirect measures derived from nucleus maturation and neonatal connectivity. They do not measure the time at which individual afferents reached the cortex.

`GfromSurf.m` weights graph edges by Euclidean mesh-edge length. Its mask removes edges only when **both** endpoints are unassigned; it retains edges between assigned and unassigned vertices. Distances are mesh-graph approximations rather than exact continuous-surface geodesics.

### 3. Test maturation-connectivity relationships

`PerformThalGramAnalysis_Fig.m` creates a partial Fig. 3 layout and computes:

| Analysis | Measure and inference |
|---|---|
| Total connectivity vs maturation | Pearson correlation between maturation index and the row sum of `norm_TCprop`; 1,000 nucleus-label permutations, two-sided p-value and Fisher-z 95% CI. |
| Connectivity spread vs maturation | Number of cortical targets with transformed connectivity `>= threshold`, at thresholds `0:0.01:1`; Pearson correlations, parametric p-values passed to `BF_FDR` at 0.05, and Fisher-z 95% CIs. |
| Selected spread threshold | Most negative correlation among thresholds passing the script's corrected-p-value mask. The printed `p-value` is the value returned by `BF_FDR`. |
| Hub connectivity vs maturation | Sum of transformed connectivity to right cortical nodes whose degree exceeds the **whole-brain 90th degree percentile**; Pearson correlation, two-sided permutation p-value and Fisher-z 95% CI. |
| Non-hub connectivity vs maturation | Equivalent analysis for the remaining right cortical nodes. |
| Participation vs maturation | `1 - sum_m (strength_to_group_m / total_strength)^2`, with hub/non-hub groups; Pearson correlation, two-sided permutation p-value and Fisher-z 95% CI. |

Permutation p-values use `(number of |null r| >= |observed r| + 1) / (1000 + 1)`. Missing maturation values are excluded from correlations.

The threshold sweep uses **parametric** p-values, while the total-strength and hub analyses use permutation p-values. The supplied `BF_FDR.m` calculates rank-based adjusted values without the cumulative-minimum step and upper clipping used in standard monotonic Benjamini–Hochberg adjusted p-values. Interpret reported values as this implementation's output. The selected threshold is chosen from the same data; its confidence interval does not account for that selection. If no threshold passes the mask, inspect the sweep before proceeding to the selected-threshold panel.

## Key outputs

### Saved analysis data

`CompileData.m` writes **`Data4Analysis.mat`**:

| Variable | Size | Contents |
|---|---|---|
| `Cconn` | 200 × 200 | Weighted, consistency-filtered and strength-thresholded group cortical connectome. |
| `Deg` | 200 × 1 | Binary node degree from `Cconn > 0`. |
| `RDeg` | 100 × 1 | Degree of nodes 101–200. |
| `norm_TCprop` | 35 × 100 | Group thalamocortical connectivity after column-wise scaled-sigmoid normalisation. |
| `ThalNucName` | 35 entries | Nucleus labels after exclusions. |
| `ThalMatEst` | 35 × 1 | Mean nucleus maturation index; unavailable estimates remain `NaN`. |
| `ThalPredAge` | 35 × 1 | Mean predicted nucleus age. |
| `Thal2Ignore` | 1 × 5 | Excluded atlas IDs: `[10 11 12 29 30]`. |
| `IsValidMatuEst` | 35 × 1 | Logical mask for non-missing maturation indices. |

### Generative-model exports

`GetData4GenMdl.m` writes the following headerless numeric files under **`Data4GenMdl/`** and **`../data/connectivity/`**:

| Filename | Size | Contents |
|---|---|---|
| `Connectivity.csv` | 200 × 200 | `Cconn`. |
| `parc.csv` | 32,492 × 1 | Surface parcel labels, including 0 for unassigned vertices. |
| `Distance_euclidean.csv` | 200 × 200 | Euclidean distance between parcel representatives, stacked left then right. |
| `Left_geo_distance.csv` | 100 × 100 | Left mesh-graph shortest-path distances between parcel representatives. |
| `Right_geo_distance.csv` | 100 × 100 | Right mesh-graph shortest-path distances between parcel representatives. |
| `cortical_innervation_PredAge.csv` | 1 × 100 | Connectivity-weighted predicted nucleus age for each right cortical parcel. |
| `cortical_innervation_MatuEst.csv` | 1 × 100 | Connectivity-weighted nucleus maturation index for each right cortical parcel. |

### Statistical and figure outputs

The main figure script prints correlation coefficients, p-values and 95% CIs to the Command Window. Statistical arrays remain in the workspace:

| Variables | Contents |
|---|---|
| `StrMatu_c`, `StrMatu_perm_p`, `StrMatu_CI` | Total transformed connectivity vs maturation. |
| `thrs`, `T_conn_thr`, `strThr_corr`, `strThr_p`, `strThr_p_pFDR`, `CI_strThr`, `BestCorr` | Threshold sweep and selected threshold. `strThr_p` is uncorrected; `strThr_p_pFDR` is returned by `BF_FDR`. |
| `HubMatu_corr`, `HubMatu_corr_perm_p`, `HubMatu_CI` | Three rows/entries: hub connectivity, non-hub connectivity, participation coefficient. |
| `Conn2Hub`, `Conn2NonHub`, `P`, `IsHub` | Nucleus-level measures and cortical hub/non-hub assignment (`1` = non-hub, `2` = hub). |

The script **does not automatically save the statistical arrays or its main figure**. To save them after running:

```matlab
save('ThalamocorticalResults.mat', ...
    'StrMatu_c','StrMatu_perm_p','StrMatu_CI', ...
    'thrs','T_conn_thr','strThr_corr','strThr_p','strThr_p_pFDR', ...
    'CI_strThr','BestCorr', ...
    'HubMatu_corr','HubMatu_corr_perm_p','HubMatu_CI', ...
    'Conn2Hub','Conn2NonHub','P','IsHub','NucPerm')

exportgraphics(gcf,'./figures/Fig3Partial.png','Resolution',300)
```

The figure layout reserves space for separate connectivity-example images. It is a partial figure, not an automatically assembled final manuscript figure.

## Optional figures and utilities

Run these in the same workspace after the main workflow:

```matlab
GetThalMatuExamples
MakeThalConnCortFigs
MakeAllTractsFigs
```

| Script | Additional requirements | Outputs |
|---|---|---|
| `GetThalMatuExamples.m` | `ThalOrd` and `CortOrd` from the main figure analysis. | `figures/HIGH_MAT_EXAMPLE.png` and `figures/LOW_MAT_EXAMPLE.png`. The high example uses the maximum maturation index; the low example uses the fifth entry after ascending sorting. |
| `MakeThalConnCortFigs.m` | `parc` in the workspace, right inflated surface, external `ExampleSurfacePlotFunction`, configured ImageMagick path. | `ThalConnCort_Figures/Nuc_01.png` … `Nuc_35.png` and `All_TCconn_Montage.png`. |
| `MakeAllTractsFigs.m` | Compiled variables, example tractogram/assignments/surface below, configured ImageMagick path. | `Thal_Streamline_Figures/Nuc_01.png` … `Nuc_35.png` and `All_Nuclei_Montage.png`. |
| `GetDemographicStats.m` | `DemographicDataQCpass.mat`. Can run independently. | Command Window summaries of birth/scan age range and median, and female count. No output file. |
| `equal_surface_area_code.m` | Right white-matter surface and `week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii` under `dHCP40wk_surface/`. | Experimental generation of 100 parcels; `new_parc` and `new_centroids` remain in the workspace. No file is saved. |

The included parcellation is sufficient for the main workflow; rerunning parcellation generation is optional. It uses random initial seeds and a dense all-pairs distance matrix, so it can be memory intensive.

### Example streamline inputs

`plot_streamlines_from_nucleus.m` hard-codes the following files:

```text
ExampleIndvTract/combined_edited_tracts_excl_subMask.tck
ExampleIndvTract/CC01236XX16_155830_thalEdits_dwi_rand100_assignments.txt
ExampleIndvTract/week-40_hemi-right_space-sub-CC01236XX16_ses-155830_dwi_wm.surf.gii
```

Edit these paths to visualise another subject. Assignments must have one row per streamline and two columns: the first is the thalamic label (`200 + atlas ID`), and the second is a right cortical label in `101:200`. The function assumes this endpoint orientation and selects all matching streamlines. The current plotting loop uses **all** matching streamlines. These figures may therefore be slow for large tractograms. Output `Nuc_XX` numbers are positions in the 35-nucleus filtered list, not original atlas IDs.

## Included helper functions

| Function | Purpose |
|---|---|
| `Get_dHCP_Indv_TC.m` | Script loading and normalising subject-level thalamocortical matrices; requires `Thal2Ignore` in the workspace. |
| `GetIndvRand100.m` | Function returning the 200 × 200 × subjects μ-scaled cortical matrices. |
| `MakeGroupAverage.m` | Consistency filtering, mean nonzero edge weights and strength thresholding; optional second output gives subject edge counts. |
| `strength_threshold.m` | Keep edges at/above a cutoff chosen by edge proportion (<1) or edge count (>=1). Requires `triu2vec`. |
| `GfromSurf.m` | Construct a weighted graph from triangular surface edges. |
| `bipartite_PC.m` | Weighted participation coefficient of each row across target groups. |
| `plot_thalamocortical_connectogram.m` | Circular bipartite display with node order, colours, sizes and edge-opacity options. |
| `plot_streamlines_from_nucleus.m`, `read_mrtrix_tracks.m` | Read and render example nucleus-to-cortex streamlines. |
| `scatterfit.m`, `find_point_on_line.m` | Scatter plots with linear fits and correlation annotations. |
| `map_values_to_colormap.m`, `coolwarm.m` | Colour mapping and a diverging colour palette. |
| `BF_NormalizeMatrix.m`, `BF_FDR.m` | Matrix normalisation and the supplied multiple-testing helper. See their headers for attribution/licensing. |
| `rank_normalise.m` | Double-sort rank transform `(rank - 1)/N + 1e-8`; ties receive distinct ranks. Included but not called in the main workflow. |

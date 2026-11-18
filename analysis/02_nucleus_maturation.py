"""
02 | Predicted maturation of thalamic nuclei at mid-gestation

The age model from step 01 is applied to microarray data from each thalamic
nucleus in each of the four prenatal specimens (15-21 PCW). Predicted age is
scaled to a maturation index reflecting the maturational state of each nucleus
relative to the rest of the thalamus.

Inputs : outputs from 01
Outputs: outputs/maturation/
           maturation_by_specimen_{tag}.csv   - nucleus x specimen predictions
           maturation_{tag}.csv               - mean over specimens, per nucleus
           maturation_dhcp40w_{tag}.nii.gz    - maturation index on 40-week template
           missing_sample_tests_{tag}.csv     - fig. S3 statistics
         matlab_analysis/data/
           thalamic_nucleus_maturation_{de}thalamic_genes_{ds}DS.csv - as maturation_{tag}.csv, read by CompileData.m
Figures: Fig. 2C, 2D (NIfTI for display), 2E (left), figs. S1, S2, S3

Results: "Gene expression indexes relative maturation of thalamic nuclei at mid-gestation"
"""

from itertools import combinations

import joblib
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.stats import pearsonr

from prenatal_thalamus import config
from prenatal_thalamus.expression import (fit_age_model, load_rnaseq_thalamus, load_thalamus_expression,
                                          model_genes, predict_maturation)
from prenatal_thalamus.utils import outdir, report, require, savefig

OUT = outdir(config.MATURATION_DIR)
TAG = config.tag()

# alternative gene selection parameters (DS, DE) for fig. S1
SUPP_PARAMETERS = [(0.3, 1.2), (0.3, 1.4), (0.5, 1.3)]


def plot_nucleus_maturation(nucleus):
    """Fig. 2C: predicted age of each nucleus, averaged over specimens."""
    nucleus = nucleus.sort_values('maturation_index')
    fig, ax = plt.subplots(figsize=(5, 6))
    ax.scatter(nucleus['predicted_age'], nucleus['thalamic_nucleus'], c=nucleus['maturation_index'],
               cmap='RdBu', s=80, edgecolor='grey')
    ax.set_xlabel('predicted age / PCW')
    ax.tick_params(axis='y', labelsize=8)
    savefig(fig, 'nucleus-maturation', panel='fig2c')


def write_maturation_image(nucleus):
    """Fig. 2D: maturation index (centred on the mean; -ve = less mature) mapped onto thalamic labels."""
    require(config.THALAMUS_LABELS_40W, config.THALAMUS_LABELS_LUT)
    img = nib.load(config.THALAMUS_LABELS_40W)
    labels = np.asarray(img.dataobj).astype(int)
    with open(config.THALAMUS_LABELS_LUT) as f:   # '<label id> <nucleus name>' per line
        lut = pd.Series({int(i): name.strip() for i, name in (line.strip().split(' ', 1) for line in f if line.strip())})
    label_to_index = lut.map(nucleus.set_index('thalamic_nucleus')['maturation_index'])

    out = np.zeros(labels.shape, dtype=np.float32)
    for label_id, value in label_to_index.dropna().items():
        out[labels == label_id] = value - nucleus['maturation_index'].mean()
    header = img.header.copy()
    header.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(out, img.affine, header), OUT / f'maturation_dhcp40w_{TAG}.nii.gz')
    print(f'  -> outputs/maturation/maturation_dhcp40w_{TAG}.nii.gz (Fig. 2D)')


def plot_maturation_by_specimen(maturation):
    """Fig. 2E (left): maturation index in each specimen."""
    order = sorted(maturation['age'].unique())
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.boxplot(x='age', y='maturation_index', hue='age', order=order, data=maturation, palette='RdBu',
                fliersize=0, saturation=0.5, legend=False, ax=ax)
    sns.stripplot(x='age', y='maturation_index', hue='age', order=order, data=maturation, palette='RdBu',
                  s=8, edgecolor='black', linewidth=1, legend=False, ax=ax)
    for _, nucleus in maturation.groupby('thalamic_nucleus'):
        nucleus = nucleus.sort_values('age')
        ax.plot([order.index(a) for a in nucleus['age']], nucleus['maturation_index'], color='grey', alpha=0.2, lw=1, zorder=0)
    ax.set_xlabel('specimen age (PCW)')
    ax.set_ylabel('maturation index')
    savefig(fig, 'maturation-by-specimen', panel='fig2e')


def specimen_correlations(maturation):
    """Fig. S2: correlation of nucleus maturation between specimens."""
    print('\nCorrelation of maturation between specimens (fig. S2)')
    wide = maturation.pivot(index='thalamic_nucleus', columns='specimen_id', values='maturation_index')
    pairs = list(combinations(config.SPECIMENS, 2))
    fig, ax = plt.subplots(figsize=(5, 5))
    for (a, b), colour in zip(pairs, sns.color_palette('Set2', len(pairs))):
        both = wide[[a, b]].dropna()
        res = pearsonr(both[a], both[b])
        ci = res.confidence_interval()   # 95% CI, Fisher z-transform
        r = res.statistic
        report(f'{a} vs {b}', f'r = {r:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}], n = {len(both)}')
        sns.regplot(x=both[a], y=both[b], ci=None, color=colour, ax=ax, label=f'{a} vs {b}: r = {r:.2f}')
    ax.set_xlabel('maturation index (specimen a)')
    ax.set_ylabel('maturation index (specimen b)')
    ax.legend(fontsize=7)
    savefig(fig, 'maturation-specimen-correlation', panel='figS2')


def missing_sample_tests(thal, nucleus):
    """Fig. S3: is maturation related to which specimens are missing for a nucleus?"""
    print('\nMaturation vs missing tissue samples (fig. S3)')
    available = thal.pivot_table(index='thalamic_nucleus', columns='specimen_id', values='expression')[config.SPECIMENS]
    df = pd.DataFrame({'maturation': nucleus.set_index('thalamic_nucleus')['maturation_index']})
    df['missing_early'] = available[['H376.IIIA.02', 'H376.IIIB.02']].isna().sum(1)   # 15/16 PCW
    df['missing_late'] = available[['H376.IV.02', 'H376.IV.03']].isna().sum(1)        # 21 PCW
    df['missing_any'] = available.isna().sum(1)
    paper_p = {'missing_early': 0.92, 'missing_late': 0.20, 'missing_any': 0.35}
    results = []
    for term in paper_p:
        fit = smf.ols(f'maturation ~ {term}', data=df).fit()
        results.append([term, fit.params[term], fit.pvalues[term]])
        report(term, f'p = {fit.pvalues[term]:.2f}', paper=f'p = {paper_p[term]}')
    pd.DataFrame(results, columns=['term', 'coef', 'p']).to_csv(OUT / f'missing_sample_tests_{TAG}.csv', index=False)

    fig, ax = plt.subplots(1, 4, figsize=(10, 4), width_ratios=(0.25, 0.25, 0.25, 0.25))
    sns.heatmap(available.notna(), cmap='Reds', cbar=False, ax=ax[0])
    ax[0].set_title('sample availability', loc='left')
    for a, term in zip(ax[1:], paper_p):
        sns.boxplot(x=term, y='maturation', data=df, hue=term, palette='Reds', legend=False, ax=a)
        sns.stripplot(x=term, y='maturation', data=df, color='black', ax=a)
        a.set_ylim(0, 1)
    fig.tight_layout()
    savefig(fig, 'missing-samples', panel='figS3')


def gene_selection_robustness(nucleus):
    """Fig. S1: predicted nucleus ages and maturation by specimen under alternative gene selection parameters."""
    print('\nAlternative gene selection parameters (fig. S1)')
    rna = load_rnaseq_thalamus()
    paper_n = {(0.3, 1.2): 459, (0.3, 1.4): 66, (0.5, 1.3): 149}
    fig, ax = plt.subplots(2, len(SUPP_PARAMETERS), figsize=(5 * len(SUPP_PARAMETERS), 10), height_ratios=(0.6, 0.4))
    for col, (ds, de) in enumerate(SUPP_PARAMETERS):
        thal = load_thalamus_expression(ds=ds, de=de)
        genes = model_genes(thal, rna)
        model, cv_r2, _, _ = fit_age_model(rna, genes)
        alt = predict_maturation(model, thal, genes)
        alt_nucleus = alt.groupby('thalamic_nucleus', as_index=False)[['predicted_age', 'maturation_index']].mean()

        # consistency with the main analysis: correlation of nucleus predicted age
        both = pd.merge(nucleus, alt_nucleus, on='thalamic_nucleus', suffixes=('_main', '_alt'))
        r = pearsonr(both['predicted_age_main'], both['predicted_age_alt'])[0]
        report(f'DS > {ds}, DE > {de}', f'{len(genes)} genes, CV R2 = {cv_r2.mean():.2f}, r with main = {r:.2f}',
               paper=f'{paper_n[(ds, de)]} genes')

        # predicted age of each nucleus (as Fig. 2C)
        alt_nucleus = alt_nucleus.sort_values('maturation_index')
        ax[0, col].scatter(alt_nucleus['predicted_age'], alt_nucleus['thalamic_nucleus'],
                           c=alt_nucleus['maturation_index'], cmap='RdBu', s=60, edgecolor='grey')
        ax[0, col].tick_params(axis='y', labelsize=7)
        ax[0, col].set_xlabel('predicted age / PCW')
        ax[0, col].set_title(f'DS = {ds}, DE = {de}, n = {len(genes)}', loc='left')

        # maturation by specimen (as Fig. 2E)
        sns.boxplot(x='age', y='maturation_index', hue='age', data=alt, palette='RdBu', fliersize=0,
                    saturation=0.5, legend=False, ax=ax[1, col])
        sns.stripplot(x='age', y='maturation_index', hue='age', data=alt, palette='RdBu', s=6,
                      edgecolor='black', linewidth=1, legend=False, ax=ax[1, col])
        ax[1, col].set_xlabel('specimen age (PCW)')
        ax[1, col].set_ylabel('maturation index')
    fig.tight_layout()
    savefig(fig, 'gene-selection-robustness', panel='figS1')


def main():
    print('\n=== 02 | Predicted maturation of thalamic nuclei ===')

    thal = load_thalamus_expression()
    genes = np.loadtxt(config.AGE_MODEL_DIR / f'thalamic_genes_{TAG}.txt', dtype=str).tolist()
    model = joblib.load(config.AGE_MODEL_DIR / f'age_model_{TAG}.joblib')

    # predicted age and maturation index per nucleus and specimen
    maturation = predict_maturation(model, thal, genes)
    maturation.to_csv(OUT / f'maturation_by_specimen_{TAG}.csv', index=False)

    nucleus = maturation.groupby('thalamic_nucleus', as_index=False)[['predicted_age', 'maturation_index']].mean()
    nucleus.to_csv(OUT / f'maturation_{TAG}.csv', index=False)
    # input to matlab_analysis/CompileData.m (first column: nucleus labels)
    nucleus.to_csv(outdir(config.MATLAB_DATA_DIR) /
                   f'thalamic_nucleus_maturation_{config.DE_CUTOFF}thalamic_genes_{config.DS_CUTOFF}DS.csv', index=False)

    print('\nNucleus maturation, mean over specimens (Fig. 2C)')
    print(nucleus.sort_values('maturation_index', ascending=False).round(2).to_string(index=False))
    plot_nucleus_maturation(nucleus)
    write_maturation_image(nucleus)

    print('\nMaturation by specimen (Fig. 2E)')
    plot_maturation_by_specimen(maturation)

    specimen_correlations(maturation)
    missing_sample_tests(thal, nucleus)
    gene_selection_robustness(nucleus)


if __name__ == '__main__':
    main()

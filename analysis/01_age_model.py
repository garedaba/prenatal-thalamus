"""
01 | Transcriptomic signature of age in the developing thalamus

1. Select genes with higher expression in the fetal thalamus than in other
   brain regions (microarray; ratio > DE_CUTOFF, DS > DS_CUTOFF).
2. Fit a ridge regression model predicting log(post-conceptional days) from
   normalised expression of these genes in PsychENCODE thalamic tissue samples
   (n = 22, 8 PCW to 1 postnatal year).
3. Evaluate with repeated split-half cross-validation (1000 repeats) and refit
   on all samples.

Inputs : config.MICROARRAY_THALAMUS, RNASEQ_LONG (preprocessed expression data)
Outputs: outputs/age_model/
           thalamic_genes_{tag}.txt        - genes used in the model (table S2)
           age_model_{tag}.joblib          - fitted model (all samples)
           age_model_coefs_{tag}.csv       - model coefficients
           age_model_cv_{tag}.csv          - split-half R2 for each repeat
Figures: Fig. 2B

Results: "Gene expression indexes relative maturation of thalamic nuclei at mid-gestation"
"""

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from prenatal_thalamus import config
from prenatal_thalamus.expression import fit_age_model, load_rnaseq_thalamus, load_thalamus_expression, model_genes
from prenatal_thalamus.utils import outdir, report, savefig

OUT = outdir(config.AGE_MODEL_DIR)
TAG = config.tag()

def main():
    print('\n=== 01 | Transcriptomic signature of age ===')

    # 1. gene selection -------------------------------------------------------
    print(f'\nGene selection (DS > {config.DS_CUTOFF}, thalamus/brain expression ratio > {config.DE_CUTOFF})')
    thal = load_thalamus_expression()
    rna = load_rnaseq_thalamus()
    genes = model_genes(thal, rna)
    report('thalamus-enriched genes', thal['gene'].nunique(), paper=207)
    report('genes also in RNA-seq', len(genes), paper=185)
    report('thalamic nuclei', thal['thalamic_nucleus'].nunique(), paper=30)
    np.savetxt(OUT / f'thalamic_genes_{TAG}.txt', genes, fmt='%s')

    # 2-3. age model ------------------------------------------------------------
    print(f'\nAge model (ridge regression, {config.CV_REPEATS} split-half repeats)')
    report('RNA-seq thalamic samples', len(rna), paper=22)
    model, cv_r2, X, y = fit_age_model(rna, genes)
    report('split-half R2 (mean ± SD)', f'{cv_r2.mean():.2f} ± {cv_r2.std():.2f}', paper='0.84 ± 0.13')

    joblib.dump(model, OUT / f'age_model_{TAG}.joblib')
    coefs = pd.Series(model.coef_, index=genes, name='coef').rename_axis('gene')
    coefs.to_csv(OUT / f'age_model_coefs_{TAG}.csv')
    pd.DataFrame({'repeat': np.arange(len(cv_r2)), 'r2': cv_r2}).to_csv(OUT / f'age_model_cv_{TAG}.csv', index=False)

    # Fig. 2B -----------------------------------------------------------------
    print()
    fig, ax = plt.subplots(2, 1, figsize=(4, 5.5), height_ratios=(0.2, 0.8))

    # tissue sample ages (top)
    ages = np.sort(rna['age'].values)
    ax[0].scatter(ages, np.zeros_like(ages), marker='|', s=200, color='black')
    ax[0].set_yticks([])
    ax[0].set_xlabel('post-conceptional days')
    ax[0].set_title('BrainSpan tissue samples', loc='left')

    # observed vs predicted age (model fit on all samples; CV performance reported above)
    ax[1].scatter(y, model.predict(X), s=80, edgecolor='grey', facecolor='white')
    ax[1].plot((4, 6.2), (4, 6.2), linestyle='--', color='grey', zorder=0)
    ax[1].set_xlim(4, 6.2)
    ax[1].set_ylim(4, 6.2)
    ax[1].set_xlabel('sample age / log(PCD)')
    ax[1].set_ylabel('predicted age / log(PCD)')
    fig.tight_layout()
    savefig(fig, 'age-model', panel='fig2b')


if __name__ == '__main__':
    main()

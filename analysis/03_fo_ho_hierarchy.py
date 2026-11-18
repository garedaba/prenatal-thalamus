"""
03 | Maturation of first-order (FO) and higher-order (HO) nuclei

1. Compare maturation of FO and HO nuclei within each specimen.
2. Correlate nucleus maturation with the mouse thalamocortical hierarchy
   (ratio of feedforward to feedback connections; Harris et al. 2019).
3. Compare age-model coefficients of genes enriched in FO vs HO neurons in
   the developing mouse thalamus (Lo Giudice et al. 2024).

Inputs : outputs from 01 and 02; config.HARRIS_HIERARCHY, LOGIUDICE_FO,
         LOGIUDICE_HO, HUMAN_MOUSE_ORTHOLOGS
Outputs: outputs/maturation/harris_hierarchy_{tag}.csv
         outputs/maturation/logiudice_coefficients_{tag}.csv
         outputs/maturation/fo_ho_hierarchy_stats_{tag}.csv   - estimates, 95% CIs, p-values
Figures: Fig. 2E (right), 2F

Results: "Gene expression indexes relative maturation of thalamic nuclei at mid-gestation"
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ptitprince as pt
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.stats import pearsonr

from prenatal_thalamus import config
from prenatal_thalamus.utils import report, require, savefig

TAG = config.tag()
FO_HO_COLOURS = {'FO': '#1f63a8', 'HO': '#a81529'}
STATS = []


def save_stats(analysis, statistic, estimate, ci_low, ci_high, p, n):
    """Collect a test result for the summary table."""
    STATS.append(dict(analysis=analysis, statistic=statistic, estimate=estimate,
                      ci95_low=ci_low, ci95_high=ci_high, p=p, n=n))

# Mouse (Allen CCF) thalamic areas -> μBrain nuclei.
# NaN: no corresponding nucleus with microarray data.
# CHECK: several mouse areas map to the same human nucleus (midline x3,
# limitans x2) and enter the correlation as separate points.
HARRIS_TO_UBRAIN = {'VPM': 'ventral-posterior-medial',
                    'VPL': 'ventral-posterior-lateral',
                    'PCN': np.nan,
                    'LGd': 'lateral-geniculate',
                    'CL': np.nan,
                    'IAD': 'midline',
                    'MG': 'medial-geniculate',
                    'AM': 'anteromedial',
                    'IMD': 'midline',
                    'MD': 'mediodorsal',
                    'PT': 'paratenial',
                    'LD': 'lateral-dorsal',
                    'PO': 'limitans',
                    'PVT': np.nan,
                    'CM': np.nan,          # central medial (not centromedian)
                    'AV': 'anteroventral',
                    'VAL': 'ventral-anterior',
                    'SMT': 'midline',
                    'LP': 'lateral-posterior',
                    'PIL': np.nan,
                    'PF': 'parafascicular',
                    'RE': 'medioventral',
                    'VM': 'ventral-medial',
                    'POL': 'limitans'}


def fo_ho_by_specimen(maturation):
    """Fig. 2E (right): mean maturation of FO and HO nuclei in each specimen."""
    print('\nMaturation of FO and HO nuclei by specimen (Fig. 2E)')
    df = maturation.assign(order=maturation['thalamic_nucleus'].map(config.NUCLEUS_ORDER)).dropna(subset=['order'])

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.barplot(x='age', y='maturation_index', hue='order', hue_order=['HO', 'FO'], data=df,
                errorbar=None, palette=FO_HO_COLOURS, ax=ax)
    sns.stripplot(x='age', y='maturation_index', hue='order', hue_order=['HO', 'FO'], data=df, dodge=True,
                  palette=FO_HO_COLOURS, edgecolor='black', linewidth=1, s=7, legend=False, ax=ax)
    ax.set_xlabel('specimen age (PCW)')
    ax.set_ylabel('maturation index')
    savefig(fig, 'fo-ho-maturation', panel='fig2e')


def harris_hierarchy(nucleus):
    """Fig. 2F (left): maturation vs mouse thalamocortical hierarchy score."""
    print('\nMaturation vs mouse thalamocortical hierarchy, Harris et al. 2019 (Fig. 2F)')
    require(config.HARRIS_HIERARCHY)
    harris = pd.read_excel(config.HARRIS_HIERARCHY, sheet_name=1)[['Area', 'Cre_conf CC+TC+CT']]
    harris = harris.rename(columns={'Cre_conf CC+TC+CT': 'hierarchy_score'})
    harris['thalamic_nucleus'] = harris['Area'].map(HARRIS_TO_UBRAIN)
    harris['maturation'] = harris['thalamic_nucleus'].map(nucleus.set_index('thalamic_nucleus')['maturation_index'])
    harris = harris.dropna()
    harris.to_csv(config.MATURATION_DIR / f'harris_hierarchy_{TAG}.csv', index=False)

    result = pearsonr(harris['hierarchy_score'], harris['maturation'])
    ci = result.confidence_interval(confidence_level=0.95)   # Fisher z-transform
    report('maturation vs hierarchy score',
           f'r = {result.statistic:.2f} [95% CI {ci.low:.2f}, {ci.high:.2f}], p = {result.pvalue:.3f}, n = {len(harris)}',
           paper='r = -0.49, p = 0.032')
    save_stats('harris_hierarchy', 'pearson_r', result.statistic, ci.low, ci.high, result.pvalue, len(harris))

    fig, ax = plt.subplots(figsize=(5, 5))
    sns.regplot(x='hierarchy_score', y='maturation', data=harris, color='black',
                scatter_kws={'s': 80, 'edgecolor': 'black', 'facecolor': '#9a2125'}, ax=ax)
    ax.set_ylim(0, 1)
    ax.set_xlabel('hierarchy score (FF>FB ← → FB>FF)')
    ax.set_ylabel('maturation index')
    ax.set_title('Harris et al. 2019', loc='right')
    savefig(fig, 'harris-hierarchy', panel='fig2f')


def logiudice_genes():
    """Fig. 2F (right): age-model coefficients of mouse FO and HO marker genes."""
    print('\nAge-model coefficients of FO and HO marker genes, Lo Giudice et al. 2024 (Fig. 2F)')
    require(config.LOGIUDICE_FO, config.LOGIUDICE_HO, config.HUMAN_MOUSE_ORTHOLOGS)
    coefs = pd.read_csv(config.AGE_MODEL_DIR / f'age_model_coefs_{TAG}.csv')
    orthologs = pd.read_csv(config.HUMAN_MOUSE_ORTHOLOGS, delimiter='\t')
    coefs['mouse'] = coefs['gene'].map(dict(zip(orthologs['Gene name'], orthologs['Mouse gene name'])))

    fo = pd.read_csv(config.LOGIUDICE_FO)['gene']
    ho = pd.read_csv(config.LOGIUDICE_HO)['gene']
    df = pd.concat([coefs[coefs['mouse'].isin(fo)].assign(order='FO'),
                    coefs[coefs['mouse'].isin(ho)].assign(order='HO')])
    df.to_csv(config.MATURATION_DIR / f'logiudice_coefficients_{TAG}.csv', index=False)

    # difference in mean coefficient, HO - FO
    fit = smf.ols('coef ~ order', data=df).fit()
    term = 'order[T.HO]'
    ci_low, ci_high = fit.conf_int(alpha=0.05).loc[term]
    report('mean coefficient, HO - FO',
           f'{fit.params[term]:.3f} [95% CI {ci_low:.3f}, {ci_high:.3f}], p = {fit.pvalues[term]:.3f} '
           f'(n FO = {(df.order == "FO").sum()}, n HO = {(df.order == "HO").sum()})')
    save_stats('logiudice_fo_ho', 'mean_difference_HO_minus_FO', fit.params[term], ci_low, ci_high,
               fit.pvalues[term], len(df))

    fig, ax = plt.subplots(figsize=(5, 5))
    palette = [FO_HO_COLOURS['FO'], FO_HO_COLOURS['HO']]
    pt.half_violinplot(x='coef', y='order', data=df, order=['FO', 'HO'], palette=palette, bw=0.3, cut=0,
                       width=0.6, inner=None, ax=ax)
    sns.stripplot(x='coef', y='order', data=df, order=['FO', 'HO'], hue='order', palette=palette, legend=False,
                  size=8, linewidth=1, edgecolor='black', zorder=0, ax=ax)
    ax.axvline(0, color='black', linestyle='--')
    ax.set_xlabel('age model coefficient (← younger | older →)')
    ax.set_ylabel('')
    ax.set_title('Lo Giudice et al. 2024', loc='right')
    savefig(fig, 'logiudice-fo-ho-genes', panel='fig2f')


def main():
    print('\n=== 03 | Maturation of FO and HO nuclei and functional hierarchy ===')
    maturation = pd.read_csv(config.MATURATION_DIR / f'maturation_by_specimen_{TAG}.csv')
    nucleus = pd.read_csv(config.MATURATION_DIR / f'maturation_{TAG}.csv')

    fo_ho_by_specimen(maturation)
    harris_hierarchy(nucleus)
    logiudice_genes()

    pd.DataFrame(STATS).to_csv(config.MATURATION_DIR / f'fo_ho_hierarchy_stats_{TAG}.csv', index=False)


if __name__ == '__main__':
    main()

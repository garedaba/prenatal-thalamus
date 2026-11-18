"""
Gene selection, transcriptomic age model and nucleus maturation (steps 01-02).
"""

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.model_selection import ShuffleSplit
from sklearn.preprocessing import MinMaxScaler

from . import config
from .utils import require


# ---------------------------------------------------------------------------
# Gene selection (microarray)
# ---------------------------------------------------------------------------
def load_thalamus_expression(ds=config.DS_CUTOFF, de=config.DE_CUTOFF):
    """
    Microarray expression for dorsal thalamic nuclei, restricted to genes with
    DS > ds, higher expression in thalamus than elsewhere (ratio > de) and data
    in at least (1 - MAX_MISSING) of nuclei in every specimen.

    The input file is already filtered at DS > 0.3, so `ds` only matters for
    stricter cutoffs (fig. S1). `thalamus_ratio` is the ratio of mean
    (log2-normalised) expression in thalamus to the rest of the brain.

    Returns long-format data (thalamic_nucleus, nuclear_complex, specimen_id,
    gene, expression, age).
    """
    require(config.MICROARRAY_THALAMUS)
    thal = pd.read_csv(config.MICROARRAY_THALAMUS)
    thal = thal[thal['DS'] > ds]
    thal = thal[~thal['thalamic_nucleus'].isin(config.EXCLUDED_NUCLEI)]

    # thalamus-enriched genes: log2(ratio) > log2(de), i.e. ratio of mean expression > de
    thal = thal[np.log2(thal['thalamus_ratio']) > np.log2(de)]

    # proportion of nuclei missing data, per gene and specimen
    wide = thal.pivot_table(index=['thalamic_nucleus', 'gene'], columns='specimen_id', values='expression')
    wide = wide[config.SPECIMENS]
    missing = wide.isna().groupby(level='gene').mean()
    keep = missing.index[missing.max(axis=1) < config.MAX_MISSING]
    thal = thal[thal['gene'].isin(keep)].copy()

    thal['age'] = thal['specimen_id'].map(config.SPECIMEN_AGE)
    return thal


# ---------------------------------------------------------------------------
# Bulk RNA-seq
# ---------------------------------------------------------------------------
def load_rnaseq_thalamus():
    """PsychENCODE thalamic samples: one row per sample, log2(RPKM) per gene."""
    require(config.RNASEQ_LONG)
    rna = pd.read_csv(config.RNASEQ_LONG)
    rna = rna[rna['region'] == 'THALAMUS'].copy()
    rna['log2rpkm'] = np.log2(rna['rpkm'])
    return rna.pivot_table(index=['sample', 'age'], columns='symbol', values='log2rpkm').reset_index()


def scale_within_sample(df, genes):
    """Min-max scale expression to [0, 1] within each sample (row)."""
    return pd.DataFrame(MinMaxScaler().fit_transform(df[genes].T).T, columns=genes, index=df.index)


# ---------------------------------------------------------------------------
# Age model
# ---------------------------------------------------------------------------
def model_genes(thal, rna):
    """Genes present in both the filtered microarray and RNA-seq data."""
    return sorted(set(thal['gene']) & set(rna.columns))


def fit_age_model(rna, genes, repeats=config.CV_REPEATS, seed=config.RANDOM_SEED):
    """
    Ridge regression of log(age in post-conceptional days) on normalised gene
    expression. Performance is estimated with repeated split-half
    cross-validation; the returned model is refit on all samples.
    """
    X = scale_within_sample(rna, genes)
    y = np.log(rna['age'].values)

    cv = ShuffleSplit(n_splits=repeats, test_size=config.CV_TEST_SIZE, random_state=seed)
    scores = []
    for train, test in cv.split(X, y):
        mdl = Ridge(alpha=config.RIDGE_ALPHA).fit(X.iloc[train], y[train])
        scores.append(r2_score(y[test], mdl.predict(X.iloc[test])))

    model = Ridge(alpha=config.RIDGE_ALPHA).fit(X, y)
    return model, np.array(scores), X, y


def predict_maturation(model, thal, genes):
    """
    Apply the age model to microarray data from each nucleus in each specimen.

    Missing values are imputed with the mean within each specimen, and
    expression is min-max scaled within each sample as for the RNA-seq data.

    Returns one row per nucleus x specimen with predicted age (PCW) and
    maturation index.
    """
    wide = thal.pivot_table(index=['thalamic_nucleus', 'specimen_id', 'age'],
                            columns='gene', values='expression').reindex(columns=genes).reset_index()

    for s in wide['specimen_id'].unique():
        rows = wide['specimen_id'] == s
        wide.loc[rows, genes] = SimpleImputer().fit_transform(wide.loc[rows, genes])

    X = scale_within_sample(wide, genes)
    predicted_age = np.exp(model.predict(X)) / 7   # days -> weeks

    out = wide[['thalamic_nucleus', 'specimen_id', 'age']].copy()
    out['predicted_age'] = predicted_age
    # maturation index: predicted age min-max scaled over all nucleus x specimen predictions
    out['maturation_index'] = MinMaxScaler().fit_transform(predicted_age.reshape(-1, 1)).ravel()
    return out

#
# Author: Miquel Anglada Girotto
# Contact: miquelangladagirotto [at] gmail [dot] com
# Last Update: 2021-02-09
#
# Script purpose
# --------------
#
#

import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn import metrics
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
from tqdm import tqdm
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs

"""
import os
ROOT = '~/repositories/target_spotter/'
model_summaries_file = os.path.join(ROOT,'data','fitted','drug_association','model_summaries.tsv.gz')
fitted_spldep_file = os.path.join(ROOT,'data','fitted','drug_association','prep_splicing_dependency.tsv.gz')
fitted_growth_rates_file = os.path.join(ROOT,'data','fitted','drug_association','growth_rates.tsv.gz')
splicing_dependency_file = os.path.join(ROOT,'workflows','splicing_dependency','mean.tsv.gz')

model_summaries = pd.read_table(model_summaries_file)
splicing_dependency = pd.read_table(splicing_dependency_file, index_col=0).iloc[:,:3]
splicing_dependency.columns = ['A','B','C']
fitted_spldep = pd.read_table(fitted_spldep_file, index_col=0)
fitted_growth_rates = pd.read_table(fitted_growth_rates_file, index_col=0)

common_samples = set(fitted_growth_rates.index).intersection(fitted_spldep.columns)
common_events = set(fitted_spldep.index).intersection(splicing_dependency.index).intersection(model_summaries["EVENT"])
splicing_dependency = splicing_dependency.loc[common_events]
fitted_spldep = fitted_spldep.loc[common_events, common_samples]
fitted_growth_rates = fitted_growth_rates.loc[common_samples]
model_summaries = model_summaries.loc[model_summaries["EVENT"].isin(common_events)]

# compute growth_rates
growth_rates = infer_growth_rates(splicing_dependency, fitted_growth_rates, fitted_spldep)

"""

##### FUNCTIONS #####
def infer_growth_rates(splicing_dependency, fitted_growth_rates, fitted_spldep):
    # combine splicing dependencies and growth rates used for model fitting
    # with new splicing dependencies
    spldep = (
        fitted_spldep.T.join(fitted_growth_rates)
        .T.add_prefix("fitted_")
        .join(splicing_dependency)
    )

    # impute growth rates
    imputer = KNNImputer()
    imputed = spldep.copy()  # create dataframe with same index and columns
    imputed.values[:, :] = imputer.fit_transform(spldep.T).T

    # prepare output
    growth_rates = pd.DataFrame(imputed.loc["growth_rate", splicing_dependency.columns])
    return growth_rates


def get_drug_pcs(drug):
    drugmat = drug.pivot_table(
        index="DRUG_ID",
        columns="ARXSPAN_ID",
        values="IC50_PUBLISHED",
        aggfunc=np.median,
    )
    drugmat = drugmat.apply(lambda x: x.fillna(np.median(x.dropna())), axis=1)
    drugmat = np.log(drugmat)
    pca = PCA(1)
    pca.fit(drugmat)
    growth_rates = pd.DataFrame(
        pca.components_.T, index=drugmat.columns, columns=["growth_rate"],
    )
    return growth_rates, pca


def get_model_lr_info(model):
    llf = model.lml()
    rank = np.linalg.matrix_rank(model.X)
    df_resid = model.nsamples - rank
    return llf, df_resid


def compare_lr_test(model_null, model_alt):
    llf_null, df_null = get_model_lr_info(model_null)
    llf_alt, df_alt = get_model_lr_info(model_alt)

    lrdf = df_null - df_alt
    lrstat = -2 * (llf_null - llf_alt)
    lr_pvalue = stats.chi2.sf(lrstat, lrdf)

    return lrstat, lr_pvalue, lrdf


def fit_limixmodel(y, X, sigma):
    event = X.columns[0]

    # fit model
    if sigma is not None:
        QS = economic_qs(sigma.loc[y.index, y.index])
    else:
        QS = None
    model = LMM(y, X, QS)
    model.fit(verbose=False)

    # get rsquared
    pred = np.dot(X, model.beta)
    rsquared = metrics.r2_score(y, pred)

    # likelihood ratio test
    model_null = LMM(y, X[["growth_rate", "intercept"]], QS)
    model_null.fit(verbose=False)
    lr_stat, lr_pvalue, lr_df = compare_lr_test(model_null, model)

    # score using test data
    pearson_coef, pearson_pvalue = stats.pearsonr(pred, y)
    spearman_coef, spearman_pvalue = stats.spearmanr(pred, y)

    # prepare output
    ## get info
    params = pd.Series(model.beta, index=X.columns)
    scanner = model.get_fast_scanner()
    bse = pd.Series(scanner.null_beta_se, index=X.columns)
    zscores = pd.Series(params / bse, index=X.columns)
    pvalues = pd.Series(stats.norm.sf(np.abs(zscores)) * 2, index=X.columns)

    # make summary
    summary = {
        "DRUG_ID": np.nan,
        "EVENT": np.nan,
        "ENSEMBL": np.nan,
        "GENE": np.nan,
        "spldep_coefficient": params[event],
        "spldep_stderr": bse[event],
        "spldep_zscore": zscores[event],
        "spldep_pvalue": pvalues[event],
        "growth_coefficient": params["growth_rate"],
        "growth_stderr": bse["growth_rate"],
        "growth_zscore": zscores["growth_rate"],
        "growth_pvalue": pvalues["growth_rate"],
        "intercept_coefficient": params["intercept"],
        "intercept_stderr": bse["intercept"],
        "intercept_zscore": zscores["intercept"],
        "intercept_pvalue": pvalues["intercept"],
        "n_obs": model.nsamples,
        "rsquared": rsquared,
        "pearson_correlation": pearson_coef,
        "pearson_pvalue": pearson_pvalue,
        "spearman_correlation": spearman_coef,
        "spearman_pvalue": spearman_pvalue,
        "lr_stat": lr_stat,
        "lr_pvalue": lr_pvalue,
        "lr_df": lr_df,
    }
    summary = pd.Series(summary)

    return summary


def fit_single_model(y, X, sigma, method):
    methods = {"limix": fit_limixmodel}
    summary = methods[method](y, X, sigma)
    return summary


def fit_model(y_drug, x_spldep, x_growth_rates, sigma, ensembl, gene, method):

    X = pd.concat([x_spldep, x_growth_rates["growth_rate"]], axis=1)
    y = y_drug

    # dropna
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan].copy()
    y = y[~is_nan].copy()

    try:
        # standardize features
        X.values[:, :] = StandardScaler().fit_transform(X)
        X["intercept"] = 1.0

        summary = fit_single_model(y, X, sigma, method)

    except:
        X["intercept"] = np.nan

        # create empty summary
        summary = pd.Series(
            np.nan,
            index=[
                "DRUG_ID",
                "EVENT",
                "ENSEMBL",
                "GENE",
                "spldep_coefficient",
                "spldep_stderr",
                "spldep_zscore",
                "spldep_pvalue",
                "growth_coefficient",  # PC1
                "growth_stderr",
                "growth_zscore",
                "growth_pvalue",
                "intercept_coefficient",
                "intercept_stderr",
                "intercept_zscore",
                "intercept_pvalue",
                "n_obs",
                "rsquared",
                "pearson_correlation",
                "pearson_pvalue",
                "spearman_correlation",
                "spearman_pvalue",
                "lr_stat",
                "lr_pvalue",
                "lr_df",
            ],
        )
    # update
    summary["DRUG_ID"] = y_drug.name
    summary["EVENT"] = x_spldep.name
    summary["ENSEMBL"] = ensembl
    summary["GENE"] = gene
    # add
    summary["spldep_mean"] = x_spldep.mean()
    summary["spldep_std"] = x_spldep.std()
    summary["growth_mean"] = x_growth_rates.mean().values[0]  # PC1
    summary["growth_std"] = x_growth_rates.std().values[0]  # PC1

    return summary


def fit_models(drug, spldep, growth_rates, mapping, n_jobs):
    sigma = spldep.cov()
    drugs_oi = drug["DRUG_ID"].unique()
    results = []
    for drug_oi in drugs_oi:
        print(drug_oi)

        # prepare drug target variable
        y_drug = drug.loc[drug["DRUG_ID"] == drug_oi]
        y_drug = pd.Series(
            np.log(y_drug["IC50_PUBLISHED"].values),  # log-normalize
            index=y_drug["ARXSPAN_ID"].values,
            name=drug_oi,
        )

        # run against all events
        res = Parallel(n_jobs=n_jobs)(
            delayed(fit_model)(
                y_drug,
                spldep.loc[event, y_drug.index],
                growth_rates.loc[y_drug.index],
                sigma,
                ensembl,
                gene,
                method="limix",
            )
            for event, ensembl, gene in tqdm(mapping.values)
        )
        res = pd.DataFrame(res)

        # compute adjusted p-values
        res["lr_padj"] = np.nan
        idx = ~res["lr_pvalue"].isnull()
        res.loc[idx, "lr_padj"] = sm.stats.multipletests(
            res.loc[idx, "lr_pvalue"], method="fdr_bh"
        )[1]
        # save
        results.append(res)

    results = pd.concat(results)
    results["DRUG_ID"] = results["DRUG_ID"].astype("int")  # formatting

    return results

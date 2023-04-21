#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#

import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn.model_selection import train_test_split
from joblib import Parallel, delayed
from tqdm import tqdm

# default variables
METHOD = "OLS"
TEST_SIZE = 0.15
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

##### FUNCTIONS #####
def get_summary_stats(df, col_oi):
    summary_stats = {
        col_oi + "_mean": np.mean(df[col_oi]),
        col_oi + "_median": np.median(df[col_oi]),
        col_oi + "_std": np.std(df[col_oi]),
        col_oi + "_q25": np.quantile(df[col_oi], 0.25),
        col_oi + "_q75": np.quantile(df[col_oi], 0.75),
    }
    return summary_stats


def fit_olsmodel(y, X, n_iterations):
    event, gene = X.columns[:2]

    summaries = []
    for i in range(n_iterations):
        # split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=TEST_SIZE, random_state=i
        )

        # fit linear model to training data
        model = sm.OLS(y_train, X_train).fit()

        # log-likelihood test
        model_null = sm.OLS(y_train, X_train[[gene, "intercept"]]).fit()
        lr_stat, lr_pvalue, lr_df = model.compare_lr_test(model_null)

        # score using test data
        prediction = model.predict(X_test)
        pearson_coef, pearson_pvalue = stats.pearsonr(prediction, y_test)
        spearman_coef, spearman_pvalue = stats.spearmanr(prediction, y_test)

        # prepare output
        summary_it = {
            "iteration": i,
            "event_coefficient": model.params[event],
            "event_stderr": model.bse[event],
            "event_zscore": model.params[event] / model.bse[event],
            "event_pvalue": model.pvalues[event],
            "gene_coefficient": model.params[gene],
            "gene_stderr": model.bse[gene],
            "gene_zscore": model.params[gene] / model.bse[gene],
            "gene_pvalue": model.pvalues[gene],
            "intercept_coefficient": model.params["intercept"],
            "intercept_stderr": model.bse["intercept"],
            "intercept_zscore": model.params["intercept"] / model.bse["intercept"],
            "intercept_pvalue": model.pvalues["intercept"],
            "n_obs": model.nobs,
            "rsquared": model.rsquared,
            "pearson_correlation": pearson_coef,
            "pearson_pvalue": pearson_pvalue,
            "spearman_correlation": spearman_coef,
            "spearman_pvalue": spearman_pvalue,
            "lr_stat": lr_stat,
            "lr_pvalue": lr_pvalue,
            "lr_df": lr_df,
        }
        summaries.append(summary_it)

    summaries = pd.DataFrame(summaries)

    # compute average likelihood-ratio test
    avg_lr_stat = np.mean(summaries["lr_stat"])
    avg_lr_df = np.round(summaries["lr_df"].mean())
    lr_pvalue = stats.chi2.sf(avg_lr_stat, avg_lr_df)

    # prepare output
    ## summary
    summary = {"EVENT": event, "ENSEMBL": gene, "GENE": y.name, "n_obs": model.nobs}
    summary.update(get_summary_stats(summaries, "event_coefficient"))
    summary.update(get_summary_stats(summaries, "gene_coefficient"))
    summary.update(get_summary_stats(summaries, "intercept_coefficient"))
    summary.update(get_summary_stats(summaries, "rsquared"))
    summary.update(get_summary_stats(summaries, "pearson_correlation"))
    summary.update(get_summary_stats(summaries, "spearman_correlation"))
    summary.update(get_summary_stats(summaries, "lr_stat"))
    summary.update(
        {"lr_df": lr_df, "lr_pvalue": lr_pvalue,}
    )
    summary = pd.Series(summary)
    ## empirical distributions of coefficients
    coefs = {
        "EVENT": summary["EVENT"],
        "GENE": summary["GENE"],
        "ENSEMBL": summary["ENSEMBL"],
        "event": summaries["event_coefficient"].values,
        "gene": summaries["gene_coefficient"].values,
        "intercept": summaries["intercept_coefficient"].values,
    }

    return summary, coefs


def fit_single_model(y, X, n_iterations, method):
    methods = {
        "OLS": fit_olsmodel,
    }
    summary, coefs = methods[method](y, X, n_iterations)
    return summary, coefs


def fit_model(x_splicing, x_genexpr, y_gene_dependency, n_iterations, method):

    X = pd.DataFrame([x_splicing, x_genexpr]).T
    y = y_gene_dependency

    # dropna
    is_nan = X.isnull().any(1) | y.isnull()
    X = X.loc[~is_nan].copy()
    y = y[~is_nan].copy()

    try:
        # standardize features
        X["intercept"] = 1.0

        summary, coefs = fit_single_model(y, X, n_iterations, method)

    except:
        X["intercept"] = np.nan

        # create empy summary
        summary = pd.Series(
            np.nan,
            index=[
                "EVENT",
                "ENSEMBL",
                "GENE",
                "n_obs",
                "event_coefficient_mean",
                "event_coefficient_median",
                "event_coefficient_std",
                "event_coefficient_q25",
                "event_coefficient_q75",
                "gene_coefficient_mean",
                "gene_coefficient_median",
                "gene_coefficient_std",
                "gene_coefficient_q25",
                "gene_coefficient_q75",
                "intercept_coefficient_mean",
                "intercept_coefficient_median",
                "intercept_coefficient_std",
                "intercept_coefficient_q25",
                "intercept_coefficient_q75",
                "rsquared_mean",
                "rsquared_median",
                "rsquared_std",
                "rsquared_q25",
                "rsquared_q75",
                "pearson_correlation_mean",
                "pearson_correlation_median",
                "pearson_correlation_std",
                "pearson_correlation_q25",
                "pearson_correlation_q75",
                "spearman_correlation_mean",
                "spearman_correlation_median",
                "spearman_correlation_std",
                "spearman_correlation_q25",
                "spearman_correlation_q75",
                "lr_stat_mean",
                "lr_stat_median",
                "lr_stat_std",
                "lr_stat_q25",
                "lr_stat_q75",
                "lr_df",
                "lr_pvalue",
            ],
        )
        summary["EVENT"] = x_splicing.name
        summary["ENSEMBL"] = x_genexpr.name
        summary["GENE"] = y_gene_dependency.name

        # create empty empirical coefficients
        coefs = {
            "EVENT": summary["EVENT"],
            "GENE": summary["GENE"],
            "ENSEMBL": summary["ENSEMBL"],
            "event": np.full(n_iterations, np.nan),
            "gene": np.full(n_iterations, np.nan),
            "intercept": np.full(n_iterations, np.nan),
        }

    # add some more info to the summary
    summary["event_mean"] = x_splicing.mean()
    summary["event_std"] = x_splicing.std()
    summary["gene_mean"] = x_genexpr.mean()
    summary["gene_std"] = x_genexpr.std()

    return summary, coefs


def get_coefs(res, coef_oi, size):
    index = ["EVENT", "GENE", "ENSEMBL"] + list(range(size))
    coefs = pd.Series(
        [res["EVENT"], res["GENE"], res["ENSEMBL"]] + list(res[coef_oi]), index=index,
    )
    return coefs


def fit_models(
    gene_dependency, splicing, genexpr, mapping, n_iterations, n_jobs, method=METHOD
):

    results = Parallel(n_jobs=n_jobs)(
        delayed(fit_model)(
            splicing.loc[event],
            genexpr.loc[ensembl],
            gene_dependency.loc[gene],
            n_iterations,
            method=method,
        )
        for event, ensembl, gene in tqdm(mapping.values)
    )

    # split results
    summaries = []
    coefs_event = []
    coefs_gene = []
    coefs_intercept = []
    for summary, coefs in results:
        summaries.append(summary)
        coefs_event.append(get_coefs(coefs, "event", n_iterations))
        coefs_gene.append(get_coefs(coefs, "gene", n_iterations))
        coefs_intercept.append(get_coefs(coefs, "intercept", n_iterations))

    summaries = pd.DataFrame(summaries)
    coefs_event = pd.DataFrame(coefs_event)
    coefs_gene = pd.DataFrame(coefs_gene)
    coefs_intercept = pd.DataFrame(coefs_intercept)

    # add FDR correction to model summaries
    summaries["lr_padj"] = np.nan
    idx = ~summaries["lr_pvalue"].isnull()
    summaries.loc[idx, "lr_padj"] = sm.stats.multipletests(
        summaries.loc[idx, "lr_pvalue"], method="fdr_bh"
    )[1]

    return summaries, coefs_event, coefs_gene, coefs_intercept

#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Compute splicing dependency using the parameters from the fitted linear
# models to predict gene dependency from splicing PSI and gene expression
# log2(TPM+1)
#
# gene_dependency = intercept + psi + genexpr + psi*genexpr
# splicing_dependency = intercept + psi + psi*genexprs


import os
import defaults
import pandas as pd
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed

GENE_LENGTHS_FILE = defaults.GENE_LENGTHS_FILE
EXAMPLE_FILES = defaults.EXAMPLE_FILES


def count_to_tpm(mrna_count, gene_lengths=None):
    if gene_lengths is None:
        gene_lengths = pd.read_table(GENE_LENGTHS_FILE, index_col=0, header=None)

    X = mrna_count / gene_lengths.loc[mrna_count.index].values
    tpm = 1e6 * X / X.sum(axis=0)
    log_tpm = np.log2(tpm + 1)

    return log_tpm


def compute_single_splicing_dependency(
    b_event, b_gene, b_interaction, b_intercept, x_splicing, x_genexpr
):

    samples = x_splicing.index
    event = x_splicing.name

    PSI = x_splicing.values.reshape(1, -1)
    TPM = x_genexpr.values.reshape(1, -1)
    PROD = PSI * TPM

    # compute
    y = b_intercept + b_event * PSI + b_gene * TPM + b_interaction * PROD

    # summarize
    mean = pd.Series(np.mean(y, axis=0), index=samples, name=event)
    median = pd.Series(np.median(y, axis=0), index=samples, name=event)
    std = pd.Series(np.std(y, axis=0), index=samples, name=event)

    summary = {"mean": mean, "median": median, "std": std}

    return summary


def compute_splicing_dependency(
    splicing,
    genexpr,
    coefs_event,
    coefs_gene,
    coefs_interaction,
    coefs_intercept,
    n_jobs,
):
    # prep coefficients
    coefs_event = coefs_event.drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    coefs_gene = coefs_gene.drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    coefs_interaction = coefs_interaction.drop(columns=["GENE"]).set_index(
        ["EVENT", "ENSEMBL"]
    )
    coefs_intercept = coefs_intercept.drop(columns=["GENE"]).set_index(
        ["EVENT", "ENSEMBL"]
    )

    # predict splicing dependency for each combination of parameters
    event_gene = coefs_event.index.to_frame()

    result = Parallel(n_jobs=n_jobs)(
        delayed(compute_single_splicing_dependency)(
            b_event=coefs_event.loc[(event, gene)].values.reshape(-1, 1),
            b_gene=coefs_gene.loc[(event, gene)].values.reshape(-1, 1),
            b_interaction=coefs_interaction.loc[(event, gene)].values.reshape(-1, 1),
            b_intercept=coefs_intercept.loc[(event, gene)].values.reshape(-1, 1),
            x_splicing=splicing.loc[event],
            x_genexpr=genexpr.loc[gene],
        )
        for event, gene in tqdm(event_gene[["EVENT", "ENSEMBL"]].values)
    )
    spldep_mean = pd.DataFrame([r["mean"] for r in result])
    spldep_median = pd.DataFrame([r["median"] for r in result])
    spldep_std = pd.DataFrame([r["std"] for r in result])

    splicing_dependency = {
        "mean": spldep_mean,
        "median": spldep_median,
        "std": spldep_std,
    }

    return splicing_dependency


def load_examples(dataset="CCLE"):
    splicing = pd.read_table(EXAMPLE_FILES[dataset]["splicing"])
    genexpr = pd.read_table(EXAMPLE_FILES[dataset]["genexpr"])

    return splicing, genexpr

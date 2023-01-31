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
import sqlalchemy
import shutil

GENE_LENGTHS_FILE = defaults.GENE_LENGTHS_FILE
EXAMPLE_FILES = defaults.EXAMPLE_FILES


def count_to_tpm(mrna_count, gene_lengths=None):
    if gene_lengths is None:
        gene_lengths = pd.read_table(GENE_LENGTHS_FILE, index_col=0, header=None)

    X = mrna_count / gene_lengths.loc[mrna_count.index].values
    X = X.replace([np.inf, -np.inf], np.nan)
    tpm = 1e6 * X / X.sum(axis=0)
    log_tpm = np.log2(tpm + 1)

    return log_tpm


def compute_single_splicing_dependency(
    b_event, b_gene, b_intercept, x_splicing, x_genexpr
):
    samples = x_splicing.index
    event = x_splicing.name
    
    PSI = x_splicing.values.reshape(1, -1)
    TPM = x_genexpr.values.reshape(1, -1)
    
    # compute
    y = b_intercept + b_event * PSI + b_gene * TPM
    
    # summarize
    mean = pd.Series(np.mean(y, axis=0), index=samples, name=event)
    median = pd.Series(np.median(y, axis=0), index=samples, name=event)
    std = pd.Series(np.std(y, axis=0), index=samples, name=event)
    q25 = pd.Series(np.quantile(y, 0.25, axis=0), index=samples, name=event)
    q75 = pd.Series(np.quantile(y, 0.75, axis=0), index=samples, name=event)

    summary = {"mean": mean, "median": median, "std": std, "q25": q25, "q75": q75}

    return summary


def compute_splicing_dependency(
    splicing, genexpr, coefs_event, coefs_gene, coefs_intercept, n_jobs,
):
    # prep coefficients
    coefs_event = coefs_event.drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    coefs_gene = coefs_gene.drop(columns=["GENE"]).set_index(["EVENT", "ENSEMBL"])
    coefs_intercept = coefs_intercept.drop(columns=["GENE"]).set_index(
        ["EVENT", "ENSEMBL"]
    )

    # predict splicing dependency for each combination of parameters
    event_gene = coefs_event.index.to_frame()

    result = Parallel(n_jobs=n_jobs)(
        delayed(compute_single_splicing_dependency)(
            b_event=coefs_event.loc[(event, gene)].values.reshape(-1, 1),
            b_gene=coefs_gene.loc[(event, gene)].values.reshape(-1, 1),
            b_intercept=coefs_intercept.loc[(event, gene)].values.reshape(-1, 1),
            x_splicing=splicing.loc[event],
            x_genexpr=genexpr.loc[gene],
        )
        for event, gene in tqdm(event_gene[["EVENT", "ENSEMBL"]].values)
    )
    spldep_mean = pd.DataFrame([r["mean"] for r in result])
    spldep_median = pd.DataFrame([r["median"] for r in result])
    spldep_std = pd.DataFrame([r["std"] for r in result])
    spldep_q25 = pd.DataFrame([r["q25"] for r in result])
    spldep_q75 = pd.DataFrame([r["q75"] for r in result])

    splicing_dependency = {
        "mean": spldep_mean,
        "median": spldep_median,
        "std": spldep_std,
        "q25": spldep_q25,
        "q75": spldep_q75,
    }

    return splicing_dependency


def compute_max_harm_score(splicing, splicing_dependency):
    """
    Delta PSI = PSI_final - PSI_initial
    Onco-event (SplDep<0): Max. DeltaPSI = 0 - PSI_initial, or
    Tumor Suppressor event (SplDep>0): Max. DeltaPSI = 100 - PSI_initial
    
    Max. Harm Score = (-1) * SplDep * DeltaPSI
    """
    # subset
    common_samples = set(splicing.columns).intersection(splicing_dependency.columns)
    common_events = set(splicing.index).intersection(splicing_dependency.index)
    splicing = splicing.loc[common_events, common_samples].copy()
    splicing_dependency = splicing_dependency.loc[common_events, common_samples].copy()

    # compute
    ## the PSI_final will depend on the SplDep sign
    psi_final = pd.DataFrame(
        np.nan, index=splicing_dependency.index, columns=splicing_dependency.columns
    )
    psi_final.values[splicing_dependency < 0] = 0  # remove onco-events
    psi_final.values[splicing_dependency > 0] = 100  # include tumor-suppressor events
    max_harm = (-1) * (psi_final - splicing) * splicing_dependency
    
    return max_harm


def load_examples(dataset="CCLE"):
    splicing = pd.read_table(EXAMPLE_FILES[dataset]["splicing"])
    genexpr = pd.read_table(EXAMPLE_FILES[dataset]["genexpr"])

    return splicing, genexpr


def save_sql(df, table_name, output_dir):
    # create SQL database
    url = os.path.join('sqlite:///'+output_dir,table_name+".sql")
    print(url)
    engine = sqlalchemy.create_engine(url, echo=False)
    df.to_sql(table_name, con=engine, if_exists='replace')
    print("Saved %s" % url)
    
            
def prep_for_webapp(
    splicing, genexpr, spldep, max_harm, pred_ic50, output_dir="webapp_inputs"
):
    # load references
    info_event_gene = pd.read_table(defaults.MAPPING_FILE)
    info_drugs = pd.read_table(defaults.INFO_DRUGS_FILE)
    
    # prepare
    ## PSI
    splicing.index.name = "EVENT"
    splicing = pd.merge(
            info_event_gene[["EVENT","GENE","ENSEMBL"]], splicing.reset_index(), 
            how="inner", left_on="EVENT", right_on="EVENT"
        )
    ## gene expression
    genexpr.index.name = "ENSEMBL"
    genexpr = pd.merge(
            info_event_gene[["GENE","ENSEMBL"]].drop_duplicates(), genexpr.reset_index(), 
            how="inner", left_on="ENSEMBL", right_on="ENSEMBL"
        )
    ## splicing dependency and maximum harm scores
    spldep.index.name = "EVENT"
    spldep = pd.merge(
        info_event_gene[["EVENT","GENE","ENSEMBL"]], spldep.reset_index(), 
        how="inner", left_on="EVENT", right_on="EVENT"
    )
    max_harm.index.name = "EVENT"
    max_harm = pd.merge(
        info_event_gene[["EVENT","GENE","ENSEMBL"]], max_harm.reset_index(), 
        how="inner", left_on="EVENT", right_on="EVENT"
    )
    ## predicted IC50s
    pred_ic50 = pd.merge(
        info_drugs[["DRUG_ID","DRUG_NAME","drug_screen","ID"]].drop_duplicates(), pred_ic50, 
        how="inner", on="ID"
    )
    
    # save
    ## create webapp directory
    os.makedirs(output_dir, exist_ok=True)
    ## save as SQL
    save_sql(splicing, "splicing", output_dir)
    save_sql(genexpr, "genexpr", output_dir)
    save_sql(spldep, "spldep", output_dir)
    save_sql(max_harm, "max_harm", output_dir)
    save_sql(pred_ic50, "pred_ic50", output_dir)
    ## archive
    shutil.make_archive(output_dir, "zip", output_dir)
    ## remove directory
    print("Compressing Webapp inputs...")
    shutil.rmtree(output_dir, ignore_errors=True)
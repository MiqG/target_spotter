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

import numpy as np
import pandas as pd
import os
import gc
import utils
from model_gene_dependency import fit_models
from make_isoform_stats import make_isoform_stats
import defaults

MAPPING_FILE = defaults.MAPPING_FILE
FITTED_SPLDEP_DIR = defaults.FITTED_SPLDEP_DIR
CCLE_STATS_FILE = defaults.CCLE_STATS_FILE
COEFS_SPLICING_FILE = defaults.COEFS_SPLICING_FILE
COEFS_GENEXPR_FILE = defaults.COEFS_GENEXPR_FILE
COEFS_INTERCEPT_FILE = defaults.COEFS_INTERCEPT_FILE
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

#### FUNCTIONS ####
class SplicingDependency:
    def __init__(
        self, normalize_counts=False, log_transform=False, n_iterations=100, n_jobs=None
    ):
        # parameters
        self.normalize_counts = normalize_counts
        self.log_transform = log_transform
        self.n_iterations = n_iterations
        self.n_jobs = n_jobs

    def _prep_fit(self):
        # unpack attributes
        gene_dependency = self.gene_dependency_
        splicing = self.splicing_
        genexpr = self.genexpr_
        isoform_stats = self.isoform_stats_
        mapping = self.mapping_

        # drop undetected & uninformative events
        splicing = splicing.dropna(thresh=2)
        splicing = splicing.loc[splicing.std(axis=1) != 0]

        # subset
        gene_annot = mapping[["ENSEMBL", "GENE"]].drop_duplicates().dropna()
        common_samples = (
            set(gene_dependency.columns)
            .intersection(splicing.columns)
            .intersection(genexpr.columns)
        )

        common_genes = set(
            gene_annot.loc[gene_annot["GENE"].isin(gene_dependency.index), "ENSEMBL"]
        ).intersection(genexpr.index)

        common_events = set(splicing.index).intersection(
            mapping.loc[mapping["ENSEMBL"].isin(common_genes), "EVENT"]
        )

        splicing = splicing.loc[common_events, common_samples]
        genexpr = genexpr.loc[common_genes, common_samples]
        gene_dependency = gene_dependency.loc[
            set(gene_annot.set_index("ENSEMBL").loc[common_genes, "GENE"]),
            common_samples,
        ]
        mapping = mapping.loc[mapping["EVENT"].isin(common_events)]

        # transform genexpr from counts to TPM?
        if self.normalize_counts:
            print("Normalizing counts to TPM...")
            genexpr = utils.count_to_tpm(genexpr)

        # log-transform TPMs?
        elif self.log_transform:
            print("Transforming TPM into log2(TPM+1)...")
            genexpr = np.log2(genexpr + 1)

        # standardize splicing and gene expression
        if isoform_stats is None:
            # compute directly from data
            isoform_stats = make_isoform_stats(splicing, genexpr, mapping)

        # keep only genes with exons
        genexpr = genexpr.loc[genexpr.index.isin(isoform_stats["ENSEMBL"])]
        
        # get summary stats
        isoform_stats = isoform_stats.set_index(["EVENT","ENSEMBL"])
        ## PSI
        splicing_mean = isoform_stats.loc[
            (splicing.index, slice(None)), "event_mean"
        ].values.reshape(-1, 1)
        splicing_std = isoform_stats.loc[
            (splicing.index, slice(None)), "event_std"
        ].values.reshape(-1, 1)
        splicing = (splicing - splicing_mean) / splicing_std
        ## TPM
        genexpr_mean = (
            isoform_stats.loc[(slice(None), genexpr.index), "gene_mean"]
            .reset_index(["ENSEMBL"])
            .drop_duplicates()["gene_mean"]
            .values.reshape(-1, 1)
        )
        genexpr_std = (
            isoform_stats.loc[(slice(None), genexpr.index), "gene_std"]
            .reset_index(["ENSEMBL"])
            .drop_duplicates()["gene_std"]
            .values.reshape(-1, 1)
        )
        genexpr = (genexpr - genexpr_mean) / genexpr_std
        
        # update attributes
        self.gene_dependency_ = gene_dependency
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.isoform_stats_ = isoform_stats
        self.mapping_ = mapping

    def fit(self, gene_dependency, splicing, genexpr, isoform_stats=None, mapping=None):
        """
        if isoform_stats is None. We compute them from data.
        """
        self.gene_dependency_ = gene_dependency
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.isoform_stats_ = isoform_stats
        self.mapping_ = mapping

        # preprocessing inputs for prediction
        print("Preprocessing inputs...")
        self._prep_fit()

        # load default mapping
        if self.mapping_ is None:
            self.mapping_ = pd.read_table(MAPPING_FILE)

        # run linear models
        (summaries, coefs_splicing, coefs_genexpr, coefs_intercept,) = fit_models(
            self.gene_dependency_,
            self.splicing_,
            self.genexpr_,
            self.mapping_,
            self.n_iterations,
            self.n_jobs,
        )

        self.summaries_ = summaries
        self.coefs_splicing_ = coefs_splicing
        self.coefs_genexpr_ = coefs_genexpr
        self.coefs_intercept_ = coefs_intercept

    def _prep_predict(self):
        # unpack attributes
        isoform_stats = self.isoform_stats_
        splicing = self.splicing_
        genexpr = self.genexpr_
        coefs_splicing = self.coefs_splicing_
        coefs_genexpr = self.coefs_genexpr_
        coefs_intercept = self.coefs_intercept_

        # subset
        ## common event - genes
        event_gene = coefs_splicing[["EVENT", "GENE", "ENSEMBL"]]
        events_avail = set(event_gene["EVENT"])
        genes_avail = set(event_gene["ENSEMBL"])

        common_events = events_avail.intersection(splicing.index).intersection(
            event_gene.loc[event_gene["ENSEMBL"].isin(genexpr.index), "EVENT"]
        )
        common_genes = genes_avail.intersection(genexpr.index).intersection(
            event_gene.loc[event_gene["EVENT"].isin(splicing.index), "ENSEMBL"]
        )
        ## common samples
        common_samples = set(splicing.columns).intersection(genexpr.columns)
        ## apply selection
        coefs_splicing = coefs_splicing.loc[
            coefs_splicing["EVENT"].isin(common_events)
            & coefs_splicing["ENSEMBL"].isin(common_genes)
        ].copy()
        coefs_genexpr = coefs_genexpr.loc[
            coefs_genexpr["EVENT"].isin(common_events)
            & coefs_genexpr["ENSEMBL"].isin(common_genes)
        ].copy()
        coefs_intercept = coefs_intercept.loc[
            coefs_intercept["EVENT"].isin(common_events)
            & coefs_intercept["ENSEMBL"].isin(common_genes)
        ].copy()
        splicing = splicing.loc[common_events, common_samples].copy()
        genexpr = genexpr.loc[common_genes, common_samples].copy()

        # transform genexpr from counts to TPM
        if self.normalize_counts:
            print("Normalizing counts to TPM...")
            genexpr = utils.count_to_tpm(genexpr)
        elif self.log_transform:
            print("Transforming TPM into log2(TPM+1)...")
            genexpr = np.log2(genexpr + 1)

        # standardize
        ## PSI
        event_mean = isoform_stats.loc[
            (splicing.index, slice(None)), "event_mean"
        ].values.reshape(-1, 1)
        event_std = isoform_stats.loc[
            (splicing.index, slice(None)), "event_std"
        ].values.reshape(-1, 1)
        splicing = (splicing - event_mean) / event_std
        ## TPM
        gene_mean = (
            isoform_stats.loc[(slice(None), genexpr.index), "gene_mean"]
            .reset_index(["ENSEMBL"])
            .drop_duplicates()["gene_mean"]
            .values.reshape(-1, 1)
        )
        gene_std = (
            isoform_stats.loc[(slice(None), genexpr.index), "gene_std"]
            .reset_index(["ENSEMBL"])
            .drop_duplicates()["gene_std"]
            .values.reshape(-1, 1)
        )
        genexpr = (genexpr - gene_mean) / gene_std

        # update attributes
        self.prep_splicing_ = splicing
        self.prep_genexpr_ = genexpr
        self.coefs_splicing_ = coefs_splicing
        self.coefs_genexpr_ = coefs_genexpr
        self.coefs_intercept_ = coefs_intercept

    def predict(
        self,
        splicing,
        genexpr,
        isoform_stats=None,
        coefs_splicing=None,
        coefs_genexpr=None,
        coefs_intercept=None,
    ):
        # prepare
        ## save inputs as attributes
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.isoform_stats_ = isoform_stats
        self.coefs_splicing_ = coefs_splicing
        self.coefs_genexpr_ = coefs_genexpr
        self.coefs_intercept_ = coefs_intercept

        ## load defaults
        print("Loading defaults...")
        if isoform_stats is None:
            self.isoform_stats_ = pd.read_table(CCLE_STATS_FILE).set_index(
                ["EVENT", "ENSEMBL"]
            )
        if coefs_splicing is None:
            self.coefs_splicing_ = pd.read_pickle(COEFS_SPLICING_FILE)
        if coefs_genexpr is None:
            self.coefs_genexpr_ = pd.read_pickle(COEFS_GENEXPR_FILE)
        if coefs_intercept is None:
            self.coefs_intercept_ = pd.read_pickle(COEFS_INTERCEPT_FILE)

        ## preprocessing inputs for prediction
        print("Preprocessing inputs...")
        self._prep_predict()

        # estimate splicing dependency
        print("Computing splicing dependencies...")
        splicing_dependency = utils.compute_splicing_dependency(
            self.prep_splicing_,
            self.prep_genexpr_,
            self.coefs_splicing_,
            self.coefs_genexpr_,
            self.coefs_intercept_,
            self.n_jobs,
        )
        self.splicing_dependency_ = splicing_dependency

        # estimate maximum harm score
        max_harm_score = {
            "mean": utils.compute_max_harm_score(
                self.splicing_, self.splicing_dependency_["mean"]
            ),
            "median": utils.compute_max_harm_score(
                self.splicing_, self.splicing_dependency_["median"]
            ),
            "q25": utils.compute_max_harm_score(
                self.splicing_, self.splicing_dependency_["q25"]
            ),
            "q75": utils.compute_max_harm_score(
                self.splicing_, self.splicing_dependency_["q75"]
            ),
        }
        self.max_harm_score_ = max_harm_score

        return self.splicing_dependency_["mean"], self.max_harm_score_["mean"]


class FitFromFiles:
    def __init__(
        self,
        gene_dependency_file,
        splicing_file,
        genexpr_file,
        isoform_stats_file=None,
        mapping_file=MAPPING_FILE,
        output_dir=FITTED_SPLDEP_DIR,
        normalize_counts=False,
        log_transform=False,
        n_iterations=100,
        n_jobs=None,
    ):

        # inputs
        self.gene_dependency_file = gene_dependency_file
        self.splicing_file = splicing_file
        self.genexpr_file = genexpr_file
        self.isoform_stats_file = isoform_stats_file
        self.mapping_file = mapping_file

        # outputs
        self.output_dir = output_dir

        # parameters
        self.normalize_counts = normalize_counts
        self.log_transform = log_transform
        self.n_iterations = n_iterations
        self.n_jobs = n_jobs

    def load_data(self):
        gene_dependency = pd.read_table(self.gene_dependency_file, index_col=0)
        splicing = pd.read_table(self.splicing_file, index_col=0)
        genexpr = pd.read_table(self.genexpr_file, index_col=0)
        mapping = pd.read_table(self.mapping_file)
        
        if self.isoform_stats_file is not None:
            isoform_stats = pd.read_table(self.isoform_stats_file).set_index(["EVENT", "ENSEMBL"])
        else:
            isoform_stats = None
        
        gc.collect()

        self.gene_dependency_ = gene_dependency
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.isoform_stats_ = isoform_stats
        self.mapping_ = mapping

    def save(self, estimator):
        isoform_stats = estimator.isoform_stats_
        summaries = estimator.summaries_
        coefs_splicing = estimator.coefs_splicing_
        coefs_genexpr = estimator.coefs_genexpr_
        coefs_intercept = estimator.coefs_intercept_

        os.makedirs(self.output_dir, exist_ok=True)

        isoform_stats.reset_index().to_csv(
            os.path.join(self.output_dir, "isoform_stats.tsv.gz"), **SAVE_PARAMS
        )
        summaries.to_csv(
            os.path.join(self.output_dir, "model_summaries.tsv.gz"), **SAVE_PARAMS
        )
        coefs_splicing.to_pickle(
            os.path.join(self.output_dir, "coefs_splicing.pickle.gz")
        )
        coefs_genexpr.to_pickle(
            os.path.join(self.output_dir, "coefs_genexpr.pickle.gz")
        )
        coefs_intercept.to_pickle(
            os.path.join(self.output_dir, "coefs_intercept.pickle.gz")
        )

    def run(self):
        print("Loading data...")
        self.load_data()
        
        print("Fitting models...")
        estimator = SplicingDependency(
            normalize_counts=self.normalize_counts,
            log_transform=self.log_transform,
            n_iterations=self.n_iterations,
            n_jobs=self.n_jobs,
        )
        estimator.fit(
            self.gene_dependency_, self.splicing_, self.genexpr_, self.isoform_stats_, self.mapping_
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)


class PredictFromFiles:
    def __init__(
        self,
        splicing_file,
        genexpr_file,
        isoform_stats_file=None,
        coefs_splicing_file=None,
        coefs_genexpr_file=None,
        coefs_intercept_file=None,
        output_dir="splicing_dependency",
        normalize_counts=False,
        log_transform=False,
        n_jobs=None,
    ):

        # inputs
        self.splicing_file = splicing_file
        self.genexpr_file = genexpr_file
        self.isoform_stats_file = isoform_stats_file
        self.coefs_splicing_file = coefs_splicing_file
        self.coefs_genexpr_file = coefs_genexpr_file
        self.coefs_intercept_file = coefs_intercept_file

        # outputs
        self.output_dir = output_dir

        # parameters
        self.normalize_counts = normalize_counts
        self.log_transform = log_transform
        self.n_jobs = n_jobs

    def load_data(self):
        splicing = pd.read_table(self.splicing_file, index_col=0)
        genexpr = pd.read_table(self.genexpr_file, index_col=0)

        # subset samples
        common_samples = set(splicing.columns).intersection(genexpr.columns)

        # take default files
        if self.isoform_stats_file is None:
            self.isoform_stats_file = CCLE_STATS_FILE
        if self.coefs_splicing_file is None:
            self.coefs_splicing_file = COEFS_SPLICING_FILE
        if self.coefs_genexpr_file is None:
            self.coefs_genexpr_file = COEFS_GENEXPR_FILE
        if self.coefs_intercept_file is None:
            self.coefs_intercept_file = COEFS_INTERCEPT_FILE

        # load coefficients
        isoform_stats = pd.read_table(self.isoform_stats_file).set_index(["EVENT", "ENSEMBL"])
        coefs_splicing = pd.read_pickle(self.coefs_splicing_file)
        coefs_genexpr = pd.read_pickle(self.coefs_genexpr_file)
        coefs_intercept = pd.read_pickle(self.coefs_intercept_file)

        gc.collect()

        # update attributes
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.isoform_stats_ = isoform_stats
        self.coefs_splicing_ = coefs_splicing
        self.coefs_genexpr_ = coefs_genexpr
        self.coefs_intercept_ = coefs_intercept

    def save(self, estimator):
        os.makedirs(self.output_dir, exist_ok=True)

        # save splicing dependency
        splicing_dependency = estimator.splicing_dependency_
        splicing_dependency["mean"].reset_index().to_csv(
            os.path.join(self.output_dir, "mean.tsv.gz"), **SAVE_PARAMS
        )
        splicing_dependency["median"].reset_index().to_csv(
            os.path.join(self.output_dir, "median.tsv.gz"), **SAVE_PARAMS
        )
        splicing_dependency["std"].reset_index().to_csv(
            os.path.join(self.output_dir, "std.tsv.gz"), **SAVE_PARAMS
        )
        splicing_dependency["q25"].reset_index().to_csv(
            os.path.join(self.output_dir, "q25.tsv.gz"), **SAVE_PARAMS
        )
        splicing_dependency["q75"].reset_index().to_csv(
            os.path.join(self.output_dir, "q75.tsv.gz"), **SAVE_PARAMS
        )

        # save max harm scores
        max_harm_score = estimator.max_harm_score_
        max_harm_score["mean"].reset_index().to_csv(
            os.path.join(self.output_dir, "max_harm_score-mean.tsv.gz"), **SAVE_PARAMS
        )
        max_harm_score["median"].reset_index().to_csv(
            os.path.join(self.output_dir, "max_harm_score-median.tsv.gz"), **SAVE_PARAMS
        )
        max_harm_score["q25"].reset_index().to_csv(
            os.path.join(self.output_dir, "max_harm_score-q25.tsv.gz"), **SAVE_PARAMS
        )
        max_harm_score["q75"].reset_index().to_csv(
            os.path.join(self.output_dir, "max_harm_score-q75.tsv.gz"), **SAVE_PARAMS
        )

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Estimating splicing dependencies...")
        estimator = SplicingDependency(
            normalize_counts=self.normalize_counts,
            log_transform=self.log_transform,
            n_jobs=self.n_jobs,
        )
        _ = estimator.predict(
            self.splicing_,
            self.genexpr_,
            self.isoform_stats_,
            self.coefs_splicing_,
            self.coefs_genexpr_,
            self.coefs_intercept_,
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)

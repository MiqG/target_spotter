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

import pandas as pd
import os
import gc
import utils
from model_gene_dependency import fit_models
import defaults

MAPPING_FILE = defaults.MAPPING_FILE
FITTED_SPLDEP_DIR = defaults.FITTED_SPLDEP_DIR
CCLE_STATS_FILE = defaults.CCLE_STATS_FILE
COEFS_SPLICING_FILE = defaults.COEFS_SPLICING_FILE
COEFS_GENEXPR_FILE = defaults.COEFS_GENEXPR_FILE
COEFS_INTERACTION_FILE = defaults.COEFS_INTERACTION_FILE
COEFS_INTERCEPT_FILE = defaults.COEFS_INTERCEPT_FILE
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

#### FUNCTIONS ####
class SplicingDependency:
    def __init__(self, normalize_counts=True, n_iterations=100, n_jobs=None):
        # parameters
        self.normalize_counts = normalize_counts
        self.n_iterations = n_iterations
        self.n_jobs = n_jobs

    def fit(self, gene_dependency, splicing, genexpr, mapping=None):
        # prepare
        ## transform genexpr from counts to TPM
        if self.normalize_counts:
            print("Normalizing counts to TPM...")
            genexpr = utils.count_to_tpm(genexpr)

        ## load default mapping
        if mapping is None:
            mapping = pd.read_table(MAPPING_FILE)

        # run linear models
        (
            summaries,
            coefs_event,
            coefs_gene,
            coefs_interaction,
            coefs_intercept,
        ) = fit_models(
            gene_dependency, splicing, genexpr, mapping, self.n_iterations, self.n_jobs,
        )

        self.summaries_ = summaries
        self.coefs_event_ = coefs_event
        self.coefs_gene_ = coefs_gene
        self.coefs_interaction_ = coefs_interaction
        self.coefs_intercept_ = coefs_intercept

    def _preprocess(self):
        # unpack attributes
        ccle_stats = self.ccle_stats_
        splicing = self.splicing_
        genexpr = self.genexpr_
        coefs_splicing = self.coefs_splicing_
        coefs_genexpr = self.coefs_genexpr_
        coefs_interaction = self.coefs_interaction_
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
        coefs_interaction = coefs_interaction.loc[
            coefs_interaction["EVENT"].isin(common_events)
            & coefs_interaction["ENSEMBL"].isin(common_genes)
        ].copy()
        coefs_interaction = coefs_interaction.loc[
            coefs_interaction["EVENT"].isin(common_events)
            & coefs_interaction["ENSEMBL"].isin(common_genes)
        ].copy()
        splicing = splicing.loc[common_events, common_samples].copy()
        genexpr = genexpr.loc[common_genes, common_samples].copy()

        # transform genexpr from counts to TPM
        if self.normalize_counts:
            print("Normalizing counts to TPM...")
            genexpr = util.count_to_tpm(genexpr)

        # standardize
        ## PSI
        event_mean = ccle_stats.loc[
            (splicing.index, slice(None)), "event_mean"
        ].values.reshape(-1, 1)
        event_std = ccle_stats.loc[
            (splicing.index, slice(None)), "event_std"
        ].values.reshape(-1, 1)
        splicing = (splicing - event_mean) / event_std
        ## TPM
        gene_mean = (
            ccle_stats.loc[(slice(None), genexpr.index), "gene_mean"]
            .reset_index(["ENSEMBL"])
            .drop_duplicates()["gene_mean"]
            .values.reshape(-1, 1)
        )
        gene_std = (
            ccle_stats.loc[(slice(None), genexpr.index), "gene_std"]
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
        self.coefs_interaction_ = coefs_interaction
        self.coefs_intercept_ = coefs_intercept

    def predict(
        self,
        splicing,
        genexpr,
        ccle_stats=None,
        coefs_splicing=None,
        coefs_genexpr=None,
        coefs_interaction=None,
        coefs_intercept=None
    ):
        # prepare
        ## save inputs as attributes
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.ccle_stats_ = ccle_stats
        self.coefs_splicing_ = coefs_splicing
        self.coefs_genexpr_ = coefs_genexpr
        self.coefs_interaction_ = coefs_interaction
        self.coefs_intercept_ = coefs_intercept

        ## load defaults
        print("Loading defaults...")
        if ccle_stats is None:
            self.ccle_stats_ = pd.read_table(CCLE_STATS_FILE).set_index(["EVENT", "ENSEMBL"])
        if coefs_splicing is None:
            self.coefs_splicing_ = pd.read_pickle(COEFS_SPLICING_FILE)
        if coefs_genexpr is None:
            self.coefs_genexpr_ = pd.read_pickle(COEFS_GENEXPR_FILE)
        if coefs_interaction is None:
            self.coefs_interaction_ = pd.read_pickle(COEFS_INTERACTION_FILE)
        if coefs_intercept is None:
            self.coefs_intercept_ = pd.read_pickle(COEFS_INTERCEPT_FILE)

        ## preprocessing inputs for prediction
        print("Preprocessing inputs...")
        self._preprocess()

        # estimate splicing dependency
        splicing_dependency = utils.compute_splicing_dependency(
            self.prep_splicing_,
            self.prep_genexpr_,
            self.coefs_splicing_,
            self.coefs_genexpr_,
            self.coefs_interaction_,
            self.coefs_intercept_,
            self.n_jobs
        )
        self.splicing_dependency_ = splicing_dependency

        return self.splicing_dependency_["median"]
    
    
class FitFromFiles:
    def __init__(
        self,
        gene_dependency_file,
        splicing_file,
        genexpr_file,
        mapping_file=MAPPING_FILE,
        output_dir=FITTED_SPLDEP_DIR,
        normalize_counts=False,
        n_iterations=100,
        n_jobs=None,
    ):

        # inputs
        self.gene_dependency_file = gene_dependency_file
        self.splicing_file = splicing_file
        self.genexpr_file = genexpr_file
        self.mapping_file = mapping_file

        # outputs
        self.output_dir = output_dir

        # parameters
        self.normalize_counts = normalize_counts
        self.n_iterations = n_iterations
        self.n_jobs = n_jobs

    def load_data(self):
        gene_dependency = pd.read_table(self.gene_dependency_file, index_col=0)
        splicing = pd.read_table(self.splicing_file, index_col=0)
        genexpr = pd.read_table(self.genexpr_file, index_col=0)
        mapping = pd.read_table(self.mapping_file).iloc[:5000]; print("DEV")

        gene_annot = mapping[["ENSEMBL", "GENE"]].drop_duplicates().dropna()

        # drop undetected & uninformative events
        splicing = splicing.dropna(thresh=2)
        splicing = splicing.loc[splicing.std(axis=1) != 0]

        # subset
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

        gc.collect()

        self.gene_dependency_ = gene_dependency
        self.splicing_ = splicing
        self.genexpr_ = genexpr
        self.mapping_ = mapping

    def save(self, estimator):
        summaries = estimator.summaries_
        coefs_event = estimator.coefs_event_
        coefs_gene = estimator.coefs_gene_
        coefs_interaction = estimator.coefs_interaction_
        coefs_intercept = estimator.coefs_intercept_

        os.makedirs(self.output_dir, exist_ok=True)

        summaries.to_csv(
            os.path.join(self.output_dir, "model_summaries.tsv.gz"), **SAVE_PARAMS
        )
        coefs_event.to_pickle(os.path.join(self.output_dir, "coefs_event.pickle.gz"))
        coefs_gene.to_pickle(os.path.join(self.output_dir, "coefs_gene.pickle.gz"))
        coefs_interaction.to_pickle(
            os.path.join(self.output_dir, "coefs_interaction.pickle.gz")
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
            n_iterations=self.n_iterations,
            n_jobs=self.n_jobs,
        )
        estimator.fit(
            self.gene_dependency_, self.splicing_, self.genexpr_, self.mapping_
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)
        
        
class PredictFromFiles:
    def __init__(
        self,
        splicing_file,
        genexpr_file,
        ccle_stats_file=None,
        coefs_splicing_file=None,
        coefs_genexpr_file=None,
        coefs_interaction_file=None,
        coefs_intercept_file=None,
        output_dir="splicing_dependency",
        normalize_counts=False,
        n_iterations=100,
        n_jobs=None,
    ):

        # inputs
        self.splicing_file = splicing_file
        self.genexpr_file = genexpr_file
        self.ccle_stats_file = ccle_stats_file
        self.coefs_splicing_file = coefs_splicing_file
        self.coefs_genexpr_file = coefs_genexpr_file
        self.coefs_interaction_file = coefs_interaction_file
        self.coefs_intercept_file = coefs_intercept_file

        # outputs
        self.output_dir = output_dir

        # parameters
        self.normalize_counts = normalize_counts
        self.n_iterations = n_iterations
        self.n_jobs = n_jobs

    def load_data(self):
        splicing = pd.read_table(self.splicing_file, index_col=0)
        genexpr = pd.read_table(self.genexpr_file, index_col=0)

        # subset samples
        common_samples = set(splicing.columns).intersection(genexpr.columns)
        
        # TODO load stats and coefs if not None
        
        gc.collect()

        self.splicing_ = splicing
        self.genexpr_ = genexpr

    def save(self, estimator):
        splicing_dependency = estimator.splicing_dependency_

        os.makedirs(self.output_dir, exist_ok=True)

        splicing_dependency["mean"].to_csv(
            os.path.join(self.output_dir, "mean.tsv.gz"), **SAVE_PARAMS
        )
        splicing_dependency["median"].to_csv(
            os.path.join(self.output_dir, "median.tsv.gz"), **SAVE_PARAMS
        )
        splicing_dependency["std"].to_csv(
            os.path.join(self.output_dir, "std.tsv.gz"), **SAVE_PARAMS
        )

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Estimating splicing dependencies...")
        estimator = SplicingDependency(
            normalize_counts=self.normalize_counts,
            n_iterations=self.n_iterations,
            n_jobs=self.n_jobs,
        )
        _ = estimator.predict(
            self.splicing_,
            self.genexpr_,
#             self.ccle_stats_,
#             self.coefs_splicing_,
#             self.coefs_genexpr_,
#             self.coefs_interaction_,
#             self.coefs_intercept_,
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)
        

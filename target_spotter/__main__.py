# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
#

import pandas as pd
import os
import argparse
import gc
import defaults
from SplicingDependency import SplicingDependency


MAPPING_FILE = defaults.MAPPING_FILE
FITTED_DIR = defaults.FITTED_DIR

SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}


"""
Development
-----------
# fit
ROOT = '/home/miquel/repositories/target_spotter'
PREP_DIR = os.path.join(ROOT,'data','prep')
splicing_file = os.path.join(PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
genexpr_file = os.path.join(PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
gene_dependency_file = os.path.join(PREP_DIR,'demeter2','CCLE.tsv.gz')
normalize_counts=False
n_jobs=10
n_iterations=2
args = {
    "cmd":{
        "gene_dependency_file": gene_dependency_file,
        "splicing_file": splicing_file,
        "genexpr_file": genexpr_file,
        "output_dir": output_dir
    }
}

# predict

"""


class fit:
    def __init__(
        self,
        gene_dependency_file,
        splicing_file,
        genexpr_file,
        mapping_file=MAPPING_FILE,
        output_dir=FITTED_DIR,
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
        mapping = pd.read_table(self.mapping_file).iloc[:5000]  ##

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

        os.makedirs(self.output_dir)

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

        print("Saving results to %s..." % self.output_dir)
        self.save(estimator)


class predict:
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

        os.makedirs(self.output_dir)

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

        print("Saving results to %s..." % self.output_dir)
        self.save(estimator)


def parse_args():
    parser = argparse.ArgumentParser(
        prog="target_spotter",
        description="Systematic prioritization of splicing targets to treat cancer.",
    )
    subparser = parser.add_subparsers(help="sub-commands", dest="cmd")

    # target_spotter fit
    fit_parser = subparser.add_parser("fit", help="model gene dependency.")
    fit_parser.add_argument("--gene_dependency_file", type=str, required=True)
    fit_parser.add_argument("--splicing_file", type=str, required=True)
    fit_parser.add_argument("--genexpr_file", type=str, required=True)
    fit_parser.add_argument("--mapping_file", type=str, default=MAPPING_FILE)
    fit_parser.add_argument("--output_dir", type=str, default=FITTED_DIR)
    fit_parser.add_argument("--normalize_counts", type=bool, default=False)
    fit_parser.add_argument("--n_iterations", type=int, default=100)
    fit_parser.add_argument("--n_jobs", type=int, default=1)

    # target_spotter predict
    pred_parser = subparser.add_parser("predict", help="estimate splicing dependency.")
    pred_parser.add_argument("--splicing_file", type=str, required=True)
    pred_parser.add_argument("--genexpr_file", type=str, required=True)
    pred_parser.add_argument("--ccle_stats_file", type=str, default=None)
    pred_parser.add_argument("--coefs_splicing_file", type=str, default=None)
    pred_parser.add_argument("--coefs_genexpr_file", type=str, default=None)
    pred_parser.add_argument("--coefs_interaction_file", type=str, default=None)
    pred_parser.add_argument("--coefs_intercept_file", type=str, default=None)
    pred_parser.add_argument("--output_dir", type=str, default="splicing_dependency")
    pred_parser.add_argument("--normalize_counts", type=bool, default=False)
    pred_parser.add_argument("--n_jobs", type=int, default=1)

    # get arguments
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.cmd == "fit":
        fit(
            gene_dependency_file=args.gene_dependency_file,
            splicing_file=args.splicing_file,
            genexpr_file=args.genexpr_file,
            mapping_file=args.mapping_file,
            output_dir=args.output_dir,
            normalize_counts=args.normalize_counts,
            n_iterations=args.n_iterations,
            n_jobs=args.n_jobs,
        ).run()

    elif args.cmd == "predict":
        predict(
            splicing_file = args.splicing_file,
            genexpr_file = args.genexpr_file,
            ccle_stats_file = args.ccle_stats_file,
            coefs_splicing_file = args.coefs_splicing_file,
            coefs_genexpr_file = args.coefs_genexpr_file,
            coefs_interaction_file = args.coefs_interaction_file,
            coefs_intercept_file = args.coefs_intercept_file,
            output_dir = args.output_dir,
            normalize_counts = args.normalize_counts,
            n_jobs = args.n_jobs
        ).run()

    # fit
    # fit(**args["cmd"]).run()


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")

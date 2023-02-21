# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
#

import pandas as pd
import os
import argparse
import defaults, SplicingDependency, OneSampleDiff, DrugAssociation

MAPPING_FILE = defaults.MAPPING_FILE
FITTED_SPLDEP_DIR = defaults.FITTED_SPLDEP_DIR
FITTED_ONEDIFF_DIR = defaults.FITTED_ONEDIFF_DIR
FITTED_DRUGASSOC_DIR = defaults.FITTED_DRUGASSOC_DIR

##### FUNCTIONS #####
def parse_args():
    parser = argparse.ArgumentParser(
        prog="target_spotter",
        description="Systematic prioritization of splicing targets to treat cancer.",
    )
    subparser = parser.add_subparsers(help="sub-commands", dest="cmd")

    # SplicingDependency
    ## target_spotter spldep_fit
    fit_parser = subparser.add_parser("spldep_fit", help="fit models of splicing dependency.")
    fit_parser.add_argument("--gene_dependency_file", type=str, required=True)
    fit_parser.add_argument("--splicing_file", type=str, required=True)
    fit_parser.add_argument("--genexpr_file", type=str, required=True)
    fit_parser.add_argument("--isoform_stats_file", type=str, default=None)
    fit_parser.add_argument("--mapping_file", type=str, default=MAPPING_FILE)
    fit_parser.add_argument("--output_dir", type=str, default=FITTED_SPLDEP_DIR)
    fit_parser.add_argument("--normalize_counts", action="store_true")
    fit_parser.add_argument("--log_transform", action="store_true")
    fit_parser.add_argument("--n_iterations", type=int, default=100)
    fit_parser.add_argument("--n_jobs", type=int, default=1)

    ## target_spotter spldep_predict
    pred_parser = subparser.add_parser("spldep_predict", help="estimate splicing dependency.")
    pred_parser.add_argument("--splicing_file", type=str, required=True)
    pred_parser.add_argument("--genexpr_file", type=str, required=True)
    pred_parser.add_argument("--isoform_stats_file", type=str, default=None)
    pred_parser.add_argument("--coefs_splicing_file", type=str, default=None)
    pred_parser.add_argument("--coefs_genexpr_file", type=str, default=None)
    pred_parser.add_argument("--coefs_intercept_file", type=str, default=None)
    pred_parser.add_argument("--output_dir", type=str, default="splicing_dependency")
    pred_parser.add_argument("--normalize_counts", action="store_true")
    pred_parser.add_argument("--log_transform", action="store_true")
    pred_parser.add_argument("--n_jobs", type=int, default=1)

    ## target_spotter drugassoc_fit
    fit_parser = subparser.add_parser("drugassoc_fit", help="estimate splicing dependency.")
    fit_parser.add_argument("--drug_response_file", type=str, required=True)
    fit_parser.add_argument("--splicing_dependency_file", type=str, required=True)
    fit_parser.add_argument("--growth_rates_file", type=str, default=None)
    fit_parser.add_argument("--mapping_file", type=str, default=MAPPING_FILE)
    fit_parser.add_argument("--selected_models_file", type=str, default=None)
    fit_parser.add_argument("--output_dir", type=str, default=FITTED_DRUGASSOC_DIR)
    fit_parser.add_argument("--n_jobs", type=int, default=1)

    ## target_spotter drugassoc_predict
    pred_parser = subparser.add_parser("drugassoc_predict", help="estimate drug response.")
    pred_parser.add_argument("--splicing_dependency_file", type=str, required=True)
    pred_parser.add_argument("--growth_rates_file", type=str, default=None)
    pred_parser.add_argument("--model_summaries_file", type=str, default=None)
    pred_parser.add_argument("--fitted_growth_rates_file", type=str, default=None)
    pred_parser.add_argument("--fitted_spldep_file", type=str, default=None)
    pred_parser.add_argument("--output_dir", type=str, default="drug_association")

    # get arguments
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print(args)
    
    # SplicingDependency
    if args.cmd == "spldep_fit":
        SplicingDependency.FitFromFiles(
            gene_dependency_file=args.gene_dependency_file,
            splicing_file=args.splicing_file,
            genexpr_file=args.genexpr_file,
            isoform_stats_file=args.isoform_stats_file,
            mapping_file=args.mapping_file,
            output_dir=args.output_dir,
            normalize_counts=args.normalize_counts,
            log_transform=args.log_transform,
            n_iterations=args.n_iterations,
            n_jobs=args.n_jobs,
        ).run()

    elif args.cmd == "spldep_predict":
        SplicingDependency.PredictFromFiles(
            splicing_file=args.splicing_file,
            genexpr_file=args.genexpr_file,
            isoform_stats_file=args.isoform_stats_file,
            coefs_splicing_file=args.coefs_splicing_file,
            coefs_genexpr_file=args.coefs_genexpr_file,
            coefs_intercept_file=args.coefs_intercept_file,
            output_dir=args.output_dir,
            normalize_counts=args.normalize_counts,
            log_transform=args.log_transform,
            n_jobs=args.n_jobs,
        ).run()

    elif args.cmd == "drugassoc_fit":
        DrugAssociation.FitFromFiles(
            drug_response_file=args.drug_response_file,
            splicing_dependency_file=args.splicing_dependency_file,
            growth_rates_file=args.growth_rates_file,
            mapping_file=args.mapping_file,
            selected_models_file=args.selected_models_file,
            n_jobs=args.n_jobs,
            output_dir=args.output_dir,
        ).run()

    elif args.cmd == "drugassoc_predict":
        DrugAssociation.PredictFromFiles(
            splicing_dependency_file=args.splicing_dependency_file,
            growth_rates_file=args.growth_rates_file,
            model_summaries_file=args.model_summaries_file,
            fitted_growth_rates_file=args.fitted_growth_rates_file,
            fitted_spldep_file=args.fitted_spldep_file,
            output_dir=args.output_dir,
        ).run()

##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")

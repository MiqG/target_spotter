# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
#

import pandas as pd
import os
import argparse
import defaults, SplicingDependency, DrugAssociation

MAPPING_FILE = defaults.MAPPING_FILE
FITTED_SPLDEP_DIR = defaults.FITTED_SPLDEP_DIR
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
    fit_parser = subparser.add_parser(
        "spldep_fit", help="fit models of splicing dependency."
    )
    fit_parser.add_argument(
        "--gene_dependency_file",
        type=str,
        required=True,
        help="Tab-separated file of gene-level dependencies. We assume the first column is the index.",
    )
    fit_parser.add_argument(
        "--splicing_file",
        type=str,
        required=True,
        help="Tab-separated file of exon inclusion PSIs. We assume the first column is the index in VastDB format.",
    )
    fit_parser.add_argument(
        "--genexpr_file",
        type=str,
        required=True,
        help="Tab-separated file of gene expression TPMs, log2(TPMs+1), or raw counts. We assume the first column is the index in ENSEMBL format and that inputs are in log2(TPMs+1).",
    )
    fit_parser.add_argument(
        "--isoform_stats_file",
        type=str,
        default=None,
        help="Tab-separated file generated running `target_spotter.make_isoform_stats`. If no file is given we compute them directly from splicing and genexpr files.",
    )
    fit_parser.add_argument(
        "--mapping_file",
        type=str,
        default=MAPPING_FILE,
        help="Tab-separated file mapping VastDB splicing event identifiers ('EVENT' column) to their corresponding ENSEMBL gene identifier ('ENSEMBL' column).",
    )
    fit_parser.add_argument(
        "--output_dir",
        type=str,
        default=FITTED_SPLDEP_DIR,
        help="Path to the output directory.",
    )
    fit_parser.add_argument(
        "--normalize_counts",
        action="store_true",
        help="If your 'genexpr_file' contains raw gene expression counts, add this argument so we will normalize and log transform these counts.",
    )
    fit_parser.add_argument(
        "--log_transform",
        action="store_true",
        help="If your 'genexpr_file' contains TPM not log-transformed, add this flag so we will log-transform the TPMs.",
    )
    fit_parser.add_argument(
        "--n_iterations",
        type=int,
        default=100,
        help="Our fitting pipeline splits the data into training and test set 'n_iterations' times to obtain the distribution of coefficients.",
    )
    fit_parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of threads to use."
    )

    ## target_spotter spldep_predict
    pred_parser = subparser.add_parser(
        "spldep_predict", help="estimate splicing dependency."
    )
    pred_parser.add_argument(
        "--splicing_file",
        type=str,
        required=True,
        help="Tab-separated file of exon inclusion PSIs. We assume the first column is the index in VastDB format.",
    )
    pred_parser.add_argument(
        "--genexpr_file",
        type=str,
        required=True,
        help="Tab-separated file of gene expression TPMs, log2(TPMs+1), or raw counts. We assume the first column is the index in ENSEMBL format and that inputs are in log2(TPMs+1).",
    )
    pred_parser.add_argument(
        "--isoform_stats_file",
        type=str,
        default=None,
        help="Tab-separated file generated running `target_spotter.make_isoform_stats`. If no file is given we retrieve `isoform_stats` used during training.",
    )
    pred_parser.add_argument(
        "--coefs_splicing_file",
        type=str,
        default=None,
        help="Splicing coefficients from our linear models outputted during training.",
    )
    pred_parser.add_argument(
        "--coefs_genexpr_file",
        type=str,
        default=None,
        help="Gene expression coefficients from our linear models outputted during training.",
    )
    pred_parser.add_argument(
        "--coefs_intercept_file",
        type=str,
        default=None,
        help="Intercept coefficients from our linear models outputted during training.",
    )
    pred_parser.add_argument(
        "--output_dir",
        type=str,
        default="splicing_dependency",
        help="Path to the output directory.",
    )
    pred_parser.add_argument(
        "--normalize_counts",
        action="store_true",
        help="If your 'genexpr_file' contains raw gene expression counts, add this argument so we will normalize and log transform these counts.",
    )
    pred_parser.add_argument(
        "--log_transform",
        action="store_true",
        help="If your 'genexpr_file' contains TPM not log-transformed, add this flag so we will log-transform the TPMs.",
    )
    pred_parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of threads to use."
    )

    ## target_spotter drugassoc_fit
    fit_parser = subparser.add_parser(
        "drugassoc_fit", help="estimate splicing dependency."
    )
    fit_parser.add_argument(
        "--drug_response_file",
        type=str,
        required=True,
        help="Tab-separated file containing drug sensitivity scores for across drugs and cell lines in long-format.",
    )
    fit_parser.add_argument(
        "--splicing_dependency_file",
        type=str,
        required=True,
        help="Tab-separated file predicted splicing dependencies across samples (columns) and exons (rows) using fitted splicing dependency models.",
    )
    fit_parser.add_argument(
        "--growth_rates_file",
        type=str,
        default=None,
        help="Tab-separated file with growth rates for each sample.",
    )
    fit_parser.add_argument(
        "--mapping_file",
        type=str,
        default=MAPPING_FILE,
        help="Tab-separated file mapping VastDB splicing event identifiers ('EVENT' column) to their corresponding ENSEMBL gene identifier ('ENSEMBL' column).",
    )
    fit_parser.add_argument(
        "--selected_models_file",
        type=str,
        default=None,
        help="List of selected cancer-driver exons as a .txt file.",
    )
    fit_parser.add_argument(
        "--output_dir",
        type=str,
        default=FITTED_DRUGASSOC_DIR,
        help="Path to the output directory.",
    )
    fit_parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of threads to use."
    )

    ## target_spotter drugassoc_predict
    pred_parser = subparser.add_parser(
        "drugassoc_predict", help="estimate drug response."
    )
    pred_parser.add_argument(
        "--splicing_dependency_file",
        type=str,
        required=True,
        help="Tab-separated file predicted splicing dependencies across samples (columns) and exons (rows) using fitted splicing dependency models.",
    )
    pred_parser.add_argument("--growth_rates_file", type=str, default=None)
    pred_parser.add_argument("--model_summaries_file", type=str, default=None)
    pred_parser.add_argument("--fitted_growth_rates_file", type=str, default=None)
    pred_parser.add_argument("--fitted_spldep_file", type=str, default=None)
    pred_parser.add_argument("--dataset", type=str, default="GDSC1")
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
            dataset=args.dataset,
            output_dir=args.output_dir,
        ).run()


##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")

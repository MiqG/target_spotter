#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------

import os
import pandas as pd
import numpy as np
import argparse
import defaults

SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}
MAPPING_FILE = defaults.MAPPING_FILE
OUTPUT_DIR = defaults.REFERENCES_DIR

"""
Development
-----------
splicing_file = os.path.join(defaults.PREP_DIR,'event_psi','CCLE-EX.tsv.gz')
genexpr_file = os.path.join(defaults.PREP_DIR,'genexpr_tpm','CCLE.tsv.gz')
mapping_file = os.path.join(defaults.REFERENCES_DIR,'mapping.tsv.gz')
output_dir = OUTPUT_DIR
"""

##### FUNCTIONS #####
def load_data(splicing_file, genexpr_file, mapping_file):
    splicing = pd.read_table(splicing_file, index_col=0)
    genexpr = pd.read_table(genexpr_file, index_col=0)
    mapping = pd.read_table(mapping_file)
    
    return splicing, genexpr, mapping


def get_summary_stats(df, label):
    summary_stats = {
        label + "_mean": df.mean(axis=1),
        label + "_median": df.median(axis=1),
        label + "_std": df.std(axis=1),
        label + "_q25": df.quantile(0.25, axis=1),
        label + "_q75": df.quantile(0.75, axis=1),
    }
    return summary_stats


def make_ccle_stats(splicing, genexpr, mapping):
    splicing_stats = pd.DataFrame(get_summary_stats(splicing, 'event'))
    genexpr_stats = pd.DataFrame(get_summary_stats(genexpr, 'gene'))
    
    ccle_stats = pd.merge(mapping, splicing_stats, how='left', left_on='EVENT', right_index=True)
    ccle_stats = pd.merge(ccle_stats, genexpr_stats, how='left', left_on='ENSEMBL', right_index=True)
    
    return ccle_stats
    
    
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--splicing_file", type=str, required=True)
    parser.add_argument("--genexpr_file", type=str, required=True)
    parser.add_argument("--mapping_file", type=str, default=MAPPING_FILE)
    parser.add_argument("--output_dir", type=str, default=OUTPUT_DIR)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    splicing_file = args.splicing_file
    genexpr_file = args.genexpr_file
    mapping_file = args.mapping_file
    output_dir = args.output_dir
    
    print("Loading data...")
    splicing, genexpr, mapping = load_data(splicing_file, genexpr_file, mapping_file)

    print("Computing CCLE statistics...")
    ccle_stats = make_ccle_stats(splicing, genexpr, mapping)

    print("Saving results...")
    ccle_stats.to_csv(os.path.join(output_dir,'ccle_stats.tsv.gz'), **SAVE_PARAMS)

    
##### SCRIPT #####
if __name__ == "__main__":
    main()
    print("Done!")
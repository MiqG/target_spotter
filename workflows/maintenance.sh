#!/bin/bash

set -o nounset
set -o errexit

# inputs
ROOT='/home/miquel/repositories/target_spotter'
GENE_DEPENDENCY_FILE=$ROOT'/data/prep/demeter2/CCLE.tsv.gz'
SPLICING_FILE=$ROOT'/data/prep/event_psi/CCLE-EX.tsv.gz'
GENEXPR_FILE=$ROOT'/data/prep/genexpr_tpm/CCLE.tsv.gz'

# make summary stats for CCLE (required to make predictions)
# python $ROOT/target_spotter/make_ccle_stats.py \
#         --splicing_file=$SPLICING_FILE \
#         --genexpr_file=$GENEXPR_FILE

# re-run model fitting
# python $ROOT/target_spotter fit \
#         --gene_dependency_file=$GENE_DEPENDENCY_FILE \
#         --splicing_file=$SPLICING_FILE \
#         --genexpr_file=$GENEXPR_FILE \
#         --n_jobs=10 \
#         --n_iterations=2

# test prediction of splicing dependencies
# python $ROOT/target_spotter predict \
#         --splicing_file=$SPLICING_FILE \
#         --genexpr_file=$GENEXPR_FILE \
#         --n_jobs=10

python $ROOT/target_spotter app

echo "Finished maintenance!"

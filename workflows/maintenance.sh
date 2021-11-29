#!/bin/bash

set -o nounset
set -o errexit

# inputs
ROOT='/home/miquel/repositories/target_spotter'
GENE_DEPENDENCY_FILE=$ROOT'/data/prep/demeter2/CCLE.tsv.gz'
SPLICING_EX_FILE=$ROOT'/data/prep/event_psi/CCLE-EX.tsv.gz'
SPLICING_ALTA_FILE=$ROOT'/data/prep/event_psi/CCLE-ALTA.tsv.gz'
SPLICING_ALTD_FILE=$ROOT'/data/prep/event_psi/CCLE-ALTD.tsv.gz'
SPLICING_INT_FILE=$ROOT'/data/prep/event_psi/CCLE-INT.tsv.gz'
SPLICING_FILES=$SPLICING_EX_FILE,$SPLICING_ALTA_FILE,$SPLICING_ALTD_FILE,$SPLICING_INT_FILE
GENEXPR_FILE=$ROOT'/data/prep/genexpr_tpm/CCLE.tsv.gz'

# make summary stats for CCLE (required to make predictions)
# python $ROOT/target_spotter/make_ccle_stats.py \
#         --splicing_file=$SPLICING_FILES \
#         --genexpr_file=$GENEXPR_FILE

# re-run model fitting
# python $ROOT/target_spotter fit \
#         --gene_dependency_file=$GENE_DEPENDENCY_FILE \
#         --splicing_file=$SPLICING_FILE \
#         --genexpr_file=$GENEXPR_FILE \
#         --n_jobs=10 \
#         --n_iterations=2

# test prediction of splicing dependencies
# python $ROOT/target_spotter spldep_predict \
#         --splicing_file=$SPLICING_FILE \
#         --genexpr_file=$GENEXPR_FILE \
#         --n_jobs=10

# python $ROOT/target_spotter app

# one-sample differential analysis
TCGA_PSI=$ROOT'/data/prep/event_psi/KICH.tsv'
TCGA_METADATA=$ROOT'/data/prep/metadata/KICH.tsv'
OUTPUT_DIR=$ROOT'/data/fitted/onesample_diff/TCGA/KICH'

# ## for each exon, create distribution of delta PSIs
# python $ROOT/target_spotter onediff_fit \
#             --data_file=$TCGA_PSI \
#             --metadata_file=$TCGA_METADATA \
#             --sample_col="sampleID" \
#             --comparison_col="sample_type" \
#             --condition_oi="Primary Tumor" \
#             --condition_ref="Solid Tissue Normal" \
#             --output_dir=$OUTPUT_DIR \
#             --n_jobs=10
            
# ## for each exon in each sample, measure the delta w.r.t. a reference population
# ## and get a p-value for this delta
# python $ROOT/target_spotter onediff_predict \
#             --data_file=$TCGA_PSI \
#             --cancer_type="KICH" \
#             --n_jobs=10

echo "Finished maintenance!"

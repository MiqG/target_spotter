import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")
PREP_DIR = os.path.join(DATA_DIR, "prep")
REFERENCES_DIR = os.path.join(DATA_DIR, "references")
EXAMPLES_DIR = os.path.join(DATA_DIR, "examples")
FITTED_DIR = os.path.join(DATA_DIR, "fitted")
FITTED_ONEDIFF_DIR = os.path.join(FITTED_DIR, "onesample_diff")
FITTED_SPLDEP_DIR = os.path.join(FITTED_DIR, "splicing_dependency")
FITTED_DRUGASSOC_DIR = os.path.join(FITTED_DIR, "drug_association")
RESULTS_DIR = "results"

# references
MAPPING_FILE = os.path.join(REFERENCES_DIR, "mapping.tsv.gz")
GENE_LENGTHS_FILE = os.path.join(REFERENCES_DIR, "gene_lengths.tsv")
CCLE_STATS_FILE = os.path.join(REFERENCES_DIR, "ccle_stats.tsv.gz")
CCLE_PCA_FILE = os.path.join(REFERENCES_DIR, "ccle_pca.pickle.gz")

# coefficients modeling gene dependency
COEFS_SPLICING_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_splicing.pickle.gz")
COEFS_GENEXPR_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_genexpr.pickle.gz")
COEFS_INTERACTION_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_interaction.pickle.gz")
COEFS_INTERCEPT_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_intercept.pickle.gz")

# default outputs drug - splicing dependency associations
MODEL_SUMMARIES_FILE = os.path.join(FITTED_DRUGASSOC_DIR, "model_summaries.tsv.gz")
FITTED_GROWTH_RATES_FILE = os.path.join(FITTED_DRUGASSOC_DIR, "growth_rates.tsv.gz")
FITTED_SPLDEP_FILE = os.path.join(FITTED_DRUGASSOC_DIR, "fitted_splicing_dependency.tsv.gz")

# app
APP_SCRIPT = os.path.join(ROOT, "target_spotter", "app.py")

# images
LOGO_FILE = os.path.join(ROOT, "images", "logo.png")

# example datasets
EXAMPLE_FILES = {
    "CCLE": {
        "splicing": os.path.join(EXAMPLES_DIR, "CCLE", "splicing_EX.tsv.gz"),
        "genexpr": os.path.join(EXAMPLES_DIR, "CCLE", "genexpr.tsv.gz"),
        "zip": os.path.join(EXAMPLES_DIR, "CCLE", "sampledata.zip"),
    }
}

# TCGA
## one-sample differential analyses
TCGA_MEDIAN_REFS_FILES = {
    "KICH": os.path.join(FITTED_ONEDIFF_DIR, "TCGA", "KICH", "median_refs.tsv.gz")
}

TCGA_POPULATION_DELTAS_FILES = {
    "KICH": os.path.join(
        FITTED_ONEDIFF_DIR, "TCGA", "KICH", "population_deltas.pickle.gz"
    )
}
##

import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")
PREP_DIR = os.path.join(DATA_DIR, "prep")
REFERENCES_DIR = os.path.join(DATA_DIR, "references")
FITTED_DIR = os.path.join(DATA_DIR, "fitted")
EXAMPLES_DIR = os.path.join(DATA_DIR, "examples")
ONESAMPLE_DIFF_DIR = os.path.join(DATA_DIR, "onesample_diff")
RESULTS_DIR = "results"

# references
MAPPING_FILE = os.path.join(REFERENCES_DIR, "mapping.tsv.gz")
GENE_LENGTHS_FILE = os.path.join(REFERENCES_DIR, "gene_lengths.tsv")
CCLE_STATS_FILE = os.path.join(REFERENCES_DIR, "ccle_stats.tsv.gz")

# coefficients
COEFS_SPLICING_FILE = os.path.join(FITTED_DIR, "coefs_event.pickle.gz")
COEFS_GENEXPR_FILE = os.path.join(FITTED_DIR, "coefs_gene.pickle.gz")
COEFS_INTERACTION_FILE = os.path.join(FITTED_DIR, "coefs_interaction.pickle.gz")
COEFS_INTERCEPT_FILE = os.path.join(FITTED_DIR, "coefs_intercept.pickle.gz")

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
    "LGG": os.path.join(ONESAMPLE_DIFF_DIR, "TCGA", "median_refs", "LGG.tsv.gz")
}

TCGA_POPULATION_DELTAS_FILES = {
    "LGG": os.path.join(
        ONESAMPLE_DIFF_DIR, "TCGA", "population_deltas", "LGG.pickle.gz"
    )
}

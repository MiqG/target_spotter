import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")
PREP_DIR = os.path.join(DATA_DIR, "prep")
REFERENCES_DIR = os.path.join(DATA_DIR, "references")
EXAMPLES_DIR = os.path.join(DATA_DIR, "examples")
FITTED_DIR = os.path.join(DATA_DIR, "fitted")
FITTED_SPLDEP_DIR = os.path.join(FITTED_DIR, "splicing_dependency")
FITTED_DRUGASSOC_DIR = os.path.join(FITTED_DIR, "drug_association")
RESULTS_DIR = "results"

# references
MAPPING_FILE = os.path.join(REFERENCES_DIR, "mapping.tsv.gz")
GENE_LENGTHS_FILE = os.path.join(REFERENCES_DIR, "gene_lengths.tsv")
CCLE_STATS_FILE = os.path.join(REFERENCES_DIR, "ccle_stats.tsv.gz")
CCLE_PCA_FILE = os.path.join(REFERENCES_DIR, "ccle_pca.pickle.gz")
INFO_DRUGS_FILE = os.path.join(REFERENCES_DIR,"info_drugs.tsv.gz")

# coefficients modeling gene dependency
COEFS_SPLICING_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_splicing.pickle.gz")
COEFS_GENEXPR_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_genexpr.pickle.gz")
COEFS_INTERCEPT_FILE = os.path.join(FITTED_SPLDEP_DIR, "coefs_intercept.pickle.gz")

# default outputs drug - splicing dependency associations
MODEL_SUMMARIES_FILES = {
    "GDSC1": os.path.join(FITTED_DRUGASSOC_DIR, "models_drug_response-GDSC1-EX", "model_summaries.tsv.gz"),
    "GDSC2": os.path.join(FITTED_DRUGASSOC_DIR, "models_drug_response-GDSC2-EX", "model_summaries.tsv.gz")
}
FITTED_GROWTH_RATES_FILES = {
    "GDSC1": os.path.join(FITTED_DRUGASSOC_DIR, "models_drug_response-GDSC1-EX", "growth_rates.tsv.gz"),
    "GDSC2": os.path.join(FITTED_DRUGASSOC_DIR, "models_drug_response-GDSC2-EX", "growth_rates.tsv.gz")
}
FITTED_SPLDEP_FILES = {
    "GDSC1": os.path.join(FITTED_DRUGASSOC_DIR, "models_drug_response-GDSC1-EX", "fitted_splicing_dependency.tsv.gz"),
    "GDSC2": os.path.join(FITTED_DRUGASSOC_DIR, "models_drug_response-GDSC2-EX", "fitted_splicing_dependency.tsv.gz")
}

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


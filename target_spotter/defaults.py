import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data")
PREP_DIR = os.path.join(DATA_DIR, "prep")
REFERENCES_DIR = os.path.join(DATA_DIR, "references")
FITTED_DIR = os.path.join(DATA_DIR, "fitted")
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

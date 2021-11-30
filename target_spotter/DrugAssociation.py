#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# 

import os
import pandas as pd
import gc
import defaults
from model_drug_screens import fit_models, get_drug_pcs

MAPPING_FILE = defaults.MAPPING_FILE
FITTED_DRUGASSOC_DIR = defaults.FITTED_DRUGASSOC_DIR
MODEL_SUMMARIES_FILE = defaults.MODEL_SUMMARIES_FILE
PCA_GROWTH_RATES_FILE = defaults.PCA_GROWTH_RATES_FILE
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}

#### FUNCTIONS ####
class DrugAssociation:
    def __init__(self, n_jobs=None):
        # parameters
        self.n_jobs = n_jobs

    def fit(self, drug_response, splicing_dependency, growth_rates=None, mapping=None):
        # prepare
        ## load default mapping
        if growth_rates is None:
            print("Inferring growth rates from Drug Response profiles ...")
            growth_rates, pca = get_drug_pcs(drug_response)
            
        if mapping is None:
            mapping = pd.read_table(MAPPING_FILE)
        mapping = mapping.loc[mapping["EVENT"].isin(splicing_dependency.index)].copy()
        
        # run linear models
        summaries = fit_models(
            drug_response, splicing_dependency, growth_rates, mapping, self.n_jobs
        )
        self.model_summaries_ = summaries
        self.pca_ = pca

    def _preprocess(self):
        # unpack attributes
        
        # standardize splicing dependency and growth rates
        
        # update attributes
        pass
    
    def predict(
        self,
        splicing_dependency,
        growth_rates=None,
        model_summaries=None
    ):
        # prepare
        ## save inputs as attributes
        self.splicing_dependency_ = splicing_dependency
        self.growth_rates_ = growth_rates
        self.model_summaries_ = model_summaries

        ## load defaults
        print("Loading defaults...")
        if self.model_summaries_ is None:
            model_summaries_file = os.path.join(FITTED_DRUGASSOC_DIR, "")
            self.model_summaries_ = pd.read_table()
        if self.growth_rates_ is None:
            # infer growth rates from fitted PCA
            infer_growth_rates()
            
        ## preprocessing inputs for prediction
        print("Preprocessing inputs...")
        self._preprocess()

        # estimate splicing dependency
        print("Estimating drug responses...")
        drug_response = estimate_drug_responses()
        self.drug_response_ = drug_response

        return self.drug_response_
    
    
class FitFromFiles:
    def __init__(
        self,
        drug_response_file,
        splicing_dependency_file,
        growth_rates_file=None,
        mapping_file=MAPPING_FILE,
        selected_models_file=None,
        output_dir=FITTED_DRUGASSOC_DIR,
        n_jobs=None,
    ):

        # inputs
        self.splicing_dependency_file = splicing_dependency_file
        self.drug_response_file = drug_response_file
        self.growth_rates_file = growth_rates_file
        self.mapping_file = mapping_file
        self.selected_models_file = selected_models_file
        
        # outputs
        self.output_dir = output_dir
        
        # parameters
        self.n_jobs = n_jobs
        
    def load_data(self):
        
        # unpack
        drug_response_file = self.drug_response_file
        splicing_dependency_file = self.splicing_dependency_file
        growth_rates_file = self.growth_rates_file
        mapping_file = self.mapping_file
        selected_models_file = self.selected_models_file
        
        # read
        splicing_dependency = pd.read_table(splicing_dependency_file, index_col=0)
        drug = pd.read_table(drug_response_file)
        mapping = pd.read_table(mapping_file)
        if growth_rates_file is not None:
            growth_rates = pd.read_table(growth_rates)
        else:
            growth_rates = None

        # drop undetected & uninformative features
        splicing_dependency = splicing_dependency.dropna(thresh=2)
        splicing_dependency = splicing_dependency.loc[splicing_dependency.std(axis=1) != 0]

        # subset
        ## events
        if selected_models_file is not None:
            selected_models = list(
                pd.read_table(selected_models_file, header=None)[0].values
            )
            idx = splicing_dependency.index.isin(selected_models)
            splicing_dependency = splicing_dependency.loc[idx].copy()

        ## samples
        common_samples = set(splicing_dependency.columns).intersection(drug["ARXSPAN_ID"])
        splicing_dependency = splicing_dependency.loc[:, common_samples].copy()
        drug = drug.loc[drug["ARXSPAN_ID"].isin(common_samples)].copy()
        mapping = mapping.loc[mapping["EVENT"].isin(splicing_dependency.index)].copy()

        gc.collect()
        
        # pack
        self.drug_response_ = drug        
        self.splicing_dependency_ = splicing_dependency
        self.growth_rates_ = growth_rates
        self.mapping_ = mapping

    def save(self, estimator):
        model_summaries = estimator.model_summaries_
        pca = estimator.pca_

        os.makedirs(self.output_dir, exist_ok=True)

        model_summaries.to_csv(
            os.path.join(self.output_dir, "model_summaries.tsv.gz"), **SAVE_PARAMS
        )
        pd.to_pickle(pca, os.path.join(self.output_dir, "pca_growth_rates.pickle.gz"))

    def run(self):
        print("Loading data...")
        self.load_data()
        
        print("Fitting models...")
        estimator = DrugAssociation(n_jobs=self.n_jobs)
        estimator.fit(self.drug_response_, self.splicing_dependency_, 
                      self.growth_rates_, self.mapping_)

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)
        
        
class PredictFromFiles:
    def __init__(
        self,
        n_jobs=None,
    ):

        # inputs

        # outputs
        self.output_dir = output_dir

        # parameters
        self.n_jobs = n_jobs

    def load_data(self):
        # read
        
        # subset samples
        
        # take default files

        # load model summaries
    
        # update attributes
        self.coefs_intercept_ = coefs_intercept
        
    def save(self, estimator):
        splicing_dependency = estimator.splicing_dependency_

        os.makedirs(self.output_dir, exist_ok=True)

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Estimating splicing dependencies...")
        estimator = DrugAssociation(
            n_jobs=self.n_jobs,
        )
        _ = estimator.predict(
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)
        

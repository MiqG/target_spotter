#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
#

import os
import pandas as pd
import numpy as np
import gc
import defaults
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from model_drug_screens import fit_models, get_drug_pcs, infer_growth_rates

MAPPING_FILE = defaults.MAPPING_FILE
FITTED_DRUGASSOC_DIR = defaults.FITTED_DRUGASSOC_DIR
MODEL_SUMMARIES_FILE = defaults.MODEL_SUMMARIES_FILE
FITTED_GROWTH_RATES_FILE = defaults.FITTED_GROWTH_RATES_FILE
FITTED_SPLDEP_FILE = defaults.FITTED_SPLDEP_FILE
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
        self.growth_rates_ = growth_rates

    def _preprocess(self):
        # unpack attributes
        splicing_dependency = self.splicing_dependency_
        growth_rates = self.growth_rates_
        model_summaries = self.model_summaries_

        # subset
        ## common samples
        common_samples = set(splicing_dependency.columns).intersection(
            growth_rates.index
        )

        ## common events
        common_events = set(model_summaries["EVENT"]).intersection(
            splicing_dependency.index
        )

        ## subset
        model_summaries = model_summaries.loc[
            model_summaries["EVENT"].isin(common_events)
        ]
        spldep = splicing_dependency.loc[model_summaries["EVENT"], common_samples]
        gr = growth_rates.loc[common_samples].T
        gr = gr.loc[gr.index.repeat(len(spldep))]

        # standardize
        ## splicing dependency
        spldep_mean = model_summaries["spldep_mean"].values.reshape(-1, 1)
        spldep_std = model_summaries["spldep_std"].values.reshape(-1, 1)
        spldep = (spldep - spldep_mean) / spldep_std

        ## growth rates
        gr_mean = model_summaries["growth_mean"].values.reshape(-1, 1)
        gr_std = model_summaries["growth_std"].values.reshape(-1, 1)
        gr = (gr - gr_mean) / gr_std

        # update attributes
        self.prep_splicing_dependency_ = spldep
        self.prep_growth_rates_ = gr
        self.prep_model_summaries_ = model_summaries

    def _estimate_drug_response(
        self, prep_splicing_dependency, prep_growth_rates, prep_model_summaries
    ):
        # prepare arrays
        spldep = prep_splicing_dependency.values
        gr = prep_growth_rates.values
        coef_spldep = prep_model_summaries["spldep_coefficient"].values.reshape(-1, 1)
        coef_gr = prep_model_summaries["growth_coefficient"].values.reshape(-1, 1)
        coef_intercept = prep_model_summaries["intercept_coefficient"].values.reshape(
            -1, 1
        )

        # log(IC50) ~ SplicingDependency + GrowthRate + Intercept
        ## inferred drug response per drug and event
        full_est = prep_model_summaries[["DRUG_ID", "EVENT", "ENSEMBL", "GENE"]].copy()
        full_est[prep_splicing_dependency.columns] = np.nan
        full_est[prep_splicing_dependency.columns] = (
            coef_spldep * spldep + coef_gr * gr + coef_intercept
        )
        ## estimated drug response per drug
        samples = prep_splicing_dependency.columns
        drug_ests = []
        for drug_id, drug_df in full_est.groupby(["DRUG_ID"]):
            # prepare weights for weighted average
            # higher pearson, higher contribution
            weights = prep_model_summaries.loc[prep_model_summaries["DRUG_ID"]==drug_id]
            weights = np.clip(weights.set_index("EVENT")["pearson_correlation"],0,1)
            weights = weights / weights.sum()
            
            # get estimations by all event
            m = drug_df.set_index("EVENT").loc[weights.index,samples].fillna(0).values
            w = weights.fillna(0).values.reshape(-1,1)
            ## (n. events x n.samples)^T dotprod. (n. events x 1) = (n.samples x 1)
            drug_est = np.dot(m.T,w)
            drug_est = pd.DataFrame({"DRUG_ID": drug_id, "sample": samples, "predicted_ic50": drug_est.ravel()})
            drug_ests.append(drug_est)
        drug_ests = pd.concat(drug_ests)
        
        return drug_ests, full_est

    def predict(
        self,
        splicing_dependency,
        growth_rates=None,
        model_summaries=None,
        fitted_growth_rates=None,
        fitted_spldep=None,
    ):
        # prepare
        ## save inputs as attributes
        self.splicing_dependency_ = splicing_dependency
        self.growth_rates_ = growth_rates
        self.model_summaries_ = model_summaries
        self.fitted_growth_rates_ = fitted_growth_rates
        self.fitted_spldep_ = fitted_spldep

        ## load defaults
        print("Loading defaults...")
        if self.model_summaries_ is None:
            self.model_summaries_ = pd.read_table(MODEL_SUMMARIES_FILE)
        if self.growth_rates_ is None:
            if self.fitted_growth_rates_ is None:
                self.fitted_growth_rates_ = pd.read_table(
                    FITTED_GROWTH_RATES_FILE, index_col=0
                )
            if self.fitted_spldep_ is None:
                self.fitted_spldep_ = pd.read_table(FITTED_SPLDEP_FILE, index_col=0)

            self.growth_rates_ = infer_growth_rates(
                self.splicing_dependency_,
                self.fitted_growth_rates_,
                self.fitted_spldep_,
            )

        ## preprocessing inputs for prediction
        print("Preprocessing inputs...")
        self._preprocess()

        # estimate splicing dependency
        print("Estimating drug responses...")
        drug_estimates, full_estimates = self._estimate_drug_response(
            self.prep_splicing_dependency_,
            self.prep_growth_rates_,
            self.prep_model_summaries_,
        )
        self.drug_estimates_ = drug_estimates
        self.full_estimates_ = full_estimates

        return self.drug_estimates_, self.full_estimates_


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
        splicing_dependency = splicing_dependency.loc[
            splicing_dependency.std(axis=1) != 0
        ]

        # subset
        ## events
        if selected_models_file is not None:
            selected_models = list(
                pd.read_table(selected_models_file, header=None)[0].values
            )
            idx = splicing_dependency.index.isin(selected_models)
            splicing_dependency = splicing_dependency.loc[idx].copy()

        ## samples
        common_samples = set(splicing_dependency.columns).intersection(
            drug["ARXSPAN_ID"]
        )
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
        prep_splicing_dependency = self.splicing_dependency_
        model_summaries = estimator.model_summaries_
        growth_rates = estimator.growth_rates_

        os.makedirs(self.output_dir, exist_ok=True)

        prep_splicing_dependency.reset_index().to_csv(
            os.path.join(self.output_dir, "fitted_splicing_dependency.tsv.gz"),
            **SAVE_PARAMS
        )
        model_summaries.to_csv(
            os.path.join(self.output_dir, "model_summaries.tsv.gz"), **SAVE_PARAMS
        )
        growth_rates.reset_index().to_csv(
            os.path.join(self.output_dir, "growth_rates.tsv.gz"), **SAVE_PARAMS
        )

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Fitting models...")
        estimator = DrugAssociation(n_jobs=self.n_jobs)
        estimator.fit(
            self.drug_response_,
            self.splicing_dependency_,
            self.growth_rates_,
            self.mapping_,
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)


class PredictFromFiles:
    def __init__(
        self,
        splicing_dependency_file,
        growth_rates_file=None,
        model_summaries_file=None,
        fitted_growth_rates_file=None,
        fitted_spldep_file=None,
        output_dir="drug_association",
    ):
        # inputs
        self.splicing_dependency_file = splicing_dependency_file
        self.growth_rates_file = growth_rates_file
        self.model_summaries_file = model_summaries_file
        self.fitted_growth_rates_file = fitted_growth_rates_file
        self.fitted_spldep_file = fitted_spldep_file

        # outputs
        self.output_dir = output_dir

    def load_data(self):
        # prep defaults
        if self.model_summaries_file is None:
            self.model_summaries_file = MODEL_SUMMARIES_FILE
        if self.fitted_growth_rates_file is None:
            self.fitted_growth_rates_file = FITTED_GROWTH_RATES_FILE
        if self.fitted_spldep_file is None:
            self.fitted_spldep_file = FITTED_SPLDEP_FILE

        # read
        splicing_dependency = pd.read_table(self.splicing_dependency_file, index_col=0)
        model_summaries = pd.read_table(self.model_summaries_file)
        if self.growth_rates_file is None:
            self.fitted_growth_rates_ = pd.read_table(
                FITTED_GROWTH_RATES_FILE, index_col=0
            )
            self.fitted_spldep_ = pd.read_table(FITTED_SPLDEP_FILE, index_col=0)
            growth_rates = infer_growth_rates(
                splicing_dependency, self.fitted_growth_rates_, self.fitted_spldep_
            )

        # subset
        common_events = set(splicing_dependency.index).intersection(
            model_summaries["EVENT"].values
        )
        common_samples = set(splicing_dependency.columns).intersection(
            growth_rates.index
        )

        # update attributes
        self.splicing_dependency_ = splicing_dependency.loc[common_events]
        self.model_summaries_ = model_summaries.loc[
            model_summaries["EVENT"].isin(common_events)
        ]
        self.growth_rates_ = growth_rates.loc[common_samples]

    def save(self, estimator):
        drug_estimates = estimator.drug_estimates_
        full_estimates = estimator.full_estimates_

        os.makedirs(self.output_dir, exist_ok=True)
        drug_estimates.to_csv(
            os.path.join(self.output_dir, "estimated_drug_response_by_drug.tsv.gz"),
            **SAVE_PARAMS
        )
        full_estimates.to_csv(
            os.path.join(self.output_dir, "estimated_drug_response_by_drug_and_event.tsv.gz"),
            **SAVE_PARAMS
        )

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Estimating splicing dependencies...")
        estimator = DrugAssociation()
        _ = estimator.predict(
            self.splicing_dependency_, self.growth_rates_, self.model_summaries_
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)

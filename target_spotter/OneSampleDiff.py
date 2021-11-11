#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Perform one-sample differential analyses based on population-level measurements.

import pandas as pd
from pandas.io.pickle import to_pickle
import numpy as np
import os
import gc
from statsmodels.distributions.empirical_distribution import ECDF
import defaults
from joblib import Parallel, delayed
from collections import ChainMap

TCGA_MEDIAN_REFS_FILES = defaults.TCGA_MEDIAN_REFS_FILES
TCGA_POPULATION_DELTAS_FILES = defaults.TCGA_POPULATION_DELTAS_FILES
FITTED_ONEDIFF_DIR = defaults.FITTED_ONEDIFF_DIR
SAVE_PARAMS = {"sep": "\t", "compression": "gzip", "index": False}


class OneSampleDiff:
    def __init__(self, n_jobs=None):
        """
        samples_oi : list of TCGA sample identifiers (independently of cancer_type).
        
        When selecting cancer type, we go to the PANCAN metadata and extract the list of samples and we continue from there.
        """
        # parameters
        self.n_jobs = n_jobs

    def get_population_deltas(self, x_oi, x_ref):
        """
        For every row in data, we return the vector of differences between x_oi and the median(x_ref).
        """
        x_ref = x_ref[~np.isnan(x_ref)]
        nobs = len(x_ref)
        if nobs > 0:
            median_ref = np.median(x_ref)
        else:
            median_ref = np.nan

        deltas = x_oi.values - median_ref
        out = {
            "index": x_oi.name,
            "median_ref": median_ref,
            "nobs": nobs,
            "deltas": deltas,
        }
        return out

    def fit(self, data, samples_oi, samples_ref):
        """
        Generates distributions of population deltas of every sample in condition_oi w.r.t. the median condition_ref.
        """
        results = Parallel(n_jobs=self.n_jobs)(
            delayed(self.get_population_deltas)(x[samples_oi], x[samples_ref])
            for idx, x in data.iterrows()
        )

        # prepare outputs
        population_deltas = []
        median_refs = []
        for res in results:
            index = res["index"]
            median_ref = res["median_ref"]
            nobs = res["nobs"]
            deltas = res["deltas"]
            median_refs.append({"index": index, "median_ref": median_ref, "nobs": nobs})
            population_deltas.append({index: deltas})

            del index, median_ref, deltas

        median_refs = pd.DataFrame(median_refs).set_index("index")
        population_deltas = dict(ChainMap(*population_deltas))

        # update attributes
        self.median_refs_ = median_refs
        self.population_deltas_ = population_deltas

    def _preprocess(self):
        # unpack
        data = self.data_
        median_refs = self.median_refs_
        population_deltas = self.population_deltas_

        # subset
        ## common index
        common_idx = set(data.index).intersection(median_refs.index)

        data = data.loc[common_idx]
        median_refs = median_refs.loc[common_idx]
        population_deltas = {idx: population_deltas[idx] for idx in common_idx}

        # update attributes
        self.data_ = data
        self.median_refs_ = median_refs
        self.population_deltas_ = population_deltas

    def _compute_onesample_pvalue(self, sample_delta, population_deltas):
        """
        as in https://github.com/comprna/SUPPA/blob/c41c6c98f962008181ac419b8e0daf5f8dd5a914/lib/diff_tools.py#L266
        """

        abs_population_deltas = np.abs(population_deltas)
        ecdf = ECDF(abs_population_deltas)
        # one-tailed test
        sample_pvalue = (1.0 - ecdf(np.abs(sample_delta))) * 0.5

        return sample_pvalue

    def predict(self, data, median_refs, population_deltas):
        """
        For every sample (column) in data, test whether every row is differentially distributed w.r.t. a reference population.
        """
        # save inputs as attributes
        self.data_ = data
        self.median_refs_ = median_refs
        self.population_deltas_ = population_deltas

        # preprocess inputs for prediction
        print("Preprocessing inputs...")
        self._preprocess()

        # compute one-sample differential analyses
        median_refs = self.median_refs_["median_ref"].values.reshape(-1, 1)
        delta_data = self.data_ - median_refs
        delta_pvalues = delta_data.copy()
        delta_pvalues.values[:, :] = np.nan
        for col in delta_data.columns:
            sample_deltas = delta_data[col]
            for idx in sample_deltas.index:
                sample_delta = sample_deltas[idx]
                population_deltas = self.population_deltas_[idx]
                pvalue = self._compute_onesample_pvalue(sample_delta, population_deltas)
                delta_pvalues.loc[idx, col] = pvalue

        # store outputs as attributes
        self.delta_data_ = delta_data
        self.delta_pvalues_ = delta_pvalues


class FitFromFiles:
    def __init__(
        self,
        data_file,
        metadata_file,
        sample_col,
        comparison_col,
        condition_oi,
        condition_ref,
        output_dir=FITTED_ONEDIFF_DIR,
        n_jobs=None,
    ):
        # inputs
        self.data_file = data_file
        self.metadata_file = metadata_file

        # outputs
        self.output_dir = output_dir

        # parameters
        self.sample_col = sample_col
        self.comparison_col = comparison_col
        self.condition_oi = condition_oi
        self.condition_ref = condition_ref
        self.n_jobs = n_jobs

    def load_data(self):
        # read
        data = pd.read_table(self.data_file, index_col=0).iloc[:500]
        print("DEV")
        metadata = pd.read_table(self.metadata_file)

        # subset
        common_samples = set(metadata[self.sample_col]).intersection(data.columns)

        # prepare samples
        idx = metadata[self.comparison_col] == self.condition_oi
        samples_oi = metadata.loc[idx, self.sample_col].values

        idx = metadata[self.comparison_col] == self.condition_ref
        samples_ref = metadata.loc[idx, self.sample_col].values

        # update attributes
        self.data_ = data[common_samples]
        self.metadata_ = metadata.loc[metadata[self.sample_col].isin(common_samples)]
        self.samples_oi_ = samples_oi
        self.samples_ref_ = samples_ref

        gc.collect()

    def save(self, estimator):
        median_refs = estimator.median_refs_
        population_deltas = estimator.population_deltas_

        os.makedirs(self.output_dir, exist_ok=True)

        median_refs.reset_index().to_csv(
            os.path.join(self.output_dir, "median_refs.tsv.gz"), **SAVE_PARAMS
        )
        to_pickle(
            population_deltas,
            os.path.join(self.output_dir, "population_deltas.pickle.gz"),
        )

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Generating empirical population of deltas...")
        estimator = OneSampleDiff(n_jobs=self.n_jobs)
        estimator.fit(
            data=self.data_, samples_oi=self.samples_oi_, samples_ref=self.samples_ref_
        )

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)


class PredictFromFiles:
    def __init__(
        self,
        data_file,
        median_refs_file=None,
        population_deltas_file=None,
        output_dir="one_sample_diff_analysis",
        cancer_type=None,
        n_jobs=None,
    ):
        # inputs
        self.data_file = data_file
        self.median_refs_file = median_refs_file
        self.population_deltas_file = population_deltas_file

        # outputs
        self.output_dir = output_dir

        # parameters
        self.cancer_type = cancer_type
        self.n_jobs = n_jobs

        # use defaults if cancer type is defined
        if self.cancer_type is not None:
            print("Using default fitted files for %s" % self.cancer_type)
            self.median_refs_file = TCGA_MEDIAN_REFS_FILES[self.cancer_type]
            self.population_deltas_file = TCGA_POPULATION_DELTAS_FILES[self.cancer_type]
        else:
            assert (self.median_refs_file is not None) & (
                self.population_deltas_file is not None
            )

    def load_data(self):
        # read
        data = pd.read_table(self.data_file, index_col=0)
        median_refs = pd.read_table(self.median_refs_file, index_col=0)
        population_deltas = pd.read_pickle(self.population_deltas_file)

        # update
        self.data_ = data
        self.median_refs_ = median_refs
        self.population_deltas_ = population_deltas

    def save(self, estimator):
        delta_data = estimator.delta_data_
        delta_pvalues = estimator.delta_pvalues_

        os.makedirs(self.output_dir, exist_ok=True)

        delta_data.reset_index().to_csv(
            os.path.join(self.output_dir, "deltas.tsv.gz"), **SAVE_PARAMS
        )
        delta_pvalues.reset_index().to_csv(
            os.path.join(self.output_dir, "pvalues.tsv.gz"), **SAVE_PARAMS
        )

    def run(self):
        print("Loading data...")
        self.load_data()

        print("Performing one-sample differential analysis...")
        estimator = OneSampleDiff(n_jobs=self.n_jobs)
        _ = estimator.predict(self.data_, self.median_refs_, self.population_deltas_)

        print("Saving results to %s ..." % self.output_dir)
        self.save(estimator)

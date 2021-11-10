#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Perform one-sample differential analyses based on population-level measurements.

import pandas as pd
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import defaults
from joblib import Parallel, delayed
from collections import ChainMap

"""
Development
-----------

"""
TCGA_MEDIAN_REFS_FILES = default.TCGA_MEDIAN_REFS_FILES
TCGA_POPULATION_DELTAS_FILES = default.TCGA_POPULATION_DELTAS_FILES
    
    
class OneSampleDiff:
    def __init__(self, cancer_type=None, samples_oi=None, n_jobs=None):
        """
        samples_oi : list of TCGA sample identifiers (independently of cancer_type).
        
        When selecting cancer type, we go to the PANCAN metadata and extract the list of samples and we continue from there.
        """
        # parameters
        self.n_jobs = n_jobs
    
    def get_population_deltas(x_oi, x_ref):
        """
        For every row in data, we return the vector of differences between x_oi and the median(x_ref).
        """
        median_ref = np.median(x_ref)
        deltas = x_oi - median_ref
        out = {'index': x_oi.name, 'median_ref': median_ref, 'deltas': deltas}
        return out
    
    def fit(self, data, samples_oi, samples_ref):
        """
        Generates distributions of population deltas of every sample in condition_oi w.r.t. the median condition_ref.
        """
        results = Parallel(n_jobs=self.n_jobs)(
            delayed(get_population_deltas)(
                x[samples_oi],
                x[samples_ref]
            )
            for idx, x in data.iterrows()
        )
        
        # prepare outputs
        population_deltas = []
        median_refs = []
        for index, median_ref, deltas in results.values():
            median_refs.append({'index': index, 'median_ref': median_ref})
            population_deltas.append({index: deltas})
        
        median_refs = pd.DataFrame(median_refs).set_index('index')
        
        # store
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
        
        # update attributes
        self.prep_data_ = data.loc[common_idx]
        self.median_refs_ = median_refs.loc[common_idx]
        self.population_deltas_ = {idx: population_deltas[idx] for idx in common_idx}
        
        
    def _compute_onesample_pvalue(self, sample_delta, population_deltas):
        """
        as in https://github.com/comprna/SUPPA/blob/c41c6c98f962008181ac419b8e0daf5f8dd5a914/lib/diff_tools.py#L266
        """

        abs_population_deltas = np.abs(population_deltas)
        ecdf = ECDF(abs_population_deltas)
        # one-tailed test
        sample_pvalue = (1.0 - ecdf(np.abs(sample_delta))) * 0.5

        return sample_pvalue
        
        
    def predict(self, data, cancer_type, median_refs=None, population_deltas=None):
        """
        For every sample (column) in data, test whether every row is differentially distributed w.r.t. a reference population.
        """
        # save inputs as attributes
        self.data_ = data
        self.cancer_type = cancer_type
        self.median_refs_ = median_refs
        self.population_deltas_ = population_deltas
        
        # load defaults
        if median_refs is None:
            self.median_refs_ = pd.read_table(TCGA_MEDIAN_REFS_FILES[cancer_type])
        if population_deltas is None:
            self.population_deltas_ = pd.read_pickle(TCGA_POPULATION_DELTAS_FILES[cancer_type])
        
        # preprocess inputs for prediction
        print("Preprocessing inputs...")
        self._preprocess()
        
        # compute one-sample differential analyses
        delta_data = self.data_ - self.median_refs_[self.data_.index].values.reshape(-1,1)
        delta_pvalues = delta_data.copy()
        delta_pvalues.values[:,:] = np.nan
        for col in delta_data.columns:
            sample_deltas = delta_data[col]
            for idx in sample_deltas.index:
                sample_delta = sample_deltas[idx]
                population_deltas = self.population_deltas_[idx]
                pvalue = self._compute_onesample_pvalue(sample_delta, population_deltas)
                delta_pvalues.loc[idx,col] = pvalue
        
        # save outputs
        self.delta_data_ = delta_data
        self.delta_pvalues_ = delta_pvalues
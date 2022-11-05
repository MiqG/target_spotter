<p align="center">
  <img src="images/logo.png" width="200">
</p>

Systematic prioritization of splicing targets to treat cancer.

[![pipy](https://img.shields.io/pypi/v/target_spotter?color=informational)](https://pypi.python.org/pypi/target_spotter)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Description
This package contains 3 modules:
- `SplicingDependency`
- `DrugAssociations`
- `PersonalizeTreatment`
- `examples`

## Requirements
### basic usage
In brackets, versions of packages used to develop `target_spotter`.
- `numpy` (1.19.2)
- `pandas` (1.1.2)
- `scipy` (1.6.2)
- `scikit-learn` (0.23.2)
- `joblib` (1.0.1)
- `tqdm` (4.59.0)
- `statsmodels`

### interactive usage
- `streamlit`

### maintenance
- `limix` (install from [here](https://github.com/limix/limix)) with:
```shell
pip install -e git+git://github.com/limix/limix.git@master#egg=limix
```


## Installation
### pip
```shell
pip install target_spotter
```
### local
```shell
git clone https://github.com/MiqG/target_spotter.git
cd target_spotter
pip install -e .
```

## Usage
### from the shell
```shell
# predict splicing dependency
target_spotter spldp_predict \
            --splicing_file='splicing.tsv.gz' \
            --genexpr_file='genexpr.tsv.gz' \
            --output_dir='results' \
            --n_jobs=4

# one-sample differential splicing analysis
target_spotter onediff_predict \
            --data_file='splicing.tsv.gz' \
            --cancer_type='KICH' \
            --output_dir='results' \
            --n_jobs=4
```
### within a script
#### Inference of splicing dependencies
```python
import target_spotter as ts
from examples import make_sampledata

splicing, genexpr = make_sampledata(n_samples=10)

estimator = ts.SplicingDependency() # already fitted
spldep_medians = estimator.predict(splicing, genexpr)
spldep_medians
```
#### One-sample differential splicing analysis
```python
import target_spotter as ts
from examples import make_sampledata

splicing, genexpr = make_sampledata(n_samples=10)

estimator = ts.OneSampleDiff(cancer_type='KIRC') # already fitted
# get differentially spliced exons with respect to solid tissue normal samples
result = estimator.predict(splicing, ref_condition='stn')

```

### through the app
```shell
# start the app and follow instructions
target_spotter app
```

## Maintenance
```shell
# fit splicing dependency models
target_spotter fit --rnai_file='rnai.tsv.gz' \
                   --splicing_file='splicing.tsv.gz' \
                   --genexpr_file='genexpr.tsv.gz' \
                   --output_dir='models'
```


## Tutorials
- [Find the splicing dependencies in your sample](https://github.com/CRG-CNAG/target_spotter/blob/main/tutorials/basics.ipynb)

## Contact
This project has been fully developed at the [Centre for Genomic Regulation](https://www.crg.eu/) within the group of [Design of Biological Systems](https://www.crg.eu/en/luis_serrano).

Please, report any issues that you experience through this repository's ["Issues"](https://github.com/CRG-CNAG/target_spotter/issues) or email:
- [Miquel Anglada-Girotto](mailto:miquel.anglada@crg.eu)
- [Luis Serrano](mailto:luis.serrano@crg.eu)

## License

`target_spotter` is distributed under a BSD 3-Clause License (see [LICENSE](https://github.com/CRG-CNAG/target_spotter/blob/main/LICENSE)).

## References
- *Himberg, J., & Hyvarinen, A.* "Icasso: software for investigating the reliability of ICA estimates by clustering and visualization". IEEE XIII Workshop on Neural Networks for Signal Processing (2003). DOI: https://doi.org/10.1109/NNSP.2003.1318025
- *Sastry, Anand V., et al.* "The Escherichia coli transcriptome mostly consists of independently regulated modules." Nature communications 10.1 (2019): 1-14. DOI: https://doi.org/10.1038/s41467-019-13483-w
- *Kairov, U., Cantini, L., Greco, A. et al.* Determining the optimal number of independent components for reproducible transcriptomic data analysis. BMC Genomics 18, 712 (2017). DOI: https://doi.org/10.1186/s12864-017-4112-9

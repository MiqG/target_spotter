<p align="center">
  <img src="images/logo.png" width="200">
</p>

Systematic prioritization of splicing targets to treat cancer.

[![pipy](https://img.shields.io/pypi/v/target_spotter?color=informational)](https://pypi.python.org/pypi/target_spotter)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Due to name availability on `PyPI`, the package name had to be `target_spotter` rather than just *spotter*.

## Description
This package contains 2 main modules:
- `SplicingDependency`
    
    Defines the class to fit and predict splicing dependencies from exon inclusion (PSI) and gene expression (TPM)

- `DrugAssociation`
    
    Defines the class to fit and predict drug sensitivities from splicing dependencies inferred with `SplicingDependency`.
    
## Requirements
In brackets, versions of packages used to develop `target_spotter`.
- `numpy` (1.19.2)
- `pandas` (1.1.2)
- `scipy` (1.6.2)
- `scikit-learn` (0.23.2)
- `joblib` (1.0.1)
- `tqdm` (4.59.0)
- `statsmodels` (0.13.1)
- `glimix-core` (3.1.12)
- `numpy-sugar` (1.5.3)

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
Here's a tutorial on how to [compute splicing dependencies and predict drug sensitivities](https://github.com/CRG-CNAG/target_spotter/blob/main/tutorials/basics.ipynb) from transcriptomic data as input.

### inputs obtained running [`vast-tools`](https://github.com/vastgroup/vast-tools)
To run our predictions you'll need:
- exon inclusion table in *percentage spliced in* (PSI); first column has exon identifiers.
```
EVENT  ACH-000415  ACH-000894  ACH-000422  ACH-000358  ACH-000468                                
HsaEX0067681        1.77        1.94        1.18        7.30        0.98   
HsaEX6078702      100.00      100.00      100.00      100.00      100.00   
HsaEX0056692        5.88       38.53       24.35        9.09       16.96   
HsaEX0056690       88.89       93.85       83.78       87.65       91.11   
HsaEX0056691      100.00      100.00       96.72      100.00       95.24
```


- gene expression table in either *transcripts per million* (TPM) or *raw counts* (Count); first column has gene identifiers.
```
ID  ACH-000415  ACH-000894  ACH-000422  ACH-000358  ACH-000468
ENSG00000000003    4.389567    7.281791    5.064366    6.165912    3.939227   
ENSG00000000005    0.000000    0.014355    0.555816    0.000000    0.000000   
ENSG00000000419    5.147714    6.803744    6.841596    5.928607    6.625417   
ENSG00000000457    1.000000    2.469886    2.931683    2.726831    1.963474   
ENSG00000000460    1.555816    3.811471    3.834913    4.347666    3.228049
```


ideally, we recommend using the outputs from running [`vast-tools`](https://github.com/vastgroup/vast-tools) on your RNA-seq data. Specifically, we run:

```
# to obtain PSI and TPM for each sample
vast-tools align \
    fw.fastq.gz \
    rv.fastq.gz \
    --sp Hs2 \
    --dbDir $DBDIR \
    --expr \
    --EEJ_counts \
    --cores $THREADS \
    --output $OUTPUTDIR


# to combine PSI and TPM into single tables
vast-tools combine \
    --cores $THREADS \
    --sp Hs2 \
    --dbDir $DBDIR \
    --keep_raw_reads \
    --keep_raw_incl \
    --output $OUTPUTDIR \
    --TPM

# to cleanup uncertain PSI quantification
vast-tools tidy \
    $COMBINED_INCLUSION_TABLE \
    -min_N 1 \
    -min_SD 0 \
    --min_ALT_use 25 \
    --noVLOW \
    --log \
    -outFile $PSI_TIDY_OUTPUT
```

### interactively
Check out our [colab notebook](https://drive.google.com/file/d/11NtO1tH5JKUvOKm7M88jUTy6nwpEG1_2/view?usp=share_link). You'll be able to run the full pipeline (i.e. obtain splicing dependencies, max. harm scores and drug sensitivities) for your samples.

### from the shell
```
# predict splicing dependency
python ./target_spotter spldep_predict \
        --splicing_file=$PSI_FILE \
        --genexpr_file=$GENEXPR_FILE \
        --output_dir="results_spldep" \
        --n_jobs=4 \
        --log_transform
        
# predict drug sensitivity
## GDSC1
python ./target_spotter drugassoc_predict \
        --splicing_dependency_file="results_spldep/mean.tsv.gz" \
        --dataset="GDSC1" \
        --output_dir="drug_response_GDSC1"
## GDSC2                    
python ./target_spotter drugassoc_predict \
        --splicing_dependency_file="results_spldep/mean.tsv.gz" \
        --dataset="GDSC2" \
        --output_dir="drug_response_GDSC2"
```

### within a script
```python
import pandas as pd
from target_spotter import SplicingDependency, DrugAssociation, defaults

# load data
splicing_file = defaults.EXAMPLE_FILES["CCLE"]["splicing"]
genexpr_file = defaults.EXAMPLE_FILES["CCLE"]["genexpr"]
splicing = pd.read_table(splicing_file).set_index("EVENT")
genexpr = pd.read_table(genexpr_file).set_index("ID")

# predict splicing dependencies
estimator = SplicingDependency(log_transform=True) # already fitted
spldep_means, max_harm_score_means = estimator.predict(splicing, genexpr)

# predict drug sensitivities
datasets = ["GDSC1","GDSC2"]
ic50_by_drugs = []
ic50_by_exons = []
for dataset in datasets:
    print(dataset)
    estimator = DrugAssociation() # already fitted
    ic50_by_drug, ic50_by_exon = estimator.predict(spldep_means, dataset=dataset)
    ic50_by_drugs.append(ic50_by_drug)
    ic50_by_exons.append(ic50_by_exon)

ic50_by_drugs = pd.concat(ic50_by_drugs)
ic50_by_exons = pd.concat(ic50_by_exons)
```
### outputs
- predicted splicing dependencies
    ```python
    print(spldep_means.iloc[:5,:5])
    ```
    ```
                  ACH-000861  ACH-000889  ACH-000036  ACH-000174  ACH-000808
    HsaEX6065058         NaN         NaN         NaN         NaN         NaN
    HsaEX6065028         NaN         NaN         NaN         NaN         NaN
    HsaEX6008208         NaN         NaN         NaN         NaN         NaN
    HsaEX1001338         NaN         NaN         NaN         NaN         NaN
    HsaEX0001886         NaN         NaN         NaN         NaN         NaN
    ```


- predicted maximum harm scores
    ```python
    print(max_harm_score_means.iloc[:5,:5])
    ```
    ```
                  ACH-000861  ACH-000036  ACH-000174  ACH-000889  ACH-000808
    HsaEX6082977  -93.496257  -93.950451  -93.425121  -95.945873  -98.724904
    HsaEX0045309   -2.716957         NaN         NaN         NaN         NaN
    HsaEX1038630         NaN         NaN         NaN         NaN         NaN
    HsaEX6002514  -96.395778 -100.157193  -84.599117  -91.125346 -100.202869
    HsaEX6066904         NaN         NaN         NaN         NaN         NaN
    ```


- predicted drug sensitivities
    ```python
    print(ic50_by_drugs.head())
    ```
    ```
      dataset           ID      sample  predicted_ic50
    0   GDSC1  1001_2000.0  ACH-000861        3.538367
    1   GDSC1  1001_2000.0  ACH-000036        3.539620
    2   GDSC1  1001_2000.0  ACH-000174        3.564975
    3   GDSC1  1001_2000.0  ACH-000889        3.581721
    4   GDSC1  1001_2000.0  ACH-000808        3.732448
    ```

## Contact
This project has been fully developed at the [Centre for Genomic Regulation](https://www.crg.eu/) within the group of [Design of Biological Systems](https://www.crg.eu/en/luis_serrano).

Please, report any issues that you experience through this repository's ["Issues"](https://github.com/MiqG/target_spotter/issues) or email:
- [Miquel Anglada-Girotto](mailto:miquel.anglada@crg.eu)
- [Luis Serrano](mailto:luis.serrano@crg.eu)

## License
`target_spotter` is distributed under a GNU License (see [LICENSE](https://github.com/MiqG/target_spotter/blob/main/LICENSE)).

## Citation
Manuscript under preparation.

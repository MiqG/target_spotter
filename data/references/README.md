# Data references

Essential files for our computational framework.

- `ccle_stats.tsv.gz`: exon-level summary statistics of the CCLE transcriptomic dataset at the splicing event and corresponding gene level. Splicing event PSI and gene TPM were quantified using `vast-tools`. We use it to standardize new input data to make new predictions with respect to our training data.
- `gene_lengths.tsv`: we obtained effective gene lengths from file “Hs2_mRNA-50-SS.eff” in the Hs2 (hg38) human annotation from VastDB. This is only used in case inputs are given as raw RNA seq counts.
- `info_drugs.tsv.gz`: drug metadata from GDSC screens that relates unique IDs to human-readable names.
- `mapping.tsv.gz`: gene exon mapping from VastDB.
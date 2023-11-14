# GeNA
`GeNA` (Genotype-Neighborhood Associations) is a tool for identifying genetic variant associations to the abundance of cell states in single-cell datasets (cell state abundance quantitative trait loci, csaQTLs). In `GeNA`, we have adapted the framework we developed for [Covarying Neighborhood Analysis](https://github.com/immunogenomics/cna) in order to enable genome-wide csaQTL surveys in single-cell data. Instead of testing associations to predefined cell types, `GeNA` identifies the granular cell states whose abundance is most associated with genetic variants. The scripts required to run `GeNA` are stored in this repo.

We have evaluated GeNA in simulation to assess calibration and statistical power and we have applied `GeNA` in a genome-wide survey to scRNA-seq profiling from a population cohort of 969 individuals. Scripts documenting our use of `GeNA` in these analyses for our manuscript are found in a separate repository, [immunogenomics/GeNA-applied](https://github.com/immunogenomics/GeNA-applied/). For more information about GeNA please refer to our preprint linked [here].

# Installation
To use `GeNA`, you can clone this repository.
Dependencies:
- Python version 3.8.10
- R version 4.1.1
- PLINK version 2.00a2.3
- CNA version 0.1.6

# Tutorial
We illustrate how to use `GeNA` in a tutorial [here](https://github.com/immunogenomics/GeNA/blob/main/tutorial/Example_csaQTL_GWAS.ipynb). First, we demonstrate how to construct the single-cell data object format `GeNA` expects, then we summarize the arguments input to and files output from a call to `GeNA`. Finally, we illustrate basic characterization of example loci.

# Citation
If you use `GeNA` in your work, you can cite our preprint as: (Coming soon!)

# Contact
If you have questions about `GeNA` or require user support, please contact Laurie Rumker (Laurie_Rumker AT hms.harvard.edu) or post an issue on this repo.

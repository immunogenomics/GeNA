# GeNA
`GeNA` (Genotype-Neighborhood Associations) is a tool for identifying genetic variant associations to the abundance of cell states in single-cell datasets (cell state abundance quantitative trait loci, csaQTLs). In `GeNA`, we have adapted the framework we developed for [Covarying Neighborhood Analysis](https://github.com/immunogenomics/cna) in order to enable genome-wide csaQTL surveys in single-cell data. Instead of testing associations to predefined cell types, `GeNA` identifies the granular cell states whose abundance is most associated with genetic variants. The scripts required to run `GeNA` are stored in this repo.

We have evaluated GeNA in simulation to assess calibration and statistical power and we have applied `GeNA` in a genome-wide survey to sc-mRNAseq profiling from a population cohort of 969 individuals. Scripts documenting our use of `GeNA` in these analyses for our manuscript are found in a separate repository, [immunogenomics/GeNA-applied](https://github.com/immunogenomics/GeNA-applied/). For more information about GeNA please refer to our preprint linked [here].

# Installation
To use `GeNA`, you can clone this repository and add it to your `PYTHONPATH`.

# Tutorial
Coming soon!

# Citation
If you use `GeNA` in your work, you can cite our preprint as: <>

# Contact
If you have questions about `GeNA` or require user support, please contact Laurie Rumker (Laurie_Rumker AT hms.harvard.edu) or post an issue on this repo.

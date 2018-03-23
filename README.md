N14N15
======

N14/N15 ratio quantitation

to install run the `N14N15`:

```r
# add path to Bioconductor repositories
source("http://bioconductor.org/biocLite.R")
options(repos=biocinstallRepos(character()))

install.packages("devtools")
library("devtools")
install_github("vladpetyuk/N14N15")
```

A brief description:
* Input: mzXML (or mzML) and mzIdentML file of LC-MS/MS dataset from N15 metabolic labeling
* Runs the analysis of N14 and N15 isotopic envelopes to determine isotopic incorporation rate
* Output: for each peptide/charge combo (LC-MS/MS feature) reports the abundance ratio of N15 specie to N14 along with a number of measurement confidence metrics.

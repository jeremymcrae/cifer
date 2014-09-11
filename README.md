## CIFER: CNV of inheritance from exome read-depth

This package attempts to predict the inheritance state of CNV calls generated by
ConVex. This package examines the log2-ratio scores in the CNV in both parents
of a proband with a CNV, and compares the scores to a model constructed from the
log2-ratio scores in a population of unrelated trios (specifically, the parents
of the trios).

### Predicting inheritance:
Using the installed package should involve calling the classification function
and accessing the inheritance prediction as:
```R
library(cifer)
prediction <- classify_exome_cnv(proband_id, maternal_id, paternal_id, chrom, start, stop)
inheritance <- prediction$inheritance
```

The classify_exome_cnv() arguments are:
* proband_id: sample ID for the proband (eg "DDDP100001")
* maternal_id: sample ID for the proband's mother
* paternal_id: sample ID for the proband's father
* chrom: chromosome that the CNV is on (eg "1", "2", ..., "X")
* start: start nucleotide of the CNV as integer
* stop: stop nucleotide of the CNV as integer

Optional arguments are:
* cnv: dataframe of CNv info, not necessary for predicting inheritance, used
    for evaluating CNVs if we plot the population distribution.
# DATAFREEZE_DIR: path to folder containing files listing study samples, and
    their relationships.

The code should be available at: https://github.com/jeremymcrae/cifer
And can be obtained with: `git clone https://github.com/jeremymcrae/cifer.git`

### Building the CIFER package:
Change to the cifer directory, and start R.

Load devtools and roxygen (requires R >= 3.0.1), construct the package 
documentation, check the package contents are correct, which should complete 
without errors, aside from a note about a non-standard python directory. 
Finally, build the package:
```R
library(devtools)
library(roxygen)
document()
check()
build()
```

# Discrete_curve_group_code

This repository contains the code used for curve group analyses of genes as performed in [Lund et al. 2016.](http://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-016-0129-z)
The code was created by the Norwegian Computing Centre.


# Repository

This method requires R and some dependencies mentioned in the scripts.

The five scripts in this project constitute the analyses performed using the breast cancer case-control data set described in the publication.
The scripts are somewhat adopted from the original analysis to analyse prospective case-control data sets in general.

# Usage

This code is adopted to analyse case-control differences in gene expression levels in cancer studies employing a prospective design.
Data on time to diagnosis must be available and the methods allow for comparisons between two strata (e.g. with and without spread).

The main script is analyseData.r which calls the other scripts.

Adjustments of reading of the data (readCancerData.r) must be performed for each data set.

There are also choices to consider with regards to choosing the number of people/
duration in time that is included in each curve group (in analyseCancerData.r)

# Citation

When used, the method should be cited with reference to [Lund et al. 2016.](http://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-016-0129-z)

Also, the code is to be cited in the acknowledgement section of the paper as 
"The code used for analysis is provided in a repository at https://github.com/theresenoest/Discrete_curve_group_code"

or using the following doi: in the method section of the paper

Lisence: MIT

# Issues 
Please report any issues regarding the contents of the repository [here]. (https://github.com/theresenoest/Discrete_curve_group_code/issues) 

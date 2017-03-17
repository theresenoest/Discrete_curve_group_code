# Discrete_curve_group_code

This repository contains the code used for curve group analyses of gene expression data with respect to time to diagnosis as described in [Lund et al. 2016.](http://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-016-0129-z)
The code was created by the Norwegian Computing Center for data in the Norwegian Women and Cancer Study and has been modified to become more general.


# Repository

This method requires R and some dependencies mentioned in the scripts.

The five scripts in this repository constitute the analyses performed using the breast cancer case-control data set described in the publication.
The scripts are adopted from the original analysis to analyse any prospective case-control data set.

# Usage

This code is adopted to analyse case-control differences in gene expression levels in cancer studies employing a prospective design.
Data on time to diagnosis must be available and the methods allow for comparisons between two strata (e.g. with and without spread).

The main script is analyseCancerData.r which calls the other scripts.

Regarding the time windows to analyse, the number of people/
duration in time that is included in each curve group must be considered and chosen (in analyseCancerData.r).
Also, adjustments of e.g. reading of the data, variable names and file paths must be performed for each data set.
Further details can be found in each script.

# Citation

This section will be updated soon.

When used, the method should be cited with reference to [Lund et al. 2016.](http://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-016-0129-z)

Also, the code is to be cited using the following doi: doi:XXX when describing the methods or using e.g. the following sentence as an acknowledgement: "The code used for analysis was provided in a repository at https://github.com/theresenoest/Discrete_curve_group_code".

# Issues 
Please report any issues regarding the contents of the repository [here.](http://github.com/theresenoest/Discrete_curve_group_code/issues)

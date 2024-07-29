# diagmeta: Meta-Analysis of Diagnostic Accuracy Studies with Several Cutpoints
Official Git repository of R package **diagmeta**

[![License: GPL (>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN Version](https://www.r-pkg.org/badges/version/diagmeta)](https://cran.r-project.org/package=diagmeta)
[![GitHub develop](https://img.shields.io/badge/develop-0.6--0-purple)](https://img.shields.io/badge/develop-0.6--0-purple)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/diagmeta)](https://cranlogs.r-pkg.org/badges/diagmeta)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/diagmeta)](https://cranlogs.r-pkg.org/badges/grand-total/diagmeta)


## Description

R package **diagmeta** implements the method by Steinhauser et
al. (2016) for meta-analysis of diagnostic accuracy studies with
several cutpoints.
 
### References

[Steinhauser S, Schumacher M, Rücker G (2016): Modelling multiple thresholds in meta-analysis of diagnostic test accuracy studies. *BMC Medical Research Methodology*, **16**, 97](https://scholar.google.com/scholar?q=Steinhauser+Schumacher+Rücker+2016+BMC)


## Installation

### Current official [![CRAN Version](https://www.r-pkg.org/badges/version/diagmeta)](https://cran.r-project.org/package=diagmeta) release:
```r
install.packages("diagmeta")
```

### Current [![GitHub develop](https://img.shields.io/badge/develop-0.6--0-purple)](https://img.shields.io/badge/develop-0.6--0-purple) release on GitHub:

Installation using R package
[**remotes**](https://cran.r-project.org/package=remotes):
```r
install.packages("remotes")
remotes::install_github("guido-s/diagmeta")
```


### Bug Reports:

You can report bugs on GitHub under
[Issues](https://github.com/guido-s/diagmeta/issues).

or using the R command

```r
bug.report(package = "diagmeta")
```

(which is not supported in RStudio).

## diagmeta, version 0.6.0 (2023-mm-dd)

### Major changes

* Data sets with negative association between biomarker values and probability
  of target condition can be used without data manipulation 

* Order of panels in plot.diagmeta() can be specified by the user

### User-visible changes

* diagmeta();
  - new argument 'direction' to specify whether the probability of the target
    condition (e.g., a disease) is increasing or decreasing with higher values
    of the biomarker

* plot.diagmeta():
  - default set of graphs (argument 'which') changed from
    c("survival", "youden", "roc", "sroc") to
    c("regression", "cdf", "sensspec","youden", "roc", "sroc")

### Internal changes

* New internal function plot.diagmeta-internal()


## diagmeta, version 0.5-1 (2022-12-21)

### User-visible changes

* Change maintainer's email address

* New branch 'release' on GitHub starting with **diagmeta**, version
  0.5-1

### Internal changes

* Rename list element 'Cov.fixed' to 'Cov.common'


## diagmeta, version 0.5-0 (2022-04-22)

### Major changes

* Behaviour of print.diagmeta() and print.summary.diagmeta() switched
  (to be in line with other print and summary functions in R)

* Do not stop with an error if optimal cut-off cannot be determined
  for logistic distribution

* Calculate area under the curve for specificity given sensitivity

### Bug fixes

* diagmeta():
  - fix for erratic confidence limits of AUC which could be outside
    the admissible range from 0 to 1 or exclude the AUC estimate

### User-visible changes

* More concise printout for summary.diagmeta()

### Internal changes

* diagmeta():
  - new list elements 'AUCSens' and 'AUCSpec' to calculate AUC for
    sensitivity given specificity or vice versa (existing list element
    'AUC' is equal to 'AUCSens')

* New internal function catch() to catch value for an argument


## diagmeta, version 0.4-1 (2021-05-11)

### Bug fixes

* plot.diagmeta():
  - print correct confidence region for specificities in SROC curves
  
* diagstats():
  - print results for requested specificity if only argument 'spec' is
    provided

### User-visible changes

* Use Markdown for NEWS

### Internal changes

* diagmeta():
  - new list element 'Cov.fixed' with covariance matrix from fixed
    effects model


## diagmeta, version 0.4-0 (2020-04-02)

### Major changes

* New default model (argument 'model') in diagmeta(), i.e., common
  random intercept and common slope ("CICS"), due to estimation
  problems with the previous default ("DIDS") after changes in R
  package **lme4**
  
### Bug fixes

* plot.diagmeta():
  - correct line types for survival functions
    
### User-visible changes

* diagmeta():
  - print a more informative error message in case of a negative
    correlation between increasing cutoffs and sensitivity

* plot.diagmeta():
  - argument 'points' considered for plots of type "regression",
    "cdf", "survival", "Youden", "ROC" and "sensspec"
    
* Help pages:
  - use the default model in all examples


## diagmeta, version 0.3-1 (2019-04-11)

### User-visible changes

* Export R functions:
  - as.data.frame.diagmeta(), plot.diagmeta(), print.diagmeta(),
    print.diagstats(), print.summary.diagmeta(), summary.diagmeta()

* plot.diagmeta():
  - argument 'col.points' can be any color defined in colours()
  - new argument 'col.ci' to specify color of curves with confidence
    limits

### Internal changes

* diagmeta():
  - check for numerical values in arguments 'TP', 'FP', 'TN', 'FN',
    and 'cutoff'


## diagmeta, version 0.3-0 (2018-12-11)

### User-visible changes

* plot.diagmeta():
  new plot type to show sensitivity and specificity curves

* New arguments 'sens' and 'spec' in diagstats()

* print.summary.diagmeta():
  - print confidence interval for optimal cutoff
    (for normal distribution)

* New function as.data.frame.diagmeta()

### Bug fixes

* plot.diagmeta():
  - correct ROC curves for datasets with decreasing cutoff values
    for individual studies (points (0, 0) and (1, 1) were connected
    with the wrong values on the ROC curve)

### Internal changes

* diagmeta():
  - calculate and return lower and upper confidence limit for
    optimal cutoff (for normal distribution)
    

## diagmeta, version 0.2-0 (2018-03-23)

### First version released on CRAN

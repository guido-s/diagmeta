## diagmeta, version 0.4-1 (2021-02-dd)

### Bug fixes

* Print correct confidence region for specificities in SROC curves
  generated with plot.diagmeta()

### User-visible changes

* Use Markdown for NEWS
    

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

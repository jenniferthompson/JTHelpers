# JTHelpers

A package of helper functions I wrote for my own frequent use. I hope they help you too.

* Related to the [rms](https://cran.r-project.org/web/packages/rms/rms.pdf) package:
    * `rms_model_results()`: Create a tidy data.frame combining results of `summary.rms()` and `anova.rms()` for a given `rms` object, plus information from the original model fit.
    * `rms_calc_comparisons()`: Calculate estimates (XB or hazard/odds ratios, depending on model
    type) for specified values of a numeric model covariate, or all levels of a factor covariate, vs
    a reference level. Stores all in a tidy data.frame. Does not take interaction terms into account.

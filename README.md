# JTHelpers

A package of helper functions I wrote for my own frequent use. I hope they help you too.

* Related to the [rms](https://cran.r-project.org/web/packages/rms/rms.pdf) package:
    * `rms_model_results()`: Create a tidy data.frame combining results of `summary.rms()` and `anova.rms()` for a given `rms` object, plus information from the original model fit.
    * `rms_calc_comparisons()`: Calculate estimates (XB or hazard/odds ratios, depending on model
    type) for specified values of a numeric model covariate, or all levels of a factor covariate, vs
    a reference level. Stores all in a tidy data.frame. Does not take interaction terms into account.
    * `rms_estimate_with_int()`: Get estimates, on original or ratio scale, for a main variable at a
    given level of another variable using `summary.rms()`. Returns a tidy data frame with columns
    for main and interacting variable names and levels and numeric columns for effects and CIs. Will
    often be used in conjunction with `lapply` or `purrr::map`/`purrr::pmap`.
    * `rms_po_assume()`: Creates figures to visually examine proportional odds assumption from an
    `lrm` model fit. Based on code and method outlined in Harrell's Regression Modeling Strategies
    (2001).
* Related to negative binomial regression:
    * `calc_nb_counts()`: For a `glm.nb` model fit, calculates predicted counts and CIs for each row in a supplied design matrix.
    * `calc_nb_ratioci()`: For a `glm.nb` model fit, calculates incidence rate ratio and CI for a
specified continuous predictor variable and comparison (eg, 75th vs 25th percentile).
* Formatting:
    * `rndformat()`: Rounds and formats a numeric value to the same number of decimal places to give
    a cleaner look.
    * `formatp()`: Formats a numeric value to be "<0.0001", "<0.001", or rounded to 3 decimal
    places, depending on value. Intended for formatting p-values.
* `as.mids.update()`: Update to `as.mids()` from the `mice` package which does not require data sets
to be sorted by .imp field. This allows for combining manually created data sets with imputed data.
Motivating example was a study in which we were calculating number of days with a given condition,
when only some days had recorded data for that condition. I created multiple imputed data sets with
the condition imputed on each missing day, then calculated durations for each patient in each 
imputed data set, combining them with this updated function in order to work with the rest of the
`mice` package.
* `calc_cat_freqs()`: Calculate frequencies and proportions of a categorical variable, optionally
including overall totals and/or using alternate denominator(s).
* `create_countprocess_data()`: Creates a data set for use in a time-dependent Cox regression model.
Majority of this code was written by Zhiguo (Alex) Zhao with edits from Cole Beck; to the best of my
knowledge it is not available in another package, so is here for my convenience.
* `lm_diagnostics()`: Creates figures to check hetereoscedasticity and normality assumptions of a
linear regression model.
* `model_lrpval()`: Get p-value for a likelihood ratio test for an original model with `removeTerm` removed.

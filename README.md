## Bayesian inference of the carbon homeostasis model using starch and sucrose time-series data

### Workflow
1. Parameter estimation of the wild-type using `fitting_wt.R`.
2. Parameter estimation of the mutant using `fitting_mt.R`.
3. Extraction of parameter samples from stan sampling file using `parameter_sample.R`.
4. Calculation of the posterior distribution of the model prediction from parameter samples using `posterior_prediction.R`.

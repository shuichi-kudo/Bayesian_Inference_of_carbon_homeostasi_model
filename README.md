## Bayesian inference of the carbon homeostasis model using starch and sucrose time-series data

### Workflow
1. Parameter estimation of the wild-type using `fitting_wt.R`
2. Parameter estimation of the mutant using `fitting_mt.R`
3. Extraction of parameter samples from stan sampling file using `parameter_sample.R`
4. Calculation of the posterior distribution of the model prediction from parameter samples using `posterior_prediction.R`
5. Creating a summary table of the posterior estimates using `posterior_summary.R`
6. Visualization

### Code description
#### Parameter estimation of the wild-type
- `fitting_wt.R`  
  This code provides the Bayesian inference of wild-type data in `dataset.csv`. It processes the rawdata,  executes a stan code `fitting_wt.stan`, and checks the convergence of the MCMC samples.
- `fitting_wt.stan`  
  This stan code defines the core processes of the Bayesian inference, such as the numerical solvor of the model ODEs and the calculation of the likelihoods.


#### Parameter estimation of the mutants
- `fitting_mt.R`  
  This code executes the Bayesian inference and the model selection for each mutant. For each mutant, the parameter estimation is done using all the stan files in `/mutant_models_stan` directory, which defines every combination of the model parameters. It returns the MCMC samples of each model and its WAIC for each mutant.
- `mutant_models_stan/fitting_mt_modelX-XX.stan`
  This stan code defines the core process of the Bayesian inference. Each code has a different combination of the parameters estimated. The file name X-XX indicates the model ID and the ID is listed in `model_list.csv`


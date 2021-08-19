# *WiBB*

*WiBB*: an integrated method for quantifying the relative importance of predictive variables

Authors: Li, Qin; Kou, Xiaojun

(Ecography；DOI: 10.1111/ecog.05651)

### Data

- full_datasets: including independent datasets for three correlation structures, each with 200 datasets. Each dataset has a sample size of 1000, with one response variable (*y*) and four predictors (*x1*, *x2*, *x3*, *x4*). The first three (genuine) predictors were set to be strongly, moderately, and weakly correlated with the response variable, respectively (denoted by large, medium, and small Pearson correlation coefficients, *r*), while the correlation between the response and the last (spurious) predictor was set to be zero.
	+ dt\_lm\_r1 (correlation structure: *r* = 0.3, 0.2, 0.1, 0.0);
	+ dt\_lm\_r2 (correlation structure: *r* = 0.6, 0.4, 0.2, 0.0);
	+ dt\_lm\_r3 (correlation structure: *r* = 0.8, 0.6, 0.3, 0.0);

- sub_datasets: resampled datasets with varying sample sizes from the full datasets for fitting LM and GLM separately;
	+ dataset_lm: R data (*.rds*) for three correlation structures $\times$ six sample sizes (25, 50, 100, 200, 500, 1000);
	+ dataset_glm: R data (*.rds*) for three correlation structures $\times$ six sample sizes (25, 50, 100, 200, 500, 1000);

- fit_result: example LM and GLM fitting results with the sample size = 1000;
	+ r1\_s6\_lm\_res.rds: an example fitting result by linear model, including estimates of variable importance based on sum of weight (sw), relative sum of weight (s_wi), standardized beta regression (b\_star), *WiBB* (wi), and model ranking of the full candidate model set (ms.rank);
	+ r1\_s6\_glm\_res.rds: an example fitting result by generalized linear model; same content as the above;

- empirical_dataset: *Mimulus* data
	+ mimulus\_occ\_var.csv: 71 *Mimulus* species occurrence coordinates and six associated climatic variables (temperature of cold- est month, *T_cold*; growing degree days above 0 Celsius degree, *GDD0*; precipitation seasonality, *P_season*; aridity of growing season, *Aridity*; isothermality, *ISO*; and temperature-precipitation synchronicity, *TP_syn*);
	+ background_pts_var.csv: background locality data (randomly sampled from the overall distribution) and associated climatic data (same climatic variables as the above);

### R scripts

- Rscripts\_0\_functions.R: custrom functions, including data simulation, data sampling, and counting the correct ranking frequency;
- Rscripts\_1\_data\_simulation.R: data simulation and sub-sampling for later LM and GLM;
- Rscripts\_2\_analysis.R: LM and GLM model fitting and calculations of variable importance by four indices: SW, SWi, beta*, WiBB, with varying sample sizes and bootstrap replicates;
- Rscripts\_3\_empirical_dataset.R: apply WiBB method to an empirical dataset of 71 *Mimulus* species;

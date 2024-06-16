# Dynamic Treatment Effect Estimation with Interactive Fixed Effects and Short Panels

[Nicholas Brown](https://sites.google.com/msu.edu/nicholasbrown)<sup>1</sup> and
[Kyle Butts](https://www.kylebutts.com/)<sup>2</sup>
<br>
<sup>1</sup>Florida State University, <sup>2</sup>University of Arkansas

#### [Paper](https://www.econ.queensu.ca/sites/econ.queensu.ca/files/wpaper/qed_wp_1495.pdf) | [Five-minute Summary](https://www.kylebutts.com/papers/generalized-imputation-estimators/)


## Abstract

We present a unifying identification strategy of dynamic average treatment effect parameters for staggered 
interventions when parallel trends are valid only after controlling for interactive fixed effects. This 
setting nests the usual parallel trends assumption, but allows treated units to have heterogeneous exposure 
to unobservable macroeconomic trends. We show that any estimator that is consistent for the unobservable trends 
up to a non-singular rotation can be used to consistently estimate heterogeneous dynamic treatment effects. 
This result can apply to data sets with either many or few pre-treatment time periods. We also demonstrate 
the robustness of two-way fixed effects imputation to certain parallel trends violations and provide a 
test for its consistency. A quasi-long differencing estimator is proposed and implemented to estimate the 
effect of Walmart openings on local economic conditions.

## Replication

### Simulations

1. `Simulation-1-twfe_vs_factor.jl` contains code to produce Table 1

2. `Simulation-2-signal_to_noise.jl` contains code to produce Figure 1

There are set of helper functions in `Simulation-helpers.jl` and `Simulation-factor_imputation.jl` that are used in both simulations. `Simulations-helpers.jl` contains our data-generating process code and our TWFE and TWFE imputation estimators. `Simulation-factor_imputation.jl` contain the code to estimate our generalized imputation estimator.

### Application

1. `Walmart-data.R` contains all the data to produce our final sample.

- It takes imputed CBP data from [Eckert et. al.](https://www.fpeckert.me/cbp/), 1967 CBP data, 1980 Census Summary files, Walmart openings from [Arcidiacono et. al. (2020)](https://www.aeaweb.org/articles?id=10.1257/app.20180047) and combines it.
- Note that `efsy_panel_naics.csv` needs to be downloaded from their website and put in the `raw-data/` folder since it is too large for github.
- The sample restrictions on the data follow from [Basker (2005)](https://direct.mit.edu/rest/article/87/1/174/57523/Job-Creation-or-Destruction-Labor-Market-Effects)

2. `Walmart-analysis-cbp_testing_p.jl` determines the correct number of factors following [Ahn, Lee, and Schmidt (2013)](https://doi.org/10.1016/j.jeconom.2012.12.002). This is run first to determine the p to use for retail (p = 2) and wholesale retail (p = 1) employment. It also produces the naive-standard error figures.

3. `Walmart-analysis-did2s_cbp_bootstrap.R` and `Walmart-analysis-factor_cbp_bootstrap.jl` produce the TWFE imputation and generalized imputation bootstrapped estiamtes and saves them in `data/` folder. The first row in each corresponds to the true point estimates.

4. `Walmart-figures-event_study_cbp.jl` produces Figures 2 and 3 using the bootstrapped simulations from the above scripts.

5. `Walmart-figures-synthetic_control_style_plot.jl` produces Figure 4 synthetic control style plots.




## Citation

```
@techreport{brown2022unified,
  title={A Unified Framework for Dynamic Treatment Effect Estimation in Interactive Fixed Effect Models},
  author={Brown, Nicholas and Butts, Kyle},
  year={2022}
}
```


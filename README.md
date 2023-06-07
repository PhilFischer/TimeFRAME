# TimeFRAME
> R package for the analysis of isotopic time series data using Bayesian hierarchical models implemented in Stan.

## Installation
Make sure R is updated and install TimeFRAME using the R console. On Windows [RTools](https://cran.r-project.org/bin/windows/Rtools/) will additionally be required.
```
if (!require("devtools")) install.packages("devtools") 
devtools::install_github("PhilFischer/TimeFRAME")
```

TimeFRAME uses Stan as sampling library in the background. If there are any issues regarding this, please refer to the [RStan Documentation](https://mc-stan.org/users/interfaces/rstan).

## Usage
For viewing in-depth guides on the usage of this package please refer to the vignettes shipped with it.
You can check the available vignettes in the folder `vignettes` or using the `browseVignettes` command in R.
```
browseVignettes("TimeFRAME")
```

### Getting Started
Once the package is installed, import it into the workspace with `library(TimeFRAME)`.

Now the contained isotopic data can be viewed as `n2o_sources` and `n2o_frac`. To specify a model using these as sources and fractionation factors use the `frame_model` constructor.
```
model <- frame_model(n2o_sources, frac = n2o_frac)
```

A data frame of isotopic measurements `df` can be fit for instance using the classical isotopic mapping technique or using the original FRAME model with independent time steps.
```
fit <- isotope_map(model, df)
print(fit)

fit2 <- fit_frame(model, df)
print(fit2)
```

### Fitting Time Series Models

Complex time series models include spline GLMs, Gaussian processes on measurements and Dirichlet-Gaussian process priors. They require a model specification `model` and a data frame of isotopic measurements `df` as well as a vector of time points `t` when the measurements were taken. 

```
# Fit independent time steps model
fit.independent <- fit_frame(model, df, t = t)

# Use M=8 degrees of freedom for splines of source contributions, and M.r=4 degrees of freedom for fractionation
fit.splines <- fit_glm(model, df, t = t, M = 8, M.r = 4)

# Use scaled correlation length of rho=0.2 for source contribution and rho.r = 0.5 for fractionation, 
# but update these values according to the posterior probability of a hierarchical model
fit.hdgp <- fit_dgp(model, df, t = t, rho = 0.2, rho.r = 0.5, estim.rho = TRUE)
```

![appl_est_th](https://github.com/PhilFischer/TimeFRAME/assets/36499405/9b94e9af-80b4-44c4-80d2-d8a5dcc75e78)

## References

Original FRAME Model
> M. P. Lewicki, D. Lewicka-Szczebak, and G. Skrzypek, “FRAME—monte carlo model for evaluation of the stable isotope mixing and fractionation,” PLOS ONE, vol. 17, no. 11, V. Kovtun, Ed., e0277204, Nov. 2022. doi: 10.1371/journal.pone.0277204.

Statistical Sampling
> Stan Development Team (2023). RStan: the R interface to Stan. R package version 2.21.8. https://mc-stan.org

N2O Isotopic Data
> L. Yu, E. Harris, D. Lewicka-Szczebak, et al., “What can we learn from N2O isotope data? – analytics, processes and modelling,” Rapid Communications in Mass Spectrometry, vol. 34, no. 20, Aug. 2020. https://doi.org/10.1002/rcm.8858

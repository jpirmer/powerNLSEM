# powerNLSEM

This is an `R` package to conduct the model-implied simulation-based power estimation (MSPE) procedures to find the minimum sample size for a given power within a nonlinear Structural Equation Model (NLSEM) for several parameters of interest (POI). 

## Install the Latest Working Version from Github
This requires the package `devtools`.

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("jpirmer/powerNLSEM", build_vignettes = TRUE)
```

Use `build_vignettes = T` to be able to see the documentation linked in "Getting Started". 


### Install the Submitted Version from GitHub
If you wish to install the version of the package as it was submitted in 2023, please use

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("jpirmer/powerNLSEM", build_vignettes = TRUE, 
                         ref = "Submitted2023")
```

`Submitted2023` is the branch name.


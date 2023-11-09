# MonotoneHazardRatio

## Overview

MonotoneeHazardRatio is a tool for nonparametric estimation and inference of a monotone non-decreasing hazard ratio, based on the [study](https://arxiv.org/abs/2205.01745#:~:text=Abstract%3A%20The%20ratio%20of%20the,time%2Dto%2Devent%20analysis.) by Y. Wu and T. Westling. 

## Dependent packages

Our packages needs the following packages to work. 
```
library(survival)
library(fdrtool)
library(KernSmooth)
library(twostageTE)
```

## Usage

It is staightforward to use this package. First you need to import the data (optional: split the data into two groups "S" and "T" such that the hazard ratio $\lambda_S/\lambda_T$ is non-decreasing). Pass the dataframes along with the evaluation grid to the function `monotoneHR()`, which takes $\alpha =0.05$ as the default confidence level, to have the hazard ratio and its confidence intervals estimated.

## Example

Here is an example with some data available [here](https://github.com/YujianWu2008/MonotoneHazardRatio/blob/main/example.csv). The estimated hazard ratio is stored in `theta$hr`, while the confidence intervals are stored in `theta$ci.lower` and `theta$ci.upper`. 

```
library(MonotoneHazardRatio)

### Use the example data in the package
data(survData)

### split it into two dataframes
f.data <- survData[survData$group == 'S']
g.data <- survData[survData$group == 'T']

### Evaluation grid
t.grid <- seq(0, 10, 1)

### Estimation and inference
theta <- monotoneHR(t.grid, f.data, g.data)
```

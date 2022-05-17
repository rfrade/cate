# cate
Conditional Average Treatment Effects of continous covariate

This package is an attempt to implement the CATE as proposed in:

Abrevaya, J., Hsu, Y. C., & Lieli, R. P. (2015). Estimating conditional average treatment effects. Journal of Business & Economic Statistics, 33(4), 485-505.

# What it does?
Given a treatment, the function __cate__ computes the Conditional ATE. More specifically:

*"First, the propensity score, the probability of treatment conditional on X, is estimated by ei- ther a kernel-based regression (the fully nonparametric case) or by a parametric model (the semiparametric case). In the second step the observed outcomes are weighted based on treatment status and the inverse of the estimated propensity score, and local averages are computed around points in the support X1"* (p. 486)


## Example

To illustrate how to use it, we will use data of a cash transfer program in Uruguay. The cash was given to pregnant women. The outcome variable is the weight of the baby. We will estimate the CATE of the weight of the babies of the mothers
who joined the program in comparison to mothers who didn't.

To follow the example, download the data from https://github.com/rfrade/cate_r/blob/main/birthdata.txt and set it up and change on the code bellow.


```
library(devtools)
install_github("rfrade/cate")

library("cate")
library(dummies)
library(ggplot2)
library(dplyr)


## DATA PREPARATION
setwd("/path/to/birthdata.csv")
data = read.csv("data_uruguai.csv", sep = ",", header = T)

data = data[data$applied == 1,]
n = dim(data)[1];

trim_score = 0.01
data = data[data$score_m > -trim_score & data$score_m <= trim_score,]
treat = dim(data[data$score_m > -trim_score & data$score_m <= 0, ])[1]
ntreat = dim(data[data$score_m > 0 & data$score_m <= trim_score,])[1]

p_hat = treat/(treat + ntreat)

y = data$weight
data$age_m = data$age_m+runif(length(data$age_m))-.5;

## CONDITIONING VARIABLE WILL BE AGE
conditioning_variable = data$age_m
d = data$treated;

data$weight = NULL
data$applied = NULL
data$score_m = NULL
data$p_low_weight = NULL
data$bajo2500 = NULL
data$treated = NULL

# Different Bandwidths for estimation
h = sd(data$age_m)*c(.25, .5, 1, 5)
kerntype='norm'

alpha=0.1
data$p_hat = p_hat

## COMPUTE CATE FOR THE 4 bandwidths
graph = ggplot() 
for (h1 in h) {
  cate_df = cate::cate(y, d, data, conditioning_variable, h1, alpha, 0.2)
  graph = graph + geom_line(aes_string(x = "grid", y = "cate_estimate"), data = cate_df)
  print(mean(cate_df$cate_estimate))
}

# The result will be the ATE conditioned on the age of the mother.
print(graph)

```

Any errors are mine. Still in development version. Based on the original matlab code developed by the authors: https://sites.google.com/site/robertplieli/research

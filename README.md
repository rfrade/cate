# cate
Conditional Average Treatment Effects of continous covariate

This package is an attempt to implement the CATE as proposed in:

Abrevaya, J., Hsu, Y. C., & Lieli, R. P. (2015). Estimating conditional average treatment effects. Journal of Business & Economic Statistics, 33(4), 485-505.

# What it does?
Given a treatment, the function __cate__ computes the Conditional ATE. More specifically:

*"First, the propensity score, the probability of treatment conditional on X, is estimated by ei- ther a kernel-based regression (the fully nonparametric case) or by a parametric model (the semiparametric case). In the second step the observed outcomes are weighted based on treatment status and the inverse of the estimated propensity score, and local averages are computed around points in the support X1"* (p. 486)

In this package, I estimate the propensity score by logit. In a future version I will add the non-parametric case.

# The main equations
Under the unconfoundedness assumption, and being Y the outcome variable, D the treatment indicator, X a group of covariates, the average treatment effect conditioned on $x_1$ is given by:

$$ \tau(x_1) = E[ E[Y|D=1,X] - E[Y|D=0,X] | X_1 = x_1 ] $$

Using the inverse probability weighted estimator of the ATE, we have:

$$ \tau(x_1) = E \left[ \frac{DY}{p(X)} - \frac{(1-D)Y}{1-p(X)}  | X_1 = x_1  \right] $$

Being $K(\cdot)$ a kernel function, $h_1$ a bandwidth, the CATE estimator was implemented using the following expression:

$$ \tau(x_1) = \frac{a}{b} $$

Where

$$ a =  \frac{1}{nh_1} \sum_{i=1}^{n} \left( \frac{DY}{p(X)} - \frac{(1-D)Y}{1-p(X)}  \right)  K \left(\frac{X_1 - x_1}{h_1} \right) $$

And

$$ b = \sum_{i=1}^{n} K \left(\frac{X_1 - x_1}{h_1} \right) $$


## Example

To illustrate how to use it, we will use data of a cash transfer program in Uruguay. The cash was given to pregnant women. The outcome variable is the weight of the baby. The continous covariate which we are condinioning to is the age of the mother. We will estimate the CATE of the weight of the babies of the mothers who joined the program in comparison to mothers who didn't. This estimator may have strong policy implications given the possibility of identifying groups in which the treatment has to be adjusted.

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

The final result will show a image like the following, which can be thought as a continuous "histogram" given the conditioning variable (age). In this case we can see older women benefit more from the program:

<img width="685" alt="image" src="https://user-images.githubusercontent.com/6012568/172379328-4f7416d3-b3e1-48be-999c-3e870f8faf25.png">


Any errors are mine. Still in development version. Based on the original matlab code developed by the authors: https://sites.google.com/site/robertplieli/research

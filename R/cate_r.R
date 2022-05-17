# CONTAINS THE MAIN CATE FUNCTION AND SOME UTILITY FUNCTIONS


#' Main function of the package. Estimate the Conditional ATE. 
#' @param y outcome variable
#' @param d treatment indicator
#' @param x_df covariates
#' @param conditioning_variable 
#' @param h bandwidth
#' @param alpha propensity score trimming parameter
#' @param grid_trim_param trimming parameter for the estimation grid
cate = function(y, d, x_df, conditioning_variable, h, alpha, grid_trim_param) {
  
  x_df$d = d
  x_df$y = y
  
  grid = get_grid(conditioning_variable, grid_trim_param)
  
  #STEP 1: Propensity Score
  if (!("p_hat" %in% colnames(x_df))) {
    x_df$p_hat = get_propensity_score(d, x_df)
  }
  
  x_df = trim_data(x_df, alpha)
  n = dim(x_df)[1]
  
  #STEP 2: CATE
  y = x_df$y
  phat = x_df$p_hat
  d = x_df$d
  a=((d*y)/phat - ((1-d)*y)/(1-phat))
  
  #lbd=0.001; %lower bound for f1 density estimate
  f1hat=kdensity(grid, conditioning_variable, h, "norm")
  
  print('Estimating CATE at gridpoints...')
  
  WeightedATEhat=numeric()
  for (j in 1:length(grid)) {
    
    z=(conditioning_variable-grid[j])/h;
    kernel_z = dnorm(z)
    WeightedATEhat[j]=(1/(n*h))*sum( a*kernel_z )
    
  }
  
  #CATE point estimate
  cate_estimate=WeightedATEhat/f1hat #NumX1Grid x 1
  
  cate_df = data.frame(cate_estimate, grid)
  
  return(cate_df)
}

#' Creates a grid based on which the CATE will be estimated
#' @param conditioning_variable the continuous variable for which CATE will
#' be estimated. Like age, salary, education...
#' @param grid_trim_param only create the grid for [param, 1-param]
get_grid = function(conditioning_variable, grid_trim_param = 0.2) {
  lower_bound = quantile(conditioning_variable, grid_trim_param)
  upper_bound = quantile(conditioning_variable, 1-grid_trim_param)
  step_size = (upper_bound-lower_bound)/50
  grid = seq(lower_bound, upper_bound, step_size)
  return(grid)
}

#' Computes the propensity score using logit 
#' @param d covariates
#' @param x binary variable indicating treatment
get_propensity_score = function(d, x) {
  logit = glm(d ~ ., data = x, family = binomial(link = "logit"))
  p_hat = predict(logit, x, type="response")
  return(p_hat)
}

#' Computes propensity score the deletes the data outside the threshold 
#' @param df data with covariates
#' @param alpha 
trim_data = function(df, alpha=0.5) {
  # TRIM OBS BASED ON ALPHA
  
  if (alpha < 0.5) {
    df$good_obs=as.integer((df$p_hat<=1-alpha) & (df$p_hat>=alpha))
    #print('Percentage of observations dropped=' + (n-sum(good_obs)) /n ))
    df = df[df$good_obs == T,]
    
  } else if (alpha == 0.5) {
    #%if alpha=.5, estimate it using the Crump et al. (2009) method
    alpha = seq(0.01, 0.49, 0.005)
    crit_value = numeric()
    for (j in 1:length(alpha)) {
      
      soma = 2*sum( (df$p_hat*(1-df$p_hat)>=alpha[j]*(1-alpha[j])) / (df$p_hat*(1-df$p_hat)) )
      divid = sum( (df$p_hat*(1-df$p_hat)>=alpha[j]*(1-alpha[j])) )
      
      crit_value[j] = soma / divid
      
    }
    
    alpha=alpha[crit_value>=1/(alpha*(1-alpha))]
    alphahat=min(alpha)
    print(paste('Estimated alpha=',alphahat))
    
    n_before = dim(df)[1]
    df = df[(df$p_hat<=1-alphahat) & (df$p_hat>=alphahat),]
    print(paste('Percentage of observations dropped=',((n_before-dim(df)[1])/(n_before))))
    
  } else {
    return("alpha has to be beteween >0 and 0.5")
  }
  return(df)
}







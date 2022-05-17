# UTILITY FUNCTIONS TO COMPUTE KERNELs

#' Uses ggplot to plot a kdensity
#' @param eval_points points at which kernel will be evalulated
#' @param pop vector with data points
#' @param h bandwidth
#' @param kernel_type "norm" for normal (default), "epa" for Epanishnikov
#' @example plot_density(seq(-1,1,by=0.1), rnorm(100), 0.2, "norm")
plot_density = function (eval_points, pop, h, kernel_type = "norm") {
  fx = kdensity(eval_points, pop, h, kernel_type)
  ggplot2::ggplot() + ggplot2::geom_line(ggplot2::aes(x = eval_points, y = fx))
}


#' Estimates a kernel density for a given vector
#' @param eval_points points at which kernel will be evalulated
#' @param pop vector with data points
#' @param h bandwidth
#' @param kernel_type "norm" for normal (default), "epa" for Epanishnikov
#' @example kdensity(-1:1, rnorm(100), 0.2, "norm") -> [0.26, 0.44, 0.30]
kdensity = function(eval_points, pop, h, kernel_type) {
  fx = double()
  i = 1
  for (x in eval_points) {
    fx[i] = kernel(x, pop, h, kernel_type = "norm")
    i = i + 1
  }
  return(fx)
}


#' This function returns exchange rates from any two currencies at the given period 
#' @param x point at which kernel will be evalulated
#' @param pop vector with data points
#' @param h bandwidth
#' @param kernel_type "norm" for normal (default), "epa" for Epanishnikov
#' @example kernel(0, rnorm(100), 0.2, "norm") - > 0.35
kernel = function(x, pop, h, kernel_type = "norm") {
  n = length(pop)
  z = (pop-x) / h
  kernel_values = c()
  if (kernel_type == "norm") {
    kernel_values = dnorm(z)
  } else {
    kernel_values = epa_kernel(z)
  }
  result = (1/(n*h))*sum(kernel_values)
  return(result)
}


#' Epanishnikov kernel function 
#' @param z point at which kernel will be evalulated
#' @example epa_kernel(0.5) - > 0.5625
epa_kernel = function(z) {
  return( 0.75*(1-z^2)*(z>-1)*(z<1) )
}





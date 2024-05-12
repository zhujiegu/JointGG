objective_function <- function(params) {
  
  
  # Numerical integration of the common part and store the results (all sample combined)
  common_part <- GH_com(Nr.cores=Nr.cores, GH_level=GH_level, dat=dat, params=params, plot_nodes=F)
  list_com <- common_part$list_com # values of common parts (list(length N) of list (dim^Q))
  list_nodes <- common_part$nodes
  
  # likelihood of each sample
  list_likl <- lapply(list_com, function(e){
    return(sum(e))
  })
  likl_log <- Reduce('+', list_likl) %>% log
  
  
  
  
  
  
  # Define the integrand as a function of x and param
  integrand <- function(x, param) {
    exp(-param * x^2)  # Example integrand
  }
  
  # Perform numerical integration
  integral_result <- integrate(integrand, lower = 0, upper = Inf, param = params)
  
  # Extract the value of the integral
  integral_value <- integral_result$value
  
  # Return the negative of integral for minimization (if maximizing the integral)
  return(-integral_value)
}


# Initial guess for the parameter
initial_guess <- c(1)

# Run optimization
result <- optim(par = initial_guess, fn = objective_function, method = "BFGS")

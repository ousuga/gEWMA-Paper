#Library
if (!require('readxl')) install.packages('readxl')
library(readxl)

# Parameter reading
parameters_normal <- read_excel("parameters.xlsx",sheet="Normal")
attach(parametros_normal)
dist <- "Normal"

# Function to obtain RL for a simulated gamma control chart   
rl_ind_normal     <- function(phi, lambda, L, k) {
  sd_x   <- sqrt(1/(1-phi^2)) 
  num_sd_z  <- (lambda*(1+phi*(1-lambda)))
  den_sd_z  <- ((1-phi^2)*(2-lambda)*(1-phi*(1-lambda)))
  sd_z      <- sqrt(num_sd_z/den_sd_z) 
  cl        <- L*sd_z
  delta     <- k*sd_x
  size      <- 5000
  X_0       <- arima.sim(list(order=c(1,0,0), ar=phi), n=size, rand.gen=rnorm)
  X         <- X_0 + delta
  Zant      <- 0
  flag      <- 0
  i         <- 0
  while (flag == 0 & i <= size) {
    i     <- i + 1
    Znew  <- lambda*X[i] + (1-lambda)*Zant
    Zant  <- Znew
    flag  <- abs(Znew) > cl
  }
  return(i)
}

# Function to obtain a summary for a RL of a 
# simulated normal control chart 
results_normal       <- function(Matriz, n, phi, lambda, L , k ) {
  rl_normal         <- replicate(Matriz[n], ar_ind_normal(Matriz[phi],
                       Matriz[lambda], Matriz[L] , Matriz[k])) 
  sum_normal        <- c(Matriz[phi], Matriz[lambda], Matriz[L],
                         Matriz[k], summary(rl_normal),
                         quantile(rl_normal, probs=c(0.01,0.99)), sd(rl_normal))
  names(sum_normal) <- c("phi","lambda","L","k", "Min", "Q1", "Median", 
                         "Mean", "Q3", "Max", "P1", "P99", "sd")
  sum_normal
}
resfinal_normal         <- apply(parameters_normal,1, results_normal, n="n", phi="phi", 
                                 lambda="lambda", L="L" , k="k" )
file.name               <- paste("Results_",dist,".csv")
write.csv(t(resfinal_normal), file=file.name) 


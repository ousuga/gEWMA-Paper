#Library
if (!require('readxl')) install.packages('readxl')
library(readxl)

# Parameter reading
parameters_gamma <- read_excel("parameters.xlsx",sheet="Gamma")
attach(parameters_gamma)
dist <- "Gamma"

# Function to obtain RL for a simulated gamma control chart   
rl_ind_gamma      <- function(phi, lambda, L, k, Shape, Scale) {
  Mean      <- (Shape*Scale)/(1-phi) 
  Variance  <- Shape*(Scale^2)
  sd_x      <- sqrt(Variance/(1-phi^2)) 
  num_sd_z  <- Variance*(lambda*(1+phi*(1-lambda)))
  den_sd_z  <- ((1-phi^2)*(2-lambda)*(1-phi*(1-lambda)))
  sd_z      <- sqrt(num_sd_z/den_sd_z) 
  cl        <- Mean + L*sd_z 
  delta     <- k*sd_x
  size      <- 5000
  X_0       <- arima.sim(list(order=c(1,0,0), ar=phi), n=size, 
                     rand.gen=rgamma, shape = Shape, scale = Scale)
  X         <- X_0  + delta
  Zant      <- Mean
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
# simulated gamma control chart 
results_gamma       <- function(matrix, n, phi, lambda, L , k ,Shape, Scale) {
  rl_gamma          <- replicate(matrix[n], rl_ind_gamma(matrix[phi], 
                      matrix[lambda], matrix[L], matrix[k] , matrix[Shape],
                      matrix[Scale])) 
  sum_gamma         <- c(matrix[phi], matrix[lambda], matrix[L], matrix[k], 
                        matrix[Shape], matrix[Scale], summary(rl_gamma),
                        quantile(rl_gamma, probs=c(0.01,0.99)), sd(rl_gamma))
  names(sum_gamma)  <- c("phi", "lambda", "L","k", "Shape", "Scale", "Min", 
                         "Q1", "Median", "Mean", "Q3", "Max.", "P1",
                         "P99", "sd")
  sum_gamma
}

resfinal_gamma    <- apply(parameters_gamma,1, results_gamma, 
                           n="n", phi="phi", lambda="lambda", L="L" , k="k",
                           Shape="Shape" , Scale="Scale")
file.name         <- paste("Results_gamma" , dist,".csv")
write.csv(t(resfinal_gamma), file=file.name) 


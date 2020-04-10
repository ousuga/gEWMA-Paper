#Library
if (!require('readxl')) install.packages('readxl')
library(readxl)

# Parameter reading
parameters_lognormal <- read_excel("parameters.xlsx",sheet="Lognormal")
attach(parameters_lognormal)
dist <- "Lognormal"

# Function to obtain RL for a simulated lognormal control chart   
rl_ind_lognormal <- function(phi, lambda, L, k, mu, sigma) {
  Mean      <- exp(mu + (sigma^2)/2 )/(1-phi) 
  Variance  <- (exp(sigma^2)-1)*exp(2*mu + (sigma^2))
  sd_x      <- sqrt(vari/(1-phi^2)) 
  num_sd_z  <- Variance*(lambda*(1+phi*(1-lambda)))
  den_sd_z  <- ((1-phi^2)*(2-lambda)*(1-phi*(1-lambda)))  
  sd_z      <- sqrt(num_sd_z/den_sd_z) 
  cl        <- Mean + L*sd_z 
  delta     <- k*sd_x
  size      <- 5000
  X_0       <- arima.sim(list(order=c(1,0,0), ar=phi), n=size,
                     rand.gen=rlnorm, meanlog = mu, sdlog = sigma)
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
results_lognormal      <- function(Matriz, n, phi, lambda, L , k , mu , sigma) {
  rl_lognormal         <- replicate(Matriz[n], rl_ind_lognormal(Matriz[phi],
                          Matriz[lambda], Matriz[L] , Matriz[k],
                          Matriz[mu], Matriz[sigma])) 
  sum_lognormal        <- c(Matriz[phi], Matriz[lambda], Matriz[L], 
                            Matriz[k], Matriz[mu], Matriz[sigma],
                            summary(rl_lognormal),quantile(rl_lognormal,
                            probs=c(0.01,0.99)), sd(rl_lognormal))
  names(sum_lognormal) <- c("phi","lambda","L","k", "mu", "sigma","Min", 
                            "Q1", "Median", "Mean", "Q3", "Max", "P1",
                            "P99", "sd")
  sum_lognormal
}
resfinal_lognormal         <- apply(parameters_lognormal,1, results_lognormal,
                                    n="n", phi="phi", lambda="lambda", L="L" ,
                                    k="k", mu="mu" , sigma="sigma" )
file.name               <- paste("Results_lognormal",dist,".csv")
write.csv(t(resfinal_lognormal), file=file.name) 

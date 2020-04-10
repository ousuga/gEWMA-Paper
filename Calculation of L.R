
# Algorithm 1 - Shiau and Hsu (2005)
# Compute an ARL estimate by Monte Carlo simulation

arl <- function(phi, lambda, L) {
  
Nsim <- 1000
size <- 5000

rl <- 0
control <- NA
#X <- 0
Z <- 0

for (j in 1:Nsim){
  Zeta <- 0
  X <- arima.sim(list(order=c(1,0,0), ar=phi), n=size, rand.gen=rnorm)
  for (i in 1:size){  
    Z[i]  <- lambda*X[i] + (1-lambda)*Zeta
    Zeta  <- Z[i]
    cl    <- L*sqrt((lambda*(1+phi*(1-lambda)))/((1-phi^2)*(2-lambda)*(1-phi*(1-lambda))))
    control <- abs(Z)>cl
  } 

  rl[j] <- which(control=="TRUE")[1]
  }
resul <- ifelse(rl==(size+1),size, rl)
mean(resul,na.rm = TRUE)

}

arl(0.0,1,3.090) 

# Algorithm 2 - Shiau and Hsu
# Compute the control limit multiple h/L/nsigmas

Lo <- 2
Hi <- 3

arl(phi= 0.1, lambda= 0.05, L=Lo) 
arl(phi= 0.1, lambda= 0.05, L=Hi) 

repeat{
  L <- (Lo + Hi)/2
  arl_l <- arl(phi= 0.1, lambda= 0.05, L=L)
  arl_l
  if(arl_l > 370.4) Hi<- L else Lo<-L
  if(abs(arl_l-370.4) <= 0.5) {break}
  L
}


if(abs(arl_l-370.4) <= 0.5) L<-L else 
  if(arl_l > 370.4) Hi<- L else Lo<-L

resul <- arl(phi= 0.1, lambda= 0.05, L=2.567462) 



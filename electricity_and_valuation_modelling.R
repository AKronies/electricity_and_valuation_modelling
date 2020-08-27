####Load Relevant Packages####
library(devtools)
library(easyGgplot2)
library(ggplot2)
library(MASS)
library(car)
library(stats4)
library(lubridate)
rm(list = ls())
set.seed(1)

####Wind Turbine#####
####...Technical Details...####
l <- 50
p <- 1.28
Cp <- 0.4
V <- 17
cut_in <- 3
cut_off <- 18
MW_max <- 3.5
A <- pi*l^2

###...Power Function of Wind Energy (approximated)...####
Power <- function(p,A,V,Cp,MW_max) {
  V <- replace(V, V > cut_off, 0)
  V <- replace(V, V < cut_in, 0)
  Power <- min(0.5*p*A*V^3*Cp*10^-6, MW_max)
  return(Power)
}

#...Plot..#
test <- matrix(NA, 1000, 2)
colnames(test) <- c("V", "Power")
test[,1] <- seq(1,25, by=((25-1)/(1000-1)))
for (i in 1:1000){
  test[i,2] <- Power(p,A,test[i,1],Cp, MW_max)
}
plot(test, type="l")


####Random Number Simulation for Valuation exercise####
extrayears <- 25 #WATCH OUT. Adjust in enddate
nPeriods <- 365*extrayears
nTrials <- 1000



####Income Modelling####
####...Wind data...####
hours_day <- 24
A_ <- 9 #7.9
k_ <- 2.5
V_ <- matrix(rweibull(nPeriods*nTrials, scale=A_, shape=k_), nrow = nPeriods, ncol = nTrials) #Weibull-Distribution Assumption

V_t <- replace(V_, V_ < 0, 0) #Wind speed
V_t_input <- sqrt(replace(V_t, V_t > cut_off, 0)*replace(V_t, V_t < cut_in, 0))
Power_t_ <- 0.5*p*A*V_t_input^3*Cp*10^-6
Power_t <- replace(Power_t_, Power_t_ > MW_max, MW_max)*hours_day

production <- sweep(Power_t,2,colMeans(Power_t),`/`)


####...Electricity Price Simulation####
#load electricity price data: System electricity prices on the Nordic regions are from Energi Data Service (https://www.energidataservice.dk/). 
#Data depict is based on hourly electricity prices through which I calculate daily averages.
df_daily <- readRDS("daily_electricity_prices.rds") #load data: System 

#Only look at data from 2014 onward until 2018
start <- 366 + 365 +1 +365  #366 for 2012  #1097 for 2014
nryears <- 4 #6 for when you start at 2012 #4 for when you start at 2014
df_daily <- df_daily[start:nrow(df_daily),]
#Remove 29/02/2012 and 29/02/2016 to have even years with 365 days a year
df_daily <- df_daily[df_daily$Date != "2012-02-29",]
df_daily <- df_daily[df_daily$Date != "2016-02-29",]


####......Calibrate seasonality f(t)###
Prices <- df_daily
PriceDates <- Prices$Date
rownames(Prices) <- Prices$Date
Prices$Date <- NULL
logPrices <- log(Prices)

PriceTimes <- seq(1/365, nryears, 1/365)

seasonMatrix <- function(t){
  matrix <- matrix(NA, length(t), ncol = 6)
  matrix[,1] <- sin(2*pi*t)
  matrix[,2] <- cos(2*pi*t) 
  matrix[,3] <- sin(4*pi*t)
  matrix[,4] <- cos(4*pi*t)
  matrix[,5] <- t
  matrix[,6] <- 1
  return(matrix)
}
C <- seasonMatrix(PriceTimes)


mldivide <- function(A, B, pinv = TRUE) {
  stopifnot(is.numeric(A) || is.complex(A),
            is.numeric(B) || is.complex(B))
  if (is.vector(A)) A <- as.matrix(A)
  if (is.vector(B)) B <- as.matrix(B)
  if (nrow(A) != nrow(B))
    stop("Matrices 'A' and 'B' must have the same number of rows.")
  if (pinv) {
    pinv(t(A) %*% A) %*% t(A) %*% B
  } else {
    qr.solve(A, B)
  }
}
mrdivide <- function(A, B, pinv = TRUE) {
  stopifnot(is.numeric(A) || is.complex(A),
            is.numeric(B) || is.complex(B))
  if (is.vector(A)) A <- t(A)
  if (is.vector(B)) B <- t(B)
  if (ncol(A) != ncol(B))
    stop("Matrices 'A' and 'B' must have the same number of columns.")
  
  t(mldivide(t(B), t(A), pinv = pinv))
}
seasonParam <- mldivide(C, logPrices[,1], pinv = FALSE)


#De-seasonalized log price
X <- logPrices-C%*%seasonParam


####......Calibrate the stochastic part X_t###
#Prices at t, X(t)
Pt <- X[2:nrow(X),]
#Prices at t-1, X(t-1)
Pt_1 <- X[1:(nrow(X)-1),]


#Discrection for daily prices 
dt <- 1/365

mrjpdf <- function(a, phi, mu_J, sigmaSq, sigmaSq_J, lambda){
  
  R <- 
    lambda * exp((-(Pt-a-phi*Pt_1-mu_J)^2)/(2*(sigmaSq+sigmaSq_J))) *
    (1/sqrt(2*pi*(sigmaSq+sigmaSq_J))) +
    
    (1-lambda) * exp((-(Pt-a-phi*Pt_1)^2)/(2*sigmaSq)) *
    (1/sqrt(2*pi*sigmaSq)) 
  
  - sum(log(R))
 
}

lb <- c(-Inf, -Inf, -Inf, 0, 0, 0)
ub <- c(Inf, 1, Inf, Inf, Inf, 1)

#Initial Values
X0 <- list(a = 0, phi = 0, mu_J = 0, sigmaSq = var(X), sigmaSq_J = var(X), lambda = 0.05)

#Solve maximum likelihood
params <- mle(mrjpdf, start = X0, method = "BFGS") #, lower = lb, upper = ub #if problems occur, try using "Nelder-Mead" or "L-BFGS-B"
#if problems occur, try different starting values in X0

alpha <- params@coef[1]/dt
kappa <- (1-params@coef[2])/dt
mu_J <- params@coef[3]
sigma <- sqrt(params@coef[4]/dt)
sigma_J <- sqrt(params@coef[5])
lambda <- params@coef[6]/dt


####......Simulate for 25 years###
n1 <- matrix(rnorm(nPeriods*nTrials, 0, 1), nPeriods, nTrials)
n2 <- matrix(rnorm(nPeriods*nTrials, 0, 1), nPeriods, nTrials)
j <- matrix(rbinom(nPeriods*nTrials, 1, lambda *dt), nPeriods, nTrials)
SimPrices <- matrix(0, nPeriods, nTrials)
SimPrices[1,] <- X[nrow(X),]
for (i in 2:nPeriods){
  SimPrices[i,] <- alpha*dt + (kappa*dt) * SimPrices[(i-1),] + #WATCH OUT: CHECK WITH BOTH: 1-kappa*dt and kappa*dt
    sigma * sqrt(dt) * n1[i,] + j[i,] * (mu_J + sigma_J*n2[i,])
}

#Add back Seasonality
enddate <- as.Date("2018/12/31") %m+% years(extrayears-1)
SimPriceDates <- seq(as.Date("2018/1/1"), enddate, "days")
SimPriceDates  <- SimPriceDates[-which(SimPriceDates=="2020-02-29")]
SimPriceDates  <- SimPriceDates[-which(SimPriceDates=="2024-02-29")]
SimPriceDates  <- SimPriceDates[-which(SimPriceDates=="2028-02-29")]
SimPriceDates  <- SimPriceDates[-which(SimPriceDates=="2032-02-29")]
SimPriceDates  <- SimPriceDates[-which(SimPriceDates=="2036-02-29")]
SimPriceDates  <- SimPriceDates[-which(SimPriceDates=="2040-02-29")]

SimPriceTimes <- seq(nryears+1/365, nryears+extrayears, 1/365) #8 because 2 years are added
CSim <- seasonMatrix(SimPriceTimes)
SimSeasonality <- CSim %*% seasonParam
SimSeasonalityMatrix <- matrix(NA, nPeriods, nTrials)
SimSeasonalityMatrix[,] <- SimSeasonality
logSimPrices <- SimPrices + SimSeasonalityMatrix - 0.045 * production #adjust prices according to specific production output dynamics of simulation: prices tend to be lower when production is high, and vice verse
PricesSim <- exp(logSimPrices)

#Plot electricity Price data
#Plot log prices and simulated log prices
plot(x = c(PriceDates, SimPriceDates)[1:k], y = c(logPrices[,1], 
                                                  logSimPrices[1:(k-nryears*365),path]), type = "l", col = "gray",
     xlab = "Date", ylab = "log(Prices)", cex.lab=1.5, cex.axis=1.5) # ,main = "log(price) and seasonality") , ylim=c(1, 4.5)
lines(x = SimPriceDates[1:(k-nryears*365)], y = logSimPrices[1:(k-nryears*365),path], col = "darkgoldenrod2", type = "l")
lines(x = c(PriceDates, SimPriceDates)[1:k], y = (seasonMatrix(c(PriceTimes,SimPriceTimes))%*%seasonParam)[1:k], col="darkgreen")
legend("bottomright", 
       legend = c("log(Prices)", "Simulated log(Prices)", "seasonality"), 
       col = c("gray","darkgoldenrod2", "darkgreen"), 
       lty = c(1), lwd=c(1), cex = 1.5)

#Plot real prices and simulates prices
PricesSim <- exp(logSimPrices)
plot(x = c(PriceDates, SimPriceDates)[1:k], y = c(Prices[,1],PricesSim[1:(k-nryears*365),path]), type = "l", col = "gray",
     xlab = "Date", ylab = "Prices", cex.lab=1.5, cex.axis=1.5) # ,main = "Actual Prices and Simulated Prices") , ylim=c(0, 100)
lines(x = SimPriceDates[1:(k-nryears*365)], y = PricesSim[1:(k-nryears*365),path], col = "darkgoldenrod2", type = "l")
lines(x = c(PriceDates, SimPriceDates)[1:k], y = (exp(seasonMatrix(c(PriceTimes,SimPriceTimes))%*%seasonParam))[1:k], col="darkgreen")
legend("topright", 
       legend = c("Prices", "Simulated Prices", "seasonality"), 
       col = c("gray","darkgoldenrod2", "darkgreen"), 
       lty = c(1), lwd=c(1), cex = 1.5)


####...Subsidy data...####
####......Old Subsidy Scheme in Denmark...####
R_0 <- 33.5
prob <- 0
delta_t <- 1/365
jump_prob <- prob * delta_t
mean_jump <- 5
sigma_jump <- 3
B_t <- matrix(3.1, nPeriods, nTrials)

set.seed(1)
z3 <- matrix(rnorm(nPeriods*nTrials), nPeriods, nTrials)
z_jump_yer_or_no <- matrix(rnorm(nPeriods*nTrials), nPeriods, nTrials)
jump <- ifelse(z_jump_yer_or_no <= qnorm(jump_prob), 1, 0)
jump_value <- -abs(mean_jump + sigma_jump*z3)


subsidy_t <- matrix(0, nPeriods, nTrials)
for (j in 1:nTrials){
  subsidy_t[1,j] <- R_0
  for (i in 2:nPeriods){
    subsidy_t[i,j] <- max(subsidy_t[i-1,j] + jump_value[i,j]*jump[i,j], 0)
  }
}
plot(subsidy_t[,7], type="l")
subsidy_t <- as.data.frame(subsidy_t)

#calculate cumulated power
cpower_t <- matrix(NA, nrow = nPeriods, ncol = nTrials)
for(j in 1:nTrials){
  cpower_t[,j] <- cumsum(Power_t[,j])
}

#power dataframe for subsidies in order to determine when subsidies run out (after 22.000 full-load hours)
power_t_subsidy <- cpower_t
power_t_subsidy[power_t_subsidy > MW_max*22000] <- MW_max*22000
for (j in 1:nTrials){
  for (i in (nPeriods-1):1){
    power_t_subsidy[i+1,j] <- power_t_subsidy[i+1,j]-power_t_subsidy[i,j]
  }
}


####......New Subsidy Scheme...####
S_0 <- 17.4 #simulated under maximal bid
prob_new <- 0
delta_t_new <- 1/365
jump_prob_new <- prob * delta_t
mean_jump_new <- 2.5
sigma_jump_new <- 2

set.seed(2)
z4 <- matrix(rnorm(nPeriods*nTrials), nPeriods, nTrials)
z_jump_yer_or_no_new <- matrix(rnorm(nPeriods*nTrials), nPeriods, nTrials)
jump_new <- ifelse(z_jump_yer_or_no_new <= qnorm(jump_prob_new), 1, 0)
jump_value_new <- -abs(mean_jump_new + sigma_jump_new*z4)


subsidy_t_new <- matrix(0, nPeriods, nTrials)
for (j in 1:nTrials){
  subsidy_t_new[1,j] <- S_0
  for (i in 2:nPeriods){
    subsidy_t_new[i,j] <- max(subsidy_t_new[i-1,j] + jump_value_new[i,j]*jump_new[i,j], 0)
  }
}
subsidy_t_new[(20*365+1):nrow(subsidy_t_new),] <- 0
plot(subsidy_t_new[,7], type="l")




####...OPEX...####
Costs <- 6000*12
C_t <- matrix(Costs/365, nPeriods, nTrials)

  
####...Total Income and Valuation......#### 
r <- 0.07
P_t <- as.data.frame(PricesSim)
B_t <- as.data.frame(B_t)
Power_t <- as.data.frame(Power_t)
subsidy_t <- as.data.frame(subsidy_t)
subsidy_t_new <- as.data.frame(subsidy_t_new)
power_t_subsidy <- as.data.frame(power_t_subsidy)
C_t <- as.data.frame(C_t)


####......Old Subsidy Scheme......####
Income <- (P_t + B_t) * Power_t + subsidy_t * power_t_subsidy - C_t
rdiscount <- matrix(NA, nPeriods, 1)
for (i in 1:nPeriods){
  rdiscount[i] <- 1/((1+r)^(i/365))
} 
dIncome <- matrix(NA, nPeriods, nTrials) #Price
for (j in 1:nTrials){
  dIncome[,j] <- Income[,j] * rdiscount
}
CAPEX <- 3500000
dcf_value <- colSums(dIncome)/CAPEX
hist(dcf_value, xlab = "PV/CAPEX", main = NA, col = "gray")

####......New Subsidy Scheme......####
Income_new <- (P_t + subsidy_t_new) * Power_t - C_t
CAPEX_new <- 3500000
dIncome_new <- matrix(NA, nPeriods, nTrials) #Price
for (j in 1:nTrials){
  dIncome_new[,j] <- Income_new[,j] * rdiscount
}
dcf_value_new <- colSums(dIncome_new)/CAPEX_new
hist(dcf_value_new, xlab = "PV/CAPEX", main = NA, col = "gray")


####......Under no Subsidy Scheme......####
Income_no_S <- P_t * Power_t - C_t
CAPEX_no_S <- 3500000
dIncome_no_S <- matrix(NA, nPeriods, nTrials) #Price
for (j in 1:nTrials){
  dIncome_no_S[,j] <- Income_no_S[,j] * rdiscount
}
dcf_value_no_S <- colSums(dIncome_no_S)/CAPEX_no_S 
hist(dcf_value_no_S, xlab = "PV/CAPEX", main = NA, col = "gray")




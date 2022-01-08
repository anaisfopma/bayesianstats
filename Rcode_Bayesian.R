########################
### DATA PREPARATION ###
########################

# Select NHANES data of the adult participants for the year 2011-2012.
library(NHANES)
library(dplyr)
data(NHANES)
dat <- data.frame(NHANES[NHANES$SurveyYr=="2011_12" & NHANES$Age>=18,])
dat <- dat %>% select(DirectChol, Pulse, BMI, PhysActiveDays, SleepHrsNight)
dat <- na.omit(dat)

# Mean center the predictors. 
dat$Pulse_mc <- dat$Pulse - mean(dat$Pulse)
dat$BMI_mc <- dat$BMI - mean(dat$BMI)
dat$PhysActiveDays_mc <- dat$PhysActiveDays - mean(dat$PhysActiveDays)
dat$SleepHrsNight_mc <- dat$SleepHrsNight - mean(dat$SleepHrsNight)

#####################
### GIBBS SAMPLER ###
#####################

# Make the function that will carry out Gibbs sampling with as input the data, the priors, the number of iterations, and the initial values.
gibber <- function(Y, X1, X2, N, mu_00, tau_00, mu_20, tau_20, alpha_0, beta_0, b0, b1, b2, var, n.iterations, iteration){ 
  
  # create storage for results with 5 columns: one each for the iteration counter, b0, b1, b2, and var
  results <- matrix(0,n.iterations,5)
  colnames(results)<-c("iteration", "b0", "b1", "b2", "var")
  results[1,] <- c(1, b0, b1, b2, var)
  
  # sample for n.iterations in a for loop
  for(iter in 2:n.iterations){
    
    # KEEP COUNT OF THE ITERATIONS
    iteration  <- iteration + 1
    
    # SAMPLE b0
    
    # step 1: get the conditional posterior mean of b0, this is mu_01
    mu_01 <- ((sum(Y-b1*X1-b2*X2)/var) + (mu_00/tau_00)) / (N/var + 1/tau_00)
    
    # step 2: get the conditional posterior variance of b0, this is tau_01
    tau_01 <- 1 / (N/var + 1/tau_00)
    
    # step 3: sample a value from the conditional posterior of b0, given
    # the mu_01 and tau_01
    b0 <- rnorm(1, mu_01, sqrt(tau_01))
    
    # SAMPLE b1
    
    # step 1: get the conditional posterior mean of b1, this is mu_11
    mu_11      <- ((sum(X1*(Y-b0-b1*X2))/var) + (mu_20/tau_20)) / (sum(X1^2)/var + 1/tau_10)
    
    # step 2: get the conditional posterior variance of b1, this is tau_11
    tau_11     <- 1 / (sum(X1^2)/var + 1/tau_10)
    
    # step 3: sample a value from the conditional posterior of b1, given the mu_11 and tau_11
    b1         <- rnorm(1, mu_11, sqrt(tau_11))
    
    # SAMPLE b2 with random walk MH sampler
    
    # step 1: get the conditional posterior mean of b2, this is mu_21
    mu_21 <- ((sum(X2*(Y-b0-b1*X1))/var) + (mu_20/tau_20)) / (sum(X2^2)/var + 1/tau_20)
    
    # step 2: get the conditional posterior variance of b2, this is tau_21
    tau_21 <- 1 / (sum(X2^2)/var + 1/tau_20)
    
    # step 3: sample the 'next' (candidate) value from the proposal distribution 
    b2_next <- rnorm(1, b2, 0.1)
    
    # step 4: sample a random value from the uniform distribution on the interval from 0 to 1
    u <- runif(1,0,1)
    
    # step 5: compute the target densities of the current and next values of b2
    target_current <- dnorm(b2, mu_21, sqrt(tau_21))
    target_next <- dnorm(b2_next, mu_21, sqrt(tau_21))
    
    # step 6: compute the acceptance ratio
    r <- (target_next/target_current)
    
    # step 7: accept the next value or retain the current value
    if (u < r)(b2 <- b2_next)
   
    # SAMPLE var
    
    # step 1: get the conditional posterior shape parameter, this is alpha
    alpha_1    <- N/2 + alpha_0
    
    # step 2: get the conditional posterior rate parameter, this is beta
    beta_1     <- sum((Y-(b0+b1*X1+b2*X2))^2)/2 + beta_0
    
    # step 3: we sample a value from the conditional posterior of s2, given the alpha and beta
    var        <- 1/rgamma(1, alpha_1, beta_1)
    
    # PUT THE RESULTS IN THE STORAGE 
    results[iter,] <- c(iteration, b0, b1, b2, var)}
  
  # make a data frame out of the results and remove the burn-in of a 1000 iterations
  as.data.frame(results[-c(1:1000), ])
}

###################
### RUN MODEL 1 ###
###################

# Specify the data.
Y = dat$DirectChol
X1 = dat$Pulse_mc
X2 = dat$BMI_mc
N = length(dat$DirectChol)

# Specify the priors.
# for b0
mu_00 = 0
tau_00 = 1000

# for b1
mu_10 = 0.01
tau_10 = 100

# for b2
mu_20 = -0.05
tau_20 = 10

# for s2
alpha_0 = 0.001
beta_0 = 0.001

# Specify the initial values.
b0 = 0
b1 = 0
b2 = 0
var = 1

# Specify the number of iterations.
n.iterations = 10000

# Start the iteration counter at 1.
iteration = 1

# Set a seed for reproducibility.
set.seed(321)

# Run the Gibbs sampler a first time: this is chain 1.
chain1 <- gibber(Y, X1, X2, N, mu_00, tau_00, mu_20, tau_20, alpha_0, beta_0, b0, b1, b2, var, n.iterations, iteration)

# Specify new initial values.
b0 = -2
b1 = 3
b2 = 0
var = 2

# Run the Gibbs sampler a second time with the new initial values: this is chain 2.
chain2 <- gibber(Y, X1, X2, N, mu_00, tau_00, mu_20, tau_20, alpha_0, beta_0, b0, b1, b2, var, n.iterations, iteration)

###########################
### CONVERGENCE MODEL 1 ###
###########################

# Assess convergence by plotting the results of each parameter with burn-in = 1000.
library(ggplot2)

# plot for intercept
p1 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = b0), color = "magenta") +
  geom_line(data = chain2, aes(x = iteration, y = b0), color = "midnightblue") +
  theme(axis.title.y = element_blank()) +
  ggtitle("b0")

# plot for slope of Pulse
p2 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = b1), color = "magenta") +
  geom_line(data = chain2, aes(x = iteration, y = b1), color = "midnightblue") +
  theme(axis.title.y = element_blank()) +
  ggtitle("b1")

# plot for slope of BMI
p3 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = b2), color = "magenta") +
  geom_line(data = chain2, aes(x = iteration, y = b2), color = "midnightblue") +
  theme(axis.title.y = element_blank()) +
  ggtitle("b2")

# plot for variance
p4 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = var), color = "magenta") +
  geom_line(data = chain2, aes(x = iteration, y = var), color = "midnightblue") +
  theme(axis.title.y = element_blank()) +
  ggtitle("var")

# Arrange the plots in one frame.
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

##################################################
### OBTAIN ESTIMATES AND INTERVALS FOR MODEL 1 ###
##################################################

# Take the results of the two chains together and get rid of the iteration column.
res <- rbind(chain1, chain2)
res <- res[ , -1]

# Make a storage for the estimates.
output <- matrix(0,4,5)
rownames(output) <- c("b0", "b1", "b2", "var")
colnames(output) <- c("Mean","SD","2.5%","97.5%", "MC Error")

# Estimate the means, standard deviations, lower and upper quantiles and put them in the storage.
output[,1] <- apply(res,2,mean)
output[,2] <- apply(res,2,sd)
output[,3] <- apply(res,2,quantile,0.025)
output[,4] <- apply(res,2,quantile,0.975)
output[,5] <- output[,2]/100

# Print the output. 
output

###################
### RUN MODEL 2 ###
###################

# Specify the data.
Y = dat$DirectChol
X1 = dat$PhysActiveDays_mc
X2 = dat$SleepHrsNight_mc 
N = length(dat$DirectChol)

# Specify the priors.
# for b0
mu_00 = 0
tau_00 = 1000

# for b1
mu_10 = 0.01
tau_10 = 10

# for b2
mu_20 = 0.02
tau_20 = 10

# for s2
alpha_0 = 0.001
beta_0 = 0.001

# Specify the initial values.
b0 = 0
b1 = 0
b2 = 0
var = 1

# Specify the number of iterations.
n.iterations = 10000

# Start the iteration counter at 1.
iteration = 1

# Set a seed for reproducibility.
set.seed(321)

# Run the Gibbs sampler a first time: this is chain 1.
chain1 <- gibber(Y, X1, X2, N, mu_00, tau_00, mu_20, tau_20, alpha_0, beta_0, b0, b1, b2, var, n.iterations, iteration)

# Specify new initial values.
b0 = -2
b1 = 3
b2 = 0
var = 2

# Run the Gibbs sampler a second time with the new initial values: this is chain 2.
chain2 <- gibber(Y, X1, X2, N, mu_00, tau_00, mu_20, tau_20, alpha_0, beta_0, b0, b1, b2, var, n.iterations, iteration)

###########################
### CONVERGENCE MODEL 2 ###
###########################

# plot for intercept
p5 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = b0), color = "chartreuse") +
  geom_line(data = chain2, aes(x = iteration, y = b0), color = "coral") +
  theme(axis.title.y = element_blank()) +
  ggtitle("b0")

# plot for slope of PhysActiveDays
p6 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = b1), color = "chartreuse") +
  geom_line(data = chain2, aes(x = iteration, y = b1), color = "coral") +
  theme(axis.title.y = element_blank()) +
  ggtitle("b1")

# plot for slope of SleepHrsNight
p7 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = b2), color = "chartreuse") +
  geom_line(data = chain2, aes(x = iteration, y = b2), color = "coral") +
  theme(axis.title.y = element_blank()) +
  ggtitle("b2")

# plot for variance
p8 = ggplot() + 
  geom_line(data = chain1, aes(x = iteration, y = var), color = "chartreuse") +
  geom_line(data = chain2, aes(x = iteration, y = var), color = "coral") +
  theme(axis.title.y = element_blank()) +
  ggtitle("var")

# Arrange the plots in one frame.
grid.arrange(p5, p6, p7, p8, ncol = 2, nrow = 2)

##################################################
### OBTAIN ESTIMATES AND INTERVALS FOR MODEL 2 ###
##################################################

# Take the results of the two chains together and get rid of the iteration column.
res <- rbind(chain1, chain2)
res <- res[ , -1]

# Make a storage for the estimates.
output <- matrix(0,4,5)
rownames(output) <- c("b0", "b1", "b2", "var")
colnames(output) <- c("Mean","SD","2.5%","97.5%", "MC Error")

# Estimate the means, standard deviations, lower and upper quantiles and put them in the storage.
output[,1] <- apply(res,2,mean)
output[,2] <- apply(res,2,sd)
output[,3] <- apply(res,2,quantile,0.025)
output[,4] <- apply(res,2,quantile,0.975)
output[,5] <- output[,2]/100

# Print the output. 
output

##################################
### POSTERIOR PREDICTIVE CHECK ###
##################################

# Create an empty storage for the simulated data sets.
datasets <- matrix(0, 1584, 18000)

# Simulate the data sets and put them in the storage.
for(x in 1:18000){
  Y_estimated <- res[x,1] + res[x,2] * X1 + res[x,3] * X2
  datasets[ ,x] <- Y_estimated
}

# Compute the test statistic for the observed data set.
TS_observed <- abs(3 * (mean(dat$DirectChol) - median(dat$DirectChol)) / sd(dat$DirectChol))

# Create an empty storage for the test statistics of the simulated data sets.
TS_simulated_all <- matrix(0, 1, 18000)
  
# Compute test statistic for each of the simulated data sets.
for(x in 1:18000){
  TS_simulated <- abs(3 * (mean(datasets[ ,x]) - median(datasets[ ,x])) / sd(datasets[ ,x]))
  TS_simulated_all[, x] <- TS_simulated
  }

# Put a 1 in the matrix if larger than the observed test statistic and a 0 if smaller.
for(x in 1:18000){ 
  if (TS_simulated_all[1, x] > TS_observed) (TS_simulated_all[1, x] <- 1) 
  else (TS_simulated_all[1, x] <- 0)
  }

# Calculate the Bayesian p-value by taking the mean.
mean(TS_simulated_all[1, ])

################################
### MODEL SELECTION WITH DIC ###
################################

# NOTE: this is done with RJAGS, therefore estimates used are different from the ones obtained with the programmed Gibbs sampler.
library(rjags)

# Specify Model 1.
mod1 <- "# Model 1
model {
# likelihood of the data
for (i in 1:length(DirectChol)) {
DirectChol[i] ~ dnorm(mu[i], tau)
mu[i] <- b[1] + b[2] * Pulse_mc[i] + b[3] * BMI_mc[i]
} 
# priors
tau ~ dgamma(0.001, 0.001)
b[1] ~ dnorm(0, 0.001)
b[2] ~ dnorm(0.01, 0.01)
b[3] ~ dnorm(-0.05, 0.1)
# calculate variance
sigma2 <- 1/tau
}"
z1 <- textConnection(mod1)

# Create a model object.
model1 <- jags.model(file = z1, data = dat, inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 321), n.chains = 2)

# Specify the burn-in period.
update(object = model1, n.iter = 1000)

# Obtain the DIC for Model 1.
dic.model1 <- dic.samples(model1, 10000, "pD")
dic.model1

# Specify Model 2.
mod2 <- "# Model 2
model {
# likelihood of the data
for (i in 1:length(DirectChol)) {
DirectChol[i] ~ dnorm(mu[i], tau)
mu[i] <- b[1] + b[2] * PhysActiveDays_mc[i] + b[3] * SleepHrsNight_mc [i]
} 
# priors
tau ~ dgamma(0.001, 0.001)
b[1] ~ dnorm(0, 0.001)
b[2] ~ dnorm(0.01, 0.1)
b[3] ~ dnorm(0.02, 0.1)
# calculate variance
sigma2 <- 1/tau
}"

z2 <- textConnection(mod2)

# Create a model object.
model2 <- jags.model(file = z2, data = dat, inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 321), n.chains = 2)

# Specify the burn-in period.
update(object = model2, n.iter = 1000)

# Obtain the DIC for Model 2. 
dic.model2 <- dic.samples(model2, 10000, "pD")
dic.model2

#########################################
### MODEL SELECTION WITH BAYES FACTOR ###
#########################################

# Load bain.
library(bain)

# Make an lm object.
mod <- lm(DirectChol ~ Pulse_mc + BMI_mc, data = dat)

# Set a seed for reproductibility. 
set.seed(321)

# Obtain the BF by specifying the hypotheses using the names in coef(mod).
results <- bain(mod,"BMI_mc < 0 & Pulse_mc > 0; BMI_mc = 0 & Pulse_mc = 0", fraction = 1, standardize = FALSE)
results

############################################
### COMPARISON WITH FREQUENTIST APPROACH ###
############################################

# Obtain estimates for Model 1.
m1 <- lm(DirectChol ~ Pulse_mc + BMI_mc , data = dat)
summary(m1)
confint(m1, level=0.95)

# Obtain estimates for Model 2. 
m2 <- lm(DirectChol ~ PhysActiveDays_mc + SleepHrsNight_mc , data = dat)
summary(m2)
confint(m2, level=0.95)
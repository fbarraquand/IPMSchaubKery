# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# 5.2 Feeding demographic data into the analysis of a matrix population model
# ===========================================================================

# 5.2.1 Using capture-recapture data in a matrix population model
# ---------------------------------------------------------------

library(IPMbook); library(jagsUI)
data(woodchat5)
str(woodchat5)
# $ ch : num [1:1902, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ age : num [1:1902] 2 2 2 2 2 2 2 2 2 2 ...
# $ repro: num [1:929, 1:3] 6 2 2 5 3 5 3 2 3 2 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:3] "Reproduction" "Year" "Age of mother"
# $ count: num [1:20] 91 119 131 88 139 145 148 116 112 106 ...

marr <- marrayAge(woodchat5$ch, woodchat5$age)

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
                  rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), mean.f=c(2.6, 3.6), T=15)
str(jags.data)
# List of 7
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ mean.f     : num [1:2] 2.6 3.6
# $ T          : num 15

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t] # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Model for population dynamics: a simple matrix population model
  # Define initial stage-specific population sizes
  N[1,1] <- 1
  N[2,1] <- 1
  # Loop over time
  for (t in 1:T){
    # Population projection
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
    # Annual growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  lambda <- ann.growth.rate[T]
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 1), mean.sa=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "lambda")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out1 <- jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na)
traceplot(out1)
print(out1, 3)
#             mean    sd    2.5%     50%   97.5% overlap0 f  Rhat n.eff
# mean.sj    0.304 0.016   0.274   0.304   0.336    FALSE 1 1.001  6000
# mean.sa    0.542 0.014   0.514   0.542   0.569    FALSE 1 1.001  2640
# mean.p     0.604 0.021   0.563   0.604   0.644    FALSE 1 1.001  6000
# lambda     1.018 0.027   0.966   1.018   1.073    FALSE 1 1.001  6000


# 5.2.2 Combining capture-recapture and productivity data in a matrix
#       population model
# -------------------------------------------------------------------

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
                  rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
                  age=woodchat5$repro[,3], T=15)
str(jags.data)
# List of 8
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ J          : num [1:929] 6 2 2 5 3 5 3 2 3 2 ...
# $ age        : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ T          : num 15

# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(mean.f[age[i]])
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }

  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                    # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Model for population dynamics: a simple matrix population model
  # Define initial stage-specific population sizes
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population projection
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
    # Annual growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  lambda <- ann.growth.rate[T]
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 1), mean.sa=runif(1, 0, 1))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "lambda")

# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000

# Call JAGS (ART <1 min), check convergence and summarize posteriors
out2 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na)
traceplot(out2)
print(out2, 3)
#               mean    sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj      0.304 0.016    0.273    0.304    0.337    FALSE 1 1.000  6000
# mean.sa      0.542 0.014    0.514    0.542    0.570    FALSE 1 1.002  1998
# mean.p       0.604 0.021    0.562    0.604    0.643    FALSE 1 1.001  6000
# mean.f[1]    2.677 0.079    2.527    2.675    2.832    FALSE 1 1.001  5055
# mean.f[2]    3.685 0.089    3.517    3.685    3.866    FALSE 1 1.001  2360
# lambda       1.030 0.029    0.973    1.030    1.089    FALSE 1 1.001  5372


# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# Run time approx. 2 mins

library(IPMbook) ; library(jagsUI)

# ~~~ requires data prepared in section 5.2 ~~~
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.3 Our first integrated population model
# =========================================

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
                  rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
                  age=woodchat5$repro[,3], C=woodchat5$count)
str(jags.data)
# List of 8
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ J          : num [1:929] 6 2 2 5 3 5 3 2 3 2 ...
# $ age        : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ C          : num [1:20] 91 119 131 88 139 145 148 116 112 106 ...

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial stage-specific population sizes: uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- N[1,t] * mean.f[1] / 2 * mean.sj + N[2,t] * mean.f[2] / 2 * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(mean.f[age[i]])
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                    # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 40000; nb <- 10000; nc <- 3; nt <- 3; na <- 3000

# Call JAGS (ART 2 min), check convergence and summarize posteriors
out3 <- jags(jags.data, inits, parameters, "model3.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             # n.thin=nt, n.adapt=na)
             n.thin=nt, n.adapt=na, parallel=TRUE)  # ~~~ for testing
traceplot(out3)
print(out3, 3)
#                         mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj                0.298  0.010    0.279    0.298    0.317    FALSE 1 1.002  1175
# mean.sa                0.539  0.012    0.515    0.539    0.564    FALSE 1 1.003   857
# mean.p                 0.608  0.019    0.571    0.608    0.645    FALSE 1 1.000  9184
# mean.f[1]              2.671  0.076    2.526    2.670    2.822    FALSE 1 1.000  7172
# mean.f[2]              3.674  0.085    3.510    3.672    3.844    FALSE 1 1.001  4281
# N[1,1]                47.751 32.297    2.676   43.470  110.762    FALSE 1 1.006   371
# N[2,1]                62.412 29.082    6.701   65.476  106.317    FALSE 1 1.006   380
# N[1,2]                53.101  4.429   44.608   53.092   61.675    FALSE 1 1.004   499
# N[2,2]                59.402  4.247   51.433   59.259   68.046    FALSE 1 1.002   897
# [ ...output truncated... ]
# N[1,20]               70.264  4.264   61.993   70.211   78.927    FALSE 1 1.000 12067
# N[2,20]               79.456  4.399   71.016   79.381   88.332    FALSE 1 1.001  2774
# sigma                 16.976  3.149   12.116   16.516   24.414    FALSE 1 1.000 10870
# ann.growth.rate[1]     1.023  0.041    0.945    1.026    1.087    FALSE 1 1.006   363
# ann.growth.rate[2]     1.016  0.007    1.003    1.015    1.029    FALSE 1 1.002   979
# [ ...output truncated... ]
# ann.growth.rate[18]    1.016  0.005    1.006    1.016    1.027    FALSE 1 1.000 18640
# ann.growth.rate[19]    1.016  0.005    1.006    1.016    1.027    FALSE 1 1.000 18640
# Ntot[1]              110.163  7.608   95.745  109.918  125.721    FALSE 1 1.001  1648
# Ntot[2]              112.503  6.546   99.570  112.493  125.646    FALSE 1 1.000 11582
# [ ...output truncated... ]
# Ntot[19]             147.325  7.114  133.524  147.263  161.763    FALSE 1 1.000 30000
# Ntot[20]             149.720  7.854  134.538  149.609  165.765    FALSE 1 1.000 30000


# ~~~~ code for Figure 5.2 ~~~~
mag <- 1.25
cex.tif <- mag * 1.25
lwd.tif <- 3 * mag
op <- par(mar=c(4, 4, 3, 0), las=1, cex=cex.tif, lwd=lwd.tif)
u <- col2rgb("grey82")
T <- length(woodchat5$count)
col.pol <- rgb(u[1], u[2], u[3], alpha=100, maxColorValue=255)
plot(out3$mean$Ntot, type="n",
     ylim=range(c(out3$q2.5$Ntot, out3$q97.5$Ntot, woodchat5$count)),
     ylab="Number", xlab="Year", las=1, cex=1.5, axes=FALSE)
axis(2, las=1, lwd=lwd.tif)
axis(2, at=c(90, 110, 130, 150), labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=lwd.tif)
polygon(c(1:T, T:1), c(out3$q2.5$Ntot, out3$q97.5$Ntot[T:1]), border=NA, col=col.pol)
points(out3$mean$Ntot, type="b", col="black", pch=16, lty=1, lwd=lwd.tif)
points(woodchat5$count, type="b", col="blue", pch=1, lty=2, lwd=lwd.tif)
legend("topleft", legend=c("Observed population counts", "Estimated population size"),
       pch=c(1, 16), lwd=c(lwd.tif, lwd.tif), col=c("blue", "black"), lty=c(2, 1), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# ~~~ requires data prepared in section 5.2 ~~~
data(woodchat5)
marr <- marrayAge(woodchat5$ch, woodchat5$age)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# 5.4 The 3-step approach to integrated population modeling
# =========================================================

# 5.4.1 Development of a model that links demographic data with population size (no code)
# 5.4.2 Formulation of the likelihood for each available data set separately (no code)
# 5.4.3 Formulation of the joint likelihood (no code)

# 5.4.4 Writing the BUGS code for the Integrated Population Model
# ---------------------------------------------------------------

# Bundle data and produce data overview
jags.data <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=dim(marr)[2],
                  rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=woodchat5$repro[,1],
                  year=woodchat5$repro[,2], age=woodchat5$repro[,3], C=woodchat5$count, pNinit=dUnif(1, 300))
str(jags.data)
# List of 10
# $ marr.j     : num [1:19, 1:20] 8 0 0 0 0 0 0 0 0 0 ...
# $ marr.a     : num [1:19, 1:20] 16 0 0 0 0 0 0 0 0 0 ...
# $ n.occasions: int 20
# $ rel.j      : num [1:19] 51 53 55 65 73 66 61 76 65 75 ...
# $ rel.a      : num [1:19] 36 39 44 61 61 50 43 61 51 53 ...
# $ J          : num [1:929] 6 2 2 5 3 5 3 2 3 2 ...
# $ year       : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ age        : num [1:929] 1 1 1 1 1 1 1 1 1 1 ...
# $ C          : num [1:20] 91 119 131 88 139 145 148 116 112 106 ...
# $ pNinit     : num [1:300] 0.00333 0.00333 0.00333 0.00333 0.00333 ...

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit.sj[t] ~ dnorm(mu.sj, tau.sj)
    sj[t] <- ilogit(logit.sj[t])          # Back-transformation from logit scale
    logit.sa[t] ~ dnorm(mu.sa, tau.sa)
    sa[t] <- ilogit(logit.sa[t])          # Back-transformation from logit scale
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    log.f[1,t] ~ dnorm(mu.f[1], tau.f[1])
    f[1,t] <- exp(log.f[1,t])             # Back-transformation from log scale
    log.f[2,t] ~ dnorm(mu.f[2], tau.f[2])
    f[2,t] <- exp(log.f[2,t])             # Back-transformation from log scale
  }

  mean.sj ~ dunif(0, 1)
  mu.sj <- logit(mean.sj)                 # Logit transformation
  mean.sa ~ dunif(0, 1)
  mu.sa <- logit(mean.sa)                 # Logit transformation
  sigma.sj ~ dunif(0, 3)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 3)
  tau.sa <- pow(sigma.sa, -2)

  for (j in 1:2){
    mean.f[j] ~ dunif(0, 10)
    mu.f[j] <- log(mean.f[j])             # Log transformation
    sigma.f[j] ~ dunif(0, 3)
    tau.f[j] <- pow(sigma.f[j], -2)
  }

  mean.p ~ dunif(0, 1)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for initial stage-spec. population sizes: discrete uniform priors
  N[1,1] ~ dcat(pNinit)
  N[2,1] ~ dcat(pNinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois(N[1,t] * f[1,t] / 2 * sj[t] + N[2,t] * f[2,t] / 2 * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i],year[i]])
  }

  # Capture-recapture data (CJS model with multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]                      # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0)
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
  }
}
")

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters <- c("mean.sj", "mean.sa", "mean.f", "mean.p", "sigma.sj", "sigma.sa", "sigma.f",
                "sigma", "sj", "sa", "f", "N", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 12000; nb <- 2000; nc <- 3; nt <- 2; na <- 1000

# Call JAGS (ART 1 min), check convergence and summarize posteriors
out4 <- jags(jags.data, inits, parameters, "model4.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
             n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out4)
print(out4, 3)
#                         mean     sd     2.5%      50%    97.5% overlap0 f  Rhat n.eff
# mean.sj                0.301  0.014    0.275    0.301    0.327    FALSE 1 1.001  2027
# mean.sa                0.542  0.016    0.511    0.542    0.573    FALSE 1 1.000  5642
# mean.f[1]              2.673  0.083    2.502    2.673    2.829    FALSE 1 1.006  1000
# mean.f[2]              3.675  0.110    3.468    3.672    3.900    FALSE 1 1.012   185
# mean.p                 0.605  0.019    0.567    0.605    0.643    FALSE 1 1.000 15000
# sigma.sj               0.114  0.077    0.007    0.103    0.291    FALSE 1 1.011   203
# sigma.sa               0.122  0.086    0.005    0.109    0.318    FALSE 1 1.001  4371
# sigma.f[1]             0.041  0.033    0.000    0.034    0.120    FALSE 1 1.001  2706
# sigma.f[2]             0.064  0.041    0.003    0.061    0.153    FALSE 1 1.002  1090
# sigma                 13.028  3.575    7.287   12.584   21.169    FALSE 1 1.001  1491
# sj[1]                  0.305  0.027    0.253    0.304    0.365    FALSE 1 1.005   980
# [ ... output truncated ... ]
# sj[19]                 0.304  0.025    0.256    0.303    0.360    FALSE 1 1.001 13466
# sa[1]                  0.554  0.036    0.490    0.548    0.641    FALSE 1 1.001  9699
# [ ... output truncated ... ]
# sa[19]                 0.538  0.032    0.468    0.539    0.604    FALSE 1 1.001  2810
# f[1,1]                 2.723  0.159    2.462    2.704    3.114    FALSE 1 1.002  5080
# f[2,1]                 4.012  0.400    3.513    3.905    5.014    FALSE 1 1.003  1260
# [ ... output truncated ... ]
# f[1,20]                2.691  0.140    2.435    2.683    3.010    FALSE 1 1.002  3466
# f[2,20]                3.814  0.229    3.442    3.783    4.342    FALSE 1 1.004   576
# N[1,1]                41.533 28.882    2.000   37.000  103.000    FALSE 1 1.002  1220
# N[2,1]                60.414 26.969    6.000   64.000  103.000    FALSE 1 1.002  1080
# [ ... output truncated ... ]
# N[1,20]               71.405  8.953   54.000   71.000   90.000    FALSE 1 1.000 15000
# N[2,20]               79.945  7.909   65.000   80.000   95.000    FALSE 1 1.001  5014
# ann.growth.rate[1]     1.130  0.106    0.936    1.124    1.351    FALSE 1 1.001  3874
# [ ... output truncated ... ]
# ann.growth.rate[19]    1.019  0.068    0.890    1.019    1.158    FALSE 1 1.000 15000
# Ntot[1]              101.947  9.925   84.000  101.000  123.000    FALSE 1 1.000 12582
# [ ... output truncated ... ]
# Ntot[20]             151.349 10.377  131.000  151.000  173.000    FALSE 1 1.000 15000


# ~~~~ code for Figure 5.4 ~~~~
mag <- 1.25
cex.tif <- mag * 1.25
lwd.tif <- 3 * mag
op <- par(mar=c(4, 4, 3, 0), las=1, cex=cex.tif, lwd=lwd.tif)
u <- col2rgb("grey82")
T <- length(woodchat5$count)
col.pol <- rgb(u[1], u[2], u[3], alpha=100, maxColorValue=255)
plot(out4$mean$Ntot, type="n",
     ylim=range(c(out4$q2.5$Ntot, out4$q97.5$Ntot, woodchat5$count)),
     ylab="Population size", xlab="Year", las=1, cex=1.5, axes=FALSE)
axis(2, las=1, lwd=lwd.tif)
axis(2, at=c(90, 110, 130, 150), labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=1:T, labels=NA, tcl=-0.25, lwd=lwd.tif)
axis(1, at=c(5, 10, 15, 20), labels=c(5, 10, 15, 20), tcl=-0.5, lwd=lwd.tif)
polygon(c(1:T, T:1), c(out4$q2.5$Ntot, out4$q97.5$Ntot[T:1]), border=NA, col=col.pol)
points(out4$mean$Ntot, type="b", col="black", pch=16, lty=1, lwd=lwd.tif)
points(woodchat5$count, type="b", col="blue", pch=1, lty=2, lwd=lwd.tif)
legend("topleft", legend=c("Observed population counts", "Estimated population size"),
       pch=c(1, 16), lwd=c(lwd.tif, lwd.tif), col=c("blue", "black"), lty=c(2, 1), bty="n")
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 5 : Introduction to integrated population models
# --------------------------------------------------------

# Run time for test script 3 mins, full run 15 hrs

library(IPMbook) ; library(jagsUI)

# 5.5 Simulation assessment of a simple IPM
# =========================================

# 5.5.1 Simulating data under an integrated population model
# ----------------------------------------------------------

# Pick values for the function arguments
T <- 20                                 # Number of years
phi <- c(0.3, 0.55)                     # Age specific survival probabilities (juv, adult)
f <- c(2.6, 3.6)                        # Age-specific productivity (1y, older)
Ni <- c(50, 50)                         # Initial pop. size for each age class (1y, older)

# Apply the function and produce data overview
set.seed(111167)                        # To initialize the RNGs at the same place
pop <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
str(pop)

# List of 14
# $ Ni         : num [1:2] 50 50
# $ phi        : num [1:3, 1:19] 0.3 0.55 0.55 0.3 0.55 0.55 0.3 0.55 ...
# $ f          : num [1:2, 1:20] 2.6 3.6 2.6 3.6 2.6 3.6 2.6 3.6 2.6 3.6 ...
# $ pBreed     : num [1:2, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ sex.ratio  : num [1:20] 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
# $ Im         : num [1:20] 0 0 0 0 0 0 0 0 0 0 ...
# $ ageOfIm    : num [1:20] 1 1 1 1 1 1 1 1 1 1 ...
# $ state      : num [1:3492, 1:20] 1 1 1 1 1 1 1 1 1 1 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ imYear     : logi [1:3492] NA NA NA NA NA NA ...
# $ reprod     : num [1:3492, 1:20, 1:3] 5 2 2 1 0 1 1 1 2 2 ...
# ..- attr(*, "dimnames")=List of 3
# .. ..$ : NULL
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# .. ..$ : chr [1:3] "F" "M" "Age"
# $ N          : num [1:6, 1:20] 50 50 100 158 336 0 49 52 101 165 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:6] "1-Year" "2-Year" "totAdults" "BornF" ...
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ breeders   : num [1:3, 1:20] 50 50 100 49 52 101 46 56 102 50 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:3] "1-Year" "2-Year" "totBreeders"
# .. ..$ : chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ totAdults  : Named num [1:20] 100 101 102 108 109 116 115 121 124 ...
# ..- attr(*, "names")= chr [1:20] "Y1" "Y2" "Y3" "Y4" ...
# $ totBreeders: Named num [1:20] 100 101 102 108 109 116 115 121 124 ...
# ..- attr(*, "names")= chr [1:20] "Y1" "Y2" "Y3" "Y4" ...

pop$state[488, 1:10]
# Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10
# NA NA  0  1  2  2 -1 NA NA  NA

pop$reprod[488,1:10,]
#      F  M Age
# Y1  NA NA  NA
# Y2  NA NA  NA
# Y3  NA NA  NA
# Y4   3  2   1
# Y5   5  3   2
# Y6   2  0   2
# Y7  NA NA  NA
# Y8  NA NA  NA
# Y9  NA NA  NA
# Y10 NA NA  NA

pop$N[,1:10]
#            Y1  Y2  Y3  Y4  Y5  Y6  Y7  Y8  Y9 Y10
# 1-Year     50  49  46  50  53  57  53  61  61  51
# 2-Year     50  52  56  58  56  59  62  60  63  62
# totAdults 100 101 102 108 109 116 115 121 124 113
# BornF     158 165 167 164 171 167 193 193 192 207
# BornT     336 346 321 337 328 356 380 369 395 389
# Im          0   0   0   0   0   0   0   0   0   0

pop$breeders[,1:10]
#              Y1  Y2  Y3  Y4  Y5  Y6  Y7  Y8  Y9 Y10
# 1-Year       50  49  46  50  53  57  53  61  61  51
# 2-Year       50  52  56  58  56  59  62  60  63  62
# totBreeders 100 101 102 108 109 116 115 121 124 113

pop1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
pop2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
pop3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)

# Pick a value of the observation error (SD) for the population survey
sigma <- 10

# Create the population survey data and produce data overview
count <- simCountNorm(N=pop1$totB, sigma=sigma)
str(count)
# List of 2
# $ sigma: num [1:20] 10 10 10 10 10 10 10 10 10 10 ...
# $ count: num [1:20] 111 84 85 94 116 132 115 89 82 61 ...

# Pick values for capture and recapture probabilities
cap <- 0.4                              # Initial capture probability (same for juv. and adults)
recap <- 0.6                            # Recapture probability

# Create the capture histories and produce data overview
ch <- simCapHist(state=pop2$state, cap=cap, recap=recap, maxAge=2)
str(ch)
# List of 5
# $ cap   : num [1:3, 1:20] 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 0.4 ...
# $ recap : num [1:2, 1:19] 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 0.6 ...
# $ maxAge: num 2
# $ ch    : num [1:2645, 1:20] 1 1 1 1 1 1 0 1 0 1 ...
# $ age   : num [1:2645] 2 2 2 2 2 2 2 2 2 2 ...

# Create m-arrays
marr <- marrayAge(ch$ch, ch$age)

# Pick probability to find a brood
pprod <- 0.3

# Create productivity data and produce data overview
pro <- simProd(reprod=pop3$reprod, pInclude=pprod)
str(pro)
# List of 4
# $ pInclude    : num [1:20] 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 ...
# $ females.only: logi FALSE
# $ prod.ind    : num [1:705, 1:3] 2 2 0 3 5 3 5 4 5 3 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:3] "Productivity" "Year" "Age of mother"
# $ prod.agg    : num [1:20, 1:2] 117 98 89 99 103 96 89 131 96 129 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:20] "1" "2" "3" "4" ...
# .. ..$ : chr [1:2] "Juveniles" "Surveyed broods"


# 5.5.2 Simulation results
# ------------------------

# ~~~~ code to run the simulations ~~~
# Define the 3 IPMs
# IPM 1: completely deterministic model

# Write JAGS model file
cat(file="model5.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- (N[1,t] * mean.f[1] / 2 + N[2,t] * mean.f[2] / 2) * mean.sj
    N[2,t+1] <- (N[1,t] + N[2,t]) * mean.sa
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  sJ1 ~ dpois(nJ1 * mean.f[1])
  sJ2 ~ dpois(nJ2 * mean.f[2])

  # Capture-recapture data (multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]           # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    resN[t] <- Ntot[t] - C[t]
  }
}
")

# IPM 2: include demographic stochasticity
# Write JAGS model file
cat(file="model6.txt", "
model {
  # Priors and linear models
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f[1] ~ dunif(0, 10)
  mean.f[2] ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pinit)
  N[2,1] ~ dcat(pinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] * mean.f[1] / 2 + N[2,t] * mean.f[2] / 2) * mean.sj)
    N[2,t+1] ~ dbin(mean.sa, N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  sJ1 ~ dpois(nJ1 * mean.f[1])
  sJ2 ~ dpois(nJ2 * mean.f[2])

  # Capture-recapture data (multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0)
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    resN[t] <- Ntot[t] - C[t]
  }
  }
")

# IPM 3: include environmental and demographic stochasticity
# Write JAGS model file
cat(file="model7.txt", "
model {
  # Priors and linear models
  for (t in 1:(n.occasions-1)){
    logit.sj[t] ~ dnorm(logit.mean.sj, tau.sj)
    sj[t] <- ilogit(logit.sj[t])      # Back-transformation from logit scale
    logit.sa[t] ~ dnorm(logit.mean.sa, tau.sa)
    sa[t] <- ilogit(logit.sa[t])      # Back-transformation from logit scale
    p[t] <- mean.p
  }

  for (t in 1:n.occasions){
    log.f[1,t] ~ dnorm(log.mean.f[1], tau.f[1])
    f[1,t] <- exp(log.f[1,t])         # Back-transformation from log scale
    log.f[2,t] ~ dnorm(log.mean.f[2], tau.f[2])
    f[2,t] <- exp(log.f[2,t])         # Back-transformation from log scale
  }

  mean.sj ~ dunif(0, 1)
  logit.mean.sj <- logit(mean.sj)      # Logit transformation
  mean.sa ~ dunif(0, 1)
  logit.mean.sa <- logit(mean.sa)      # Logit transformation
  sigma.sj ~ dunif(0, 3)
  tau.sj <- pow(sigma.sj, -2)
  sigma.sa ~ dunif(0, 3)
  tau.sa <- pow(sigma.sa, -2)

  for (j in 1:2){
    mean.f[j] ~ dunif(0, 10)
    log.mean.f[j] <- log(mean.f[j])   # Log transformation
    sigma.f[j] ~ dunif(0, 3)
    tau.f[j] <- pow(sigma.f[j], -2)
  }

  mean.p ~ dunif(0, 1)

  sigma ~ dunif(0.5, 100)
  tau <- pow(sigma, -2)

  # Population count data (state-space model)
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dcat(pinit)
  N[2,1] ~ dcat(pinit)

  # Process model over time: our model of population dynamics
  for (t in 1:(n.occasions-1)){
    N[1,t+1] ~ dpois((N[1,t] * f[1,t] / 2 + N[2,t] * f[2,t] / 2) * sj[t])
    N[2,t+1] ~ dbin(sa[t], N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
    C[t] ~ dnorm(N[1,t] + N[2,t], tau)
  }

  # Productivity data (Poisson regression model)
  for (i in 1:length(J)){
    J[i] ~ dpois(f[age[i],year[i]])
  }

  # Capture-recapture data (multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1 - p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t] * p[t]
    pr.a[t,t] <- sa[t] * p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t] * prod(sa[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
      pr.a[t,j] <- prod(sa[t:j]) * prod(q[t:(j-1)]) * p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
      pr.j[t,j] <- 0
      pr.a[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occasions-1)){
    pr.j[t,n.occasions] <- 1-sum(pr.j[t,1:(n.occasions-1)])
    pr.a[t,n.occasions] <- 1-sum(pr.a[t,1:(n.occasions-1)])
  }

  # Derived parameters
  # Annual population growth rate (added 0.001 to avoid possible division by 0)
  for (t in 1:(n.occasions-1)){
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t] + 0.001)
  }
  # Total population size
  for (t in 1:n.occasions){
    Ntot[t] <- N[1,t] + N[2,t]
    resN[t] <- Ntot[t] - C[t]
  }
}
")


# Simulation parameters

# Age specific survival probabilities (juv, adult)
phi <- c(0.3, 0.55)

# Age-specific productivity (1y, older)
f <- c(2.6, 3.6)

# Initial population size per age class
Ni <- c(50, 50)

# Number of years
T <- 20

# Observation error for the population survey
sigma <- 10

# Capture and recapture probabilities
cap <- 0.4             # Initial capture probability
recap <- 0.6           # Recapture probability

# Probability to find a brood whose reproductive output is recorded
pprod <- 0.3

# Number of simulations
# nsim <- 1200
nsim <- 4  # ~~~~ for testing

# Define matrices to store the result
res1 <- res2 <- array(NA, dim=c(106, 11, nsim))
res3 <- array(NA, dim=c(188, 11, nsim))
lam <- matrix(NA, nrow=T-1, ncol=nsim)

# Initial values
inits <- function(){list(mean.sj=runif(1, 0, 0.5))}

# Parameters monitored
parameters1 <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma",
                 "ann.growth.rate", "Ntot", "resN")

parameters2 <- c("mean.sj", "mean.sa", "mean.f", "mean.p", "sigma.sj", "sigma.sa",
                 "sigma.f", "sj", "sa", "f", "N", "sigma", "ann.growth.rate", "Ntot", "resN")

# MCMC settings
ni <- 10000; nb <- 5000; nc <- 3; nt <- 1; na <- 1000


# Start simulations
system.time(
  for (s in 1:nsim){
    set.seed(s)
    
    # Simulate 3 populations such that independent data sets are analyzed
    pop1 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
    pop2 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
    pop3 <- simPop(Ni=Ni, phi=phi, f=f, nYears=T)
    # Calculate the annual population growth rates
    lam[,s] <- pop1$totA[-1] / pop1$totA[-T]
    
    # Simulate the population survey data
    count <- simCountNorm(N=pop1$totB, sigma=sigma)$count
    
    # Simulate the capture histories and the corresponding m-arrays
    ch <- simCapHist(state=pop2$state, cap=cap, recap=recap, maxAge=2)
    marr <- marrayAge(ch$ch, ch$age)
    
    # Simulate productivity data
    P <- simProd(reprod=pop3$reprod, pInclude=pprod)
    
    # Aggregate productivity data to make the models 1 and 2 run faster
    J1 <- P$prod.ind[P$prod.ind[,3]==1,1]
    J2 <- P$prod.ind[P$prod.ind[,3]==2,1]
    sJ1 <- sum(J1)
    nJ1 <- length(J1)
    sJ2 <- sum(J2)
    nJ2 <- length(J2)
    
    # Bundle data
    jags.data1 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
                       rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ1=sJ1, nJ1=nJ1,
                       sJ2=sJ2, nJ2=nJ2, C=count)
    
    jags.data2 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
                       rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ1=sJ1, nJ1=nJ1,
                       sJ2=sJ2, nJ2=nJ2, C=count, pinit=dUnif(1, 300))
    
    jags.data3 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
                       rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), J=P$prod.ind[,1],
                       year=P$prod.ind[,2], age=P$prod.ind[,3], C=count, pinit=dUnif(1, 300))
    
    # Call JAGS
    out5 <- try(jags(jags.data1, inits, parameters1, "model5.txt",
                     n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
    if(!inherits(out5, "try-error"))
      res1[,,s] <- out5$summary
    
    out6 <- try(jags(jags.data2, inits, parameters1, "model6.txt",
                     n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
    if(!inherits(out6, "try-error"))
      res2[,,s] <- out6$summary
    
    out7 <- try(jags(jags.data3, inits, parameters2, "model7.txt",
                     n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, n.adapt=na, parallel=TRUE))
    if(!inherits(out7, "try-error"))
      res3[,,s] <- out7$summary
    
    print(s)
  } ) #s

save(res1, res2, res3, lam, phi, f, sigma, recap,
     file="IPM Simulation chapter 5.5.Rdata")

load("IPM Simulation chapter 5.5.Rdata")
library(RColorBrewer)
co <- brewer.pal(n=8, name='Blues')[c(7,5,3)]

# Function to ensure that only converged estimates are included
r.incl <- function(res, param, r.crit=1.1){
  h <- which(res[param, 8, ] < r.crit, arr.ind=TRUE)
  u <- as.numeric(which(table(h[,2])==length(param)))
  return(u)
}

# Number of simulations of converged runs that are included
sample.size <- 1000

i1 <- r.incl(res1, param=c(1:5, 46:65))[1:sample.size]
i2 <- r.incl(res2, param=c(1:5, 46:65))[1:sample.size]
i3 <- r.incl(res3, param=c(1:5, 128:147))[1:sample.size]

# Calculate the deviations of the realized population growth rates from population 1
#   and the estimated annual growth rates from the IPM
diff.ipm1 <- matrix(NA, ncol=length(i1), nrow=nrow(lam))
diff.ipm2 <- matrix(NA, ncol=length(i2), nrow=nrow(lam))
diff.ipm3 <- matrix(NA, ncol=length(i3), nrow=nrow(lam))

for (s in 1:length(i1)){
  diff.ipm1[,s] <- res1[47:65,1,i1[s]] - lam[,i1[s]]
}
for (s in 1:length(i2)){
  diff.ipm2[,s] <- res2[47:65,1,i2[s]] - lam[,i2[s]]
}
for (s in 1:length(i3)){
  diff.ipm3[,s] <- res3[129:147,1,i3[s]] - lam[,i3[s]]
}


# Code for figure 5.5
mag <- 1
cex.tif <- mag * 1.25
lwd.tif <- 1.25*mag
op <- par(mar=c(3, 4.5, 4, 0.5), las=1, cex=cex.tif, "lwd", "mfrow")
layout(matrix(1:6, 2, 3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1.1), TRUE)

labx <- c(expression(IPM[1]), expression(IPM[2]), expression(IPM[3]))

laby <- expression('Juvenile survival ('*phi[j]*')')
boxplot(cbind(res1[1,1,i1], res2[1,1,i2], res3[1,1,i3]), outline=FALSE,
        ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=NA, lwd=lwd.tif)
abline(h=phi[1], col="red", lwd=lwd.tif)
mtext("A", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Adult survival ('*phi[a]*')')
boxplot(cbind(res1[2,1,i1], res2[2,1,i2], res3[2,1,i3]), outline=FALSE,
        ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=NA, lwd=lwd.tif)
abline(h=phi[2], col="red", lwd=lwd.tif)
mtext("B", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Fecundity (f'[1]*')')
boxplot(cbind(res1[4,1,i1], res2[4,1,i2], res3[3,1,i3]), outline=FALSE,
        ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=NA, lwd=lwd.tif)
abline(h=f[1], col="red", lwd=lwd.tif)
mtext("C", at=0.5, line=0.5, cex=1.5 * cex.tif)

par(mar=c(4, 4.5, 3, 0.5), las=1, cex=cex.tif, lwd=4)
laby <- expression('Fecundity (f'[a]*')')
boxplot(cbind(res1[5,1,i1], res2[5,1,i2], res3[4,1,i3]), outline=FALSE,
        ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=labx, lwd=lwd.tif)
abline(h=f[2], col="red", lwd=lwd.tif)
mtext("D", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Bias in population growth rate ('*lambda*')')
k1 <- c(as.vector(diff.ipm1), as.vector(diff.ipm2), as.vector(diff.ipm3))
k2 <- c(rep(1, length(as.vector(diff.ipm1))), rep(2, length(as.vector(diff.ipm2))),
        rep(3, length(as.vector(diff.ipm3))))
boxplot(k1 ~ k2, outline=FALSE, ylab=laby, xlab=NA, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=labx, lwd=lwd.tif)
abline(h=0, col="red", lwd=lwd.tif)
mtext("E", at=0.5, line=0.5, cex=1.5 * cex.tif)

laby <- expression('Residual error ('*sigma*')')
boxplot(cbind(res1[46,1,i1], res2[46,1,i2], res3[128,1,i3]), outline=FALSE,
        ylab=laby, axes=FALSE, col=co, lwd=lwd.tif)
axis(2, las=1, lwd=lwd.tif)
axis(1, at=c(1, 2, 3), labels=labx, lwd=lwd.tif)
abline(h=sigma, col="red", lwd=lwd.tif)
mtext("F", at=0.5, line=0.5, cex=1.5 * cex.tif)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

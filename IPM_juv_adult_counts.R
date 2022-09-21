list.files("/home/imb/tbeauchard/Bureau/StageCNRS/CodeJAGS/IPMSchaubKery")
library(IPMbook) ; library(jagsUI)
list.files(full.names = TRUE)
dirname("ipm4.txt")
dirname("ipm1.txt")

# 6.3 Estimation of demographic parameters for which there is no explicit data
# ============================================================================

# Additional code for the simulations for Figure 6.4

# Models
# ------

# 1. IPM all data sets

# Specify the model in BUGS language
cat(file="ipm1.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }
  sigma.obs ~ dunif(0.5, 100)
  tau.obs <- pow(sigma.obs, -2)
  # State-space model for count data
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)
  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f/2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }
  # Observation model
  for (t in 1:n.occasions){
    C_j[t] ~ dnorm(N[1,t], tau.obs)
    C_ad[t] ~ dnorm(N[2,t], tau.obs)
  }
  # Poisson regression model for productivity data
  sJ ~ dpois(nJ * mean.f)
  # Capture-recapture model (multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
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
  } #t
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

# 2. IPM: count + CMR

# Specify the model in BUGS language
cat(file="ipm2.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }
  sigma.obs ~ dunif(0.5, 100)
  tau.obs <- pow(sigma.obs, -2)
  # State-space model for count data
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)
  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f/2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }
  # Observation model
  for (t in 1:n.occasions){
    C_j[t] ~ dnorm(N[1,t], tau.obs)
    C_ad[t] ~ dnorm(N[2,t], tau.obs)
  }
  # Capture-recapture model (multinomial likelihood)
  # Define the multinomial likelihood
  for (t in 1:(n.occasions-1)){
    marr.j[t,1:n.occasions] ~ dmulti(pr.j[t,], rel.j[t])
    marr.a[t,1:n.occasions] ~ dmulti(pr.a[t,], rel.a[t])
  }
  # Define the cell probabilities of the m-arrays
  for (t in 1:(n.occasions-1)){
    # Main diagonal
    q[t] <- 1-p[t]   # Probability of non-recapture
    pr.j[t,t] <- sj[t]*p[t]
    pr.a[t,t] <- sa[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(n.occasions-1)){
      pr.j[t,j] <- sj[t]*prod(sa[(t+1):j])*prod(q[t:(j-1)])*p[j]
      pr.a[t,j] <- prod(sa[t:j])*prod(q[t:(j-1)])*p[j]
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
  } #t
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

# 3. IPM: counts + productivity

# Specify the model in BUGS language
cat(file="ipm3.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)
  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }
  sigma.obs ~ dunif(0.5, 100)
  tau.obs <- pow(sigma.obs, -2)
  # State-space model for count data
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)
  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f/2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }
  # Observation model
  for (t in 1:n.occasions){
    C_j[t] ~ dnorm(N[1,t], tau.obs)
    C_ad[t] ~ dnorm(N[2,t], tau.obs)
  }
  # Poisson regression model for productivity data
  sJ ~ dpois(nJ * mean.f)
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


# 4. IPM: counts only

# Specify the model in BUGS language
cat(file="ipm4.txt", "
model {
  # Priors and constraints
  mean.sj ~ dunif(0, 1)
  mean.sa ~ dunif(0, 1)
  mean.p ~ dunif(0, 1)
  mean.f ~ dunif(0, 10)

  for (t in 1:(n.occasions-1)){
    sj[t] <- mean.sj
    sa[t] <- mean.sa
    p[t] <- mean.p
  }

  sigma.obs ~ dunif(0.5, 100)
  tau.obs <- pow(sigma.obs, -2)

  # State-space model for count data
  # Model for the initial population size: discrete uniform priors
  N[1,1] ~ dunif(1, 300)
  N[2,1] ~ dunif(1, 300)

  # Process model over time
  for (t in 1:(n.occasions-1)){
    N[1,t+1] <- mean.f/2 * mean.sj * (N[1,t] + N[2,t])
    N[2,t+1] <- mean.sa * (N[1,t] + N[2,t])
  }

  # Observation model
  for (t in 1:n.occasions){
   C_j[t] ~ dnorm(N[1,t], tau.obs)
    C_ad[t] ~ dnorm(N[2,t], tau.obs)
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

# Simulation parameters
# ---------------------
# Number of simulations
 nsim <- 1500
#nsim <- 3# ~~~ for testing

# Age specific survival probabilities (juv, adult)
sj <- 0.3
sa <- 0.55

# Fecundity rate (females)
fl1 <- 3.1            # productivity of one year old females
fl2 <- 3.1            # productivity of females older than one year

# Initial population size per age class
Ni <- c(50, 50)

# Number of years
T <- 10

# Observation error for the population survey
sigma <- 10

# Capture and recapture probabilities
cap <- 0.4                     # initial capture probability
prec <- 0.6                    # recapture probability

# Probability to find a brood whose reproductive ouput is recorded
pprod <- 0.5

# Initial values
inits.ipm <- function(){list(mean.sj=runif(1, 0.2, 0.4), mean.sa=runif(1, 0.45, 0.65))}

# Parameters monitored
parameters.ipm <- c("mean.sj", "mean.sa", "mean.p", "mean.f", "N", "sigma.obs",
                    "ann.growth.rate", "Ntot")

# Define matrices to store the result
res1DD <- res2DD <- res3DD <- res4DD <- array(NA, dim = c(45, 11, nsim))
#res4D  <- array(NA, dim = c(64, 11, nsim))

# Start simulations
system.time(
  for (s in 1:nsim){
    
    set.seed(s)
    
    ind1 <- simPop(phi=c(sj, sa), f=c(fl1, fl2), nYears=T, sex.ratio=0.5, Im=0, Ni=Ni)
    ind2 <- simPop(phi=c(sj, sa), f=c(fl1, fl2), nYears=T, sex.ratio=0.5, Im=0, Ni=Ni)
    ind3 <- simPop(phi=c(sj, sa), f=c(fl1, fl2), nYears=T, sex.ratio=0.5, Im=0, Ni=Ni)
    
    # Create the population survey data
   countJ <- simCountNorm(ind1$breeders[1,], sigma)$count
   countA <- simCountNorm(ind1$breeders[2,], sigma)$count
    
    
    
    # Create the capture histories and the corresponding m-arrays
    ch <- simCapHist(ind2$state, cap=cap, recap=prec, maxAge=2)
    marr <- marrayAge(ch$ch, ch$age)
    
    # Create productivity data
    P <- simProd(ind3$reprod, pprod)
    # Aggregate productivity data to make the model run faster
    sJ <- colSums(P$prod.agg)[1]
    nJ <- colSums(P$prod.agg)[2]
    
    # Bundle data
    jags.data.ipm1 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
                           rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), sJ=sJ, nJ=nJ, C_j=countJ, C_ad=countA)
    jags.data.ipm2 <- list(marr.j=marr[,,1], marr.a=marr[,,2], n.occasions=T,
                           rel.j=rowSums(marr[,,1]), rel.a=rowSums(marr[,,2]), C_j=countJ, C_ad=countA)
    jags.data.ipm3 <- list(n.occasions=T, sJ=sJ, nJ=nJ, C_j=countJ, C_ad=countA)
    jags.data.ipm4 <- list(n.occasions=T, C_j=countJ, C_ad=countA)
    
    # Call JAGS from R (jagsUI)
    # MCMC settings
    ni <- 10000; nt <- 1; nb <- 5000; nc <- 3; na <- 1000
    
    m1 <- try(jags(jags.data.ipm1, inits.ipm, parameters.ipm, "ipm1.txt",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE))
    if (!inherits(m1, "try-error"))
      res1D[,,s] <- m1$summary
    
    m2 <- try(jags(jags.data.ipm2, inits.ipm, parameters.ipm, "ipm2.txt",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE))
    if (!inherits(m2, "try-error"))
      res2D[,,s] <- m2$summary
    
    ni <- 200000; nt <- 1; nb <- 50000; nc <- 3
    
    m3 <- try(jags(jags.data.ipm3, inits.ipm, parameters.ipm, "ipm3.txt",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE))
    if (!inherits(m3,  "try-error"))
      res3D[,,s] <- m3$summary
  
    m4 <- try(jags(jags.data.ipm4, inits.ipm, parameters.ipm, "ipm4.txt",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na, parallel = TRUE))
    if (!inherits(m4, "try-error"))
      res4D[,,s] <- m4$summary
  
    print(s)
  }  )  # 3 sims took 90 secsm 1500 took 12 hrs
#exemple de sortie
print(m1)
#sauvegarde des matrices de résumés
save(res1D, res2D, res3D, res4D, sj, sa, fl1, fl2, sigma, prec, file="Data Fig 6.4.1500simJAD.20092022.Rdata")

#vérification de la convergence des chaînes MCMC
#traceplot(m1)
#traceplot(m2)
#traceplot(m3)
#traceplot(m4)

load("Data Fig 6.4.1500simJAD.20092022.Rdata")
m1$summary
m2$summary
m3$summary
m4$summary
# Select only the 1500 simulations that have converged
#si jamais res1D ne fonctionne pas, faire res1DD
inclD <- which(res1D[1,8,]<1.05 & res1D[2,8,]<1.05 & res1D[4,8,]<1.05 & res1D[34,8,]<1.05 &
                res2D[1,8,]<1.05 & res2D[2,8,]<1.05 & res2D[4,8,]<1.05 & res2D[34,8,]<1.05 &
                res3D[1,8,]<1.05 & res3D[2,8,]<1.05 & res3D[4,8,]<1.05 & res3D[34,8,]<1.05 &
                res4D[1,8,]<1.05 & res4D[2,8,]<1.05 & res4D[4,8,]<1.05 & res4D[34,8,]<1.05)
if(length(inclD) > 1500)
  inclD <- inclD[1:1500]

op <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name <- c("CR & P & C", "CR & C", "P & C", "C")

lab <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1DD[1,1,inclD], res2D[1,1,inclD], res3D[1,1,inclD], res4D[1,1,inclD]),
        ylab=lab, outline=FALSE, names=NA,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1DD[1,2,inclD], res2D[1,2,inclD], res3D[1,2,inclD], res4D[1,2,inclD]),
        ylab=lab, outline=FALSE, names=NA,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1DD[2,1,inclD], res2D[2,1,inclD], res3D[2,1,inclD], res4D[2,1,inclD]),
        ylab=lab, outline=FALSE, names=NA,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1DD[2,2,inclD], res2D[2,2,inclD], res3D[2,2,inclD], res4D[2,2,inclD]),
        ylab=lab, outline=FALSE, names=NA,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab <- expression('Productivity ('*italic('f')*')')
boxplot(cbind(res1DD[4,1,inclD], res2D[4,1,inclD], res3D[4,1,inclD], res4D[4,1,inclD]),
        ylab=lab, outline=FALSE, names=NA,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab <- expression('SD ('*italic('f')*')')
boxplot(cbind(res1DD[4,2,inclD], res2D[4,2,inclD], res3D[4,2,inclD], res4D[4,2,inclD]),
        ylab=lab, outline=FALSE, names=NA,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab <- expression('Population growth rate ('*lambda*')')
boxplot(cbind(res1D[34,1,inclD], res2D[34,1,inclD], res3D[34,1,inclD], res4D[34,1,inclD]),
        ylab=lab,  outline= FALSE, names=name,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
axis(1, labels=name, at=1:4)
axis(2)

lab <- expression('SD ('*lambda*')')
boxplot(cbind(res1D[34,2,inclD], res2D[34,2,inclD], res3D[34,2,inclD], res4D[34,2,inclD]),
        ylab=lab,  outline=FALSE, names=name,
        col=c("red", "dodgerblue", "darkolivegreen", "orange"), axes=FALSE)
axis(1, labels=name, at=1:4)
axis(2)
par(op)


#graphique isolant C & CR

op3 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name3 <- c("CR & C")

lab3 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(res2D[1,1,inclD],
        ylab=lab3, outline=FALSE, names=NA,
        col=c("dodgerblue"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab3 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(res2D[1,2,inclD],
        ylab=lab3, outline=FALSE, names=NA,
        col=c("dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab3 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(res2D[2,1,inclD],
        ylab=lab3, outline=FALSE, names=NA,
        col=c("dodgerblue"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab3 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(res2D[2,2,inclD],
        ylab=lab3, outline=FALSE, names=NA,
        col=c("dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab3 <- expression('Productivity ('*italic('f')*')')
boxplot(res2D[4,1,inclD],
        ylab=lab3, outline=FALSE, names=NA,
        col=c("dodgerblue"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab3 <- expression('SD ('*italic('f')*')')
boxplot(res2D[4,2,inclD],
        ylab=lab3, outline=FALSE, names=NA,
        col=c("dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab3 <- expression('Population growth rate ('*lambda*')')
boxplot(res2D[34,1,inclD],
        ylab=lab3,  outline= FALSE, names=name3,
        col=c("dodgerblue"), axes=FALSE)
axis(1, labels=name3, at=1)
axis(2)

lab3 <- expression('SD ('*lambda*')')
boxplot(res2D[34,2,inclD],
        ylab=lab3,  outline=FALSE, names=name3,
        col=c("dodgerblue"), axes=FALSE)
axis(1, labels=name3, at=1)
axis(2)
par(op3)


#graphique isolant C & P

op4 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name4 <- c("C & P")

lab4 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(res3D[1,1,inclD],
        ylab=lab4, outline=FALSE, names=NA,
        col=c("darkolivegreen"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab4 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(res3D[1,2,inclD],
        ylab=lab4, outline=FALSE, names=NA,
        col=c("darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab4 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(res3D[2,1,inclD],
        ylab=lab4, outline=FALSE, names=NA,
        col=c("darkolivegreen"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab4 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(res3D[2,2,inclD],
        ylab=lab4, outline=FALSE, names=NA,
        col=c("darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab4 <- expression('Productivity ('*italic('f')*')')
boxplot(res3D[4,1,inclD],
        ylab=lab4, outline=FALSE, names=NA,
        col=c("darkolivegreen"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab4 <- expression('SD ('*italic('f')*')')
boxplot(res3D[4,2,inclD],
        ylab=lab4, outline=FALSE, names=NA,
        col=c("darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab4 <- expression('Population growth rate ('*lambda*')')
boxplot(res3D[34,1,inclD],
        ylab=lab4,  outline= FALSE, names=name4,
        col=c("darkolivegreen"), axes=FALSE)
axis(1, labels=name4, at=1)
axis(2)

lab4 <- expression('SD ('*lambda*')')
boxplot(res3D[34,2,inclD],
        ylab=lab4,  outline=FALSE, names=name4,
        col=c("darkolivegreen"), axes=FALSE)
axis(1, labels=name4, at=1)
axis(2)
par(op4)

#graphique isolant CR & P & C

op6 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name6 <- c("CR & P & C")

lab6 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(res1D[1,1,inclD],
        ylab=lab6, outline=FALSE, names=NA,
        col=c("red"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab6 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(res1D[1,2,inclD], 
        ylab=lab6, outline=FALSE, names=NA,
        col=c("red"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab6 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(res1D[2,1,inclD],
        ylab=lab6, outline=FALSE, names=NA,
        col=c("red"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab6 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(res1D[2,2,inclD],
        ylab=lab6, outline=FALSE, names=NA,
        col=c("red"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab6 <- expression('Productivity ('*italic('f')*')')
boxplot(res1D[4,1,inclD],
        ylab=lab6, outline=FALSE, names=NA,
        col=c("red"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab6 <- expression('SD ('*italic('f')*')')
boxplot(res1D[4,2,inclD], 
        ylab=lab6, outline=FALSE, names=NA,
        col=c("red"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab6 <- expression('Population growth rate ('*lambda*')')
boxplot(res1D[34,1,inclD],
        ylab=lab6,  outline= FALSE, names=name6,
        col=c("red"), axes=FALSE)
axis(1, labels=name6, at=1)
axis(2)

lab6 <- expression('SD ('*lambda*')')
boxplot(res1D[34,2,inclD],
        ylab=lab6,  outline=FALSE, names=name6,
        col=c("red"), axes=FALSE)
axis(1, labels=name6, at=1)
axis(2)

par(op6)


#test combinaison de boxplots ("CR & P & C" sans dissocier,"C & CR" en dissociant)
op7 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name7 <- c("CR & P & C","C & CR")


lab7 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,1,incl],res2D[1,1,inclD]),
        ylab=lab7, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab7 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,2,incl],res2D[1,2,inclD]),
        ylab=lab7, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab7 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,1,incl],res2D[2,1,inclD]),
        ylab=lab7, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab7 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,2,incl],res2D[2,2,inclD]),
        ylab=lab7, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab7 <- expression('Productivity ('*italic('f')*')')
boxplot(cbind(res1[4,1,incl],res2D[4,1,inclD]),
        ylab=lab7, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab7 <- expression('SD ('*italic('f')*')')
boxplot(cbind(res1[4,2,incl],res2D[4,2,inclD]),
        ylab=lab7, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab7 <- expression('Population growth rate ('*lambda*')')
boxplot(cbind(res1[34,1,incl],res2D[34,1,inclD]),
        ylab=lab7,  outline= FALSE, names=name7,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=name7, at=1:2)
axis(2)

lab7 <- expression('SD ('*lambda*')')
boxplot(cbind(res1[34,2,incl],res2D[34,2,inclD]),
        ylab=lab7,  outline=FALSE, names=name7,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=name7, at=1:2)
axis(2)
par(op7)



#test combinaison de boxplots ("CR & P & C" sans dissocier,"C & P" en dissociant)
op8 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name8 <- c("CR & P & C","C & P")


lab8 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,1,incl],res3D[1,1,inclD]),
        ylab=lab8, outline=FALSE, names=NA,
        col=c("red","darkolivegreen"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab8 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,2,incl],res3D[1,2,inclD]),
        ylab=lab8, outline=FALSE, names=NA,
        col=c("red","darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab8 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,1,incl],res3D[2,1,inclD]),
        ylab=lab8, outline=FALSE, names=NA,
        col=c("red","darkolivegreen"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab8 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,2,incl],res3D[2,2,inclD]),
        ylab=lab8, outline=FALSE, names=NA,
        col=c("red","darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab8 <- expression('Productivity ('*italic('f')*')')
boxplot(cbind(res1[4,1,incl],res3D[4,1,inclD]),
        ylab=lab8, outline=FALSE, names=NA,
        col=c("red","darkolivegreen"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab8 <- expression('SD ('*italic('f')*')')
boxplot(cbind(res1[4,2,incl],res3D[4,2,inclD]),
        ylab=lab8, outline=FALSE, names=NA,
        col=c("red","darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab8 <- expression('Population growth rate ('*lambda*')')
boxplot(cbind(res1[34,1,incl],res3D[34,1,inclD]),
        ylab=lab8,  outline= FALSE, names=name8,
        col=c("red","darkolivegreen"), axes=FALSE)
axis(1, labels=name8, at=1:2)
axis(2)

lab8 <- expression('SD ('*lambda*')')
boxplot(cbind(res1[34,2,incl],res3D[34,2,inclD]),
        ylab=lab8,  outline=FALSE, names=name8,
        col=c("red","darkolivegreen"), axes=FALSE)
axis(1, labels=name8, at=1:2)
axis(2)
par(op8)

#test combinaison de boxplots ("CR & P & C" en dissociant,"C & CR" en dissociant)
op9 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name9 <- c("CR & P & C","C & CR")


lab9 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1D[1,1,inclD],res2D[1,1,inclD]),
        ylab=lab9, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab9 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1D[1,2,incl],res2D[1,2,inclD]),
        ylab=lab9, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab9 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1D[2,1,incl],res2D[2,1,inclD]),
        ylab=lab9, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab9 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1D[2,2,incl],res2D[2,2,inclD]),
        ylab=lab9, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab9 <- expression('Productivity ('*italic('f')*')')
boxplot(cbind(res1D[4,1,incl],res2D[4,1,inclD]),
        ylab=lab9, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab9 <- expression('SD ('*italic('f')*')')
boxplot(cbind(res1D[4,2,incl],res2D[4,2,inclD]),
        ylab=lab9, outline=FALSE, names=NA,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab9 <- expression('Population growth rate ('*lambda*')')
boxplot(cbind(res1D[34,1,incl],res2D[34,1,inclD]),
        ylab=lab9,  outline= FALSE, names=name9,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=name9, at=1:2)
axis(2)

lab9 <- expression('SD ('*lambda*')')
boxplot(cbind(res1D[34,2,incl],res2D[34,2,inclD]),
        ylab=lab9,  outline=FALSE, names=name9,
        col=c("red","dodgerblue"), axes=FALSE)
axis(1, labels=name9, at=1:2)
axis(2)
par(op9)


#test combinaison de boxplots ("CR & P & C" sans dissocier,"C & P" en dissociant, "C & CR" en dissociant)
op11 <- par(las=1, mar=c(2.5, 4.2, 1, 1), mfrow=c(4, 2))
name11 <- c("CR & P & C", "C & CR", "C & P")


lab11 <- expression('Juvenile survival ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,1,incl],res2D[1,1,inclD],res3D[1,1,inclD]),
        ylab=lab11, outline=FALSE, names=NA,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
abline(h=sj, lty=2)
axis(1, labels=NA)
axis(2)

lab11 <- expression('SD ('*italic('s')[italic(j)]*')')
boxplot(cbind(res1[1,2,incl],res2D[1,2,inclD],res3D[1,2,inclD]),
        ylab=lab11, outline=FALSE, names=NA,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab11 <- expression('Adult survival ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,1,incl],res2D[2,1,inclD],res3D[2,1,inclD]),
        ylab=lab11, outline=FALSE, names=NA,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
abline(h=sa, lty=2)
axis(1, labels=NA)
axis(2)

lab11 <- expression('SD ('*italic('s')[italic(a)]*')')
boxplot(cbind(res1[2,2,incl],res2D[2,2,inclD],res3D[2,2,inclD]),
        ylab=lab11, outline=FALSE, names=NA,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab11 <- expression('Productivity ('*italic('f')*')')
boxplot(cbind(res1[4,1,incl],res2D[4,1,inclD],res3D[4,1,inclD]),
        ylab=lab11, outline=FALSE, names=NA,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
abline(h=fl1, lty=2)
axis(1, labels=NA)
axis(2)

lab11 <- expression('SD ('*italic('f')*')')
boxplot(cbind(res1[4,2,incl],res2D[4,2,inclD],res3D[4,2,inclD]),
        ylab=lab11, outline=FALSE, names=NA,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
axis(1, labels=NA)
axis(2)

lab11 <- expression('Population growth rate ('*lambda*')')
boxplot(cbind(res1[34,1,incl],res2D[34,1,inclD],res3D[34,1,inclD]),
        ylab=lab11,  outline= FALSE, names=name11,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
axis(1, labels=name11, at=1:3)
axis(2)

lab11 <- expression('SD ('*lambda*')')
boxplot(cbind(res1[34,2,incl],res2D[34,2,inclD],res3D[34,2,inclD]),
        ylab=lab11,  outline=FALSE, names=name11,
        col=c("red","dodgerblue","darkolivegreen"), axes=FALSE)
axis(1, labels=name11, at=1:3)
axis(2)
par(op11)


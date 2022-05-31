# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------
# Code from final MS.
library(IPMbook)

# 3.1 Introduction
# ================

# ~~~~ Create Figure 3.1 ~~~~
T <- 9                              # Number of years
N <- matrix(NA, nrow=T+1, ncol=3)   # Matrix to hold size of 3 populations
lambda <- c(0.8, 1, 1.2)            # Pop. growth rate of 3 pops
N[1,] <- 10                         # Initial population size is 10
for (t in 1:T){                     # Loop over years
  for (k in 1:3){                   # Loop over 3 populations
    N[t+1,k] <- lambda[k] * N[t,k]  # Project pop. one time step forward
  } #k
} #t

op <- par(cex=1.4)
plot(N[,1], type="l", ylim=range(N), ylab="Population size", xlab="Year", axes=FALSE, lwd=2.5)
lines(N[,2], lwd=2.5)
lines(N[,3], lwd=2.5)
axis(1, at=1:(T+1), labels=1:(T+1))
axis(2, las=1)
text(x=T, y=N[5,1], expression(paste(lambda, "=0.8")))
text(x=T, y=N[10,2] + 3, expression(paste(lambda, "=1.0")))
text(x=T, y=N[8,3], expression(paste(lambda, "=1.2")))
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------
# Code from final MS.

# 3.2 Age- and stage-structured population models
# ===============================================

# 3.2.1 From a life-cycle graph to population equations
# -----------------------------------------------------

# ~~~~ Code for Fig. 3.3 ~~~~
sj <- 0.3      # Juvenile survival
sa <- 0.55     # Adult survival
f1 <- 1.3      # Number of female fledglings per 1-year old female and year
fa <- 1.8      # Number of female fledglings per adult female and year

N <- numeric(14)
N[1] <- 10
for (i in 2:14){
  if (i %% 2 == 0){   # even numbers
    N[i] <- N[i-1] / 2 * f1 + N[i-1] / 2 * fa
  } else {               # odd numbers
    N[i] <- N[i-2] * sa + N[i-1] * sj
  }
}
x <- c(rbind(1:7, 1:7 + 0.1))

plot(x=x, y=N, type="l", ylab="Total population size", ylim=c(7, 18),
     col="blue", lwd=2, axes=FALSE, xlab="Year")
axis(2, las=1)
axis(1)
points(x, N, pch=c(22, 1), cex=2, col=c("red", "darkgreen"), lwd=2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3.2.2 Age-structured, pre-birth-pulse model (no code)

# 3.2.3 Stage-structured, pre-birth-pulse model (no code)

# 3.2.4 Age-structured, post-birth-pulse model (no code)

# 3.2.5 Stage-structured, post-birth-pulse model (no code)
# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.1 Analysis of a matrix population model without stochasticity
#       and parameter uncertainty
# -----------------------------------------------------------------

# Define the demographic rates
sj <- 0.3                                           # Juvenile survival
sa <- 0.55                                          # Adult survival
f1 <- 1.3                                           # Number of female fledglings per 1-year old female
fa <- 1.8                                           # Number of female fledglings per adult female

# Define the transition matrix A (inspired by woodchat shrike demography)
A <- matrix(c(f1 * sj, fa * sj, sa, sa), ncol=2, byrow=TRUE)
dimnames(A) <- list(c("1y", "ad"), c("1y", "ad"))
A                                                   # print the matrix
#      1y   ad
# 1y 0.39 0.54
# ad 0.55 0.55

# Asymptotic population growth rate:
# ''''''''''''''''''''''''''''''''''

N1 <- c(10, 1)

T <- 7
N <- matrix(NA, nrow=2, ncol=T)                     # Population sizes
gr <- matrix(NA, nrow=3, ncol=T-1)                  # Growth rates
N[,1] <- N1
for (t in 2:T){
  N[,t] <- A %*% N[,t-1]
  gr[1,t-1] <- N[1,t] / N[1,t-1]
  gr[2,t-1] <- N[2,t] / N[2,t-1]
  gr[3,t-1] <- (N[1,t] + N[2,t]) / (N[1,t-1] + N[2,t-1])
}
sr <- N[1, ] / colSums(N)                           # State distribution (proportion of 1y)

# ~~~~ extra code for Figure 3.8 ~~~~
op <- par(las=1, "mfrow")
layout(matrix(1:3, 1, 3,byrow=TRUE), widths=c(1, 1, 1), heights=1, TRUE)
# Plot stage-specific population size
plot(N[1,], type="b", pch=16, ylim=range(c(min(N), colSums(N))), axes=FALSE,
     ylab="Population size", xlab="Time", col="red")
points(N[2,], type="b", pch=16, col="orange")
points(colSums(N), type="b", pch=16, col="blue")
axis(1)
axis(2)
mtext("A", at=1, line=0.5, cex=1.5)
legend("bottomright", legend=c("1y", "ad", "total"), bty="n", pch=rep(16, 3),
       lty=rep(1, 3), col=c("red", "orange", "blue"))

# Plot stage-specific growth rates
plot(x=(1:(T-1)) + 0.5, y=gr[1,], type="b", pch=16, ylim=range(gr),
     xlim=c(1, T), axes=FALSE, ylab="Population growth rate", xlab="Time", col="red")
points(x=(1:(T-1)) + 0.5, y=gr[2,], type="b", pch=16, col="orange")
points(x=(1:(T-1)) + 0.5, y=gr[3,], type="b", pch=16, col="blue")
axis(1, at=1:T, labels=1:T)
axis(2)
abline(h=max(Re(eigen(A)$values)), lty=2)
mtext("B", at=1, line=0.5, cex=1.5)

# Plot actual stage distribution
plot(sr, type="b", pch=16, ylab="Proportion of stage classes", axes=FALSE,
     col="red", xlab="Time", ylim=range(c(sr, 1-sr)))
points(1-sr, type="b", pch=16, col="orange")
axis(1)
axis(2)
u <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,u])
abline(h=revec[1]/sum(revec), col="red", lty=2)
abline(h=revec[2]/sum(revec), col="orange", lty=2)
mtext("C", at=1, line=0.5, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute asymptotic population growth rate (lambda)
lambda <- max(Re(eigen(A)$values))
# [1] 1.020818

# ~~~~ extra code for Figure 3.9 ~~~~
# Define the demographic rates
s1 <- 0.5     # Juvenile survival
s2 <- 0.8     # 1y survival
s3 <- 0.9     # Ad survival
rho1 <- 0.5   # Productivity of 2-years old females
rho2 <- 0.9   # Productivity of adult females

# Define the transition matrix A
Abat <- matrix(c(0, rho1*s1/2, rho2*s1/2,
                 s2, 0, 0,
                 0, s3, s3), ncol=3, byrow=TRUE)

# Define R structures and do computations
T <- 9
N1 <- c(10, 5, 5)                   # Initial population sizes
N <- matrix(NA, nrow=3, ncol=T)     # Population size
gr <- matrix(NA, nrow=4, ncol=T-1)  # Growth rate
N[,1] <- N1
for (t in 2:T){
  N[,t] <- Abat %*% N[,t-1]  # stage-specific population sizes
  gr[1,t-1] <- N[1,t] / N[1,t-1]
  gr[2,t-1] <- N[2,t] / N[2,t-1]
  gr[3,t-1] <- N[3,t] / N[3,t-1]
  gr[4,t-1] <- (N[1,t] + N[2,t] + N[3,t]) / (N[1,t-1] + N[2,t-1] + N[3,t-1])
}
sr <- sweep(N, 2, colSums(N), "/")   # Stage distribution

# Plots
co <- c("red", "orange", "tomato2")
op <- par(las=1, cex=1.1, "mfrow")
layout(matrix(1:3, 1, 3, byrow=TRUE), widths=c(1, 1, 1), heights=1, TRUE)
plot(N[1,], type="b", pch=16, ylim=range(c(min(N), colSums(N))), axes=FALSE,
     ylab="Population size", xlab="Time", col=co[1])
points(N[2,], type="b", pch=16, col=co[2])
points(N[3,], type="b", pch=16, col=co[3])
points(colSums(N), type="b", pch=16, col="blue")
axis(1, at=1:T)
axis(2)
mtext("A", at=1, line=0.5, cex=1.5)

plot(x=(1:(T-1)) + 0.5, y=gr[1,], type="b", pch=16, ylim=range(gr), xlim=c(1, T),
     axes=FALSE, ylab="Population growth rate", xlab="Time", col=co[1])
points(x=(1:(T-1)) + 0.5, y=gr[2,], type="b", pch=16, col=co[2])
points(x=(1:(T-1)) + 0.5, y=gr[3,], type="b", pch=16, col=co[3])
points(x=(1:(T-1)) + 0.5, y=gr[4,], type="b", pch=16, col="blue")
axis(1, at=1:T)
axis(2)
abline(h=max(Re(eigen(Abat)$values)), lty=2)
mtext("B", at=1, line=0.5, cex=1.5)
legend("topright", legend=c("1y", "2y", "ad", "total"), bty="n", pch=rep(16,4),
       lty=rep(1, 4), col=c(co, "blue"))

plot(sr[1,], type="b", pch=16, ylab="Proportion of stage classes", axes=FALSE,
     col=co[1], ylim=c(0, 0.8), xlab="Time")
points(sr[2,], type="b", pch=16, col=co[2])
points(sr[3,], type="b", pch=16, col=co[3])
axis(1, at=1:T)
axis(2)
u <- which.max(Re(eigen(Abat)$values))
revec <- Re(eigen(Abat)$vectors[,u])
abline(h=revec[1]/sum(revec), col=co[1], lty=2)
abline(h=revec[2]/sum(revec), col=co[2], lty=2)
abline(h=revec[3]/sum(revec), col=co[3], lty=2)
mtext("C", at=1, line=0.5, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Stable stage distribution:
# ''''''''''''''''''''''''''

u <- which.max(Re(eigen(A)$values))
# [1] 1
revec <- Re(eigen(A)$vectors[,u])
# [1] -0.6503047 -0.7596734

revec/sum(revec)
# [1] 0.4612162 0.5387838

# Stage-specific reproductive values:
# '''''''''''''''''''''''''''''''''''

u <- which.max(Re(eigen(A)$values))
levec <- Re((eigen(A)$vectors)[u,])
levec / sum(levec)
# [1] 0.4631655 0.5368345

# Net reproductive rate:
# ''''''''''''''''''''''

i <- 1:100                                          # 100 as our approximation to infinity
R0 <- sj * f1 + sj * fa * sum(sa^i)
# [1] 1.05

# Generation time
# '''''''''''''''

Q <- 100                                            # Our approximation to infinity
i <- 2:Q
G <- sj * f1 / lambda + sj * fa * sum(i * sa^(i-1) * lambda^(-i))
# [1] 2.339834

GT <- log(R0) / log(lambda)
# [1] 2.368012


# Sensitivity and elasticity (perturbation analysis):
# '''''''''''''''''''''''''''''''''''''''''''''''''''

senmat <- levec %*% t(revec)
#           [,1]      [,2]
# [1,] 0.4228963 0.4940192
# [2,] 0.4901602 0.5725957

elasmat <- senmat * A / lambda
#           1y        ad
# 1y 0.1615661 0.2613301
# ad 0.2640904 0.3085053

derivmat <- matrix(c(f1, fa, 0, 0), ncol=2, byrow=TRUE)
ssj <- sum(senmat * derivmat)                       # Sensitivity
esj <- sj / lambda * ssj                            # Elasticity
ssj; esj
# [1] 1.439
# [1] 0.4228963

# ~~~~ extra code for Figure 3.10 ~~~~
# Define the demographic rates
sj <- 0.3      # Juvenile survival
sa <- 0.55     # Adult survival
f1 <- 1.3      # Number of female fledglings per 1-year old female and year
fa <- 1.8      # Number of female fledglings per adult female and year

# Sensitivity
ds <- seq(-0.1, 0.1, length.out=501)
lams <- matrix(NA, ncol=4, nrow=501)
for (i in 1:501){
  A <- matrix(c((f1+ds[i])*sj, fa*sj, sa, sa), ncol=2, byrow=TRUE)
  lams[i,1] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, (fa+ds[i])*sj, sa, sa), ncol=2, byrow=TRUE)
  lams[i,2] <- max(eigen(A)$values)
  A <- matrix(c(f1*(sj+ds[i]), fa*(sj+ds[i]), sa, sa), ncol=2, byrow=TRUE)
  lams[i,3] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, fa*sj, sa+ds[i], sa+ds[i]), ncol=2, byrow=TRUE)
  lams[i,4] <- max(eigen(A)$values)
}

# Elasticity
es <- seq(0.9, 1.1, length.out=501)
lame <- matrix(NA, ncol=4, nrow=501)
for (i in 1:501){
  A <- matrix(c(f1*es[i]*sj, fa*sj, sa, sa), ncol=2, byrow=TRUE)
  lame[i,1] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, fa*es[i]*sj, sa, sa), ncol=2, byrow=TRUE)
  lame[i,2] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj*es[i], fa*sj*es[i], sa, sa), ncol=2, byrow=TRUE)
  lame[i,3] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, fa*sj, sa*es[i], sa*es[i]), ncol=2, byrow=TRUE)
  lame[i,4] <- max(eigen(A)$values)
}

# Plots
op <- par("mfrow", "las")
layout(matrix(1:2, 1, 2, byrow=TRUE), widths=c(1, 1), heights=1, TRUE)
co <- c("red", "orange", "darkgreen", "blue")
lw <- 2
plot(y=lams[,1], x=ds, type="l", ylim=range(cbind(lams, lame)), axes=FALSE,
     ylab="Population growth rate", xlab="Absolute change in demographic rate",
     col=co[1], lwd=lw)
axis(1)
axis(1, at=c(-0.075, -0.025, 0.025, 0.075), labels=NA, tcl=-0.25)
axis(2, las=1)
lines(y=lams[,2], x=ds, col=co[2], lwd=lw)
lines(y=lams[,3], x=ds, col=co[3], lwd=lw)
lines(y=lams[,4], x=ds, col=co[4], lwd=lw)
mtext("A", at=-0.1, line=0.5, cex=1.5)
legend("topleft", lwd=rep(lw, 4), col=co, legend=c(expression(italic(f)[1]),
                                                   expression(italic(f)[italic(a)]), expression(italic(s)[italic(j)]),
                                                   expression(italic(s)[italic(a)])), bty="n")

plot(y=lame[,1], x=es, type="l", ylim=range(cbind(lams, lame)), axes=FALSE,
     ylab="", xlab="Relative change in demographic rate (%)", col=co[1], lwd=lw)
axis(1, at=c(0.925, 0.975, 1.025, 1.075), labels=NA, tcl=-0.25)
axis(1, at=c(0.9, 0.95, 1, 1.05, 1.1), labels=c(-10, -5, 0, 5, 10))
axis(2, las=1)
lines(y=lame[,2], x=es, col=co[2], lwd=lw)
lines(y=lame[,3], x=es, col=co[3], lwd=lw)
lines(y=lame[,4], x=es, col=co[4], lwd=lw)
mtext("B", at=0.9, line=0.5, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook)

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.2 Analysis of a matrix population model with parameter uncertainty
# ----------------------------------------------------------------------

# Define mean and SD of the demographic parameters
mean.sj <- 0.3                  # Point estimate of juv. survival
se.sj.e <- 0.03                 # Uncertainty of juv. survival as SE on natural scale
mean.sa <- 0.55                 # Point estimate of ad. survival
se.sa.e <- 0.015                # Uncertainty of ad. survival as SE on natural scale
mean.f1 <- 1.3                  # Point estimate of productivity of 1y females
se.f1.e <- 0.3                  # Uncertainty of 1y productivity as SE on natural scale
mean.fa <- 1.8                  # Point estimate of productivity of ad. females
se.fa.e <- 0.1                  # Uncertainty of ad. productivity as SE on natural scale

# Define number of simulations, vectors and matrices to store results
nsim <- 100000
lambda <- R0 <- GT <- numeric(nsim)
stable.stage <- matrix(NA, ncol=2, nrow=nsim)
sensitivity <- matrix(NA, ncol=4, nrow=nsim)

# Generate demographic values from beta and normal distributions
library(IPMbook)
sj.sim <- rbeta2(nsim, mean.sj, se.sj.e)
sa.sim <- rbeta2(nsim, mean.sa, se.sa.e)
f1.sim <- rnorm(nsim, mean.f1, se.f1.e)
fa.sim <- rnorm(nsim, mean.fa, se.fa.e)

# Perform Monte Carlo simulations
for (s in 1:nsim){
  if(s %% 1000 == 0) {cat(paste("*** Simrep", s, "***\n"))}   # Counter
  
  # Transition matrix
  A <- matrix(c(sj.sim[s] * f1.sim[s], sj.sim[s] * fa.sim[s], sa.sim[s], sa.sim[s]),
              ncol=2, byrow=TRUE)
  eigenA <- eigen(A)
  
  # Asymptotic population growth rate
  lambda[s] <- max(Re(eigenA$values))
  
  # Stable stage distribution
  u <- which.max(Re(eigenA$values))
  revec <- Re(eigenA$vectors[,u])
  stable.stage[s,] <- revec / sum(revec)
  
  # Lower level sensitivities
  levec <- Re(solve(eigenA$vectors)[u,])
  senmat <- levec %*% t(revec)
  derivmat <- matrix(c(f1.sim[s], fa.sim[s], 0, 0), ncol=2, byrow=TRUE)
  sensitivity[s,1] <- sum(senmat * derivmat)
  derivmat <- matrix(c(sj.sim[s], 0, 0, 0), ncol=2, byrow=TRUE)
  sensitivity[s,2] <- sum(senmat * derivmat)
  derivmat <- matrix(c(0, sj.sim[s], 0, 0), ncol=2, byrow=TRUE)
  sensitivity[s,3] <- sum(senmat * derivmat)
  derivmat <- matrix(c(0, 0, 1, 1), ncol=2, byrow=TRUE)
  sensitivity[s,4] <- sum(senmat * derivmat)
  
  # Net reproductive rate
  i <- 1:100
  R0[s] <- sj.sim[s] * f1.sim[s] + sj.sim[s] * fa.sim[s] * sum(sa.sim[s]^i)
  
  # Generation time
  GT[s] <- log(R0[s]) / log(lambda[s])
}

# ~~~ save the results for use in section 3.4.2 ~~~
save(lambda, stable.stage, sensitivity, GT, file="IPM_03.3.2_output.RData")

# ~~~~ code for Table 3.1 ~~~~
# Summaries of the quantities of interest
summary_table <- rbind(
  "Population growth rate"=c(mean(lambda),sd(lambda), quantile(lambda, c(0.025, 0.975))),
  "Stable stage distribution (1y)"=c(mean(stable.stage[,1]), sd(stable.stage[,1]),
                                     quantile(stable.stage[,1], c(0.025, 0.975))),
  "Stable stage distribution (adults)"=c(mean(stable.stage[,2]), sd(stable.stage[,2]),
                                         quantile(stable.stage[,2], c(0.025, 0.975))),
  "Growth rate sensitivity of sj"=c(mean(sensitivity[,1]), sd(sensitivity[,1]),
                                    quantile(sensitivity[,1], c(0.025, 0.975))),
  "Growth rate sensitivity of f1"=c(mean(sensitivity[,2]), sd(sensitivity[,2]),
                                    quantile(sensitivity[,2], c(0.025, 0.975))),
  "Growth rate sensitivity of fa"=c(mean(sensitivity[,3]), sd(sensitivity[,3]),
                                    quantile(sensitivity[,3], c(0.025, 0.975))),
  "Growth rate sensitivity of sa"=c(mean(sensitivity[,4]), sd(sensitivity[,4]),
                                    quantile(sensitivity[,4], c(0.025, 0.975))),
  "Net reproductive rate"=c(mean(R0), sd(R0), quantile(R0, c(0.025, 0.975))),
  "Generation time"=c(mean(GT), sd(GT), quantile(GT, c(0.025, 0.975))))
colnames(summary_table)[1:2] <- c("Mean", "Monte Carlo SE")
round(summary_table, 3)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code to produce figure 3.11 ~~~~
op <- par(mar=c(4, 4, 3, 0), las=1, cex=1.1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.1, 1), heights=c(1, 1), TRUE)
a <- hist(lambda, nclass=50, col="dodgerblue", main="",
          xlab=expression(paste("Asymptotic population growth rate (", lambda, ")")), prob=TRUE)
mtext("A", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(stable.stage[,1], nclass=50, col="dodgerblue", main="",
          xlab="Proportion of 1-year old individuals", prob=TRUE)
mtext("B", at=a$mids[2], cex=1.5)
par(mar=c(4, 4, 3, 0))
a <- hist(sensitivity[,1], nclass=50, col="dodgerblue", main="",
          xlab=expression('Sensitivity of '*lambda*' to variation of '*italic(s)[italic(j)]),
          prob=TRUE)
mtext("C", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(GT, nclass=30, col="dodgerblue", main="", xlab="Generation time", prob=TRUE)
mtext("D", at=a$mids[2], cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.3 Analysis of a matrix population model with environmental stochasticity
# ----------------------------------------------------------------------------

# Define mean and temporal variability (SD) of the demographic parameters
mean.sj <- 0.3          # Mean juvenile survival (probability scale)
sd.sj.t <- 0.25         # Temporal variability on the logit scale
mean.sa <- 0.55         # Mean adult survival (probability scale)
sd.sa.t <- 0.07         # Temporal variability on the logit scale
mean.f1 <- 1.3          # Mean productivity of 1y old females
sd.f1.t <- 0.3          # Temporal variability on the natural scale
mean.fa <- 1.8          # Mean productivity of adult females
sd.fa.t <- 0.3          # Temporal variability on the natural scale

# Define the number of years with predictions and burn-in length
T <- 100000             # Length of Markov chain
u <- 1000               # Length of burn-in period

# Generate demographic values from normal distributions
sj <- plogis(rnorm(T, qlogis(mean.sj), sd.sj.t))
sa <- plogis(rnorm(T, qlogis(mean.sa), sd.sa.t))
f1 <- rnorm(T, mean.f1, sd.f1.t)
fa <- rnorm(T, mean.fa, sd.fa.t)

# Define population matrix and initial stage-specific population sizes
N <- matrix(NA, nrow=2, ncol=T+1)
N[,1] <- c(1, 1)

# Project population forwards
r <- numeric(T)
for (t in 1:T){
  if(t %% 1000 == 0) {cat(paste("*** Simrep", t, "***\n")) }  # Counter
  A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol=2, byrow=TRUE)
  N[,t+1] <- A %*% N[,t]
  r[t] <- log(sum(N[,t+1])) - log(sum(N[,t]))                 # Annual population growth rate
  N[,t+1] <- N[,t+1] / sum(N[,t+1])                           # Scale N to avoid numerical overflow
}

mean(r[u:T])
# [1] 0.01944997

exp(mean(r[u:T]))
# [1] 1.01964

# Generate demographic values from normal distributions
sj <- plogis(rnorm(T, qlogis(mean.sj), sd.sj.t))
sa <- plogis(rnorm(T, qlogis(mean.sa), sd.sa.t))
f1 <- rnorm(T, mean.f1, sd.f1.t)
fa <- rnorm(T, mean.fa, sd.fa.t)

# Define population matrix and initial stage-specific population sizes
N <- N.star <- matrix(NA, nrow=2, ncol=T+1)
N[,1] <- N.star[,1] <- c(1, 1)

# Project population and calculate stochastic population growth rate
delta <- 0.001                                                # Magnitude of perturbation (= "small change")
r <- r.star <- numeric(T)
for (t in 1:T){
  if(t %% 1000 == 0) {cat(paste("*** Simrep", t, "***\n")) }  # Counter
  # Projection using sj
  A <- matrix(c(sj[t] * f1[t], sj[t] * fa[t], sa[t], sa[t]), ncol=2, byrow=TRUE)
  N[,t+1] <- A %*% N[,t]
  r[t] <- log(sum(N[,t+1])) - log(sum(N[,t]))                 # Annual population growth rate
  N[,t+1] <- N[,t+1] / sum(N[,t+1])                           # Scale N to avoid numerical overflow
  # Projection using sj.star
  A.star <- matrix(c((sj[t] + delta) * f1[t], (sj[t] + delta) * fa[t], sa[t], sa[t]), ncol=2,
                   byrow=TRUE)
  N.star[,t+1] <- A.star %*% N.star[,t]
  r.star[t] <- log(sum(N.star[,t+1])) - log(sum(N.star[,t]))  # Annual population growth rate
  N.star[,t+1] <- N.star[,t+1] / sum(N.star[,t+1])            # Scale N.star to avoid numerical overflow
}

# Compute stochastic sensitivity (juvenile survival)
(exp(mean(r.star[u:T])) - exp(mean(r[u:T]))) / delta
# [1] 1.450347

# Compute stochastic elasticity (juvenile survival)
(exp(mean(r.star[u:T])) - exp(mean(r[u:T]))) / delta * mean.sj / exp(mean(r[u:T]))
# [1] 0.4263838
# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 3 mins

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.4 Analysis of a matrix population model with demographic stochasticity
# --------------------------------------------------------------------------

rbinom(n=1, size=10, prob=0.55)
# [1] 4

rbinom(n=20, size=10, prob=0.55)
# [1] 4 6 7 6 7 4 4 7 9 3 5 5 6 6 1 6 2 6 1 6

size <- c(10, 25, 50, 75, 100, 250, 500, 750, 1000, 5000, 10000)
s <- matrix(NA, nrow=1000, ncol=length(size))
for (i in 1:length(size)){
  s[,i] <- rbinom(n=1000, size=size[i], prob=0.55) / size[i]
}
boxplot(s, ylab="Realized proportion of surviving birds",
        xlab="Population size (number of females)", names=size, outline=FALSE,
        col="dodgerblue", las=1, frame=FALSE)
abline(h=0.55, col="red", lwd=1)

# Define mean of the demographic parameters
mean.sj <- 0.3
mean.sa <- 0.55
mean.f1 <- 1.3
mean.fa <- 1.8

# Define the number of years over which we let the population evolve
T <- 200

# Define population matrix and initial stage-specific population sizes
N <- matrix(NA, nrow=2, ncol=T+1)
N[,1] <- c(10, 10)

# Project population
r <- numeric(T)
for (t in 1:T){
  N[1,t+1] <- rpois(1, mean.sj * mean.f1 * N[1,t] + mean.sj * mean.fa * N[2,t])
  N[2,t+1] <- rbinom(1, sum(N[,t]), mean.sa)
  if (sum(N[,t+1]) == 0) break                                # Stop calculation if pop. extinct
  r[t] <- log(sum(N[,t+1])) - log(sum(N[,t]))
}

mean(r)                                                       # Compute stochastic growth rate
# [1] 0.02552

# Plot graphs with estimated population growth rates and size
op <- par("mfrow", "cex")
layout(matrix(1:2, 1, 2, byrow=TRUE), widths=c(1, 1), heights=1, TRUE)
plot(r, type="l", lwd=1.5, ylab="Annual population growth rate", xlab="Time", axes=FALSE)
axis(1)
axis(2, las=1)
abline(h=mean(r))
mtext("A", at=1, cex=1.5)
plot(N[1,], type="l", lwd=1.5, ylab="Population size", xlab="Time", axes=FALSE,
     ylim=range(N, na.rm=TRUE))
axis(1)
axis(2, las=1)
mtext("B", at=3, cex=1.5)
lines(N[2,], lwd=1.5, col="red")
legend("topleft", lwd=c(1.5, 1.5), col=c("black", "red"), legend=c("1-year-old", ">1-year-old"), bty="n")
par(op)

# Define mean of the demographic parameters
mean.sj <- 0.3
mean.sa <- 0.55
mean.f1 <- 1.3
mean.fa <- 1.8

# Define the number of years with predictions and the Monte Carlo setting
T <- 200
nsim <- 100000

# Define population matrix and initial stage-specific population sizes
N <- array(NA, dim=c(2, T+1, nsim))
N[,1,] <- c(10, 10)
r <- matrix(NA, nrow=T, ncol=nsim)
alive <- matrix(NA, nrow=T, ncol=nsim)
mean.r <- numeric(nsim)

# Project population
for (s in 1:nsim){
  if(s %% 1000 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  for (t in 1:T){
    N[1,t+1,s] <- rpois(1, mean.sj * mean.f1 * N[1,t,s] + mean.sj * mean.fa * N[2,t,s])
    N[2,t+1,s] <- rbinom(1, sum(N[,t,s]), mean.sa)
    if (sum(N[,t+1,s]) == 0) break
    r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s]))
    alive[t,s] <- t
  } #t
  mean.r[s] <- mean(r[min(alive[,s], na.rm=TRUE):max(alive[,s], na.rm=TRUE),s])
} #s

# Mean and SD of the population growth rate
mean(mean.r)
# [1] -0.009097617
sd(mean.r)
# [1] 0.05689681

# Mean and SD of the population growth rate for populations that did not go extinct
not.extinct <- which(!is.na(alive[T,]))
mean(mean.r[not.extinct])
# [1] 0.01966304
sd(mean.r[not.extinct])
# [1] 0.006095664

# Extinction probability (after T years)
sum(is.na(alive[T,])) / nsim
# [1] 0.28917

# Plot graph with population growth rates (Fig. 3.14 and 3.15)
op <- par(mar=c(4, 4, 2, 1), las=1, cex=1.1, "mfrow")
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow=TRUE), widths=c(1, 1), heights=c(1, 1), TRUE)
# par(mar=c(4, 4, 2, 1), las=1, cex=1.1)
plot(r[,1], type="l", lwd=0.5, ylab="Annual population growth rate", xlab="Time",
     ylim=range(r[which(!is.na(alive))]), col="lightgrey", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(r[!is.na(alive[,s]),s], lwd=0.5, col="lightgrey")
}
lines(apply(r, 1, mean, na.rm=TRUE), lwd=1.5)
mtext("A", at=0, cex=1.5)
a <- hist(mean.r, nclass=100, col="dodgerblue", main="", xlab="Population growth rate",
          xlim=c(-0.5, 0.1), axes=FALSE)
axis(1)
axis(2, at=c(0, 10000, 20000, 30000, 40000), labels=c(0, 10, 20, 30, 40))
mtext("B", at=a$mids[26], cex=1.5)
a <- hist(mean.r[not.extinct], nclass=25, col="dodgerblue", main="",
          xlab="Population growth rate", axes=FALSE)
axis(1)
axis(2, at=c(0, 2000, 4000, 6000, 8000, 10000), labels=c(0, 2, 4, 6, 8, 10))
mtext("C", at=a$mids[1], cex=1.5)
par(op)

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 3 mins

library(IPMbook)

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.5 Analysis of a matrix population model with different sources of
#       stochasticity and parameter uncertainty
# ---------------------------------------------------------------------

# Define mean, measurement error and temporal variability of the demographic parameters
mean.sj <- 0.3          # Mean value of juv. survival
sd.sj.e <- 0.005        # Uncertainty of mean juv. survival as SD on natural scale
sd.sj.t <- 0.25         # Temporal variability of juv. survival as SD on logit scale
mean.sa <- 0.55         # Mean value of ad. survival
sd.sa.e <- 0.005        # Uncertainty of mean ad. survival as SD on natural scale
sd.sa.t <- 0.07         # Temporal variability of ad. survival as SD on logit scale
mean.f1 <- 1.3          # Mean value of productivity of 1y females
sd.f1.e <- 0.05         # Uncertainty of mean productivity as SD on natural scale
sd.f1.t <- 0.3          # Temporal variability of productivity as SD on natural scale
mean.fa <- 1.8          # Mean value of productivity of adult females
sd.fa.e <- 0.03         # Uncertainty of mean productivity as SD on natural scale
sd.fa.t <- 0.3          # Temporal variability of productivity as SD on natural scale

# Define the number of years with predictions and the Monte Carlo setting
T <- 200                # Number of years (projection time frame)
nsim <- 100000          # Number of replicate populations simulated

# Define population matrix and initial stage-specific population sizes
N <- array(NA, dim=c(2, T+1, nsim))
N[,1,] <- c(10, 10)
r <- matrix(NA, nrow=T, ncol=nsim)
alive <- matrix(NA, nrow=T, ncol=nsim)
mean.r <- numeric(nsim)

# Project population
for (s in 1:nsim){ # Loop over replicate populations
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) }   # Counter
  # Generate a mean of the demographic rates (subject to measurement error)
  msj <- rbeta2(1, mean.sj, sd.sj.e)
  msa <- rbeta2(1, mean.sa, sd.sa.e)
  mf1 <- rnorm(1, mean.f1, sd.f1.e)
  mfa <- rnorm(1, mean.fa, sd.fa.e)
  
  # Generate annual demographic rates (subject to temporal variability)
  sj <- plogis(rnorm(T, qlogis(msj), sd.sj.t))
  sa <- plogis(rnorm(T, qlogis(msa), sd.sa.t))
  f1 <- pmax(0, rnorm(T, mf1, sd.f1.t))                       # Avoids negative values
  fa <- pmax(0, rnorm(T, mfa, sd.fa.t))                       # Dito
  
  # Project population (include demographic stochasticity)
  for (t in 1:T){                                             # Loop over years
    N[1,t+1,s] <- rpois(1, sj[t] * (f1[t] * N[1,t,s] + fa[t] * N[2,t,s]))
    N[2,t+1,s] <- rbinom(1, sum(N[,t,s]), sa[t])
    if (sum(N[,t+1,s]) == 0) break
    r[t,s] <- log(sum(N[,t+1,s])) - log(sum(N[,t,s]))
    alive[t,s] <- t
  } #t
  mean.r[s] <- mean(r[min(alive[,s], na.rm=TRUE):max(alive[,s], na.rm=TRUE),s])
} #s

mean(mean.r)
# [1] -0.01562592
sd(mean.r)
# [1] 0.06659048

not.extinct <- which(!is.na(alive[T,]))
mean(mean.r[not.extinct])
# [1] 0.02353411
sd(mean.r[not.extinct])
# [1] 0.01242122

# Extinction probability (after T years)
sum(is.na(alive[T,])) / nsim
# [1] 0.38023

# ~~~~ Plot graph with population growth rates (Fig. 3.15)  ~~~~
op <- par(mar=c(4, 4, 2, 1), las=1, cex=1.1, "mfrow")
layout(matrix(c(1, 1, 2, 3), 2, 2, byrow=TRUE), widths=c(1, 1), heights=c(1, 1), TRUE)
plot(r[,1], type="l", lwd=0.5, ylab="Annual population growth rate", xlab="Time",
     ylim=range(r[which(!is.na(alive))]), col="lightgrey", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(r[!is.na(alive[,s]),s], lwd=0.5, col="lightgrey")
}
lines(apply(r, 1, mean, na.rm=TRUE), lwd=1.5)
mtext("A", at=0, cex=1.5)
a <- hist(mean.r, nclass=100, col="dodgerblue", main="",
          xlab="Population growth rate", xlim=c(-0.5, 0.1), axes=FALSE)
axis(1)
axis(2, at=c(0, 5000, 10000, 15000, 20000, 25000, 30000), labels=c(0, 5, 10, 15, 20, 25, 30))
mtext("B", at=a$mids[51], cex=1.5)
a <- hist(mean.r[not.extinct], nclass=25, col="dodgerblue", main="",
          xlab="Population growth rate", axes=FALSE)
axis(1)
axis(2, at=c(0, 2000, 4000, 6000, 8000, 10000), labels=c(0, 2, 4, 6, 8, 10))
mtext("C", at=a$mids[1], cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.6 Matrix population models with density-dependence and
#       demographic stochasticity
# -----------------------------------------------------------

# Define means of the demographic rates and strength of density dependence
mean.sj <- 0.3
mean.sa <- 0.55
f1.int <- 2.3         # Productivity of 1y females when population size is 0
f1.beta <- -0.02      # Strength of density dependence on productivity of 1y
# This is simply the slope of the regression of f1 on N
fa.int <- 2.3         # Productivity of adult females when population size is 0
fa.beta <- -0.01      # Strength of density dependence on productivity of adults

# Define the number of years with predictions
T <- 200
nsim <- 1000

# Define population matrix and initial stage-specific population sizes
N <- array(NA, dim=c(2, T+1, nsim))
N[,1,] <- c(10, 10)
r <- f1 <- fa <- matrix(NA, nrow=T, ncol=nsim)

# Project population
for (s in 1:nsim){
  if(s %% 100 == 0) {cat(paste("*** Simrep", s, "***\n")) } # Counter
  for (t in 1:T){
    f1[t,s] <- f1.int + f1.beta * sum(N[,t,s])
    fa[t,s] <- fa.int + fa.beta * sum(N[,t,s])
    N[1,t+1,s] <- rpois(1, mean.sj * (f1[t,s] * N[1,t,s] + fa[t,s] * N[2,t,s]))
    N[2,t+1,s] <- rbinom(1, sum(N[,t,s]), mean.sa)
  } #t
} #s

# ~~~~ code for Fig. 3.16 ~~~~
op <- par(mar=c(2, 4, 2, 1), las=1, cex=1.1, "mfrow")
layout(matrix(1:3, 3, 1, byrow=TRUE), widths=2, heights=c(1, 1, 1.1), TRUE)
plot(colSums(N[,,1]), type="l", col="lightgrey", ylim=c(0, 120),
     ylab="Population size", xlab="", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(colSums(N[,,s]), col="lightgrey")
}
Nall <- apply(N, c(2,3), sum)
lines(apply(Nall, 1, mean), col="black", lwd=1.5)
mtext("A", at=1, cex=1.5)

f1[f1==f1.int] <- NA
plot(f1[,1], type="l", col="lightgrey", ylim=c(0.3, 2.3),
     ylab=expression(italic(f)[1]), xlab=NA, axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(f1[,s], col="lightgrey")
}
lines(apply(f1, 1, mean, na.rm=TRUE), col="black", lwd=1.5)
mtext("B", at=1, cex=1.5)

par(mar=c(4,4,2,1))
fa[fa==fa.int] <- NA
plot(fa[,1], type="l", col="lightgrey", ylim=c(0.3, 2.3),
     ylab=expression(italic(f)[italic(a)]), xlab="Time", axes=FALSE)
axis(1); axis(2)
for (s in 2:nsim){
  lines(fa[,s], col="lightgrey")
}
lines(apply(fa, 1, mean, na.rm=TRUE), col="black", lwd=1.5)
mtext("C", at=1, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mean(N[1,T+1,] + N[2,T+1,])
# [1] 52.687

mean(f1[T,], na.rm=TRUE)
# [1] 1.237938
mean(fa[T,], na.rm=TRUE)
# [1] 1.768969

mean(N[1,T+1,] + N[2,T+1,] == 0)
# [1] 0.001

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.1 Analysis of a matrix population model without stochasticity
#       and parameter uncertainty
# ------------------------------------------------------------------

# Bundle data
sj <- 0.3
sa <- 0.55
f1 <- 1.3
fa <- 1.8
jags.data <- list(sj=sj, sa=sa, f1=f1, fa=fa, T=5)

# Write JAGS model file
cat(file="model1.txt", "
model {
  # Model for initial state
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population projection
    N[1,t+1] <- sj * (f1 * N[1,t] + fa * N[2,t])
    N[2,t+1] <- sa * (N[1,t] + N[2,t])
    # Calculation of population quantities
    # Annual (realized) population growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])
  }
  lambda <- ann.growth.rate[T] # gr in final interval is our estimate of lambda
}
")

# Parameters monitored
parameters <- c("ann.growth.rate", "lambda", "N")

# MCMC settings
ni <- 1; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out1 <- jags(jags.data, NULL, parameters, "model1.txt", n.adapt=na, n.chains=nc, n.thin=nt,
             n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out1, 3)
#                     mean sd  2.5%   50% 97.5% overlap0 f
# ann.growth.rate[1] 1.015 NA 1.015 1.015 1.015    FALSE 1
# ann.growth.rate[2] 1.021 NA 1.021 1.021 1.021    FALSE 1
# ann.growth.rate[3] 1.021 NA 1.021 1.021 1.021    FALSE 1
# ann.growth.rate[4] 1.021 NA 1.021 1.021 1.021    FALSE 1
# ann.growth.rate[5] 1.021 NA 1.021 1.021 1.021    FALSE 1
# lambda             1.021 NA 1.021 1.021 1.021    FALSE 1
# N[1,1]             1.000 NA 1.000 1.000 1.000    FALSE 1
# N[2,1]             1.000 NA 1.000 1.000 1.000    FALSE 1
# N[1,2]             0.930 NA 0.930 0.930 0.930    FALSE 1
# N[2,2]             1.100 NA 1.100 1.100 1.100    FALSE 1
# [ ...output truncated... ]
# N[1,6]             1.017 NA 1.017 1.017 1.017    FALSE 1
# N[2,6]             1.188 NA 1.188 1.188 1.188    FALSE 1


# Write JAGS model file
cat(file="model2.txt", "
model {
  # Model for initial state
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population projection
    N[1,t+1] <- sj * (f1 * N[1,t] + fa * N[2,t])
    N[2,t+1] <- sa * (N[1,t] + N[2,t])

    # Calculation of population quantities
    # Annual (realized) population growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])

    # Scaled annual stage distributions
    stage.distr[1,t] <- N[1,t+1] / (N[1,t+1] + N[2,t+1])
    stage.distr[2,t] <- N[2,t+1] / (N[1,t+1] + N[2,t+1])
  }
  lambda <- ann.growth.rate[T]
  stable.stage.distr <- stage.distr[,T]

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                                              # Size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj + delta) * (f1 * N.star[1,t] + fa * N.star[2,t])
    N.star[2,t+1] <- sa * (N.star[1,t] + N.star[2,t])
    ann.growth.rate.star[t] <- (N.star[1,t+1] + N.star[2,t+1]) / (N.star[1,t] + N.star[2,t])
  }
  s.sj <- (ann.growth.rate.star[T] - ann.growth.rate[T]) / delta
  e.sj <- s.sj * sj / lambda

  # Calculation of net reproductive rate (R0)
  for (i in 1:100){
    u[i] <- pow(sa, i)
  }
  R0 <- sj * f1 + sj * fa * sum(u[])

  # Calculation of generation time (GT)
  GT <- log(R0) / log(lambda)
}
")

# Parameters monitored
parameters <- c("lambda", "stable.stage.distr", "s.sj", "e.sj", "R0", "GT")

# MCMC settings
ni <- 1; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out2 <- jags(jags.data, NULL, parameters, "model2.txt", n.adapt=na, n.chains=nc, n.thin=nt,
             n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out2, 3)

#                        mean sd  2.5%   50% 97.5% overlap0 f
# lambda                1.021 NA 1.021 1.021 1.021    FALSE 1
# stable.stage.distr[1] 0.461 NA 0.461 0.461 0.461    FALSE 1
# stable.stage.distr[2] 0.539 NA 0.539 0.539 0.539    FALSE 1
# s.sj                  1.454 NA 1.454 1.454 1.454    FALSE 1
# e.sj                  0.427 NA 0.427 0.427 0.427    FALSE 1
# R0                    1.050 NA 1.050 1.050 1.050    FALSE 1
# GT                    2.368 NA 2.368 2.368 2.368    FALSE 1

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.2 Analysis of a matrix population model with parameter uncertainty
# ----------------------------------------------------------------------

# Define mean and SD of the demographic parameters
mean.sj <- 0.3        # Point estimate of juv. survival
se.sj.e <- 0.03       # Uncertainty of juv. survival as SE on natural scale
mean.sa <- 0.55       # Point estimate of ad. survival
se.sa.e <- 0.015      # Uncertainty of ad. survival as SE on natural scale
mean.f1 <- 1.3        # Point estimate of productivity of 1y females
se.f1.e <- 0.3        # Uncertainty of productivity as SE on natural scale
mean.fa <- 1.8        # Point estimate of productivity of adult females
se.fa.e <- 0.1        # Uncertainty of productivity as SE on natural scale

# Bundle data
jags.data <- list(alpha.sj=getBeta2Par(mean.sj, se.sj.e)[1], beta.sj=getBeta2Par(mean.sj,
                                                                                 se.sj.e)[2], alpha.sa=getBeta2Par(mean.sa, se.sa.e)[1], beta.sa=getBeta2Par(mean.sa, se.sa.e)[2],
                  mean.f1=mean.f1, tau.f1=1/se.f1.e^2, mean.fa=mean.fa, tau.fa=1/se.fa.e^2, T=50)

# Write JAGS model file
cat(file="model3.txt", "
model {
  # Random number generators (RNGs)
  sj ~ dbeta(alpha.sj, beta.sj)           # These only *look* like priors
  sa ~ dbeta(alpha.sa, beta.sa)           # ... but they are not
  f1 ~ dnorm(mean.f1, tau.f1)             # ... as there is no estimation in this model
  fa ~ dnorm(mean.fa, tau.fa)

  # Initialize the population size nodes
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] <- sj * (f1 * N[1,t] + fa * N[2,t])
    N[2,t+1] <- sa * (N[1,t] + N[2,t])

    # Annual (realized) growth rate
    ann.growth.rate[t] <- (N[1,t+1] + N[2,t+1]) / (N[1,t] + N[2,t])

    # Scaled stage distributions
    stage.distr[1,t] <- N[1,t+1] / (N[1,t+1] + N[2,t+1])
    stage.distr[2,t] <- N[2,t+1] / (N[1,t+1] + N[2,t+1])
  }
  lambda <- ann.growth.rate[T]
  stable.stage.distr <- stage.distr[,T]

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                          # size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj + delta) * (f1 * N.star[1,t] + fa * N.star[2,t])
    N.star[2,t+1] <- sa * (N.star[1,t] + N.star[2,t])
    ann.growth.rate.star[t] <- (N.star[1,t+1] + N.star[2,t+1]) / (N.star[1,t] + N.star[2,t])
  }
  s.sj <- (ann.growth.rate.star[T] - ann.growth.rate[T]) / delta
  e.sj <- s.sj * sj / lambda

  # Calculation of net reproductive rate (R0)
  for (i in 1:100){
    u[i] <- pow(sa, i)
  }
  R0 <- sj * f1 + sj * fa * sum(u[])

  # Calculation of generation time (GT)
  GT <- log(R0) / log(lambda)
}
")

# Parameters monitored
parameters <- c("lambda", "stable.stage.distr", "s.sj", "e.sj", "R0", "GT")

# MCMC settings
ni <- 100000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out3 <- jags(jags.data, NULL, parameters, "model3.txt", n.adapt=na, n.chains=nc, n.thin=nt,
             n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out3, 4)
#                         mean     sd   2.5%    50%  97.5% overlap0 f
# lambda                1.0223 0.0627 0.9102 1.0188 1.1553    FALSE 1
# stable.stage.distr[1] 0.4603 0.0323 0.3979 0.4600 0.5241    FALSE 1
# stable.stage.distr[2] 0.5397 0.0323 0.4759 0.5400 0.6021    FALSE 1
# s.sj                  1.4658 0.1935 1.1135 1.4567 1.8693    FALSE 1
# e.sj                  0.4273 0.0452 0.3438 0.4258 0.5203    FALSE 1
# R0                    1.0516 0.1490 0.7798 1.0453 1.3619    FALSE 1
# GT                    2.3842 0.1910 2.0426 2.3727 2.7911    FALSE 1

# ~~~~ Produce figure 3.17 ~~~~
op <- par(mar=c(4, 4, 3, 0), las=1, cex=1.1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.1, 1), heights=c(1, 1), TRUE)
a <- hist(out3$sims.list$lambda, nclass=50, col="dodgerblue", main="",
          xlab=expression(paste("Asymptotic population growth rate (", lambda, ")")), prob=TRUE)
mtext("A", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(out3$sims.list$stable.stage[,1], nclass=50, col="dodgerblue", main="",
          xlab="Proportion of 1-year old individuals", prob=TRUE)
mtext("B", at=a$mids[2], cex=1.5)
par(mar=c(4, 4, 3, 0))
a <- hist(out3$sims.list$s.sj, nclass=50, col="dodgerblue", main="",
          xlab=expression('Sensitivity of '*lambda*' to variation of '*italic(s)[italic(j)]),
          prob=TRUE)
mtext("C", at=a$mids[2], cex=1.5)
par(mar=c(4, 2, 3, 2))
a <- hist(out3$sims.list$GT, nclass=40, col="dodgerblue", main="",
          xlab="Generation time", prob=TRUE)
mtext("D", at=a$mids[2], cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~ Produce figure 3.18 ~~~~
# Load output from section 3.3.2
load("IPM_03.3.2_output.RData")
op <- par(mar=c(4, 4, 3, 0), las=1, cex=1.1, "mfrow")
layout(matrix(1:4, 2, 2, byrow=TRUE), widths=c(1.1, 1), heights=c(1, 1), TRUE)
plot(density(lambda), main="",
     xlab=expression(paste("Asymptotic population growth rate (", lambda, ")")),
     col="blue", lwd=2, axes=FALSE)
lines(density(out3$sims.list$lambda), col="red", lwd=2, lty=2)
axis(1); axis(2)
# legend("topright", legend=c("classical MC", "Bayesian MCMC"),
#    lwd=c(2,2), col=c("blue","red"), bty="n")
mtext("A", at=0.8, cex=1.5)

par(mar=c(4, 2, 3, 2))
plot(density(stable.stage[,1]), main="", xlab="Proportion of 1-year old individuals",
     ylab="", col="blue", lwd=2, axes=FALSE)
axis(1); axis(2)
lines(density(out3$sims.list$stable.stage[,1]), col="red", lwd=2, lty=2)
mtext("B", at=0.35, cex=1.5)

par(mar=c(4, 4, 3, 0))
plot(density(sensitivity[,1]), main="",
     xlab=expression('Sensitivity of '*lambda*' to variation of '*italic(s)[italic(j)]),
     col="blue", lwd=2, axes=FALSE)
axis(1); axis(2)
lines(density(out3$sims.list$s.sj), col="red", lwd=2, lty=2)
mtext("C", at=0.83, cex=1.5)

par(mar=c(4, 2, 3, 2))
plot(density(GT), main="", xlab="Generation time", ylab="", col="blue",
     lwd=2, axes=FALSE)
axis(1); axis(2)
lines(density(out3$sims.list$GT), col="red", lwd=2, lty=2)
mtext("D", at=1.9, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.3 Analysis of a matrix population model with environmental stochasticity
# ----------------------------------------------------------------------------

# Define mean and temporal variability (SD) of the demographic parameters
mean.sj <- 0.3            # Mean juv. survival on probability scale
sd.sj.t <- 0.25           # Temporal variability of sj on the logit scale
mean.sa <- 0.55           # Mean ad. survival on probability scale
sd.sa.t <- 0.07           # Temporal variability of sa on the logit scale
mean.f1 <- 1.3            # Mean productivity of 1y females on natural scale
sd.f1.t <- 0.3            # Temporal variability of f1 on the natural scale
mean.fa <- 1.8            # Mean productivity of adult females on natural scale
sd.fa.t <- 0.3            # Temporal variability of fa on the natural scale

# Bundle data
jags.data <- list(mean.sj=mean.sj, mean.sa=mean.sa, mean.f1=mean.f1, mean.fa=mean.fa,
                  sd.sj.t=sd.sj.t, sd.sa.t=sd.sa.t, sd.f1.t=sd.f1.t, sd.fa.t=sd.fa.t, T=21000, u=1000)

# Write JAGS model file
cat(file="model4.txt", "
model {
  # Calculate precision for temporal variability of demographic rates
  tau.logit.sj <- pow(sd.sj.t, -2)
  tau.logit.sa <- pow(sd.sa.t, -2)
  tau.f1 <- pow(sd.f1.t, -2)
  tau.fa <- pow(sd.fa.t, -2)

  # Use of RNG to accomodate temporal variability of demographic rates (process variability)
  for (t in 1:T){
    sj[t] <- ilogit(logit.sj[t])                      # Backt. from logit to natural scale
    logit.sj[t] ~ dnorm(logit(mean.sj), tau.logit.sj)
    sa[t] <- ilogit(logit.sa[t])                      # Backt. from logit to natural scale
    logit.sa[t] ~ dnorm(logit(mean.sa), tau.logit.sa)
    f1[t] ~ dnorm(mean.f1, tau.f1)
    fa[t] ~ dnorm(mean.fa, tau.fa)
  }

  # Model for initial state
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] <- sj[t] * (f1[t] * N[1,t] + fa[t] * N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])

    # Annual growth rate on log scale
    r.annual[t] <- log(N[1,t+1] + N[2,t+1]) - log(N[1,t] + N[2,t])
  }
  r <- mean(r.annual[u:T])
  lambda <- exp(r)

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                                      # Size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj[t] + delta) * (f1[t] * N.star[1,t] + fa[t] * N.star[2,t])
    N.star[2,t+1] <- sa[t] * (N.star[1,t] + N.star[2,t])
    r.annual.star[t] <- log(N.star[1,t+1] + N.star[2,t+1]) - log(N.star[1,t] + N.star[2,t])
  }
  r.star <- mean(r.annual.star[u:T])
  s.sj <- (exp(r.star) - lambda) / delta
  e.sj <- s.sj * mean.sj / lambda
}
")

# Parameters monitored
parameters <- c("r", "lambda", "s.sj", "e.sj")

# MCMC settings
ni <- 1; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out4 <- jags(jags.data, NULL, parameters, "model4.txt", n.adapt=na, n.chains=nc, n.thin=nt, n.iter=ni,
             n.burnin=nb, DIC=FALSE)
print(out4, 4)

#          mean sd   2.5%    50%  97.5% overlap0 f
# r      0.0191 NA 0.0191 0.0191 0.0191    FALSE 1
# lambda 1.0193 NA 1.0193 1.0193 1.0193    FALSE 1
# s.sj   1.4486 NA 1.4486 1.4486 1.4486    FALSE 1
# e.sj   0.4263 NA 0.4263 0.4263 0.4263    FALSE 1


sj.t <- c(0.303, 0.300, 0.288, 0.310, 0.306, 0.291, 0.303, 0.281, 0.314, 0.303)
sa.t <- c(0.546, 0.545, 0.547, 0.545, 0.547, 0.556, 0.540, 0.551, 0.548, 0.550)
f1.t <- c(1.56, 1.15, 1.38, 1.89, 1.06, 0.69, 1.16, 1.10, 1.48, 1.28)
fa.t <- c(2.03, 1.69, 1.95, 1.98, 1.63, 1.43, 1.63, 1.45, 1.79, 1.91)

T <- 21000                                            # Number of predicted years
t <- 1:10                                             # Number of years with actual data
year <- sample(t, T, replace=TRUE)
sj <- sj.t[year]
sa <- sa.t[year]
f1 <- f1.t[year]
fa <- fa.t[year]

# Bundle data
jags.data <- list(sj=sj, sa=sa, f1=f1, fa=fa, T=21000, u=1000)

# Write JAGS model file
cat(file="model5.txt", "
model {
  # Model for initial state
  N[1,1] <- 1
  N[2,1] <- 1

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] <- sj[t] * (f1[t] * N[1,t] + fa[t] * N[2,t])
    N[2,t+1] <- sa[t] * (N[1,t] + N[2,t])

    # Annual (realized) growth rate on log scale
    r.annual[t] <- log(N[1,t+1] + N[2,t+1]) - log(N[1,t] + N[2,t])
  }
  r <- mean(r.annual[u:T])
  lambda <- exp(r)

  # Sensitivity and elasticity of lambda to changes in sj
  delta <- 0.001                                                # size of perturbation
  N.star[1,1] <- 1
  N.star[2,1] <- 1
  for (t in 1:T){
    N.star[1,t+1] <- (sj[t] + delta) * (f1[t] * N.star[1,t] + fa[t] * N.star[2,t])
    N.star[2,t+1] <- sa[t] * (N.star[1,t] + N.star[2,t])
    r.annual.star[t] <- log(N.star[1,t+1] + N.star[2,t+1]) - log(N.star[1,t] + N.star[2,t])
  }
  r.star <- mean(r.annual.star[u:T])
  s.sj <- (exp(r.star) - lambda) / delta
  e.sj <- s.sj * mean(sj) / lambda
}
")

# Parameters monitored
parameters <- c("r", "lambda", "s.sj", "e.sj")

# MCMC settings
ni <- 1; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out5 <- jags(jags.data, NULL, parameters, "model5.txt", n.adapt=na, n.chains=nc, n.thin=nt,
             n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out5, 4)

#          mean sd   2.5%    50%  97.5% overlap0 f
# r      0.0061 NA 0.0061 0.0061 0.0061    FALSE 1
# lambda 1.0061 NA 1.0061 1.0061 1.0061    FALSE 1
# s.sj   1.4112 NA 1.4112 1.4112 1.4112    FALSE 1
# e.sj   0.4207 NA 0.4207 0.4207 0.4207    FALSE 1

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.4 Analysis of a matrix population model with demographic stochasticity
# --------------------------------------------------------------------------

# Define mean of the demographic parameters
mean.sj <- 0.3
mean.sa <- 0.55
mean.f1 <- 1.3
mean.fa <- 1.8

# Define the number of years with predictions and the Monte Carlo setting
T <- 200

# Define initial stage-specific population sizes
N1 <- 10
N2 <- 10

# Bundle data
jags.data <- list(sj=mean.sj, sa=mean.sa, f1=mean.f1, fa=mean.fa, T=T, N1=N1, N2=N2)

# Write JAGS model file
cat(file="model6.txt", "
model {
  # Model for initial state
  N[1,1] <- N1
  N[2,1] <- N2

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] ~ dpois(sj * (f1 * N[1,t] + fa * N[2,t]))
    N[2,t+1] ~ dbin(sa, (N[1,t] + N[2,t]))
    extinct[t] <- equals(N[1,t+1] + N[2,t+1], 0)      # Determines whether
        # population is still thriving (extinct = 0) or went extinct (extinct = 1)
  }
}
")

# Initial values
Ninit <- matrix(10, nrow=2, ncol=T+1)
Ninit[,1] <- NA
inits <- function(){list(N=Ninit)}

# Parameters monitored
parameters <- c("N", "extinct")

# MCMC settings
ni <- 50000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out6 <- jags(jags.data, inits, parameters, "model6.txt", n.chains=nc, n.thin=nt, n.iter=ni,
             n.burnin=nb, DIC=FALSE)
print(out6, 4)

#                  mean       sd 2.5%   50%    97.5% overlap0 f
# N[1,1]        10.0000   0.0000   10  10.0   10.000    FALSE 1
# N[2,1]        10.0000   0.0000   10  10.0   10.000    FALSE 1
# N[1,2]         9.3110   3.0491    4   9.0   16.000    FALSE 1
# N[2,2]        10.9922   2.2313    7  11.0   15.000    FALSE 1
# [ ...output truncated... ]
# extinct[1]     0.0000   0.0000    0   0.0    0.000    FALSE 1
# extinct[2]     0.0000   0.0000    0   0.0    0.000    FALSE 1
# extinct[3]     0.0000   0.0000    0   0.0    0.000    FALSE 1
# [output truncated... ]
# extinct[198]   0.2876   0.4526    0   0.0    1.000     TRUE 1
# extinct[199]   0.2876   0.4527    0   0.0    1.000     TRUE 1
# extinct[200]   0.2878   0.4527    0   0.0    1.000     TRUE 1

# Calculation of the stochastic population growth rate and of the extinction probability
# outside JAGS
dimensions <- dim(out6$sims.list$extinct)
r.annual <- matrix(NA, nrow=dimensions[2], ncol=dimensions[1])
r <- lambda <- numeric(dimensions[1])
for (s in 1:dimensions[1]){
  for (t in 1:dimensions[2]){
    # Calculate annual growth rate on log scale
    if (out6$sims.list$extinct[s,t] == 1) break
    r.annual[t,s] <- log(out6$sims.list$N[s,1,t+1] + out6$sims.list$N[s,2,t+1]) -
      log(out6$sims.list$N[s,1,t] + out6$sims.list$N[s,2,t])
  } #t
  r[s] <- mean(r.annual[which(out6$sims.list$extinct[s,] == 0),s])
  lambda[s] <- exp(r[s])
} #s

mean(r)
# [1] -0.008533999
sd(r)
# [1] 0.05592566

mean(r[out6$sims.list$extinct[,T]==0])
# [1] 0.01973898
sd(r[out6$sims.list$extinct[,T]==0])
# [1] 0.006009628

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.5 Analysis of a matrix population model with multiple sources of
#       stochasticity and parameter uncertainty
# ---------------------------------------------------------------------

# Define mean, measurement error and temporal variability of the demographic parameters
mean.sj <- 0.3          # Mean juv. survival on probability scale
se.sj.e <- 0.005        # Uncertainty of mean sj as SE on natural scale
sd.sj.t <- 0.25         # Temporal variability of sj as SD on logit scale
mean.sa <- 0.55         # Mean ad. survival on probability scale
se.sa.e <- 0.005        # Uncertainty of mean saas SE on natural scale
sd.sa.t <- 0.07         # Temporal variability of sa as SD on logit scale
mean.f1 <- 1.3          # Mean productivity of 1y females on natural scale
se.f1.e <- 0.05         # Uncertainty of f1 as SE on natural scale
sd.f1.t <- 0.3          # Temporal variability of f1 as SD on natural scale
mean.fa <- 1.8          # Mean productivity of adult females on natural scale
se.fa.e <- 0.03         # Uncertainty of fa as SE on natural scale
sd.fa.t <- 0.3          # Temporal variability of fa as SD on natural scale

# Define the number of years with predictions and the Monte Carlo setting
T <- 200

# Bundle data
N1 <- 10
N2 <- 10
jags.data <- list(alpha.sj=getBeta2Par(mean.sj, se.sj.e)[1], beta.sj=getBeta2Par(mean.sj,
                                                                                 se.sj.e)[2], alpha.sa=getBeta2Par(mean.sa, se.sa.e)[1], beta.sa=getBeta2Par(mean.sa,
                                                                                                                                                             se.sa.e)[2], mean.f1=mean.f1, mean.fa=mean.fa, se.f1.e=se.f1.e, se.fa.e=se.fa.e,
                  sd.sj.t=sd.sj.t, sd.sa.t=sd.sa.t, sd.f1.t=sd.f1.t, sd.fa.t=sd.fa.t, T=T, N1=N1, N2=N2)

# Write JAGS model file
cat(file="model7.txt", "
model {
  # Use of RNG to account for measurement error (uncertainty)
  # in the means of the demographic rates on the natural scale
  mean.sj.e ~ dbeta(alpha.sj, beta.sj)
  mean.sa.e ~ dbeta(alpha.sa, beta.sa)
  mean.f1.e ~ dnorm(mean.f1, tau.f1.e)
  mean.fa.e ~ dnorm(mean.fa, tau.fa.e)

  # Calculate precision of the measurement error of productivity
  tau.f1.e <- pow(se.f1.e, -2)
  tau.fa.e <- pow(se.fa.e, -2)

  # Calculate precision for the temporal variability of demographic rates
  tau.logit.sj.t <- pow(sd.sj.t, -2)
  tau.logit.sa.t <- pow(sd.sa.t, -2)
  tau.f1.t <- pow(sd.f1.t, -2)
  tau.fa.t <- pow(sd.fa.t, -2)

  # Use of RNG to accomodate temporal variability of demographic rates
  # (process variability)
  for (t in 1:T){
    sj[t] <- ilogit(logit.sj[t])                    # Backtransformation from logit scale
    logit.sj[t] ~ dnorm(logit(mean.sj.e), tau.logit.sj.t)
    sa[t] <- ilogit(logit.sa[t])                    # Backtransformation from logit scale
    logit.sa[t] ~ dnorm(logit(mean.sa.e), tau.logit.sa.t)
    f1[t] ~ dnorm(mean.f1.e, tau.f1.t)T(0,)         # Truncation to positive values
    fa[t] ~ dnorm(mean.fa.e, tau.fa.t)T(0,)         # Truncation to positive values
  }

  # Model for initial state
  N[1,1] <- N1
  N[2,1] <- N2

  # Loop over time
  for (t in 1:T){
    # Population model
    N[1,t+1] ~ dpois(sj[t] * (f1[t] * N[1,t] + fa[t] * N[2,t]))
    N[2,t+1] ~ dbin(sa[t], (N[1,t] + N[2,t]))
    extinct[t] <- equals(N[1,t+1] + N[2,t+1], 0)    # Determines whether
        # population is still thriving (extinct = 0) or went extinct (extinct = 1)
  }
}
")

# Parameters monitored
parameters <- c("N", "extinct")

# MCMC settings
ni <- 50000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out7 <- jags(jags.data, NULL, parameters, "model7.txt", n.adapt=na, n.chains=nc, n.thin=nt,
             n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out7, 4)

#                    mean          sd 2.5%   50%      97.5% overlap0 f
# N[1,1]          10.0000      0.0000   10  10.0     10.000    FALSE 1
# N[2,1]          10.0000      0.0000   10  10.0     10.000    FALSE 1
# N[1,2]           9.3800      3.7034    3   9.0     17.000    FALSE 1
# N[2,2]          10.9935      2.2543    7  11.0     15.000    FALSE 1
# [ ... output truncated ... ]
# extinct[1]       0.0000      0.0000    0   0.0      0.000    FALSE 1
# extinct[2]       0.0000      0.0000    0   0.0      0.000    FALSE 1
# extinct[3]       0.0000      0.0000    0   0.0      0.000    FALSE 1
# [ ... output truncated ... ]
# extinct[198]     0.3838      0.4863    0   0.0      1.000     TRUE 1
# extinct[199]     0.3842      0.4864    0   0.0      1.000     TRUE 1
# extinct[200]     0.3846      0.4865    0   0.0      1.000     TRUE 1

# Calculation of the stochastic population growth rate outside JAGS
dimensions <- dim(out7$sims.list$extinct)
r.annual <- matrix(NA, nrow=dimensions[2], ncol=dimensions[1])
r <- lambda <- numeric(dimensions[1])
for (s in 1:dimensions[1]){
  for (t in 1:dimensions[2]){
    # Calculate annual growth rate on log scale
    if (out7$sims.list$extinct[s,t] == 1) break
    r.annual[t,s] <- log(out7$sims.list$N[s,1,t+1] + out7$sims.list$N[s,2,t+1]) -
      log(out7$sims.list$N[s,1,t] + out7$sims.list$N[s,2,t])
  } #t
  r[s] <- mean(r.annual[which(out7$sims.list$extinct[s,] == 0),s])
  lambda[s] <- exp(r[s])
} #s

mean(r)
# [1] -0.01570117
sd(r)
# [1] 0.06555304

# Population growth rate of populations that did not go extinct
mean(r[out7$sims.list$extinct[,T]==0])
# [1] 0.02353183
sd(r[out7$sims.list$extinct[,T]==0])
# [1] 0.01234406

# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.4 Analysis of matrix population models with Markov Chain
#     Monte Carlo (MCMC) software
# ==========================================================

# 3.4.6 Matrix population models with density-dependence and demographic stochasticity
# ------------------------------------------------------------------------------------

# Define means of the demographic rates and strength of density dependence
mean.sj <- 0.3
mean.sa <- 0.55
f1.int <- 2.3           # Productivity of 1y females when population size is 0
f1.beta <- -0.02        # Strength of density dependence on 1y productivity
fa.int <- 2.3           # Productivity of adult females when population size is 0
fa.beta <- -0.01        # Strength of density dependence on ad productivity

# Define the number of years with predictions
T <- 200

# Bundle data
jags.data <- list(sj=mean.sj, sa=mean.sa, f1.int=f1.int, f1.beta=f1.beta, fa.int=fa.int,
                  fa.beta=fa.beta, T=T)

# Write JAGS model file
cat(file="model8.txt", "
model {
  # Model for initial state
  N[1,1] <- 10
  N[2,1] <- 10

  # Loop over time
  for (t in 1:T){
    # Calculate actual productivity
    f1[t] <- f1.int + f1.beta * (N[1,t] + N[2,t])
    fa[t] <- fa.int + fa.beta * (N[1,t] + N[2,t])

    # Population model
    N[1,t+1] ~ dpois(sj * (f1[t] * N[1,t] + fa[t] * N[2,t]))
    N[2,t+1] ~ dbin(sa, (N[1,t] + N[2,t]))
    extinct[t] <- equals(N[1,t+1] + N[2,t+1], 0)          # Determines whether
        # population is still thriving (extinct = 0) or went extinct (extinct = 1)
  }
}
")

# Initial values
N <- matrix(NA, nrow=2, ncol=T+1)
N[,2:(T+1)] <- 10
inits <- function(){list(N=N)}

# Parameters monitored
parameters <- c("N", "f1", "fa", "extinct")

# MCMC settings
ni <- 1000; nt <- 1; nb <- 0; nc <- 1; na <- 0

# Call JAGS (ART <1 min) and summarize results
out8 <- jags(jags.data, inits, parameters, "model8.txt", n.adapt=na, n.chains=nc, n.thin=nt,
             n.iter=ni, n.burnin=nb, DIC=FALSE)
print(out8, 4)

#                 mean     sd    2.5%    50%   97.5% overlap0 f
# N[1,1]       10.0000 0.0000 10.0000 10.000 10.0000    FALSE 1
# N[2,1]       10.0000 0.0000 10.0000 10.000 10.0000    FALSE 1
# N[1,2]       11.9210 3.3838  6.0000 12.000 19.0000    FALSE 1
# N[2,2]       11.0280 2.2804  7.0000 11.000 15.0000    FALSE 1
# N[1,3]       13.4220 4.2102  6.0000 13.000 22.0000    FALSE 1
# N[2,3]       12.6680 3.3665  6.0000 12.000 19.0250    FALSE 1
# [ ... output truncated ... ]
# f1[1]         1.9000 0.0000  1.9000  1.900  1.9000    FALSE 1
# f1[2]         1.8410 0.0814  1.6795  1.840  2.0000    FALSE 1
# f1[3]         1.7782 0.1257  1.5195  1.780  2.0200    FALSE 1
# f1[4]         1.7197 0.1594  1.4000  1.720  2.0200    FALSE 1
# [ ... output truncated ... ]
# fa[1]         2.1000 0.0000  2.1000  2.100  2.1000    FALSE 1
# fa[2]         2.0705 0.0407  1.9897  2.070  2.1500    FALSE 1
# fa[3]         2.0391 0.0628  1.9097  2.040  2.1600    FALSE 1
# fa[4]         2.0099 0.0797  1.8500  2.010  2.1600    FALSE 1
# [ ... output truncated ... ]


# Carrying capacity
mean(rowSums(out8$sims.list$N[,,T+1]))
# [1] 53.382

# Productivity at carrying capacity
mean(out8$sims.list$f1[,200])
# [1] 1.22968
mean(out8$sims.list$fa[,200])
# [1] 1.76484

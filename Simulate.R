library(JMbayes)

n <- 400 # number of subjects
K <- 10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################

# parameters for the linear mixed effects model 1
betas <- c("Group0" = 3.0317, "Group1" = 3.1749, "Time1" =  0.1578)
sigma.y <- 0.6920  # measurement error standard deviation


# parameters for the survival model
gammas <- c("(Intercept)" =  -7.8526 , "Group" = -0.015) # coefficients for baseline covariates
alpha <- c(0.1912, 0.35, 0.5, 0.9, 1.3, 1.9, 2.2, 2.5912,2.9, 3.1912) # association parameter - value

phi <- 1.2 # shape for the Weibull baseline hazard
mean.Cens <-30 # mean of the exponential distribution for the censoring mechanism

D <- diag(c(0.9337, 0.1560)^2)

################################################

Bkn <- c(0, 19.5)
kn <- c(2.1, 5.5)

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
age <- rnorm(n, 46.39673, 13.69578)
DF <- data.frame(year = times, drug = factor(rep(group, each = K)), age = rep(age, each = K))


X <- model.matrix(~ 0 + drug + year, data = DF)
Z <- model.matrix(~ year, data = DF)

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group)

################################################

#simulate random effects
library(MASS)

b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)

# simulate event times

tpower <- function(x, t, p)
  # Function for truncated p-th power function
  (x - t) ^ p * (x > t)



bbase <- function(x, xl, xr, ndx, deg){
  # Function for B-spline basis
  dx <- (xr - xl) / ndx
  knots <- c(-5.0, -2.5,  0.0,  2.5,  5.0,  7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B 
}


thi = floor(19.5) +1#floor(t.max) +1
tlo = 0
nseg = 8
bdeg = 2

# simulate event times
eta.t <- as.vector(as.matrix(W) %*% gammas)
invS <- function (t, u, i) {
  h <- function (s) {
    group0 <- 1 - group[i]
    group1 <- group[i]
    ages <- age[i]
    NS <- s
    XX <- cbind(group0, group1, (NS))
    ZZ <- cbind(1, NS)
    
    fval <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
 
    Lam <- bbase(s, tlo, thi, nseg, bdeg)
    
    fvalT <- fval * Lam           
    
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + fvalT %*% alpha)
 
  }
  integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
  Up <- 5000
  tries <- 5
  Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  while(inherits(Root, "try-error") && tries > 0) {
    tries <- tries - 1
    Up <- Up + 5000
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  }
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}
na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from an exponential distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- runif(n, 0, 2 * mean.Cens)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################


# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))



dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
dat.id <- data.frame(IDnr = unique(id), years = Time, status = event, group = W[, 2])
names(dat) <- c("echotime", "group", "age","IDnr", "y0", "years", "status")


#############################
ind2 <- unique(id)[Time > 19.5]
'%!in%' <- function(x,y)!('%in%'(x,y))
dat$years[dat$IDnr %in% ind2] <- 19.5
dat$status[dat$IDnr %in% ind2] <- 0
dat.id <- dat[!duplicated(dat$IDnr), ]
#############################

data <- dat
data.id <- dat.id

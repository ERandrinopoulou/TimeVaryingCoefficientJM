## specify longitudinal and survival submodel

fm1 <- lme(y0 ~ -1 + group + echotime, data = data, na.action = na.exclude,
           random = ~ echotime | IDnr)

lmeObject <- fm1

timeVar <- "echotime"
lag <- 0
survMod <- "spline-PH"

#data.id$Time <- data.id$years
Time <- data.id$years
event <- data.id$status

#data.id$group <- as.factor(data.id$group)
#data$group <- as.factor(data$group)

W <- model.matrix(~ -1 + group, data.id)[, c(2)]
W <- as.matrix(W)

id <- data$IDnr #as.vector(unclass(lmeObject$groups[[1]]))

offset <- as.vector(c(1, 1 + cumsum(tapply(data$IDnr, data$IDnr, length))))

#################
## Design matrices
formYx <- formula(lmeObject)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = data)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = data)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

###
data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)


mfX.id <- model.frame(TermsX, data = data.id)  
mfZ.id <- model.frame(TermsZ, data = data.id)  
Xtime <- model.matrix(formYx, mfX.id)
Ztime <- model.matrix(formYz, mfZ.id)


## 15-point Gauss-Kronrod rule and design matrices
nT <- length(Time)
zeros <- numeric(nT)

y.long <- model.response(mfX, "numeric")
y <- list(y = y.long, offset = offset, logT = log(Time),
          event = event, zeros = zeros, lag = lag)

gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)


P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))


##
mfX <- model.frame(TermsX, data = data.id2)   
mfZ <- model.frame(TermsZ, data = data.id2)    
Xs <- model.matrix(formYx, mfX)
Zs <- model.matrix(formYz, mfZ)

###########################


## Details of MCMC
con <- list(program = "JAGS", n.chains = 1, n.iter = 55000,
            n.burnin = 35000, n.thin = 2, n.adapt = 500, K = 100,
            C = 5000, working.directory = getwd(), 
            openbugs.directory = NULL, clearWD = TRUE, over.relax = TRUE,
            knots = NULL, 
            bugs.seed = 1, quiet = FALSE)


x <- list(X = X, Z = Z, W = if (survMod == "weibull-PH") {
  if (is.null(W)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                           W)
} else {
  if (is.null(W)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(W) == 1) cbind(W, rep(0, nT)) else W
  }
})



# P-slines for baseline hazard

tpower <- function(x, t, p)
  # Function for truncated p-th power function
  (x - t) ^ p * (x > t)

bbase1 <- function(x, xl, xr, ndx, deg){
  # Function for B-spline basis
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B 
}

thi = floor(max(Time)) +1
tlo = 0
nseg = 8
bdeg = 2
#con$knots <- rr <- c(-5.0, -2.5,  0.0,  2.5,  5.0,  7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0)

B0 = bbase1(Time, tlo, thi, nseg, bdeg)
WBH <- B0


DD <- diag(ncol(WBH))
priorTau.Bs.gammas <- crossprod(diff(DD, diff = 2)) + 1e-06 * DD

B0 = bbase1(c(t(st)), tlo, thi, nseg, bdeg)
WBHs <- B0


#

x <- c(x, list(WBH = WBH, WBHs = WBHs))
##################


# P-slines for the coefficients that link the longitudinal and the survival outcome 


B0 = bbase1(Time, tlo, thi, nseg, bdeg)
Lam <- B0


B0 = bbase1(c(t(st)), tlo, thi, nseg, bdeg)
Lam.s <- B0


DDal <- diag(ncol(Lam))
priorTau.alphas <- crossprod(diff(DDal, diff = 2)) + 1e-06 * DDal


##################

ncX <- ncol(X)
ncZ <- ncol(Z)
ncW <- ncol(x$W)
ncWBH <- ncol(x$WBH)
ncWBHs <- ncol(x$WBHs)
ncF <- ncol(Lam)
ncFs <- ncol(Lam.s)

C <- con$C

nb <- ncZ 




mu0 <- rep(0,(ncZ))


betas <- rep(0, ncX)
var.betas <- rep(con$K, ncX)


alphas <- rep(0, (ncF))
var.alphas <- rep(con$K/10, (ncF))

Dalphas <- rep(0, (ncF))
var.alphas <- rep(con$K/10, (ncF))

gammas <- rep(0,(ncW))
var.gammas <- rep(con$K, (ncW))

Bs.gammas <- rep(0, (ncWBH))
var.Bs.gammas <- rep(con$K/10, (ncWBH))


b <- cbind(data.matrix(ranef(lmeObject)))


nY <- nrow(b)
sigma2 <- lmeObject$sigma^2


############################################

Data <- list(N = nY, K = K, offset = offset, X = X, Xtime = Xtime, 
             y = y$y, 
             Xs = Xs, 
             Z = Z, Ztime = Ztime,  
             Zs = Zs,
             
             event = event, zeros = zeros, 
             ncZ = ncol(Z), 
             ncX = ncol(X), 
             W = x$W,
             ncW = ncol(x$W),  
             
           
             WBH = x$WBH,  WBHs = x$WBHs,
             ncWBH = ncWBH, ncWBHs = ncWBHs , C = C, P = P,
             
             Lam = Lam, Lam.s = Lam.s,
             ncF = ncF, ncFs = ncFs,
             
             wk = wk, nb = nb, 
             mu0 = mu0, 
             priorMean.betas = betas, 
             priorTau.betas = diag(1/var.betas),
             
             priorA.tau = (1/sigma2)^2/10,
             priorB.tau = (1/sigma2)/10, 
             
             priorMean.gammas = gammas,
             priorTau.gammas = diag(1/var.gammas),
             
             priorMean.alphas = alphas,
             priorTau.alphas = priorTau.alphas, 
             
             priorMean.Bs.gammas = Bs.gammas,
             priorTau.Bs.gammas = priorTau.Bs.gammas,
             
             priorR.D = diag(1,(ncZ)), priorK.D = (ncZ),
             priorA.tausBs = 1,
             priorB.tausBs = 0.005
)

parms <- c("betas", "tau", "inv.D","b", "gammas", "alphas", "Bs.gammas")

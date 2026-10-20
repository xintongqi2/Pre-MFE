###############################################################################
#  Statistics for Financial Engineering -- Pre-MFE Seminar
#  R-LAB 1:  Probability distributions, LLN/CLT, and the multivariate normal
#            -- checked by simulation and on real stock-index returns
#
#  Companion to Lecture 1 (normal, Gamma, chi-square, F, t; LLN & CLT;
#  sample mean & variance; multivariate normal, Cov(AX) = A Cov(X) A^T,
#  Wishart, multivariate t).
#  Textbook: Ruppert, "Statistics and Data Analysis for Financial
#  Engineering", Chapters 4-5.
#
#  HOW TO USE THIS FILE (in RStudio):
#    * Put the cursor on a line and press  Ctrl+Enter  (Mac: Cmd+Return)
#      to run that line.  Highlight several lines to run them together.
#    * Lines starting with  #  are comments -- R ignores them.
#    * Plots appear in the "Plots" pane (bottom-right).  Use the arrows
#      there to page back through earlier plots.
###############################################################################


## ===========================================================================
## PART 0.  Setup  (run once per R session)
## ===========================================================================

# Run the next line ONCE on your computer (it downloads packages from CRAN).
# After that you can leave it commented out.
# install.packages(c("fGarch", "mvtnorm", "moments", "Ecdat"))

library(fGarch)    # standardized t distribution: dstd(), stdFit()
library(mvtnorm)   # multivariate normal and t: rmvnorm(), rmvt()
library(moments)   # skewness(), kurtosis()

set.seed(2026)     # makes the random numbers reproducible:
                   # everyone in class gets the same simulations


## ===========================================================================
## PART 1.  R in five minutes + the d / p / q / r family
## ===========================================================================

x <- c(2.1, -0.4, 1.7, 0.3)    # "<-" assigns; c() combines into a vector
mean(x); sd(x); length(x)
x^2                            # operations act element-by-element
?dnorm                         # opens the help page (works for any function)

# Every distribution in R comes with four functions:
#   d...  density (pdf)          p...  cdf  P(X <= q)
#   q...  quantile (inverse cdf) r...  random draws
dnorm(0)                 # pdf of N(0,1) at 0  = 1/sqrt(2*pi)
pnorm(1.96)              # P(Z <= 1.96)
qnorm(0.975)             # the 97.5% quantile
rnorm(5)                 # five N(0,1) draws
pnorm(110, mean = 100, sd = 15)   # X ~ N(100, 15^2): note R uses sd, not variance

# Names for the Lecture-1 distributions:
#   Normal  norm(mean, sd)         Gamma  gamma(shape = alpha, rate = beta)
#   Chi-sq  chisq(df)              F      f(df1, df2)
#   t       t(df)
# CAUTION: the lecture's Gamma(alpha, beta) has pdf
#   beta^alpha / Gamma(alpha) * x^(alpha-1) * exp(-beta*x),
# so beta is a RATE.  Always write  rate = beta  explicitly in R.

# Facts from lecture, checked numerically
gamma(5); factorial(4)                         # Gamma(n) = (n-1)!
gamma(0.5)^2; pi                               # Gamma(1/2) = sqrt(pi)
integrate(function(z) z^2 * dnorm(z), -Inf, Inf)   # E[Z^2] = 1
integrate(function(z) z   * dnorm(z), -Inf, Inf)   # E[Z]   = 0

# Picture: t_n densities approach N(0,1) as n grows
curve(dnorm(x), -5, 5, lwd = 3, ylab = "density",
      main = "Student t_n versus N(0,1)")
cols <- c("firebrick", "darkorange", "forestgreen")
dfs  <- c(1, 3, 10)
for (j in 1:3) curve(dt(x, df = dfs[j]), add = TRUE, lwd = 2, col = cols[j])
legend("topright", c("N(0,1)", "t_1 (Cauchy)", "t_3", "t_10"),
       col = c("black", cols), lwd = c(3, 2, 2, 2))

# Picture: chi-square densities
curve(dchisq(x, df = 1), 0.05, 20, ylim = c(0, 0.5), lwd = 2, ylab = "density",
      main = "Chi-square densities")
for (k in c(2, 4, 8)) curve(dchisq(x, df = k), add = TRUE, lwd = 2, lty = k / 2 + 1)
legend("topright", paste("df =", c(1, 2, 4, 8)), lwd = 2, lty = c(1, 2, 3, 5))

# Critical values used constantly in finance / statistics
qnorm(c(0.95, 0.975, 0.99))       # z-values
qt(0.975, df = c(5, 10, 30, 100)) # t-values -> 1.96 as df grows
qchisq(c(0.025, 0.975), df = 9)   # used for confidence intervals for sigma^2
qf(0.95, df1 = 4, df2 = 100)


## ===========================================================================
## PART 2.  Simulation: checking the theorems of Lecture 1
## ===========================================================================
M <- 100000   # number of Monte Carlo replications

# --- 2a.  Z^2 ~ chi^2_1 = Gamma(1/2, 1/2) ----------------------------------
Z  <- rnorm(M)
hist(Z^2, breaks = 200, freq = FALSE, xlim = c(0, 6), ylim = c(0, 1.2),
     col = "grey85", border = "white", main = "Z^2 with Z ~ N(0,1)")
curve(dgamma(x, shape = 1/2, rate = 1/2), add = TRUE, lwd = 3, col = "navy")
curve(dchisq(x, df = 1), add = TRUE, lwd = 2, lty = 2, col = "orange")
legend("topright", c("Gamma(1/2, 1/2)", "chi^2_1"),
       col = c("navy", "orange"), lwd = c(3, 2), lty = c(1, 2))
dchisq(2.3, df = 1); dgamma(2.3, shape = 0.5, rate = 0.5)   # identical

# --- 2b.  Sum of independent Gammas with the same beta (FACT 4) -------------
S <- rgamma(M, shape = 2, rate = 3) + rgamma(M, shape = 1.5, rate = 3)
hist(S, breaks = 100, freq = FALSE, col = "grey85", border = "white",
     main = "Gamma(2,3) + Gamma(1.5,3)")
curve(dgamma(x, shape = 3.5, rate = 3), add = TRUE, lwd = 3, col = "navy")
c(sample.mean = mean(S), theory = 3.5 / 3)      # E = alpha/beta
c(sample.var  = var(S),  theory = 3.5 / 9)      # Var = alpha/beta^2

# --- 2c.  F_{k,m} = (X/k)/(Y/m) ----------------------------------------------
k <- 5; m <- 20
Fsim <- (rchisq(M, k) / k) / (rchisq(M, m) / m)
hist(Fsim[Fsim < 6], breaks = 100, freq = FALSE, col = "grey85",
     border = "white", main = "(chi2_5/5) / (chi2_20/20)")
curve(df(x, df1 = k, df2 = m), add = TRUE, lwd = 3, col = "navy")

# --- 2d.  t_n = X / sqrt(Y/n),  and  T^2 ~ F_{1,n} --------------------------
n  <- 4
Tsim <- rnorm(M) / sqrt(rchisq(M, n) / n)
hist(Tsim[abs(Tsim) < 8], breaks = 150, freq = FALSE, col = "grey85",
     border = "white", main = "N(0,1) / sqrt(chi2_4 / 4)")
curve(dt(x, df = n), add = TRUE, lwd = 3, col = "navy")
curve(dnorm(x), add = TRUE, lwd = 2, lty = 2, col = "firebrick")
legend("topright", c("t_4", "N(0,1)"), col = c("navy", "firebrick"),
       lwd = c(3, 2), lty = c(1, 2))

# QQ-plot: if T^2 ~ F_{1,4} the points lie on the 45-degree line
p <- ppoints(2000)
qqplot(qf(p, 1, n), quantile(Tsim^2, p), log = "xy",
       xlab = "F(1,4) quantiles", ylab = "quantiles of simulated T^2",
       main = "T^2 versus F(1,4)")
abline(0, 1, col = "red", lwd = 2)

# Homework T1 sanity check: how stable are sample moments of t_n?
# (Run several times.  What happens for df = 1 and df = 2?)
for (d in c(1, 2, 3, 10)) {
  draws <- rt(1e6, df = d)
  cat("df =", d, "  sample mean =", round(mean(draws), 3),
      "  sample variance =", round(var(draws), 3), "\n")
}

# --- 2e.  Sample mean and sample variance of a NORMAL sample ----------------
# Theorem:  Ybar ~ N(mu, sigma^2/n),  (n-1)s^2/sigma^2 ~ chi^2_{n-1},
#           and Ybar, s^2 are INDEPENDENT.
mu <- 0.05; sigma <- 0.2; n <- 10; R <- 20000
samples <- matrix(rnorm(R * n, mean = mu, sd = sigma), nrow = R)  # R samples of size n
ybar <- rowMeans(samples)
s2   <- apply(samples, 1, var)                                    # var() uses n-1

c(sd.of.ybar = sd(ybar), theory = sigma / sqrt(n))
hist((n - 1) * s2 / sigma^2, breaks = 80, freq = FALSE, col = "grey85",
     border = "white", main = "(n-1) s^2 / sigma^2  with n = 10")
curve(dchisq(x, df = n - 1), add = TRUE, lwd = 3, col = "navy")

# Independence: Normal versus Exponential samples
expo  <- matrix(rexp(R * n, rate = 1), nrow = R)
par(mfrow = c(1, 2))
plot(ybar, s2, pch = ".", main = "Normal samples",
     xlab = "sample mean", ylab = "sample variance")
plot(rowMeans(expo), apply(expo, 1, var), pch = ".",
     main = "Exponential samples", xlab = "sample mean", ylab = "sample variance")
par(mfrow = c(1, 1))
cor(ybar, s2)                               # approx 0
cor(rowMeans(expo), apply(expo, 1, var))    # clearly positive


## ===========================================================================
## PART 3.  Law of Large Numbers and Central Limit Theorem
## ===========================================================================

# --- 3a.  LLN: running averages -------------------------------------------
N <- 10000
run.exp    <- cumsum(rexp(N, rate = 1)) / (1:N)    # E[X] = 1 exists
run.cauchy <- cumsum(rt(N, df = 1))      / (1:N)   # t_1 = Cauchy: no mean!
plot(run.exp, type = "l", log = "x", ylim = c(-3, 3), lwd = 2, col = "navy",
     xlab = "n (log scale)", ylab = "running average",
     main = "LLN: Exponential(1) converges, Cauchy does not")
lines(run.cauchy, col = "firebrick", lwd = 2)
abline(h = 1, lty = 2)
legend("topright", c("Exponential(1)", "Cauchy = t_1"),
       col = c("navy", "firebrick"), lwd = 2)

# --- 3b.  CLT: standardized means of a skewed distribution ----------------
par(mfrow = c(1, 3))
for (n in c(1, 5, 30)) {
  zbar <- sqrt(n) * (rowMeans(matrix(rexp(R * n), nrow = R)) - 1) / 1
  hist(zbar, breaks = 80, freq = FALSE, xlim = c(-4, 4), ylim = c(0, 1), col = "grey85",
       border = "white", main = paste("Exponential, n =", n), xlab = "")
  curve(dnorm(x), add = TRUE, lwd = 2, col = "navy")
}
par(mfrow = c(1, 1))


## ===========================================================================
## PART 4.  Real data: four European stock indices (Ruppert, Ch. 4 lab)
## ===========================================================================
data(EuStockMarkets)        # built into R: daily closes 1991-1998
class(EuStockMarkets)       # "mts" = multivariate time series
head(EuStockMarkets)
plot(EuStockMarkets, main = "DAX, SMI, CAC, FTSE: daily closing levels")

# Log returns  r_t = log(P_t) - log(P_{t-1})
logR <- diff(log(EuStockMarkets))
plot(logR, main = "Daily log returns")
dax  <- as.numeric(logR[, "DAX"])
length(dax)

# --- 4a.  Summary statistics -----------------------------------------------
stats <- rbind(mean     = colMeans(logR),
               sd       = apply(logR, 2, sd),
               skewness = apply(logR, 2, skewness),
               kurtosis = apply(logR, 2, kurtosis))   # normal => 3
round(stats, 4)

# Annualizing (sum of ~252 daily returns; FACT 1 / variance of a sum)
c(annual.mean = 252 * mean(dax), annual.vol = sqrt(252) * sd(dax))

# --- 4b.  Is the DAX return normal? ---------------------------------------
hist(dax, breaks = 80, freq = FALSE, col = "grey85", border = "white",
     main = "DAX daily log returns", xlab = "log return")
curve(dnorm(x, mean = mean(dax), sd = sd(dax)), add = TRUE, lwd = 3, col = "navy")

qqnorm(dax, main = "Normal QQ-plot: DAX"); qqline(dax, col = "red", lwd = 2)
shapiro.test(dax)

# Tail counts: how many "4-sigma" days did we see vs. the normal model?
z <- (dax - mean(dax)) / sd(dax)
c(observed       = sum(abs(z) > 4),
  normal.predicts = length(z) * 2 * pnorm(-4))

# --- 4c.  Fit a (standardized) t distribution by maximum likelihood --------
fit <- stdFit(dax)     # fGarch: estimates mean, sd and degrees of freedom nu
fit$par

xx <- seq(-0.1, 0.1, by = 0.001)
plot(density(dax), lwd = 3, main = "DAX: KDE vs fitted t vs normal",
     xlim = c(-0.06, 0.06), xlab = "log return")
lines(xx, dstd(xx, mean = fit$par["mean"], sd = fit$par["sd"],
               nu = fit$par["nu"]), lwd = 2, lty = 2, col = "navy")
lines(xx, dnorm(xx, mean(dax), sd(dax)), lwd = 2, lty = 3, col = "firebrick")
legend("topleft", c("KDE (data)", "fitted t", "normal"),
       lwd = c(3, 2, 2), lty = c(1, 2, 3), col = c("black", "navy", "firebrick"))

# t QQ-plot with the fitted degrees of freedom
nu <- fit$par["nu"]
qqplot(qt(ppoints(length(dax)), df = nu), dax,
       xlab = paste0("t quantiles (df = ", round(nu, 1), ")"),
       ylab = "DAX returns", main = "t QQ-plot: DAX")
qqline(dax, distribution = function(p) qt(p, df = nu), col = "red", lwd = 2)

# --- 4d.  CLT in the market: aggregate daily returns -----------------------
# Non-overlapping sums of h daily log returns = h-day log returns
agg <- function(r, h) colSums(matrix(r[1:(h * floor(length(r) / h))], nrow = h))
sapply(c(day = 1, week = 5, month = 20), function(h) kurtosis(agg(dax, h)))
# Kurtosis falls toward 3 as h grows: "aggregational Gaussianity".


## ===========================================================================
## PART 5.  Multivariate: covariance matrices, MVN, Wishart
## ===========================================================================

# --- 5a.  Sample mean vector, covariance and correlation -------------------
mu.hat    <- colMeans(logR)
Sigma.hat <- cov(logR)            # 4 x 4 sample covariance (divides by n-1)
round(Sigma.hat * 1e4, 3)         # in units of 10^-4
round(cor(logR), 2)
pairs(as.data.frame(logR), pch = ".", main = "Pairwise scatter: daily log returns")

# --- 5b.  Cov(AX) = A Cov(X) A^T :  a portfolio --------------------------
w    <- c(DAX = 0.4, SMI = 0.2, CAC = 0.2, FTSE = 0.2)   # weights, A = w^T (1 x 4)
port <- logR %*% w                                     # portfolio returns
c(direct  = var(as.numeric(port)),
  formula = as.numeric(t(w) %*% Sigma.hat %*% w))      # identical
sqrt(252 * t(w) %*% Sigma.hat %*% w)                   # annual portfolio vol

# --- 5c.  Simulating N(mu, Sigma) as  mu + A Z  with  A A^T = Sigma -----
A <- t(chol(Sigma.hat))          # Cholesky factor (lower triangular)
A %*% t(A) - Sigma.hat           # ~ zero matrix
nsim <- nrow(logR)
Zmat <- matrix(rnorm(4 * nsim), nrow = 4)          # 4 x nsim iid N(0,1)
Ysim <- t(mu.hat + A %*% Zmat)                      # nsim x 4 MVN draws
round(cov(Ysim) * 1e4, 3)                           # ~ Sigma.hat * 1e4
# Same thing with a package: rmvnorm(nsim, mean = mu.hat, sigma = Sigma.hat)

# Same mean & covariance, but do simulated returns LOOK like real ones?
rbind(real.kurtosis = apply(logR, 2, kurtosis),
      mvn.kurtosis  = apply(Ysim, 2, kurtosis))

# Multivariate t: heavier tails with the same correlation structure
# (Cov of multivariate t_nu with scale matrix S is S * nu/(nu-2), so rescale.)
nu  <- 5                                   # t_5: theoretical kurtosis = 9
Tsim <- rmvt(nsim, sigma = Sigma.hat * (nu - 2) / nu, df = nu)  # Cov = Sigma.hat
round(apply(Tsim, 2, kurtosis), 1)         # noisy: 4th moments are hard to estimate

# --- 5d.  Wishart: W = sum_{i=1}^n Y_i Y_i^T,   E[W] = n Sigma ------------
Sig2 <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
Wdraws <- rWishart(5000, df = 10, Sigma = Sig2)   # 2 x 2 x 5000 array
apply(Wdraws, c(1, 2), mean)                      # ~ 10 * Sig2
10 * Sig2
hist(Wdraws[1, 1, ] / Sig2[1, 1], breaks = 80, freq = FALSE, col = "grey85",
     border = "white", main = "W[1,1] / Sigma[1,1]  ~  chi^2_10")
curve(dchisq(x, df = 10), add = TRUE, lwd = 3, col = "navy")

## ===========================================================================
## PART 6.  For homework (Problems P1, P2 = Ruppert Ch. 4 and Ch. 5 R labs)
## ===========================================================================
# Notes on running the textbook code:
#  * The book uses  windows()  to open a new plot window.  That works only
#    on Windows.  On a Mac/Linux, or in RStudio, replace it with  dev.new()
#    or simply delete the line.
#  * The Ch. 5 lab needs  library(Ecdat)  and  library(fGarch)  and  MASS
#    (MASS comes with R).
#  * To save a plot:  in RStudio use  Plots -> Export,  or wrap the plot in
#      pdf("myplot.pdf", width = 6, height = 5);  <plot code>;  dev.off()

## OPTIONAL (needs internet): download current prices from Yahoo Finance
# install.packages("quantmod")
# library(quantmod)
# getSymbols("SPY", from = "2015-01-01")       # S&P 500 ETF
# spy <- as.numeric(na.omit(diff(log(Ad(SPY)))))
# c(kurtosis = kurtosis(spy)); qqnorm(spy); qqline(spy)

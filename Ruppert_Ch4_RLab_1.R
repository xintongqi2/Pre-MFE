###############################################################################
#  Ruppert, "Statistics and Data Analysis for Financial Engineering"
#  Section 4.11  R Lab: European Stock Indices   (Homework 1, Problem P1)
#
#  The textbook's code, run as printed, with two changes so it works on
#  every computer and in RStudio:
#    * windows()  (Windows-only)  is replaced by  dev.new()
#    * in the t-plot loop, only the DAX is drawn by default (set ALL <- TRUE
#      to get all four indices, as in the book)
###############################################################################

library(fGarch)     # for dstd(): the standardized t density

## ---------------------------------------------------------------------------
## Problem 1: the four price series
## ---------------------------------------------------------------------------
data(EuStockMarkets)
mode(EuStockMarkets)       # "numeric"
class(EuStockMarkets)      # "mts" "ts" "matrix": a multivariate time series
plot(EuStockMarkets)

# Save the plot to a PDF in the working directory (optional)
pdf("EuStocks.pdf", width = 6, height = 5)
plot(EuStockMarkets)
graphics.off()

## ---------------------------------------------------------------------------
## Problem 2: log returns
## ---------------------------------------------------------------------------
logR <- diff(log(EuStockMarkets))
plot(logR)

# The same data as a data frame: plot() now draws a scatterplot matrix
plot(as.data.frame(logR))

## ---------------------------------------------------------------------------
## Problem 3: normal plots and Shapiro-Wilk tests for all four indices
## ---------------------------------------------------------------------------
index.names <- dimnames(logR)[[2]]
par(mfrow = c(2, 2))
for (i in 1:4) {
  qqnorm(logR[, i], datax = TRUE, main = index.names[i])
  qqline(logR[, i], datax = TRUE)
  print(shapiro.test(logR[, i]))
}
par(mfrow = c(1, 1))

## ---------------------------------------------------------------------------
## Problems 4-5: t-plots with 1, 4, 6, 10, 20, 30 degrees of freedom
## ---------------------------------------------------------------------------
n      <- dim(logR)[1]
q.grid <- (1:n) / (n + 1)           # probabilities 1/(n+1), ..., n/(n+1)
df     <- c(1, 4, 6, 10, 20, 30)
ALL    <- FALSE                     # TRUE = all four indices, as in the book
for (i in if (ALL) 1:4 else 1) {
  dev.new()                         # the book has windows() here
  par(mfrow = c(3, 2))
  for (j in 1:6) {
    qqplot(logR[, i], qt(q.grid, df = df[j]),
           main = paste(index.names[i], ", df=", df[j]))
    abline(lm(qt(c(.25, .75), df = df[j]) ~ quantile(logR[, i], c(.25, .75))))
  }
}
par(mfrow = c(1, 1))

## extra (not in the book): how straight is each DAX t-plot?
## correlation of the QQ-plot points; the largest value = straightest plot
sapply(c(1, 3, 4, 6, 10, 30), function(d) cor(sort(logR[, 1]), qt(q.grid, d)))
## 0.556 0.990 0.995 0.991 0.985 0.978   -> df = 4 is best

## ---------------------------------------------------------------------------
## Problem 6: kernel density estimate vs t (df = 5) vs normal, DAX
## ---------------------------------------------------------------------------
x <- seq(-0.1, 0.1, by = 0.001)
plot(density(logR[, 1]), lwd = 2, ylim = c(0, 60))
lines(x, dstd(x, mean = median(logR[, 1]), sd = mad(logR[, 1]), nu = 5),
      lty = 5, lwd = 2)
lines(x, dnorm(x, mean = mean(logR[, 1]), sd = sd(logR[, 1])),
      lty = 3, lwd = 4)
legend("topleft", c("KDE", "t: df=5", "normal"), lwd = c(2, 2, 4),
       lty = c(1, 5, 3))

# Zoom in on the LEFT tail (adjust = 2 smooths the KDE more)
plot(density(logR[, 1], adjust = 2), lwd = 2,
     xlim = c(-0.08, -0.025), ylim = c(0, 3), main = "DAX: left tail")
lines(x, dstd(x, mean = median(logR[, 1]), sd = mad(logR[, 1]), nu = 5),
      lty = 5, lwd = 2)
lines(x, dnorm(x, mean = mean(logR[, 1]), sd = sd(logR[, 1])), lty = 3, lwd = 4)
legend("topleft", c("KDE", "t: df=5", "normal"), lwd = c(2, 2, 4), lty = c(1, 5, 3))

# Zoom in on the RIGHT tail
plot(density(logR[, 1], adjust = 2), lwd = 2,
     xlim = c(0.025, 0.06), ylim = c(0, 3), main = "DAX: right tail")
lines(x, dstd(x, mean = median(logR[, 1]), sd = mad(logR[, 1]), nu = 5),
      lty = 5, lwd = 2)
lines(x, dnorm(x, mean = mean(logR[, 1]), sd = sd(logR[, 1])), lty = 3, lwd = 4)
legend("topright", c("KDE", "t: df=5", "normal"), lwd = c(2, 2, 4), lty = c(1, 5, 3))

## ---------------------------------------------------------------------------
## Problem 7: what does density() use by default?
## ---------------------------------------------------------------------------
?density                 # read "bw" and "kernel" in the help page
density(logR[, 1])       # printout shows the bandwidth it chose
bw.nrd0(logR[, 1])       # the default bandwidth rule, computed directly

## check_setup.R -- run this BEFORE each R-lab (it checks everything Labs 1-2 need).
## In RStudio: File -> Open File -> check_setup.R, then click "Source" (top-right of the editor).

pkgs <- c("fGarch", "mvtnorm", "moments", "Ecdat")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing) > 0) {
  cat("Installing:", missing, "\n")
  install.packages(missing)
}

suppressPackageStartupMessages({
  library(fGarch); library(mvtnorm); library(moments); library(Ecdat); library(MASS)
})
data(CPSch3)                                   # earnings data used in R-Lab 2

data(EuStockMarkets)
logR <- diff(log(EuStockMarkets))
plot(logR, main = "If you can see four return series, you are ready")

cat("\nR version:", R.version.string, "\n")
cat("DAX kurtosis:", round(kurtosis(as.numeric(logR[, "DAX"])), 2), "(should print 9.28)\n")
cat("CPSch3 rows:", nrow(CPSch3), "(should print 11130)\n")
cat("\nAll set -- see you in the lab!\n")

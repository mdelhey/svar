# Test on simulated data.
# Orcale: vars::VAR(p) model

context("Simulated data")

### Generate n = 500 x p = 3 multivariate time-series
N <- 300
w <- c(0, 0, 0)              # omega: intercept term
C <- rbind(c(1.0, 0.5, 0.25) # sigma: noise covariance
         , c(0.5, 1.5, 0.25)
         , c(0.5, 1.0, 1.5))
A <- rbind(c(0.2, 0.5, 0.1)  # matrix of AR coefficients
         , c(0.0, 0.5, 0.1)
         , c(0.0, 0.5, 0))
ts <- as.ts(mAr::mAr.sim(w, A, C, N))

# Fit oracle & svar
fit.var  <- vars::VAR(ts, p = 5, type = "both")
fit.svar <- svar(ts, p = 5, method = "lm", type = "both")

# Compare predictions
pred.svar <- predict(fit.svar, n.ahead = 10)
pred.var  <- predict(fit.var,  n.ahead = 10)

preds.var  <- as.numeric(unlist(lapply(pred.var$fcst, function(i) i[, "fcst"])))
preds.svar <- as.numeric(pred.svar)

test_that("predictions are the same for: p = 5 and n.timeseries = 3", {
    expect_equivalent(preds.var, preds.svar)
})

### Generate n = 500 x p = 5 multivariate time-series
ts.larger <- cbind(ts, ts[, 1] * abs(ts[, 2]) + ts[, 3], ts[, 1]^3 - ts[, 2] + log(abs(ts[, 3])))

fit.var.larger  <- vars::VAR(ts.larger, p = 5, type = "both")
fit.svar.larger <- svar(ts.larger, p = 5, method = "lm", type = "both")

pred.svar.larger <- predict(fit.svar, n.ahead = 2)
pred.var.larger  <- predict(fit.var, n.ahead = 2)

preds.var.larger  <- as.numeric(unlist(lapply(pred.var.larger$fcst, function(i) i[, "fcst"])))
preds.svar.larger <- as.numeric(pred.svar.larger)

test_that("predictions are the same for: p = 5 and n.timeseries = 5", {
    expect_equivalent(preds.var.larger, preds.svar.larger)
})

# More tests.
#fit.svar <- svar(ts, p = 1, method = "l1.ls", type = "none", cv.lambda = TRUE, epsilon.min = 5)
#fit.svar <- svar(ts, p = 1, method = "lm", type = "const", cv.lambda = TRUE, epsilon.min = 5)

#' Test documentation
svar <- function(ts, p, lambda = NULL, cv.lambda = TRUE, verbose = TRUE,
                 method = c("lm", "ols", "olsl2", "l1.ls", "l1.ll"), 
                 type = c("none", "const", "trend", "both"),
                 epsilon.min = 0.01
                 ) {
    ### Sparse vector estimation with a few different methods.
    ### Created usign the VAR package.
    # Init house-keeping variables.
    alpha <- 1 # We only care about lasso in this experiment.
    epsilon <- Inf # init
    n.obs <- dim(ts)[1]
    n.predictors <- dim(ts)[2]
    type <- match.arg(type)
    method <- match.arg(method)
    n.constants <- switch(type, none = 0, const = 1, trend = 1, both = 2)
    # Generate lagged matrix
    ts.lag <- embed(ts, dimension = p + 1)[, -(1:n.predictors)]
    # Constant term: column of ones. Append to end.
    if (any(c("const", "both") %in% type))
        ts.lag <- cbind(ts.lag, rep(1, n.obs - p))
    # Trend term: integer sequence. Append to end.
    if (any(c("trend", "both") %in% type))
        ts.lag <- cbind(ts.lag, seq(p + 1, length = n.obs - p))    
    # Collect coefficients and residuals. Regress each variable on all lagged others (OLS estimates)
    fits.coef <- matrix(
        rep(NA, n.predictors * (n.predictors * p + n.constants)) # O(n.predictors^2 * p)
      , nrow = n.predictors
      , ncol = (n.predictors * p + n.constants)
    )
    fits.resi <- matrix(
        rep(NA, n.predictors * (n.obs - p))
      , nrow = n.obs - p
      , ncol = n.predictors
    )
    # Fix names.
    rownames(fits.coef) <- colnames(ts)
    cnames <- rep(colnames(ts), p)
    for (j in 1:p) cnames <- paste0(cnames[], j)
    if (n.constants == 1 && type == "trend") cnames <- append(cnames, "trend")
    if (n.constants == 1 && type == "const") cnames <- append(cnames, "const")
    if (n.constants == 2 && type == "both")  cnames <- append(cnames, c("const", "trend"))
    colnames(fits.coef) <- cnames
    # These methods need their own loop because they must converge.
    if (method == "l1.ll") {
        if (verbose) cat("Initializing sigma.hat and beta.hat using l1.ls \n")
        # Need to init Beta.hat and Sigma.hat for l1-ll.
        init <- svar(ts = ts, p = p, method = "l1.ls", type = type)
        sigma.hat <- init$covariance
        beta.hat <- init$coefficients
    } 
    # Loop through each predictor, estimate coefficients invidually.
    while (epsilon > epsilon.min) {
        for (j in 1:n.predictors) {
            if (verbose && (j %% 10 == 0)) cat("Estimating predictor:", j, "\n")
            # Must remove first p rows from non-lagged matrix to use OLS.
            # No bias term. Already taken care of in design matrix.
            if (method == "lm") {            
                fit.formula <- as.formula(ts[(1+p):n.obs, j] ~ -1 + .)
                fit <- lm(fit.formula, data = as.data.frame(ts.lag))
                # Save coefficients.
                if (any(is.na(coef(fit)))) {
                    warning(sprintf("[%i] coefficients are NA! Likely due to p > n.",
                                    length(which(is.na(coef(fit))))))
                    fit$coefficients[is.na(fit$coefficients)] <- 0
                }                
                fits.coef[j, ] <- coef(fit)
                fits.resi[, j] <- fit$residuals
                epsilon <- -Inf # no convergence
            }
            if (method == "ols") {
                X <- as.matrix(ts.lag)
                Y <- as.matrix(ts[(1+p):n.obs, j]) 
                beta.hat <- solve(t(X) %*% X) %*% t(X) %*% Y
                y.hat <- X %*% beta.hat
                resid <- Y - y.hat
                # Save coefficients.
                fits.coef[j, ] <- t(beta.hat)
                fits.resi[, j] <- resid
                epsilon <- -Inf # no convergence
            }
            if (method == "ols.l2") {
                X <- as.matrix(ts.lag)
                Y <- as.matrix(ts[(1+p):n.obs, j])
                I <- diag(rep(1, dim(X)[2]))
                beta.hat <- solve(t(X) %*% X %*% + (lambda * I)) %*% t(X) %*% Y
                y.hat <- X %*% beta.hat
                resid <- Y - y.hat
                # Save coefficients.
                fits.coef[j, ] <- t(beta.hat)
                fits.resi[, j] <- resid
                epsilon <- -Inf # no convergence
            }
            if (method == "l1.ls") {
                X <- as.matrix(ts.lag)
                # Subset to only 1 of p time-series
                Y <- as.matrix(ts[(1+p):n.obs, j])
                if (all(Y == Y[1]))
                    stop(sprintf("Time-series at column %i is constant.", j))
                if (cv.lambda) 
                    lambda.star <- cv.glmnet(x = X, y = Y, intercept = FALSE, family = "gaussian",
                                             nfolds = 8, alpha = alpha, parallel = TRUE,
                                             nlambda = 100, lambda.min.ratio = 1e-10,
                                             type.gaussian = "naive")$lambda.min
                else 
                    lambda.star <- lambda
                print(lambda.star)
                fit <- glmnet(x = X, y = Y, family = "gaussian", lambda = lambda.star,
                              alpha = alpha, intercept = FALSE)
                fits.coef[j, ] <- as.numeric(coef(fit))[-1] # Ignore intercept term.
                fits.resi[, j] <- Y - predict(fit, newx = X)
                epsilon <- -Inf # no convergence
            }            
            if (method == "l1.ll") {                
                X <- as.matrix(ts.lag)
                # Used for estimating B_j
                Y_j <- as.matrix(ts[(1+p):n.obs, j])
                B_j <- beta.hat[j, ]   
                # Used for estimating resid
                Y_i <- as.matrix(ts[(1+p):n.obs, -j])
                B_i <- t(beta.hat[-j, ])          
                # Useful sigma aliases
                sigma.hat.ji <- sigma.hat[upper.tri(sigma.hat)[-j, -j]]
                sigma.hat.jj <- sqrt(sigma.hat[j, j])
                # Calculate resid
                Y <- as.matrix(ts[(1+p):n.obs, ])
                B <- beta.hat
                summand <- sigma.hat %*% t((Y - X %*% B))
                r_j <- (1/2) * sigma.hat.jj * sum(summand[, -j])
                #r_j <- (1/2) * sigma.hat.jj * sum(sigma.hat.ji * (Y_i - X %*% B_i))
                if (all(Y == Y[1]))
                    stop(sprintf("Time-series at column %i is constant.", j))
                if (cv.lambda) 
                    lambda.star <- cv.glmnet(x = X, y = Y + r_j, intercept = FALSE, family = "gaussian",
                                             nfolds = 8, alpha = alpha, parallel = TRUE,
                                             nlambda = 200, type.gaussian = "covariance")$lambda.min
                else 
                    lambda.star <- lambda
                # Update beta hat
                fit <- glmnet(x = X, y = Y_j + r_j, family = "gaussian",
                              alpha = alpha, intercept = FALSE,
                              lambda = lambda.star)
                B_j.new <- as.numeric(coef(fit))[-1] # ignore intercept term
                # Check for convergance
                epsilon <- mean(abs(B_j.new - B_j))
                beta.hat[j, ] <- B_j.new
                # Save coefficients
                fits.coef[j, ] <- beta.hat[j, ]
                fits.resi[, j] <- Y_j - predict(fit, newx = X)                
            }
            if (verbose && epsilon != -Inf) cat("Epsilon = ", epsilon, "\n")
        }
    }
    # Estimate covariance & correlation
    est.var <- list(
        coefficients  = fits.coef
      , covariance    = suppressWarnings(cov(fits.resi))
      , correlation   = suppressWarnings(cor(fits.resi))
      , p             = p
      , n.constants   = n.constants
      , type          = type
      , method        = method
      , ts            = ts
      , lambda        = ifelse(exists("lambda.star"), lambda.star, NA)
    )
    # Assign class. This lets us use generics with svar methods.
    class(est.var) <- "svar"
    return(est.var)
}

print.svar <- function(est.var, ...) {
    cat("VAR(p):\t p =", est.var$p,      "\n")
    cat("Type:\t",       est.var$type,   "\n")
    cat("Method:\t",     est.var$method, "\n")
    cat("Lambda:\t",     est.var$lambda, "\n")
    
    cat("\nCoefficients:\n")
    printCoefmat(est.var$coefficients, P.value = FALSE, has.Pvalue = FALSE)
}

predict.svar <- function(est.var, n.ahead, ...) {
    ### predict generic for svar class.
    ### currently only does n.ahead prediction, no s.e.'s
    if (!is.numeric(n.ahead) || n.ahead <= 0) stop("n.ahead must be positive integer")
    # Init house-keeping variables
    ts <- est.var$ts
    n.obs <- dim(ts)[1]
    n.predictors <- dim(ts)[2]
    p <- est.var$p
    n.constants <- est.var$n.constants
    type <- est.var$type
    A <- est.var$coefficients
    # We need to edit the TS for recurssion.
    ts.pred <- ts
    for (n in 1:n.ahead) {
        summands <- matrix(rep(NA, n.predictors * p), nrow = p, ncol = n.predictors)
        for (j in 1:p) {
            # A contains p (n.predictors x n.predictors) submatricies
            A_p <- A[, ((n.predictors * (j-1)) + 1):(n.predictors * (j))]
            Y_p <- ts.pred[(nrow(ts.pred) + 1) - j, ]
            summands[j, ] <- t(A_p %*% Y_p)
        }
        # Constant and Trend doesn't depend on the TS, so also not on the order P.
        constant <- switch(
            type
          , none  = 0,
          , const = A[, "const"]
          , trend = A[, "trend"] * seq(nrow(ts) + 1, length = n.ahead)[n]
          , both  = A[, "const"] + (A[, "trend"] * seq(nrow(ts) + 1, length = n.ahead)[n])
        )
        # Generate prediction.
        new.row <- constant + apply(summands, 2, sum)
        # Append to ts for recursion.
        ts.pred <- rbind(ts.pred, new.row)
    }
    # Fix names
    n.ahead.prediction <- tail(ts.pred, n.ahead)
    rownames(n.ahead.prediction) <- paste0("n.ahead.", 1:n.ahead)
    # Assign a new class. This lets us define new plot/print methods.
    class(n.ahead.prediction) <- "svar.pred"
    return(n.ahead.prediction)
}

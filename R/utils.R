rmse <- function(y.pred, y.known, na.rm = FALSE, ..., type = c("vec", "row", "col")) {
    type <- match.arg(type)
    func <- switch(type
      , vec = "mean"
      , row = "rowMeans"
      , col = "colMeans"
    )
    cmd <- sprintf("sqrt( %s( (y.pred - y.known)^2, na.rm = na.rm, ...) )", func)
    eval(parse(text = cmd))
}

lambda.seq <- function(lambda.min, lambda.max, n.lambda) {
    seq(from = lambda.min, to = lambda.max, length = n.lambda)
}

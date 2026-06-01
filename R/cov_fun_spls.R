cov_fun_spls <- function(y, x) {

  y <- as.numeric(y)
  x <- as.matrix(x)

  yc <- y - mean(y)
  xc <- scale(x, center = TRUE, scale = FALSE)

  out <- crossprod(xc, yc) / max(1, (nrow(x) - 1))
  out <- matrix(out, ncol = 1)

  return(out)
}

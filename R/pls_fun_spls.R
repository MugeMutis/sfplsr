#' @importFrom MASS ginv
#' @importFrom stats coef lm
pls_fun_spls <- function(y, x, h) {

  x0 <- as.matrix(x)
  y0 <- as.numeric(y)

  n <- nrow(x0)
  p <- ncol(x0)

  if (h < 1) {
    stop("'h' must be a positive integer.")
  }

  h <- min(h, p, n - 1)

  x_center <- colMeans(x0)
  y_center <- mean(y0)

  x_res <- scale(x0, center = x_center, scale = FALSE)
  y_res <- matrix(y0 - y_center, ncol = 1)

  xw <- matrix(0, p, h)
  yw <- matrix(0, 1, h)
  xs <- matrix(0, n, h)
  ys <- matrix(0, n, h)
  xl <- matrix(0, p, h)
  yl <- numeric(h)

  h_eff <- 0

  for (i in 1:h) {

    cv <- cov_fun_spls(y_res, x_res)
    sv <- svd(cv)

    xw_i <- sv$u[, 1]
    yw_i <- sv$v[, 1]

    buff <- svd_flip(xw_i, yw_i)
    xw_i <- matrix(buff$x, ncol = 1)
    yw_i <- matrix(buff$y, ncol = 1)

    xs_i <- x_res %*% xw_i
    ys_i <- y_res %*% yw_i

    denom <- drop(crossprod(xs_i))
    if (!is.finite(denom) || denom < 1e-12) {
      break
    }

    xl_i <- crossprod(x_res, xs_i) / denom
    yl_i <- as.numeric(crossprod(xs_i, y_res) / denom)

    x_res <- x_res - xs_i %*% t(xl_i)
    y_res <- y_res - xs_i %*% matrix(yl_i, nrow = 1)

    xw[, i] <- xw_i[, 1]
    yw[, i] <- yw_i[, 1]
    xs[, i] <- xs_i[, 1]
    ys[, i] <- ys_i[, 1]
    xl[, i] <- xl_i[, 1]
    yl[i] <- yl_i

    h_eff <- i
  }

  if (h_eff == 0) {
    stop("PLS extraction failed because no valid latent component was obtained.")
  }

  xw <- xw[, 1:h_eff, drop = FALSE]
  yw <- yw[, 1:h_eff, drop = FALSE]
  xs <- xs[, 1:h_eff, drop = FALSE]
  ys <- ys[, 1:h_eff, drop = FALSE]
  xl <- xl[, 1:h_eff, drop = FALSE]
  yl <- yl[1:h_eff]

  xr <- xw %*% ginv(t(xl) %*% xw)
  fin.model <- lm(y0 ~ xs)
  beta_T <- as.matrix(coef(fin.model))
  fin.cf <- xr %*% beta_T[-1, , drop = FALSE]

  fits <- as.vector(cbind(1, xs) %*% beta_T)
  resds <- y0 - fits

  return(list(
    x = x0,
    y = y0,
    T = xs,
    R = xr,
    P = xl,
    W = xw,
    YW = yw,
    Yload = yl,
    d.coef = fin.cf,
    pqr.coef = beta_T,
    fitted.values = fits,
    residuals = resds,
    xcenter = x_center,
    ycenter = y_center,
    h_eff = h_eff
  ))
}

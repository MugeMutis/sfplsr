fit_spls_fixed_rho <- function(y, A, h, W) {

  n <- length(y)
  x_center <- colMeans(A)
  A_centered <- scale(A, center = x_center, scale = FALSE)

  pls_fit <- pls_fun_spls(y = y, x = A, h = h)

  alpha_hat <- as.numeric(pls_fit$pqr.coef[1, 1])
  coef_hat  <- as.matrix(pls_fit$d.coef)

  z_hat <- as.vector(alpha_hat + A_centered %*% coef_hat)
  res_z <- y - z_hat

  k_eff <- pls_fit$h_eff + 1
  mse <- mean(res_z^2)
  bic <- n * log(max(mean(res_z^2), 1e-12)) + log(n) * k_eff
  aic <- n * log(max(mean(res_z^2), 1e-12)) + 2 * k_eff

  return(list(
    alpha = alpha_hat,
    coef = coef_hat,
    fitted.values = z_hat,
    residuals = res_z,
    mse = mse,
    bic = bic,
    aic = aic,
    pls_fit = pls_fit,
    xcenter = x_center,
    h_eff = pls_fit$h_eff
  ))
}

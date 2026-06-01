#' @importFrom MASS ginv
fsac_spls <- function(y,
                      x,
                      h,
                      nbasis = NULL,
                      gpx = NULL,
                      wei_mat,
                      rho_grid = seq(-0.95, 0.95, by = 0.01),
                      criterion = c("train_mse", "bic", "aic"),
                      verbose = TRUE) {

  criterion <- match.arg(criterion)

  y <- as.numeric(y)
  x <- as.matrix(x)
  W <- normalize_weights(wei_mat)

  n <- length(y)
  px <- ncol(x)

  if (nrow(W) != n || ncol(W) != n) {
    stop("The spatial weight matrix has incompatible dimensions.")
  }

  if (is.null(gpx)) {
    gpx <- seq(1 / px, 1 - 1 / px, length.out = px)
  }

  if (is.null(nbasis)) {
    nbasis <- round(min(n / 4, 40))
  }

  BS.sol <- getAmat(data = x, nbasis = nbasis, gp = gpx)
  A <- BS.sol$Amat
  evalbase <- BS.sol$evalbase
  sinp_mat <- BS.sol$sinp_mat

  score_vec <- numeric(length(rho_grid))
  mse_vec   <- numeric(length(rho_grid))
  bic_vec   <- numeric(length(rho_grid))
  aic_vec   <- numeric(length(rho_grid))
  fit_list  <- vector("list", length(rho_grid))

  I_n <- diag(n)

  for (i in seq_along(rho_grid)) {

    rho_i <- rho_grid[i]
    z_i <- as.vector((I_n - rho_i * W) %*% y)

    fit_i <- fit_spls_fixed_rho(y = z_i, A = A, h = h, W = W)

    zhat_i <- fit_i$fitted.values
    yhat_i <- as.vector(ginv(I_n - rho_i * W) %*% zhat_i)
    res_i  <- y - yhat_i

    mse_y <- mean(res_i^2)
    bic_y <- n * log(max(mse_y, 1e-12)) + log(n) * (fit_i$h_eff + 2)
    aic_y <- n * log(max(mse_y, 1e-12)) + 2 * (fit_i$h_eff + 2)

    fit_i$rho <- rho_i
    fit_i$z <- z_i
    fit_i$zhat <- zhat_i
    fit_i$yhat <- yhat_i
    fit_i$res_y <- res_i
    fit_i$mse_y <- mse_y
    fit_i$bic_y <- bic_y
    fit_i$aic_y <- aic_y

    mse_vec[i] <- mse_y
    bic_vec[i] <- bic_y
    aic_vec[i] <- aic_y

    if (criterion == "train_mse") score_vec[i] <- mse_y
    if (criterion == "bic")       score_vec[i] <- bic_y
    if (criterion == "aic")       score_vec[i] <- aic_y

    fit_list[[i]] <- fit_i
  }

  best_id <- which.min(score_vec)
  best_fit <- fit_list[[best_id]]

  beta_basis <- solve(sinp_mat, best_fit$coef)
  beta_hat_t <- as.vector(evalbase %*% beta_basis)

  sig_hat <- sqrt(sum(best_fit$res_y^2) / max(1, n - best_fit$h_eff - 2))

  search_table <- data.frame(
    rho = rho_grid,
    train_mse = mse_vec,
    bic = bic_vec,
    aic = aic_vec,
    criterion_value = score_vec
  )

  if (verbose) {
    message("Selected rho = ", round(best_fit$rho, 4),
            "; effective PLS components = ", best_fit$h_eff,
            "; criterion = ", criterion)
  }

  fd_details <- list(
    nbasis = nbasis,
    gpx = gpx,
    m.tr = best_fit$xcenter,
    fin.cf = best_fit$coef,
    alpha = best_fit$alpha,
    bs_basis = BS.sol$bs_basis,
    inp_mat = BS.sol$inp_mat,
    sinp_mat = BS.sol$sinp_mat,
    evalbase = BS.sol$evalbase,
    search_table = search_table
  )

  return(list(
    alpha = best_fit$alpha,
    b = best_fit$coef,
    bhat = beta_hat_t,
    beta_basis = beta_basis,
    rho = best_fit$rho,
    sig = sig_hat,
    fitted.values = best_fit$yhat,
    residuals = best_fit$res_y,
    transformed.response = best_fit$z,
    transformed.fitted = best_fit$zhat,
    criterion = criterion,
    criterion.value = score_vec[best_id],
    h = h,
    h_eff = best_fit$h_eff,
    fdd = fd_details,
    A = A,
    x = x,
    y = y,
    w = W,
    pls_fit = best_fit$pls_fit
  ))
}


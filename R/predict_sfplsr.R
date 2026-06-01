#' Predict Method for Spatial FPLS Regression
#'
#' @param object An object obtained from the \code{splsr} function.
#' @param xnew A numeric matrix of new functional predictors (new observations).
#' @param wnew A square spatial weights matrix for the new observations or test validation block.
#'
#' @return A numeric vector of predicted values for the scalar response.
#' @importFrom MASS ginv
#' @export
predict_sfplsr <- function(object, xnew, wnew) {

  xnew <- as.matrix(xnew)
  wnew <- normalize_weights(wnew)

  n <- nrow(xnew)
  rho <- object$rho
  alpha <- object$alpha
  nbasis <- object$fdd$nbasis
  gpx <- object$fdd$gpx
  m.tr <- object$fdd$m.tr
  fin.cf <- object$fdd$fin.cf

  if (nrow(wnew) != n || ncol(wnew) != n) {
    stop("'wnew' has incompatible dimensions.")
  }

  BS.sol.test <- getAmat(data = xnew, nbasis = nbasis, gp = gpx)
  Anew <- BS.sol.test$Amat
  Anew_centered <- scale(Anew, center = m.tr, scale = FALSE)

  zhat_new <- as.vector(alpha + Anew_centered %*% fin.cf)
  yhat_new <- as.vector(ginv(diag(n) - rho * wnew) %*% zhat_new)

  return(yhat_new)
}


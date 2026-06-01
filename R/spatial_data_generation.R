#' Generate Spatial Functional Data
#'
#' @param n An integer specifying the sample size (number of spatial locations).
#' @param j An integer specifying the number of grid points for the functional predictor.
#' @param rho A numeric value representing the spatial autoregressive parameter.
#' @param sig.e A numeric value specifying the standard deviation of the error terms.
#'
#' @return A list containing the generated spatial functional dataset components.
#' @importFrom fda.usc fdata inprod.fdata
#' @importFrom stats dist rnorm
#' @export
#'
spatial_data_generation = function(n, j, rho, sig.e) {

  ##wei_mat_______________________
  side = ceiling(sqrt(n))
  coords = expand.grid(x = 1:side, y = 1:side)[1:n,]
  dist_mat = as.matrix(dist(coords))

  W = matrix(0, n, n)
  for(i in 1:n) {
    knn = order(dist_mat[i,])[2:6]
    W[i, knn] = 1
    W[i,] = W[i,] / sum(W[i,])
  }
  ##_____________________________

  ##functional_predictor_________
  s = seq(0, 1, length.out = j)

  ksi_mat = matrix(0, n, 5)

  for(k in 1:5) {
    base_signal = rnorm(n, 0, sd = 10 * k^(-1.5))
    spatial_noise = as.numeric(solve(diag(n) - rho * W) %*% rnorm(n, 0, sig.e))
    ksi_mat[,k] = base_signal + spatial_noise
  }


  phi = list()
  for(k in 1:5) phi[[k]] = sqrt(2) * sin(k * pi * s)

  fX_data = ksi_mat %*% do.call(rbind, phi)
  fX = fda.usc::fdata(fX_data, argvals = s)
  ##_____________________________

  vBeta_val = 0.2 * phi[[1]] + 0.4 * phi[[2]] + 0.6 * phi[[3]] + 0.8 * phi[[4]] + 1 * phi[[5]]
  vBeta = fda.usc::fdata(vBeta_val, argvals = s)

  argx = as.numeric(fda.usc::inprod.fdata(fX, vBeta))
  ##_____________________________

  volatility = seq(0.5, 1.5, length.out = n)
  err = rnorm(n, mean = 0, sd = sig.e * volatility)
  ##_____________________________

  ST = solve(diag(n) - rho*W)
  y = as.numeric(ST %*% (argx + err))

  return(list("y" = y, "x" = fX$data, "w" = W, "tcoefs" = vBeta))
}

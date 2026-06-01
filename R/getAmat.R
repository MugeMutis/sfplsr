#' @importFrom expm sqrtm
#' @importFrom fda create.bspline.basis eval.basis fdPar inprod int2Lfd smooth.basisPar
getAmat <- function(data, nbasis, gp){

  # data      : functional data
  # nbasis    : number of basis functions
  # gp        : grid points

  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  bs_basis <- create.bspline.basis(rangeval = c(gp[1], gp[p]), nbasis = nbasis)
  inp_mat <- inprod(bs_basis, bs_basis)
  sinp_mat <- sqrtm(inp_mat)
  evalbase = eval.basis(gp, bs_basis)
  fdobj <- fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  Amat <- t(pcaobj$coefs) %*% sinp_mat

  return(list(Amat = Amat, bs_basis = bs_basis, inp_mat = inp_mat,
              sinp_mat = sinp_mat, evalbase = evalbase))
}

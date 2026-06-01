normalize_weights <- function(W, tol = 1e-8) {
  row_sums <- rowSums(W)

  if (all(abs(row_sums - 1) < tol)) {
    message("The weight matrix is already row-normalized.")
    return(W)
  } else {
    message("The weight matrix is not row-normalized. Normalizing now.")
    # Normalize the weight matrix row-wise
    W_normalized <- sweep(W, 1, row_sums, FUN = "/")
    return(W_normalized)
  }
}

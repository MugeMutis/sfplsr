#' Spatial Functional Partial Least Squares Regression
#'
#' @param y A numeric vector of the scalar response variable.
#' @param x A numeric matrix of functional explanatory variable.
#' @param h An integer or vector specifying the number of latent components.
#' @param nbasis An integer or vector specifying the number of basis functions.
#' @param gpx A numeric vector defining the evaluation grid points for \code{x}.
#' @param wei_mat A square spatial weights matrix.
#' @param rho_grid A numeric vector defining the search grid for spatial parameter \eqn{\rho}.
#' @param criterion Optimization metric: \code{"train_mse"}, \code{"bic"}, or \code{"aic"}.
#' @param scv Logical; if \code{TRUE}, performs Spatial Cross-Validation for \code{nbasis} and \code{h}.
#' @param block_ids A vector indicating the spatial block assignment for each observation (required if \code{scv = TRUE}).
#' @param verbose Logical; if \code{TRUE}, prints progress and tuning results.
#'
#' @importFrom MASS ginv
#' @export
sfplsr <- function(y, x, h, nbasis = NULL, gpx = NULL, wei_mat,
                  rho_grid = seq(-0.95, 0.95, by = 0.01),
                  criterion = c("train_mse", "bic", "aic"),
                  scv = FALSE, block_ids = NULL, verbose = TRUE) {

  criterion <- match.arg(criterion)
  n_total <- length(as.numeric(y))

  y <- as.numeric(y)
  x <- as.matrix(x)
  W <- normalize_weights(wei_mat)

  n <- length(y)
  px <- ncol(x)

  # ---- CASE 1: Spatial Cross-Validation (scv == TRUE) ----
  if (scv == TRUE) {
    if (is.null(block_ids)) {
      stop("When scv = TRUE, the 'block_ids' argument indicating spatial blocks cannot be left empty!")
    }

    # Default nbasis assignment if it is left NULL by the user
    if (is.null(nbasis)) {
      nbasis <- round(min(n_total / 4, 40))
    }

    unique_blocks <- unique(block_ids)
    num_blocks <- length(unique_blocks)

    # Create the candidate grid matrix
    candidate_grid <- expand.grid(K = nbasis, h = h)
    scv_scores <- rep(NA, nrow(candidate_grid))

    if (verbose) {
      cat("Starting Spatial Cross-Validation (SCV) grid search...\n")
    }

    # Loop over all (K, h) combinations
    for (g in 1:nrow(candidate_grid)) {
      current_K <- candidate_grid$K[g]
      current_h <- candidate_grid$h[g]

      # The number of components (h) cannot exceed the number of basis functions (K)
      if (current_h > current_K) next

      block_squared_errors <- numeric(num_blocks)
      is_valid_comb <- TRUE

      # Loop over all spatial blocks
      for (v_idx in seq_along(unique_blocks)) {
        v <- unique_blocks[v_idx]
        test_idx <- which(block_ids == v)
        n_v <- length(test_idx)

        # Split the data to get the training set without block v (-v)
        y_train <- y[-test_idx]
        x_train <- x[-test_idx, , drop = FALSE]
        w_train <- W[-test_idx, -test_idx, drop = FALSE]

        # Fit the model without block v (called with scv = FALSE internally)
        fit_val <- tryCatch({
          fsac_spls(y = y_train, x = x_train, h = current_h, nbasis = current_K,
                    gpx = gpx, wei_mat = w_train, rho_grid = rho_grid,
                    criterion = criterion, verbose = FALSE)
        }, error = function(e) { NULL })

        if (is.null(fit_val)) {
          is_valid_comb <- FALSE
          break
        }

        # --- OUT-OF-SAMPLE PREDICTION STEP ---
        # 1. Project the test data onto the basis functions estimated by the model
        BS.test <- getAmat(data = x[test_idx, , drop = FALSE], nbasis = current_K, gp = fit_val$fdd$gpx)
        A_test <- BS.test$Amat

        # 2. Compute predictions in the transformed space (Z_hat_v)
        # Note: Centering is performed using the training set mean (m.tr)
        x_test_centered <- t(t(A_test) - fit_val$fdd$m.tr)
        zhat_pred_v <- as.vector(x_test_centered %*% fit_val$b + fit_val$alpha)

        # 3. Re-introduce the spatial dependency structure (rho and W) to return to the original Y scale
        # Y_hat_v = (I - rho * W_v)^(-1) %*% Z_hat_v
        W_v <- W[test_idx, test_idx, drop = FALSE]
        if (nrow(W_v) > 0) {
          # Row-standardize the spatial weight matrix for the test block
          row_sums_v <- rowSums(W_v)
          row_sums_v[row_sums_v == 0] <- 1
          W_v <- W_v / row_sums_v
        }
        I_v <- diag(n_v)

        yhat_pred_v <- tryCatch({
          as.vector(MASS::ginv(I_v - fit_val$rho * W_v) %*% zhat_pred_v)
        }, error = function(e) { zhat_pred_v }) # Safe fallback if the matrix cannot be inverted

        # Squared error term ||Y_v - Y_hat_v||^2
        sq_err_v <- sum((y[test_idx] - yhat_pred_v)^2)

        # Divide by n_v: (1 / n_v) * ||...||^2
        block_squared_errors[v_idx] <- sq_err_v / n_v
      }

      if (is_valid_comb) {
        # Calculate the final SCV score for the current grid combination: (1 / |V|) * \sum ...
        scv_scores[g] <- mean(block_squared_errors)
      }
    }

    # Select the optimal combination that minimizes the SCV score
    if (all(is.na(scv_scores))) {
      stop("All candidate combinations failed during SCV loop. Check data dimensions or boundaries.")
    }

    best_idx <- which.min(scv_scores)
    best_K   <- candidate_grid$K[best_idx]
    best_h   <- candidate_grid$h[best_idx]

    if (verbose) {
      cat("Optimum parameters selected via SCV -> nbasis (K):", best_K, "| h:", best_h, "\n")
    }

    # Train the final model on the full dataset using the optimal K and h parameters
    final_model <- fsac_spls(y = y, x = x, h = best_h, nbasis = best_K,
                             gpx = gpx, wei_mat = W, rho_grid = rho_grid,
                             criterion = criterion, verbose = verbose)

    # Store SCV history details inside the final output list
    final_model$scv_history <- data.frame(
      nbasis = candidate_grid$K,
      h = candidate_grid$h,
      scv_score = scv_scores
    )
    final_model$best_nbasis <- best_K
    final_model$best_h      <- best_h

    return(final_model)
  }

  # ---- CASE 2: Standard Estimation (scv == FALSE) ----
  if (scv == FALSE) {
    # If h or nbasis are passed as vectors but scv is FALSE, select the first element and throw a warning
    if (length(h) > 1) {
      warning("scv=FALSE implies a single model fit. 'h' cannot be a vector; only the first element will be used.")
      h <- h[1]
    }
    if (!is.null(nbasis) && length(nbasis) > 1) {
      warning("scv=FALSE implies a single model fit. 'nbasis' cannot be a vector; only the first element will be used.")
      nbasis <- nbasis[1]
    }

    # Default nbasis handling for the standard fit case
    if (is.null(nbasis)) {
      nbasis <- round(min(n_total / 4, 40))
    }

    # Forward parameters directly to your original fsac_spls execution block
    return(fsac_spls(y = y, x = x, h = h, nbasis = nbasis, gpx = gpx,
                     wei_mat = W, rho_grid = rho_grid,
                     criterion = criterion, verbose = verbose))
  }
}

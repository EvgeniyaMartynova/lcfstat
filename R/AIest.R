
AIest <- function(pp, correction="Ripley", dim=NULL, dim_lims=NULL, dr=NULL, ...) {

  # For now override the correction argument and allow only a single char
  # value due to performance reasons

  # Check the value
  if (!is.character(correction)) {
    rlang::abort("'correction' argument have to have character type")
  }

  if (length(correction) > 1) {
    rlang::abort("Due to performance reasons only a single option is allowed
                 in the 'correction' argument")
  }

  if (correction == "all") {
    rlang::abort("Won't use all edge correction methods due to performance reasons,
                 choose a single option")
  }

  k_est <- spatstat.explore::Kest(pp, correction = correction, ...)
  k_est_df <- as.data.frame(k_est)

  intensity <- pp$n / spatstat.geom::area(pp$window)
  k_est_df$pn <- intensity * k_est_df[[3]]

  pn_est_df <- k_est_df[c("r", "pn")]

  # Estimated number of points is treated as a piece-wise function
  # that equals to 0 until the closest distance between points is reached
  # and then increases monotonously.
  # Only monotonously increasing part is smoothed
  first_non_zero_ind <- which(pn_est_df$pn != 0)[1]
  pn_est_non_zero <- pn_est_df[first_non_zero_ind:nrow(pn_est_df), ]

  if (is.null(dim)) {
    dim <- choose_basis_dim(pp$n, dim_lims=dim_lims)
  }

  if (is.null(dr)) {
    dr <- k_est_df$r[2] / 5
  }

  model <- scam::scam(pn ~ s(r, k=dim, bs='mpi'), data=pn_est_non_zero)
  deriv_func <- approx_deriv(pn_est_non_zero$r, model, dr)

  # Smooth a non zero part of the function
  r_arg <- pn_est_non_zero$r[2:nrow(pn_est_non_zero)]
  pn_smoothed <- predict(model, newdata=list(r=r_arg))
  pn_deriv <- deriv_func(r_arg)
  # Workaround when getting a small negative derivative
  pn_deriv <- dplyr::if_else(pn_deriv >= 0, pn_deriv, 0)

  # Concatenate part of the function where it is 0
  # and the smoothed part
  pn <- c(rep(0, first_non_zero_ind), pn_smoothed)
  pn_deriv <- c(rep(0, first_non_zero_ind), pn_deriv)

  ai_tib <- tibble::tibble(pn=pn,
                           pn_deriv=pn_deriv,
                           r=k_est_df$r,
                           theo=0)

  # TODO: change `ai` to border correction method name
  correction_name <- colnames(k_est_df)[3]
  ai_tib <- dplyr::mutate(ai_tib, ai = dplyr::if_else(pn > 0 & pn_deriv >= 0,
                                                      exp(-log(2) / 2 * r * pn_deriv / pn) * 2 - 1,
                                                      -1))

  # TODO: make an object of a class `fv`
  ai_tib <- ai_tib[c("r", "theo", "ai")]

  return(ai_tib)
}


approx_deriv <- function(x, y_model, dx) {
  # For now compute derivative with finite diff method, central diff
  # Then perhaps improve
  x_breaks <- seq(min(x) + dx / 2, max(x) + dx / 2, dx)
  y_smoothed <- predict(y_model, newdata=list(r=x_breaks))

  y_deriv_num <- diff(y_smoothed) / dx
  x_data <- x_breaks + dx / 2
  deriv_formula <- approxfun(x_data[1:(length(x_data) - 1)], y_deriv_num)

  deriv_formula
}


choose_basis_dim <- function(sample_size, dim_lims=NULL) {

  dim <- round(sqrt(sample_size))

  if (is.null(dim_lims)) {
    return(dim)
  }

  if (length(dim_lims) != 2 || !is.integer(dim_lims)) {
    rlang::abort("Incorrect dim_lims format, should a the vector of 2 integer values")
  }

  if (dim_lims[1] > dim_lims[2]) {
    rlang::abort("Incorrect dim_lims: the first limit must be lower than the second")
  }

  dim <- min(dim_lims[2], max(dim, dim_lims[1]))
  dim
}



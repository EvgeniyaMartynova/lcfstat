
#' Attraction Index
#'
#' @param pp The observed point pattern, from which an estimate of AI(r) will be computed.
#' An object of class "ppp", or data in any format acceptable to as.ppp().
#' @param correction Optional. A character vector containing one and only one of
#' the options "none", "border", "bord.modif", "isotropic", "Ripley", "translate",
#' "translation", "rigid", "none", "periodic", "good" or "best".
#' It specifies the edge correction to be applied. Note that the option "all"
#' or providing multiple edge correction methods is not supported due to
#' performance reasons. Defaults to "Ripley".
#' @param dim Optional. The dimension of the basis used to represent the smooth term within
#' the scam model formula. If not provided, the rule of thumb is used, i.e.
#' dim = sqrt(number of points).
#' @param dim_lims Optional. The integer vector with 2 values: the lower and
#' upper limits of the possible value of dim, c(lower, upper). lower > upper
#' is not allowed.
#' @param dr Optional. The delta r to compute the derivative of the estimated
#' number of points. If not provided a rule of thumb is supplied, i.e. 1/5 of
#' the distance between subsequent r values in the K-function estimate is used.
#' @param rmax Optional. Maximum desired value of the argument.
#' @param nlarge Optional. Efficiency threshold. If the number of points exceeds
#' nlarge, then only the border correction will be computed (by default),
#' using a fast algorithm.
#'
#' @return An object of class "fv", see fv.object, which can be plotted directly using plot.fv.
#' @export
#'
#' @examples
#' #TBD
AIest <- function(pp,
                  correction="Ripley",
                  dim=NULL,
                  dim_lims=NULL,
                  dr=NULL,
                  rmax=NULL,
                  nlarge=3000) {

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

  k_est <- spatstat.explore::Kest(pp,
                                  correction = correction,
                                  rmax=rmax,
                                  nlarge=nlarge)

  intensity <- pp$n / spatstat.geom::area(pp$window)
  pn <- intensity * k_est[[3]]

  pn_est_df <- data.frame(r=k_est$r,
                          pn=pn)

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
    dr <- k_est$r[2] / 5
  }

  model <- scam::scam(pn ~ s(r, k=dim, bs='mpi'), data=pn_est_non_zero)
  deriv_func <- approx_deriv(pn_est_non_zero$r, model, dr)

  # Smooth a non zero part of the function
  r_arg <- pn_est_non_zero$r[2:nrow(pn_est_non_zero)]
  pn_smoothed <- stats::predict(model, newdata=list(r=r_arg))
  pn_deriv <- deriv_func(r_arg)
  # Workaround when getting a small negative derivative
  pn_deriv <- dplyr::if_else(pn_deriv >= 0, pn_deriv, 0)

  # Concatenate part of the function where it is 0
  # and the smoothed part
  pn <- c(rep(0, first_non_zero_ind), pn_smoothed)
  pn_deriv <- c(rep(0, first_non_zero_ind), pn_deriv)
  r <- k_est$r
  ai <- dplyr::if_else(pn > 0 & pn_deriv >= 0,
                       compute_ai(r, pn, pn_deriv),
                       -1)

  ai_df <- data.frame(r=k_est$r,
                       theo=0)
  # Name the column with AI estimate after the used border correction method
  correction_name <- colnames(k_est)[3]
  ai_df[correction_name] <- ai

  ai_fv <- spatstat.explore::fv(ai_df, valu=correction_name, fname="AI", fmla = ".~r",
                                ylab=quote(AI(r)), yexp=quote(AI(r)),
                                labl=attr(k_est, "labl"),
                                desc = attr(k_est, "desc"))

  ai_fv
}


approx_deriv <- function(x, y_model, dx) {
  # For now compute derivative with finite diff method, central diff
  # Then perhaps improve
  x_breaks <- seq(min(x) + dx / 2, max(x) + dx / 2, dx)
  y_smoothed <- stats::predict(y_model, newdata=list(r=x_breaks))

  y_deriv_num <- diff(y_smoothed) / dx
  x_data <- x_breaks + dx / 2
  deriv_formula <- stats::approxfun(x_data[1:(length(x_data) - 1)], y_deriv_num)

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


compute_ai <- function(r, pn, pn_deriv, low=-1, high=1) {
  scale <- high - low
  shift <- low
  ai <- exp(-log(2) / 2 * r * pn_deriv / pn) * scale + shift
  ai
}

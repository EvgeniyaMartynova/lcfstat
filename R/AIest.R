#' Attraction Index
#'
#'Estimates the Attraction Index function from a point pattern in a window of arbitrary shape.
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
#' @param rmax Optional. Maximum desired value of the argument.
#' @param nlarge Optional. Efficiency threshold. If the number of points exceeds
#' nlarge, then only the border correction will be computed (by default),
#' using a fast algorithm.
#'
#' @return An object of class "aifv", inherited from fv.object, which can be
#' plotted directly using plot.aifv.
#' @export
#'
#' @examples
#'
#' library(spatstat.random)
#'
#' # AI for a random point pattern
#' rpp <- rpoispp(500)
#'
#' ai_rand <- AIest(rpp)
#' plot(ai_rand, main = "AI for a random point pattern")
#'
#' ai <- AIest(rpp, "border")
#' plot(ai, main = "AI for a random point pattern")
#'
#' # AI for a clustered point pattern
#' clust_pp <- rMatClust(20, 0.05, 25)
#'
#' ai_clust <- AIest(clust_pp)
#' plot(ai_clust, main = "AI for a clustered point pattern")
#'
#' # AI for a point pattern with dispersion
#' hardcore_pp <- rHardcore(300, R=0.05)
#' ai_disp <- AIest(hardcore_pp)
#' plot(ai_disp, main = "AI for a point pattern with dispersion")
#'
#' # Plot AI for three different point pattern together
#' plot(ai_rand$r, ai_rand$iso, type="l", ylim=c(-1, 1), col=4)
#' lines(ai_rand$r, rep(0, nrow(ai_rand)), lty=2)
#' lines(ai_clust$r, ai_clust$iso, col=2)
#' lines(ai_disp$r, ai_disp$iso, col=7)
#' legend("bottomright",
#'        c("theoretical", "random", "clustered", "dispersed"),
#'        col=c(1, 4, 2, 7),
#'        lty=c(2, 1, 1, 1))
#'
AIest <- function(pp,
                  correction="Ripley",
                  dim=NULL,
                  dim_lims=NULL,
                  rmax=NULL,
                  nlarge=3000) {

  # For now override the correction argument and allow only a single char
  # value due to performance reasons

  # Check the value
  if (!is.character(correction) || length(correction) > 1) {
    rlang::abort(class = "ai_error_bad_correction",
                 message="'correction' argument have to be a character vector
                 with a single value")
  }

  if (correction == "all") {
    rlang::abort(class = "ai_error_inefficient_correction",
                 message="Won't use all edge correction methods due to performance reasons,
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
  # Handling the situation when rmax is too large and Kest returns NAs
  last_ind <- min(nrow(pn_est_df), which(is.na(pn_est_df$pn))[1] - 1, na.rm = TRUE)
  pn_est_defined <- pn_est_df[first_non_zero_ind:last_ind, ]

  if (is.null(dim)) {
    dim <- choose_basis_dim(pp$n, dim_lims=dim_lims)
  }

  model <- scam::scam(pn ~ s(r, k=dim, bs='mpi'), data=pn_est_defined)

  # Smooth a non zero part of the function
  r_arg <- pn_est_defined$r
  pn <- stats::predict(model, newdata=list(r=r_arg))
  pn_deriv <- single_mpi_derivative(model, data=list(r=r_arg))

  # Workaround when getting a small negative derivative (TODO: is it needed now when derivative is analytical?)
  pn_deriv <- dplyr::if_else(pn_deriv >= 0, pn_deriv, 0)

  ai <- dplyr::if_else(pn > 0 & pn_deriv >= 0,
                       compute_ai(r_arg, pn, pn_deriv),
                       -1)

  # Pad the AI value with -1 at distances where there is no neigbours
  # and with NA when the K function is undefined
  num_ll_pad <- first_non_zero_ind - 1
  num_na_pad <- nrow(k_est) - last_ind
  ai <- c(rep(-1, num_ll_pad), ai, rep(NA_real_, num_na_pad))

  ai_df <- data.frame(r=k_est$r,
                      theo=0)

  # Name the column with AI estimate after the used border correction method
  correction_name <- colnames(k_est)[3]
  ai_df[correction_name] <- ai

  ai_fv <- spatstat.explore::fv(ai_df, valu=correction_name, fname="AI", fmla = ".~r",
                                ylab=quote(AI(r)), yexp=quote(AI(r)),
                                labl=attr(k_est, "labl"),
                                desc = attr(k_est, "desc"))

  class(ai_fv) <- c("aifv", class(ai_fv))

  ai_fv
}


#' Derivative of a scam model with a single MPI term
#'
#' @param model A scam model with a single monotone increasing P-spline (MPI)
#' term
#' @param data A data frame containing the values of the named covariates
#' at which the smooth term is to be evaluated.
#'
#' @return A vector that contains the derivative of the given scam model at
#' the valuse provided with the data parameter
single_mpi_derivative <- function(model, data) {

  # check that it is a scam objects
  if (!inherits(model, "scam")) {
    rlang::abort(class = "ai_derivative_error_unsupported_model",
                 message = "Only SCAM objects are supported")
  }

  # check that there is only one smooth and it is an mpi.smooth
  if (length(model$smooth) != 1 ||
      !inherits(model$smooth[[1]], "mpi.smooth") ||
      length(model$assign) > 1) {
    rlang::abort(class = "ai_derivative_error_unsupported_model_formula",
                 message ="Only SCAMs with a slingle mpi smooth terms are supported")
  }

  # taken from Predict.matrix.mpi.smooth with a modification to extract the derivative
  smooth <- model$smooth[[1]]

  if (!(smooth$term %in% names(data))) {
    rlang::abort(class = "ai_derivative_error_missing_argument",
                 message = "A required argument is not provided in the data")
  }

  order <- smooth$m + 2
  q <- smooth$df + 1
  Sig <- matrix(1, q, q)
  Sig[upper.tri(Sig)] <- 0
  ll <- smooth$knots[order]
  ul <- smooth$knots[length(smooth$knots) - order + 1]
  x <- data[[smooth$term]]
  n <- length(x)
  ind <- x <= ul & x >= ll

  if (sum(ind) != n) {
    rlang::abort(class = "ai_derivative_error_invalid_argument",
                 message = "Won't evaluate the derivative outside of the function domain")
  }

  Xdm <- splines::splineDesign(smooth$knots, x, order, derivs = c(1))
  gamma_coef <- Sig %*% model$coefficients.t

  Xd <- Xdm %*% gamma_coef
  c(Xd)
}


#' The number of B-splines to use in the scam model
#'
#' @param sample_size Number of points in a point pattern
#' @param dim_lims Optional. The integer vector with 2 values: the lower and
#' upper limits of the possible value of dim, c(lower, upper). lower > upper
#' is not allowed.
#'
#' @return An integer value that represents the number of B-splines to use
#' in a scam model. The rule of thumb dim = sqrt(sample_size) is used and
#' if dim_lims is provided it is used to clip the value
choose_basis_dim <- function(sample_size, dim_lims=NULL) {

  if (length(sample_size) > 1 || sample_size < 0 || sample_size %% 1 != 0) {
    rlang::abort(class = "ai_dim_error_invalid_sample_size",
                 message = "Sample size should contain a single positive integer value")
  }

  dim <- round(sqrt(sample_size))

  if (is.null(dim_lims)) {
    return(dim)
  }

  if (length(dim_lims) != 2 || sum(dim_lims %% 1 != 0) > 0 || sum(dim_lims < 0) > 0) {
    rlang::abort(class = "ai_dim_error_invalid_dim_lims",
                 message = "Incorrect dim_lims format, should be a vector of 2 integer values")
  }

  if (dim_lims[1] > dim_lims[2]) {
    rlang::abort(class = "ai_dim_error_unordered_dim_lims",
                 message = "Incorrect dim_lims, the first limit must be lower than the second")
  }

  dim <- min(dim_lims[2], max(dim, dim_lims[1]))
  dim
}

#' The AI value computation
#'
#' @param r A vector of distances
#' @param pn A vector of estimated expected number of points within a distance r
#' @param pn_deriv A vector of estimated derivative of the expected number of points
#' within a distance r
#' @param ai_lims Optional. The double vector with 2 values: the lower and
#' upper limits of the AI, c(lower, upper). The lower value corresponds to the AI value
#' for maximal dispersion, the upper to the AI value for maximal clustering.
#' lower > upper is not allowed.
#'
#' @returns A vector with AI estimate at r
compute_ai <- function(r, pn, pn_deriv, ai_lims=c(-1, 1)) {

  if (length(ai_lims) != 2 || !is.numeric(ai_lims)) {
    rlang::abort(class = "ai_compute_error_invalid_ai_lims",
                 message = "Incorrect ai_lims format, should be a vector of 2 double values")
  }

  if (ai_lims[1] > ai_lims[2]) {
    rlang::abort(class = "ai_compute_error_unordered_ai_lims",
                 message = "Incorrect ai_lims, the first limit must be lower than the second")
  }

  if (!is.numeric(r) || !is.numeric(pn) || !is.numeric(pn_deriv)) {
    rlang::abort(class = "ai_compute_error_invalid_arg",
                 message = "r, pn and pn_deriv have to be double vectors")
  }

  if (length(r) != length(pn) || length(pn) != length(pn_deriv)) {
    rlang::abort(class = "ai_compute_error_arg_length_mismatch",
                 message = "r, pn and pn_deriv vectors should have the same length")
  }

  scale <- ai_lims[2] - ai_lims[1]
  shift <- ai_lims[1]
  ai <- exp(-log(2) / 2 * r * pn_deriv / pn) * scale + shift
  ai
}

#' Plot Function Value for AI
#'
#' @param x An object of the class "aifv" that contains the variables to be
#' plotted.
#' @param main A title of the plot.
#' @param ... Extra arguments passed to the spatstat.explore::plot.fv
#'
#' @return Invisible: either NULL, or a data frame giving the meaning of the
#' different line types and colours.
#' @export
#' @export plot.aifv
#'
#' @examples
#' library(spatstat.random)
#'
#' rpp <- rpoispp(500)
#' ai_rand <- AIest(rpp)
#' plot(ai_rand, main = "AI for a random point pattern")
#'
plot.aifv <- function(x, main="", ...) {
  spatstat.explore::plot.fv(x, ylim = c(-1,1), main = main, ...)
}

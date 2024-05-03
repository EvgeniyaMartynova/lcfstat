#' Local Correlation Function
#'
#'Estimates the local correlation function function from a point pattern in a window of arbitrary shape.
#'
#' @param pp The observed point pattern, from which an estimate of LCF(r) will be computed
#' or N(r) estimate. In case of the point pattern, it should be an object of class "ppp",
#' or data in any format acceptable to as.ppp(). For N(r), a data frame with columns
#' "r" (distance) and "pn" (estimated number of points)
#' @param correction Optional. A character vector containing one and only one of
#' the options "none", "border", "bord.modif", "isotropic", "Ripley", "translate",
#' "translation", "rigid", "none", "periodic", "good" or "best".
#' It specifies the edge correction to be applied. Note that the option "all"
#' or providing multiple edge correction methods is not supported due to
#' performance reasons. Defaults to "Ripley".
#' @param r Optional. Vector of values for the argument r at which LCF(r)
#' should be evaluated. The values must be in increasing order. Advanced use only.
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
#' @return An object of class "lcffv", inherited from fv.object, which can be
#' plotted directly using plot.lcffv.
#' @export
#'
#' @examples
#'
#' library(spatstat.random)
#'
#' # LCF for a random point pattern
#' rpp <- rpoispp(500)
#'
#' lcf_rand <- LCFest(rpp)
#' plot(lcf_rand, main = "LCF for a random point pattern")
#'
#' lcf <- LCFest(rpp, "border")
#' plot(lcf, main = "LCF for a random point pattern")
#'
#' # LCF for a clustered point pattern
#' clust_pp <- rMatClust(20, 0.05, 25)
#'
#' lcf_clust <- LCFest(clust_pp)
#' plot(lcf_clust, main = "LCF for a clustered point pattern")
#'
#' # LCF for a point pattern with dispersion
#' hardcore_pp <- rHardcore(300, R=0.05)
#' lcf_disp <- LCFest(hardcore_pp)
#' plot(lcf_disp, main = "LCF for a point pattern with dispersion")
#'
#' # Plot LCF for three different point pattern together
#' plot(lcf_rand$r, lcf_rand$iso, type="l", ylim=c(-1, 1), col=4)
#' lines(lcf_rand$r, rep(0, nrow(lcf_rand)), lty=2)
#' lines(lcf_clust$r, lcf_clust$iso, col=2)
#' lines(lcf_disp$r, lcf_disp$iso, col=7)
#' legend("bottomright",
#'        c("theoretical", "random", "clustered", "dispersed"),
#'        col=c(1, 4, 2, 7),
#'        lty=c(2, 1, 1, 1))
#'
LCFest <- function(pp,
                   correction="Ripley",
                   r=NULL,
                   dim=NULL,
                   dim_lims=NULL,
                   rmax=NULL,
                   nlarge=3000) {

  # For now override the correction argument and allow only a single char
  # value due to performance reasons

  if (inherits(pp, "ppp")) {
    # Check the value
    if (!is.character(correction) || length(correction) > 1) {
      rlang::abort(class = "lcf_error_bad_correction",
                   message="'correction' argument has to be a character vector
                 with a single value")
    }

    if (correction == "all") {
      rlang::abort(class = "lcf_error_inefficient_correction",
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

    if (is.null(dim)) {
      dim <- choose_basis_dim(pp$n, dim_lims=dim_lims)
    }
  } else if (inherits(pp, "data.frame") && ("r" %in% colnames(pp)) && ("pn" %in% colnames(pp))) {
    pn_est_df <- pp

    if (is.null(dim)) {
      rlang::abort(class = "lcf_error_no_dim",
                   message="When point number estimate is supplied, \"dim\" argument is mandatory.")
    }
  } else {
    rlang::abort(class = "lcf_error_invalid_arg",
                 message="pp shoud be either an object of class \"ppp\" or a data
                 frame with the point number estimate and colums \"r\" and \"pn\"")
  }

  # Name the column with LCF estimate after the used border correction method
  correction_name <- if (exists("k_est")) colnames(k_est)[3] else "empirical"
  lcf_df <- LCF(pn_est_df, r, dim, correction_name)

  lcf_labl <- if (exists("k_est")) attr(k_est, "labl") else c("r", "%s[pois](r)", "hat(%s)[empirical](r)")
  lcf_desc <- if (exists("k_est")) attr(k_est, "desc") else c("distance argument r", "theoretical Poisson %s", "empirical estimate of %s")

  lcf_fv <- spatstat.explore::fv(lcf_df, valu=correction_name, fname="LCF", fmla = ".~r",
                                 ylab=quote(LCF(r)), yexp=quote(LCF(r)),
                                 labl=lcf_labl,
                                 desc = lcf_desc)

  class(lcf_fv) <- c("lcffv", class(lcf_fv))

  lcf_fv
}

LCFcross <- function(pp,
                     i,
                     j,
                     correction="Ripley",
                     r=NULL,
                     dim=NULL,
                     dim_lims=NULL,
                     rmax=NULL) {

  if (!inherits(pp, "ppp") && "marks" %in% attributes(pp)$names) {
    rlang::abort(class = "lcf_cross_error_invalid_arg",
                 message="'pp' argument has to be a multitype point pattern")
  }

  # Check the value
  if (!is.character(correction) || length(correction) > 1) {
    rlang::abort(class = "lcf_error_bad_correction",
    message="'correction' argument has to be a character vector
             with a single value")
  }

  if (correction == "all") {
    rlang::abort(class = "lcf_error_inefficient_correction",
                 message="Won't use all edge correction methods due to performance reasons,
                          choose a single option")
  }

  # TODO: try catch for as.character(i)

  k_ij <- spatstat.explore::Kcross(pp,
                                   i, j,
                                   correction = correction,
                                   rmax=rmax)

  n_j <- pp[pp$marks == j]$n
  intensity_j <- n_j / spatstat.geom::area(pp$window)
  pn <- intensity_j * k_ij[[3]]

  pn_est_df <- data.frame(r=k_ij$r,
                          pn=pn)

  if (is.null(dim)) {
    n_i <- pp[pp$marks == i]$n
    n_geom <- round(sqrt(n_i * n_j))
    dim <- choose_basis_dim(n_geom, dim_lims=dim_lims)
  }

  # Name the column with LCF estimate after the used border correction method
  correction_name <- names(k_ij)[3]

  lcf_df <- LCF(pn_est_df, r, dim, correction_name)

  lcf_labl <- attr(k_ij, "labl")
  lcf_desc <- attr(k_ij, "desc")
  lcf_exp <- substitute("LCF"[list(i, j)](r), list(i=i, j=j))
  lcf_lab <- substitute("LCF"[i, j](r), list(i=i, j=j))

  lcf_fv <- spatstat.explore::fv(lcf_df, argu = "r", valu=correction_name,
                                 fname=c("LCF", paste0("list(", i, ",", j, ")")),
                                 fmla = ".~r", ylab=lcf_lab, yexp=lcf_exp,
                                 labl=lcf_labl, desc=lcf_desc)

  class(lcf_fv) <- c("lcffv", class(lcf_fv))
  lcf_fv
}

LCF <- function(pn_est_df, r, dim, lcf_name) {

  # Estimated number of points is treated as a piece-wise function
  # that equals to 0 until the closest distance between points is reached
  # and then increases monotonously.
  # Only monotonously increasing part is smoothed
  first_non_zero_ind <- which(pn_est_df$pn != 0)[1]
  # Handling the situation when rmax is too large and Kest returns NAs
  last_ind <- min(nrow(pn_est_df), which(is.na(pn_est_df$pn))[1] - 1, na.rm = TRUE)
  pn_est_defined <- pn_est_df[first_non_zero_ind:last_ind, ]

  model <- scam::scam(pn ~ s(r, k=dim, bs='mpi'), data=pn_est_defined)

  # Smooth a non zero part of the function
  if (is.null(r)) {
    r <- pn_est_df$r
  }

  r_non_increasing <- any(diff(r) <= 0)
  if (r_non_increasing) {
    rlang::abort(class = "lcf_error_bad_r",
                 message = "r values should be increasing")
  }

  r_li <- which(r >= pn_est_defined$r[1])[1]
  r_hi <- min(length(r), which(r > pn_est_defined$r[nrow(pn_est_defined)])[1] - 1, na.rm = TRUE)

  if (!is.na(r_li) & r_li <= r_hi) {
    r_def <- r[r_li:r_hi]

    pn <- stats::predict(model, newdata=list(r=r_def))
    pn_deriv <- single_mpi_derivative(model, data=list(r=r_def))

    # Workaround when getting a small negative derivative (TODO: is it needed now when derivative is analytical?)
    pn_deriv <- dplyr::if_else(pn_deriv >= 0, pn_deriv, 0)

    lcf <- dplyr::if_else(pn > 0 & pn_deriv >= 0,
                          compute_lcf(r_def, pn, pn_deriv),
                          -1)

    num_ll_pad <- r_li - 1
    num_na_pad <- length(r) - r_hi
  } else {
    lcf <- NULL
    if (is.na(r_li)) {
      num_ll_pad <- length(r)
      num_na_pad <- 0
    } else {
      num_ll_pad <- 0
      num_na_pad <- length(r)
    }
  }

  # Pad the LCF value with -1 at distances where there is no neighbors
  # and with NA when the K function is undefined
  lcf <- c(rep(-1, num_ll_pad), lcf, rep(NA_real_, num_na_pad))
  lcf_df <- data.frame(r=r,
                       theo=0)
  lcf_df[lcf_name] <- lcf

  lcf_df
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
    rlang::abort(class = "lcf_derivative_error_unsupported_model",
                 message = "Only SCAM objects are supported")
  }

  # check that there is only one smooth and it is an mpi.smooth
  if (length(model$smooth) != 1 ||
      !inherits(model$smooth[[1]], "mpi.smooth") ||
      length(model$assign) > 1) {
    rlang::abort(class = "lcf_derivative_error_unsupported_model_formula",
                 message ="Only SCAMs with a slingle mpi smooth terms are supported")
  }

  # taken from Predict.matrix.mpi.smooth with a modification to extract the derivative
  smooth <- model$smooth[[1]]

  if (!(smooth$term %in% names(data))) {
    rlang::abort(class = "lcf_derivative_error_missing_argument",
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
    rlang::abort(class = "lcf_derivative_error_invalid_argument",
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
    rlang::abort(class = "lcf_dim_error_invalid_sample_size",
                 message = "Sample size should contain a single positive integer value")
  }

  dim <- round(sqrt(sample_size))

  if (is.null(dim_lims)) {
    return(dim)
  }

  if (length(dim_lims) != 2 || sum(dim_lims %% 1 != 0) > 0 || sum(dim_lims < 0) > 0) {
    rlang::abort(class = "lcf_dim_error_invalid_dim_lims",
                 message = "Incorrect dim_lims format, should be a vector of 2 integer values")
  }

  if (dim_lims[1] > dim_lims[2]) {
    rlang::abort(class = "lcf_dim_error_unordered_dim_lims",
                 message = "Incorrect dim_lims, the first limit must be lower than the second")
  }

  dim <- min(dim_lims[2], max(dim, dim_lims[1]))
  dim
}

#' The LCF value computation
#'
#' @param r A vector of distances
#' @param pn A vector of estimated expected number of points within a distance r
#' @param pn_deriv A vector of estimated derivative of the expected number of points
#' within a distance r
#' @param lcf_lims Optional. The double vector with 2 values: the lower and
#' upper limits of the LCF, c(lower, upper). The lower value corresponds to the LCF value
#' for maximal dispersion, the upper to the LCF value for maximal clustering.
#' lower > upper is not allowed.
#'
#' @returns A vector with LCF estimate at r
compute_lcf <- function(r, pn, pn_deriv, lcf_lims=c(-1, 1)) {

  if (length(lcf_lims) != 2 || !is.numeric(lcf_lims)) {
    rlang::abort(class = "lcf_compute_error_invalid_lcf_lims",
                 message = "Incorrect lcf_lims format, should be a vector of 2 double values")
  }

  if (lcf_lims[1] > lcf_lims[2]) {
    rlang::abort(class = "lcf_compute_error_unordered_lcf_lims",
                 message = "Incorrect lcf_lims, the first limit must be lower than the second")
  }

  if (!is.numeric(r) || !is.numeric(pn) || !is.numeric(pn_deriv)) {
    rlang::abort(class = "lcf_compute_error_invalid_arg",
                 message = "r, pn and pn_deriv have to be double vectors")
  }

  if (length(r) != length(pn) || length(pn) != length(pn_deriv)) {
    rlang::abort(class = "lcf_compute_error_arg_length_mismatch",
                 message = "r, pn and pn_deriv vectors should have the same length")
  }

  scale <- lcf_lims[2] - lcf_lims[1]
  shift <- lcf_lims[1]
  lcf <- exp(-log(2) / 2 * r * pn_deriv / pn) * scale + shift
  lcf
}

#' Plot Function Value for LCF
#'
#' @param x An object of the class "lcffv" that contains the variables to be
#' plotted.
#' @param ylim (optional) range of y axis. Default is set to the lower and
#' upper limits of the LCF.
#' @param main (optional) A title of the plot.
#' @param ... Extra arguments passed to the spatstat.explore::plot.fv
#'
#' @return Invisible: either NULL, or a data frame giving the meaning of the
#' different line types and colours.
#' @export
#' @export plot.lcffv
#'
#' @examples
#' library(spatstat.random)
#'
#' rpp <- rpoispp(500)
#' lcf_rand <- LCFest(rpp)
#' plot(lcf_rand, main = "LCF for a random point pattern")
#'
plot.lcffv <- function(x, ylim = c(-1,1), main="", ...) {
  spatstat.explore::plot.fv(x, ylim = ylim, main = main, ...)
}

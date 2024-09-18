#' Local Correlation Function
#'
#'Estimates the local correlation function from a point pattern in a window of arbitrary shape.
#'
#' The local correlation function, \eqn{LCF(r)}, is a summary function of spatial statistics which estimates
#' a degree of clustering or dispersion at distance \eqn{r} in a point process (Martynova, 2024). It is based on
#' the influential Ripley's \eqn{K}-function, \code{\link[spatstat.explore]{Kest}}, or more precisely, on the
#' underlying summary function \eqn{N(r)}, the expected number of random points within a distance \eqn{r}
#' of a typical random point.
#' The following properties make it possible to interpret LCF as a degree of clustering or dispersion:
#'
#' \enumerate{
#'   \item \eqn{LCF(r)} is limited to the range \eqn{[-1,1]},
#'   \item \eqn{LCF(r)} is asymptotically \eqn{0} for CSR,
#'   \item \eqn{LCF(r) = 1} for maximal clustering (e.g. for \eqn{r > r_d}  if all points are concentrated inside the disc with the diameter \eqn{r_d}),
#'   \item \eqn{LCF(r) = -1} for maximal dispersion (e.g. lattice at \eqn{r} smaller than the distance between nearest neighbor),
#' }
#'
#' The scale of LCF is inspired by the linear correlation coefficient, and highlights the
#' opposite meaning of clustering and dispersion.
#' The advantage of the LCF in comparison to other summary functions of spatial statistics is
#' that it is uniformly bounded, which means that its extreme values \eqn{-1} and \eqn{1} are identical
#' and attainable for all scales \eqn{r}. Altogether, this makes LCF an attractive metric
#' for researches who are interested in extracting spatial information in an interpretable way.
#' As with other \eqn{K}-based spatial statistics, researchers can analyze the entire
#' function graph of LCF for a single or a few patterns or construct summary measures for joint analysis
#' of many point patterns, e.g. functon value at a certain distance or the area under the curve.
#'
#' Formally, maximal clustering and dispersion can be defined using the the expected number of additional
#' random points within a distance \eqn{r} of a typical random point, \eqn{N(r)}.
#'
#' \itemize{
#'   \item A point process with \eqn{N(r) = 0} is \emph{maximally dispersed} at distance \eqn{r}.
#'   \item If \eqn{N(r) > 0} and \eqn{N(r) = N(hr)} for some \eqn{h > 1}, then a point process
#'   is \emph{maximally clustered at distance} \eqn{r}.
#' }
#'
#' LCF can be estimated as
#'
#' \deqn{LCF(r) =  \begin{cases} 2\exp{ \left(  - \frac{\ln{2}}{2} \frac{r\, N'(r)}{N(r)} \right)} - 1, \quad N(r) > 0  \\ -1, \quad N(r) = 0 \end{cases}}
#'
#' In this package, we extract the estimate of \eqn{N(r)} using the \code{\link[spatstat.explore]{Kest}}
#' function from the \code{spatstat} package and obtains its smooth monotonically increasing spline
#' approximation with \code{\link[scam]{scam}} package (Pya, 2015) to calculate the derivative, \eqn{N'(r)}.
#'
#' Importantly, LCF is bounded under the assumption that the estimated \eqn{N(r)} is monotonically
#' non-decreasing. Although, this is typically the case, there might be some edge cases
#' when the appropriate edge correction method should be chosen to preserve this property.
#'

#' @param pp The observed point pattern, from which an estimate of \eqn{LCF(r)} will be computed
#' or the estimated average number of neighbors of a point \eqn{N(r)}. In case of the point pattern, it should be an object of class "ppp".
#' For \eqn{N(r)}, a data frame with columns \code{r} (distance) and \code{pn} (estimated number of points) is expected.
#' @param correction Optional. A string containing one of
#' the options "\code{none}", "\code{border}", "\code{bord.modif}", "\code{isotropic}",
#' "\code{Ripley}", "\code{translate}", "\code{translation}", "\code{rigid}",
#' "\code{periodic}", "\code{good}" or "\code{best}".
#' It specifies the edge correction to be applied. Note that the option "all"
#' or providing multiple edge correction methods is not supported due to
#' performance reasons. Defaults to "\code{Ripley}"/"\code{isotropic}".
#' @param r Optional. Vector of values for the argument \eqn{r} at which \eqn{LCF(r)}
#' should be evaluated. The values must be in increasing order. Advanced use only.
#' @param dim Optional. The dimension of the basis used to represent the smooth term within
#' the scam model formula. If not provided, the rule of thumb is used, i.e.
#' \eqn{dim = \sqrt{n}} where \eqn{n} is the number of points. Must be provided when the estimate
#' of \eqn{N(r)} is passed to \code{pp}.
#' @param dim_lims Optional. The integer vector of 2 values: the lower and
#' upper limits of the possible value of dim, \code{c(lower, upper)}. \code{lower} > \code{upper}
#' is not allowed. Only applied when \code{dim} is not provided and computed with the rule of thumb.
#' Clips the computed \eqn{dim} to \eqn{[lower, upper]} interval.
#' @param rmax Optional. Maximum desired value of the argument.
#' @param nlarge Optional. Efficiency threshold. If the number of points exceeds
#' \code{nlarge}, then only the border correction will be computed (by default),
#' using a fast algorithm. The default is 3000.
#'
#' @return An object of class "lcffv", inherited from \code{\link[spatstat.explore]{fv.object}}, which can be
#' plotted directly using \code{\link{plot.lcffv}}. Essentially a data frame that contains
#' columns:
#'
#' \code{r} --     the value of the argument \eqn{r} at which \eqn{LCF(r)} is estimated
#'
#' \code{theo} --   the theoretical value of \eqn{LCF(r)} for Poisson process, \eqn{0}
#'
#' and a column named after the used edge correction method, e.g. "\code{iso}" or "\code{border}",
#' with the empirical estimate of \eqn{LCF(r)} obtained from the point pattern.
#'
#' @author Evgenia Martynova \email{evg.martynova@@gmail.com}
#'
#' @references Martynova, E. and Textor, J., 2024, August.
#' A Uniformly Bounded Correlation Function for Spatial Point Patterns.
#' In \emph{Proceedings of the 30th ACM SIGKDD Conference on Knowledge Discovery and Data Mining} (pp. 2177-2188).
#'
#' @references Pya, N. and Wood, S.N., 2015.
#' Shape constrained additive models. In \emph{Statistics and Computing, 25(3)} (pp. 543-559).
#'
#'@seealso \code{\link{LCFcross}} to estimate the colocalization between two types of objects
#' in a multitype point pattern.
#
#'
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
#' # Plot LCF for three different point patterns together
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

  if (inherits(pp, "ppp")) {
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

    r_arg <- get_r_arg(r)
    k_est <- spatstat.explore::Kest(pp,
                                    correction = correction,
                                    rmax=rmax,
                                    r=r_arg,
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

  # Name the column with LCF estimate after the used edge correction method
  correction_name <- if (exists("k_est")) colnames(k_est)[3] else "empirical"
  lcf_df <- LCF(pn_est_df, r, dim, correction_name)

  lcf_labl <- if (exists("k_est")) attr(k_est, "labl") else c("r", "%s[pois](r)", "hat(%s)[empirical](r)")
  lcf_desc <- if (exists("k_est")) attr(k_est, "desc") else c("distance argument r", "theoretical Poisson %s", "empirical estimate of %s")

  lcf_fv <- spatstat.explore::fv(lcf_df, valu=correction_name, fname="LCF", fmla = ".~r",
                                 ylab=quote(LCF(r)), yexp=quote(LCF(r)),
                                 labl=lcf_labl, desc = lcf_desc)

  class(lcf_fv) <- c("lcffv", class(lcf_fv))

  lcf_fv
}

#' Multitype LCF (Cross-type)
#'
#' For a multitype point pattern, estimate the cross-LCF which estimates the spatial
#' distribution of objects of type \eqn{j} with respect to objects of type \eqn{i}.
#'
#' This function is variant of the function \code{\link{LCFest}} extended to estimate
#' the colocalisation of two types of points in multitype point patterns.
#'
#' In a multitype point pattern, points can be classified into a finite number types.
#' We follow the conventions of the \code{spatstat} package, which represents a multitype pattern
#' as a single pattern of points with marks that determine the type of points.
#'
#' The argument pp must be a marked point pattern (object of class "\code{ppp}") with the mark vector
#' \code{pp$marks} of a \code{factor} type.
#'
#' The "cross-LCF" can be estimated by substituting \eqn{N(r)} (the expected number
#' of neighbors within distance \eqn{r} of a typical point) in the formula of LCF with
#' \eqn{N_{ij}(r)} - the expected number of points of type \eqn{j} within distance \eqn{r}
#' of a typical point of type \eqn{i}.
#'
#' \deqn{LCF_{ij}(r) =  \begin{cases} 2\exp{ \left(  - \frac{\ln{2}}{2} \frac{r\, N_{ij}'(r)}{N_{ij}(r)} \right)} - 1, \quad N_{ij}(r) > 0  \\ -1, \quad N_{ij}(r) = 0 \end{cases}}
#'
#' It estimates the degree of clustering of the points of type \eqn{j} around the
#' points of type \eqn{i}. If the process that generate points \eqn{i} and \eqn{j} are
#' independent, \eqn{LCF_{ij}} is asymptotically \eqn{0}. \eqn{LCF_{ij} > 0} suggest clustering of
#' points of type \eqn{j} around the points of type \eqn{i}, while \eqn{LCF_{ij} < 0} suggests dispersion.
#'
#' @param pp The observed point pattern, from which an estimate of the cross-LCF, \eqn{LCF_{ij}(r)},
#' will be computed. It must be a multitype point pattern (a marked point pattern whose marks are a factor).
#' See the details of \code{\link[spatstat.explore]{Kcross}}.
#' @param i The type (mark value) of the points in \code{pp} from which distances are measured.
#' Must be a character string. Defaults to the first level of \code{marks(pp)}.
#' @param j The type (mark value) of the points in \code{pp} to which distances are measured.
#' Must be a character string. Defaults to the second level of \code{marks(pp)}.
#' @param correction Optional. A string containing one of
#' the options "\code{none}", "\code{border}", "\code{bord.modif}", "\code{isotropic}",
#' "\code{Ripley}", "\code{translate}", "\code{translation}", "\code{rigid}",
#' "\code{periodic}", "\code{good}" or "\code{best}".
#' It specifies the edge correction to be applied. Note that the option "all"
#' or providing multiple edge correction methods is not supported due to
#' performance reasons. Defaults to "\code{Ripley}"/"\code{isotropic}".
#' @param r Optional. Vector of values for the argument \eqn{r} at which cross-LCF(r)
#' should be evaluated. The values must be in increasing order. Advanced use only.
#' @param dim Optional. The dimension of the basis used to represent the smooth term within
#' the scam model formula. If not provided, the rule of thumb is used, i.e.
#' \eqn{dim = \sqrt{n}} where \eqn{n} is the total number of points in a point pattern.
#' @param dim_lims Optional. The integer vector of 2 values: the lower and
#' upper limits of the possible value of dim, \code{c(lower, upper)}. \code{lower} > \code{upper}
#' is not allowed. Only applied when \code{dim} is not provided and computed with the rule of thumb.
#' Clips the computed \eqn{dim} to \eqn{[lower, upper]} interval.
#' @param rmax Optional. Maximum desired value of the argument.
#'
#' @return An object of class "lcffv", inherited from \code{\link[spatstat.explore]{fv.object}}, which can be
#' plotted directly using \code{\link{plot.lcffv}}. Essentially a data frame that contains
#' columns:
#'
#' \code{r} --     the value of the argument \eqn{r} at which \eqn{LCF(r)} is estimated
#'
#' \code{theo} --   the theoretical value of \eqn{LCF_{ij}(r)} for Poisson process, \eqn{0}
#'
#' and a column named after the used edge correction method, e.g. "\code{iso}" or "\code{border}",
#' with the empirical estimate of \eqn{LCF_{ij}(r)} obtained from the point pattern.
#'
#' @author Evgenia Martynova \email{evg.martynova@@gmail.com}
#'
#' @references Martynova, E. and Textor, J., 2024, August.
#' A Uniformly Bounded Correlation Function for Spatial Point Patterns.
#' In \emph{Proceedings of the 30th ACM SIGKDD Conference on Knowledge Discovery and Data Mining} (pp. 2177-2188).
#'
#'
#' @export
#'
#' @examples
#'
#' library(spatstat.random)
#' library(spatstat.geom)
#'
#' # This example compares cross-LCF for point patterns
#' # with two types of objects and different colocalization
#' # of these objects.
#'
#' # No dependence between positions of points of two types
#' # Draw a point pattern from CSR
#' pp_rand <- rpoispp(1000)
#' # Randomly assign two types of points
#' pp_rand <- pp_rand %mark% factor(sample(0:1, npoints(pp_rand), replace=TRUE))
#' plot(pp_rand)
#'
#' lcf_rand <- LCFcross(pp_rand, "0", "1")
#'
#' # Points of type 1 are clustered around points of type 0
#' # Simulated using the Matern cluster process
#' # and assigning type 0 to parent points
#' # and type 1 to offspring.
#' clust_rad <- 0.025
#' mat_clust <- rMatClust(25, clust_rad, 50,
#'                        saveparents = TRUE)
#'
#' parents <- attr(mat_clust, "parents")
#'
#' x <- c(parents$x, mat_clust$x)
#' y <- c(parents$y, mat_clust$y)
#' marks <- factor(c(rep(0, length(parents$x)), rep(1, mat_clust$n)))
#'
#' pp_attr <- ppp(x, y, marks = marks)
#' plot(pp_attr)
#'
#' lcf_attr <- LCFcross(pp_attr, "0", "1")
#'
#' # Points of different types do not mix.
#' # The points with type 0 occupy the left side of the window
#' win1 <- owin(c(0, 0.49), c(0, 1))
#' pp1 <- rpoispp(500, win=win1)
#'
#' # The points with type 1 occupy the right side of the window
#' win2 <- owin(c(0.51, 1), c(0, 1))
#' pp2 <- rpoispp(500, win=win2)
#'
#' # Combine them into a mutitype point pattern
#' xs <- c(pp1$x, pp2$x)
#' ys <- c(pp1$y, pp2$y)
#' marks <- factor(c(rep("0", pp1$n), rep("1", pp2$n)))
#'
#' pp_comp <- ppp(xs, ys, marks = marks, win=owin())
#' plot(pp_comp)
#'
#' lcf_comp <- LCFcross(pp_comp, "0", "1")
#'
#' # Plot LCF for three point patterns together
#' plot(lcf_rand$r, lcf_rand$iso, xlab = "r", ylab = "LCF", type="l", ylim=c(-1, 1), col=4)
#' lines(lcf_rand$r, rep(0, nrow(lcf_rand)), lty=2)
#' lines(lcf_attr$r, lcf_attr$iso, col=2)
#' lines(lcf_comp$r, lcf_comp$iso, col=7)
#' legend("bottomright",
#'        c("theoretical", "random", "co-occurence", "separation"),
#'        col=c(1, 4, 2, 7),
#'        lty=c(2, 1, 1, 1))
#'
LCFcross <- function(pp,
                     i,
                     j,
                     correction="Ripley",
                     r=NULL,
                     dim=NULL,
                     dim_lims=NULL,
                     rmax=NULL) {

  if (!inherits(pp, "ppp") || !("marks" %in% attributes(pp)$names)) {
    rlang::abort(class = "lcf_cross_error_invalid_arg",
                 message="'pp' argument has to be a multitype point pattern")
  }

  mark_levels <- attr(pp$marks, "levels")

  if (is.null(mark_levels) || length(mark_levels) < 2) {
    rlang::abort(class = "lcf_cross_error_bad_marks",
                 message="'marks' should be a factor with at least two levels")
  }

  if (!is.character(correction) || length(correction) > 1) {
    rlang::abort(class = "lcf_cross_error_bad_correction",
    message="'correction' argument has to be a character vector
             with a single value")
  }

  if (correction == "all") {
    rlang::abort(class = "lcf_cross_error_inefficient_correction",
                 message="Won't use all edge correction methods due to performance reasons,
                          choose a single option")
  }

  if (!missing(i)) {
    if (!(i %in% mark_levels)) {
      rlang::abort(class = "lcf_cross_error_wrong_type",
                   message="'i' should be one of the factor levels
                   of marks assigned to a point pattern")
    }
  } else {
    i <- mark_levels[1]
  }

  if (!missing(j)) {
    if (!(j %in% mark_levels)) {
      rlang::abort(class = "lcf_cross_error_wrong_type",
                   message="'j' should be one of the factor levels
                   of marks assigned to a point pattern")
    }
  } else {
    j <- mark_levels[2]
  }

  r_arg <- get_r_arg(r)
  k_ij <- spatstat.explore::Kcross(pp,
                                   i, j,
                                   r=r_arg,
                                   correction = correction,
                                   rmax=rmax)

  n_j <- pp[pp$marks == j]$n
  intensity_j <- n_j / spatstat.geom::area(pp$window)
  pn <- intensity_j * k_ij[[3]]

  pn_est_df <- data.frame(r=k_ij$r,
                          pn=pn)

  if (is.null(dim)) {
    dim <- choose_basis_dim(pp$n, dim_lims=dim_lims)
  }

  # Name the column with LCF estimate after the used border correction method
  correction_name <- colnames(k_ij)[3]

  lcf_df <- LCF(pn_est_df, r, dim, correction_name)

  lcf_labl <- attr(k_ij, "labl")
  lcf_desc <- attr(k_ij, "desc")
  lcf_exp <- substitute("LCF"[list(i, j)](r), list(i=i, j=j))
  lcf_lab <- substitute("LCF"[i, j](r), list(i=i, j=j))

  lcf_fv <- spatstat.explore::fv(lcf_df, valu=correction_name,
                                 fname=c("LCF", paste0("list(", i, ",", j, ")")),
                                 fmla = ".~r", ylab=lcf_lab, yexp=lcf_exp,
                                 labl=lcf_labl, desc=lcf_desc)

  class(lcf_fv) <- c("lcffv", class(lcf_fv))
  lcf_fv
}

#' LCF estimate
#'
#' Estimates LCF by using smooth approximation of the empirical estimate of
#' the number of points, \eqn{N(r)}, and computing its derivative.
#'
#' @param pn_est_df The data frame with estimated number of points at distance \eqn{r}
#' that will be used to estimate LCF. Should have two columns
#' "\code{r}" (distance) and "\code{pn}" (estimated number of points)
#' @param r Optional. Optional. Vector of values for the argument \eqn{r} at which \eqn{LCF(r)}
#' should be evaluated. The values must be in increasing order. Advanced use only.
#' @param dim The dimension of the basis used to represent the smooth term within
#' the scam model formula. If not provided, the rule of thumb is used, i.e.
#' \eqn{dim = \sqrt{n}} where \eqn{n} is the number of points.
#' @param est_name The name of the column that contains the LCF's empirical estimate.
#'
#' @return A data frame with 3 columns:
#'
#' \code{r} --     the value of the argument \eqn{r} at which \eqn{LCF(r)} is estimated
#'
#' \code{theo} --   the theoretical value of \eqn{LCF_{ij}(r)} for Poisson process, \eqn{0}
#'
#' and a column named after the used edge correction method, e.g. "\code{iso}" or "\code{border}",
#' with the empirical estimate of \eqn{LCF_{ij}(r)} obtained from the point pattern.
#'
#' @noRd
LCF <- function(pn_est_df, r=NULL, dim, est_name) {
  # If r is not specified, return the result for all r in pn_est_df
  if (is.null(r)) {
    r <- pn_est_df$r
  }

  first_non_zero_ind <- which(pn_est_df$pn != 0)[1]
  # If there is no non-zero estimate of point number, return LCF = -1 everywhere
  if (is.na(first_non_zero_ind)) {
    lcf_df <- data.frame(r=r,
                         theo=0)
    lcf_df[est_name] <- -1
    return(lcf_df)
  }

  # Handling the situation when rmax is too large and Kest returns NAs
  last_ind <- min(nrow(pn_est_df), which(is.na(pn_est_df$pn))[1] - 1, na.rm = TRUE)
  pn_est_defined <- pn_est_df[first_non_zero_ind:last_ind, ]

  model <- scam::scam(pn ~ s(r, k=dim, bs='mpi'), data=pn_est_defined)

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

    # Workaround when getting a small negative derivative
    pn_deriv[pn_deriv < 0] <- 0

    lcf <- ifelse(pn > 0 & pn_deriv >= 0,
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
  lcf_df[est_name] <- lcf

  lcf_df
}


#' Derivative of a scam model with a single MPI term
#'
#' @param model A scam model with a single monotone increasing P-spline (MPI)
#' term.
#' @param data A data frame containing the values of the named covariates
#' at which the smooth term is to be evaluated.
#'
#' @return A vector that contains the derivative of the given scam model at
#' the value provided with the data parameter
#'
#' @noRd
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
#' @param dim_lims Optional. The integer vector of 2 values: the lower and
#' upper limits of the possible value of dim, \code{c(lower, upper)}. \code{lower} > \code{upper}
#' is not allowed.
#'
#' @return An integer value that represents the number of B-splines to use
#' in a scam model. The rule of thumb \enq{dim = \sqrt(sample\_size)} is used and
#' if \code{dim\_lims} is provided it is used to clip the value
#'
#' @noRd
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
#' @param r A vector of distances.
#' @param pn A vector of estimated expected number of points within a distance \eqn{r}.
#' @param pn_deriv A vector of estimated derivative of the expected number of points
#' within a distance \eqn{r}.
#' @param lcf_lims Optional. The double vector with 2 values: the lower and
#' upper limits of the LCF, \code{c(lower, upper)}. The lower value corresponds to the LCF value
#' for maximal dispersion, the upper to the LCF value for maximal clustering.
#' \code{lower} > \code{upper} is not allowed.
#'
#' @returns A vector with LCF estimate at r
#'
#' @noRd
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

get_r_arg <- function(r) {
  r_arg <- NULL

  if (!is.null(r) && length(r) > 513) {
    r_arg <- r
  }

  r_arg
}

#' Plot Function Value for LCF
#'
#' Plot method for the class "lcffv".
#'
#' Calls \code{\link[spatstat.explore]{plot.fv}} from the \code{spatstat.explore}
#' package and sets y-axis to the range of LCF, \eqn{[-1,1]}, to aid visual interpretation.
#' See \code{\link[spatstat.explore]{plot.fv}} for the information about the plotting parameters.
#'
#' @param x An object of the class "lcffv" that contains the variables to be
#' plotted.
#' @param ylim (optional) range of y axis. Default is set to the lower and
#' upper limits of the LCF.
#' @param main (optional) A title of the plot.
#' @param ... Extra arguments passed to the \code{\link[spatstat.explore]{plot.fv}}.
#'
#' @return Invisible: either NULL, or a data frame giving the meaning of the
#' different line types and colours.
#'
#' @author Evgenia Martynova \email{evg.martynova@@gmail.com}
#'
#' @seealso \code{\link{LCFest}}, \code{\link{LCFcross}}
#'
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

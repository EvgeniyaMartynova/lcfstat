
# n - number of datapoints
gen_poisson_data <- function(n) {
  r <- seq(0.01, 1, length.out=n)
  pn <- r^2
  pn_der <- 2*r

  data.frame(r=r, pn=pn, pn_der=pn_der)
}

# TODO: perhaps change to return the funtions instead
gen_sigmoid_data <- function(n) {
  r <- seq(-1, 1, length.out=n)
  pn <- exp(4 * r) / (1 + exp(4 * r))
  pn_der <- 4 * exp(4 * r) / (1 + exp(4 * r))^2
  lcf <- exp(-2 * log(2)  * r / (1 + exp(4 * r))) * 2 - 1

  data.frame(r=r, pn=pn, pn_der=pn_der, lcf=lcf)
}


gen_rpp_in_circle <- function(rad, point_num, centre=c(0, 0)) {
  circle <- spatstat.geom::disc(radius = rad, centre=centre)
  circle_a <- spatstat.geom::area(circle)
  intensity <- point_num / circle_a
  rpp <- spatstat.random::rpoispp(intensity, win=circle)

  rpp
}


gen_max_clustered_pattern <- function(circle_rad, point_num) {
  pp_loc <- gen_rpp_in_circle(circle_rad, point_num)
  wind <- spatstat.geom::owin(c(-0.5, 0.5), c(-0.5, 0.5))

  # Change window to be the target domain
  pp_mc <- spatstat.geom::ppp(pp_loc$x,
                              pp_loc$y,
                              window = wind)
  pp_mc
}


gen_three_clusters_pattern <- function(circle_rad, circle_point_num, cluster_dist) {

  column_distance <- sqrt(cluster_dist^2 - (cluster_dist / 2)^2)
  centers <- data.frame(x=c(-cluster_dist / 2, 0, cluster_dist / 2),
                        y=c(-column_distance/2, column_distance/2, -column_distance/2))

  xs <- NULL
  ys <- NULL
  for (i in 1:nrow(centers)) {
    centre <- c(centers$x[i], centers$y[i])
    pp <- gen_rpp_in_circle(circle_rad, circle_point_num, centre)
    xs <- c(xs, pp[["x"]])
    ys <- c(ys, pp[["y"]])
  }

  wind <- spatstat.geom::owin(c(-0.5, 0.5), c(-0.5, 0.5))
  pp_three <- spatstat.geom::ppp(xs, ys, window = wind)

  pp_three
}


# n - number of datapoints
gen_poisson_data <- function(n) {
  r <- seq(0.01, 1, length.out=n)
  pn <- r^2
  pn_der <- 2*r

  data.frame(r=r, pn=pn, pn_der=pn_der)
}

# TODO: perhaps change to return the funtion instead
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


two_type_pp_random <- function(intensity=500) {

  pp <- spatstat.random::rpoispp(intensity)
  pp$marks <- factor(sample(0:1, spatstat.geom::npoints(pp), replace=TRUE))

  pp
}

two_type_attraction <- function(intensity=500) {

  # Generate a random point pattern with half of intensity
  pn <- intensity / 2
  pp <- spatstat.random::rpoispp(pn)

  x_displ <- 5e-05
  y_displ <- 5e-05

  # Add a very small displacement to original coordinates and combine them
  xs <- c(pp$x, pp$x + x_displ)
  ys <- c(pp$y, pp$y + y_displ)
  marks <- factor(c(rep("0", pp$n), rep("1", pp$n)))

  # Make a mutlitype point pattern
  pp_mult <- spatstat.geom::ppp(xs, ys, marks = marks)

  pp_mult
}

two_type_repulsion <- function(intensity=2000) {

  main_win <- spatstat.geom::owin()

  # Generate point random point pattern in the opposite halves of the windo
  # with a gap of 0.02
  pn <- intensity / 2

  win1 <- spatstat.geom::owin(c(0, 0.49), c(0, 1))
  pp1 <- spatstat.random::rpoispp(pn, win=win1)

  win2 <- spatstat.geom::owin(c(0.51, 1), c(0, 1))
  pp2 <- spatstat.random::rpoispp(pn, win=win2)

  # Combine coordinates and make marks
  xs <- c(pp1$x, pp2$x)
  ys <- c(pp1$y, pp2$y)
  marks <- factor(c(rep("0", pp1$n), rep("1", pp2$n)))

  # Make a mutlitype point pattern in the main window
  pp_mult <- spatstat.geom::ppp(xs, ys, marks = marks, win=main_win)
  pp_mult
}

multitype_pp_random <- function(intensity=500, type_num=2) {

  pp <- spatstat.random::rpoispp(intensity)
  pp$marks <- factor(sample(0:(type_num-1), spatstat.geom::npoints(pp), replace=TRUE))

  pp
}

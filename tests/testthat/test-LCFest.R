
# choose_basis_dim() tests
test_that("choose_basis_dim() returns rounded sqrt of the sample size as dimension", {
  expect_equal(choose_basis_dim(100), 10)
  expect_equal(choose_basis_dim(500), 22)
  expect_equal(choose_basis_dim(15), 4)
})


test_that("choose_basis_dim() clips output to provided limits", {
  expect_equal(choose_basis_dim(100, dim_lims=c(15L, 50L)), 15)
  expect_equal(choose_basis_dim(100, dim_lims=c(5, 9)), 9)

  expect_equal(choose_basis_dim(101, dim_lims=c(10, 15)), 10)
  expect_equal(choose_basis_dim(99, dim_lims=c(5L, 10L)), 10)
})


test_that("choose_basis_dim() handles illegal input sample size", {
  expect_error(choose_basis_dim(-100), class = "lcf_dim_error_invalid_sample_size")
  expect_error(choose_basis_dim(c(30, 400)), class = "lcf_dim_error_invalid_sample_size")
  expect_error(choose_basis_dim(15.4), class = "lcf_dim_error_invalid_sample_size")
})


test_that("choose_basis_dim() handles illegal input dimension limits", {
  dim <- 49
  expect_error(choose_basis_dim(dim, c(-1, 10)), class = "lcf_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(1, -10)), class = "lcf_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(-1, -10)), class = "lcf_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(9, 5)), class = "lcf_dim_error_unordered_dim_lims")
  expect_error(choose_basis_dim(dim, c(5, 9, 13)), class = "lcf_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(5.2, 9.3)), class = "lcf_dim_error_invalid_dim_lims")
})


# compute_lcf() tests
test_that("compute_lcf() returns expected values for edge cases", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  # Theoretical Poisson curve
  expect_equal(compute_lcf(r, pn, pn_der), rep(0, n))
  # Zero derivative should give maximum LCF value, default: 1
  expect_equal(compute_lcf(r, pn, rep(0, n)), rep(1, n))
  # Point number value close to zero should give minimum LCF value, default: -1
  expect_equal(compute_lcf(r, rep(1e-13, n), pn_der), rep(-1, n))
})


test_that("compute_lcf() returns expected values for a custom sigmoid function", {
  n <- 100
  sigmoid_data <- gen_sigmoid_data(n)
  r <- sigmoid_data$r
  pn <- sigmoid_data$pn
  pn_der <- sigmoid_data$pn_der
  expected_lcf <- sigmoid_data$lcf

  expect_equal(compute_lcf(r, pn, pn_der), expected_lcf)
})


test_that("compute_lcf() applies scaling according to the lcf_lims parameter", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  lcf_lims <- c(0, 1)
  expect_equal(compute_lcf(r, pn, pn_der, lcf_lims), rep(0.5, n))
  expect_equal(compute_lcf(r, pn, rep(0, n), lcf_lims), rep(1, n))
  expect_equal(compute_lcf(r, rep(1e-13, n), pn_der, lcf_lims), rep(0, n))

  lcf_lims <- c(-10, 0)
  expect_equal(compute_lcf(r, pn, pn_der, lcf_lims), rep(-5, n))
  expect_equal(compute_lcf(r, pn, rep(0, n), lcf_lims), rep(0, n))
  expect_equal(compute_lcf(r, rep(1e-13, n), pn_der, lcf_lims), rep(-10, n))
})


test_that("compute_lcf() fails when type of main arguments is illegal", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  expect_error(compute_lcf(rep("a", n), pn, pn_der), class = "lcf_compute_error_invalid_arg")
  expect_error(compute_lcf(r, rep("a", n), pn_der), class = "lcf_compute_error_invalid_arg")
  expect_error(compute_lcf(r, pn, rep("a", n)), class = "lcf_compute_error_invalid_arg")
})


test_that("compute_lcf() fails when length of main arguments is different", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  expect_error(compute_lcf(r[2:n], pn, pn_der), class = "lcf_compute_error_arg_length_mismatch")
  expect_error(compute_lcf(r, pn[2:n], pn_der), class = "lcf_compute_error_arg_length_mismatch")
  expect_error(compute_lcf(r, pn, pn_der[2:n]), class = "lcf_compute_error_arg_length_mismatch")
})


test_that("compute_lcf() fails when illegal lcf_lims argument is provided", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  expect_error(compute_lcf(r, pn, pn_der, c(-1, 1, 2)), class = "lcf_compute_error_invalid_lcf_lims")
  expect_error(compute_lcf(r, pn, pn_der, c("a", "b")), class = "lcf_compute_error_invalid_lcf_lims")
  expect_error(compute_lcf(r, pn, pn_der, c(1, -1)), class = "lcf_compute_error_unordered_lcf_lims")
})


# single_mpi_derivative() tests
test_that("single_mpi_derivative() returns output in expected format", {
  n <- 100
  sigmoid_data <- gen_sigmoid_data(n)
  r <- sigmoid_data$r
  pn <- sigmoid_data$pn

  pn_noise <- pn + rnorm(n) * 0.1

  model <- scam::scam(pn_noise ~ s(r, bs="mpi"))
  r_new <- seq(-1, 1, length.out=2 * n)
  derivative <- single_mpi_derivative(model, data=list(r=r_new))

  expect_type(derivative, "double")
  expect_null(attributes(derivative))
  expect_length(derivative, 2 * n)
})


test_that("single_mpi_derivative() computes derivative close to the true derivative", {
  n <- 1000
  sigmoid_data <- gen_sigmoid_data(n)
  r <- sigmoid_data$r
  pn <- sigmoid_data$pn

  local_seed <- 1511
  withr::with_seed(local_seed, {
    pn_noise <- pn + rnorm(n) * 0.05
  })

  model <- scam::scam(pn_noise ~ s(r, k=20, bs="mpi"))

  r_new <- seq(-1, 1, length.out=2 * n)
  deriv_analytical <- 4 * exp(4 * r_new) / (1 + exp(4 * r_new))^2
  deriv_approx <- single_mpi_derivative(model, data=list(r=r_new))

  max_dev <- max(abs(deriv_analytical - deriv_approx))
  expect_lt(max_dev, 0.1)
})


test_that("single_mpi_derivative() throws error when upsupported model is provided", {
  n <- 100
  sigmoid_data <- gen_sigmoid_data(n)
  x <- sigmoid_data$r
  f <- sigmoid_data$pn
  y <- f + rnorm(n) * 0.1

  model <- mgcv::gam(y ~ s(x))
  expect_error(single_mpi_derivative(model, data=list(x=x)), class = "lcf_derivative_error_unsupported_model")

  model <- scam::scam(y ~ s(x))
  expect_error(single_mpi_derivative(model, data=list(x=x)), class = "lcf_derivative_error_unsupported_model_formula")

  x1 <- sort(runif(n)*3-1)
  f1 <- exp(-1.3 * x1)
  y <- f + f1 + + rnorm(n) * 0.2
  model <- scam::scam(y ~ s(x, k=20, bs="mpi") + s(x1, k=15, bs="mpd"))
  expect_error(single_mpi_derivative(model, data=list(x=x, x1=x1)), class = "lcf_derivative_error_unsupported_model_formula")

  model <- scam::scam(y ~ s(x, bs="mpi") + x1)
  expect_error(single_mpi_derivative(model, data=list(x=x, x1=x1)), class = "lcf_derivative_error_unsupported_model_formula")
})


test_that("single_mpi_derivative() throws error when illegal data is provided", {
  n <- 100
  sigmoid_data <- gen_sigmoid_data(n)
  x <- sigmoid_data$r
  f <- sigmoid_data$pn
  y <- f + rnorm(n) * 0.1

  model <- scam::scam(y ~ s(x, bs="mpi"))
  x_out <- seq(-2, 0, length.out= 2 * n)
  expect_error(single_mpi_derivative(model, data=list(x=x_out)), class = "lcf_derivative_error_invalid_argument")

  model <- scam::scam(y ~ s(x, bs="mpi"))
  expect_error(single_mpi_derivative(model, data=list(r=x)), class = "lcf_derivative_error_missing_argument")
})


# LCFest() tests
test_that("LCFest() returns an lcffv object with expected attribute values", {
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  lcf <- LCFest(rpp)
  expect_equal(ncol(lcf), 3)
  expect_named(lcf, c("r", "theo", "iso"))
  expect_s3_class(lcf, "lcffv")
  expect_s3_class(lcf, "fv")
  expect_equal(attr(lcf, "fname"), "LCF")
  expect_equal(attr(lcf, "fmla"), ".~r")
  expect_equal(attr(lcf, "ylab"), quote(LCF(r)))
  expect_equal(attr(lcf, "yexp"), quote(LCF(r)))
  expect_gte(min(lcf[[3]]), -1)
  expect_lte(max(lcf[[3]]), 1)
})


test_that("LCFest() applies provided correction attribute", {
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  lcf <- LCFest(rpp, "border")
  expect_named(lcf, c("r", "theo", "border"))

  lcf <- LCFest(rpp, "translate")
  expect_named(lcf, c("r", "theo", "trans"))
})


test_that("LCFest() applies provided rmax attribute", {
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  lcf <- LCFest(rpp, rmax=0.1)
  expect_gte(min(lcf$r), 0)
  expect_lte(max(lcf$r), 0.1)

  lcf <- LCFest(rpp, rmax=0.5)
  expect_gte(min(lcf$r), 0)
  expect_lte(max(lcf$r), 0.5)

  # Passed rmax is too large
  lcf <- LCFest(rpp, rmax=1)
  first_na_ind <- which(is.na(lcf$iso))[1]
  expect_equal(sum(is.na(lcf$iso)), nrow(lcf) - first_na_ind + 1)
})


test_that("LCFest() applies provided r_arg attribute", {
  seed <- 1511
  point_num <- 100
  withr::with_seed(seed, {
    rpp <- spatstat.random::rpoispp(point_num)
  })

  r_arg <- seq(0.1, 0.2, by=0.01)
  lcf <- LCFest(rpp, r = r_arg)
  expect_equal(nrow(lcf), length(r_arg))

  r_arg <- c(0.05, 0.1, 0.2)
  lcf <- LCFest(rpp, r = r_arg)
  expect_equal(nrow(lcf), 3)

  r_arg <- c(0.1)
  lcf <- LCFest(rpp, r = r_arg)
  expect_equal(nrow(lcf), 1)

  r_arg <- c(0.001)
  lcf <- LCFest(rpp, r = r_arg)
  expect_equal(nrow(lcf), 1)

  r_arg <- c(0.699)
  lcf <- LCFest(rpp, r = r_arg)
  expect_equal(nrow(lcf), 1)

  # All r in LCF == -1 region
  r_arg <- c(0.001, 0.0011, 0.0012)
  lcf <- LCFest(rpp, r = r_arg)
  expect_equal(lcf$iso, c(-1, -1, -1))

  # All r are too large
  r_arg <- c(0.8, 0.85, 0.9)
  lcf <- LCFest(rpp, r = r_arg)
  expect_true(all(is.na(lcf$iso)))
})


test_that("LCFest() value for clustered pattern is greater than 0 at domain size", {
  parent_num <- 20
  clust_pn <- 25
  rad <- 0.05
  mat_clust_pp <- spatstat.random::rMatClust(parent_num,
                                             rad,
                                             clust_pn)

  lcf <- LCFest(mat_clust_pp)

  domain_size <- rad * 2
  domain_size_ind <- which.min(abs(lcf$r - domain_size))
  expect_gt(lcf[[3]][domain_size_ind], 0)
})


test_that("LCFest() value for the Hardcore pattern equals to -1 until the inhibition distance", {
  id <- 0.05
  point_num <- 300
  hardcore_pp <- spatstat.random::rHardcore(point_num, R=id)

  lcf <- LCFest(hardcore_pp)
  id_ind <- which.min(abs(lcf$r - id))
  expect_equal(max(lcf[[3]][1:id_ind]), -1)
})


test_that("LCFest() value reaches 1 and plateaues for maximum clustering pattern", {
  rad <- 0.05
  point_num <- 500

  seed <- 1511
  withr::with_seed(seed, {
    pp_mc <- gen_max_clustered_pattern(rad, point_num)
  })

  lcf <- LCFest(pp_mc)

  domain_size <- rad * 2
  domain_size_ind <- which.min(abs(lcf$r - domain_size))
  expect_gt(min(lcf[[3]][domain_size_ind:length(lcf$r)]), 0.975)
  expect_equal(max(lcf[[3]][domain_size_ind:length(lcf$r)]), 1)
})


test_that("LCFest() reaches expected values for the pattern with three clusters", {
  rad <- 0.03
  cluster_dist <- 0.15
  point_num <- 200

  seed <- 1510
  withr::with_seed(seed, {
    pp_three <- gen_three_clusters_pattern(rad, point_num, cluster_dist)
  })

  lcf <- LCFest(pp_three)

  # Find indices of the r vector that corresponds to distances related to
  # spatial arrangement of the point pattern
  domain_size <- rad * 2
  domain_size_ind <- which.min(abs(lcf$r - domain_size))

  min_point_dist <- cluster_dist - domain_size
  min_point_dist_ind <- which.min(abs(lcf$r - min_point_dist))

  cluster_dist_ind <- which.min(abs(lcf$r - cluster_dist))

  max_point_dist <- cluster_dist + domain_size
  max_point_dist_ind <- which.min(abs(lcf$r - max_point_dist))

  # TODO: perhaps check that function grows and decreases at the
  # corresponding intervals
  expect_gt(min(lcf[[3]][domain_size_ind:min_point_dist_ind]), 0.95)
  expect_equal(max(lcf[[3]][domain_size_ind:min_point_dist_ind]), 1)

  expect_lt(lcf[[3]][cluster_dist_ind], 0)

  expect_gt(min(lcf[[3]][max_point_dist_ind:length(lcf$r)]), 0.95)
  expect_equal(max(lcf[[3]][max_point_dist_ind:length(lcf$r)]), 1)
})


test_that("LCFest() throws error when the first argument is not of spatstat.geom::ppp type", {
  # First object has other type than spatstat.geom::ppp
  n <- 100
  df <- data.frame(x = runif(n), y = runif(n))
  expect_error(LCFest(df))
})


test_that("LCFest() throws error when illegal correctiion argument is provided", {
  # First object has other type than spatstat.geom::ppp
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  # Illegal correction parameter
  expect_error(LCFest(rpp, 1), class = "lcf_error_bad_correction")
  expect_error(LCFest(rpp, c("Ripley", "border")), class = "lcf_error_bad_correction")
  expect_error(LCFest(rpp, c("all")), class = "lcf_error_inefficient_correction")
  expect_error(LCFest(rpp, c("qwerty")), "unrecognised correction")
})



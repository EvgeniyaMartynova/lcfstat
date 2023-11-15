
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
  expect_error(choose_basis_dim(-100), class = "ai_dim_error_invalid_sample_size")
  expect_error(choose_basis_dim(c(30, 400)), class = "ai_dim_error_invalid_sample_size")
  expect_error(choose_basis_dim(15.4), class = "ai_dim_error_invalid_sample_size")
})


test_that("choose_basis_dim() handles illegal input dimension limits", {
  dim <- 49
  expect_error(choose_basis_dim(dim, c(-1, 10)), class = "ai_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(1, -10)), class = "ai_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(-1, -10)), class = "ai_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(9, 5)), class = "ai_dim_error_unordered_dim_lims")
  expect_error(choose_basis_dim(dim, c(5, 9, 13)), class = "ai_dim_error_invalid_dim_lims")
  expect_error(choose_basis_dim(dim, c(5.2, 9.3)), class = "ai_dim_error_invalid_dim_lims")
})


# compute_ai() tests
test_that("compute_ai() returns expected values for edge cases", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  # Theoretical Poisson curve
  expect_equal(compute_ai(r, pn, pn_der), rep(0, n))
  # Zero derivative should give maximum AI value, default: 1
  expect_equal(compute_ai(r, pn, rep(0, n)), rep(1, n))
  # Point number value close to zero should give minimum AI value, default: -1
  expect_equal(compute_ai(r, rep(1e-13, n), pn_der), rep(-1, n))
})


test_that("compute_ai() returns expected values for a custom sigmoid function", {
  n <- 100
  sigmoid_data <- gen_sigmoid_data(n)
  r <- sigmoid_data$r
  pn <- sigmoid_data$pn
  pn_der <- sigmoid_data$pn_der
  expected_ai <- sigmoid_data$ai

  expect_equal(compute_ai(r, pn, pn_der), expected_ai)
})


test_that("compute_ai() applies scaling according to the ai_lims parameter", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  ai_lims <- c(0, 1)
  expect_equal(compute_ai(r, pn, pn_der, ai_lims), rep(0.5, n))
  expect_equal(compute_ai(r, pn, rep(0, n), ai_lims), rep(1, n))
  expect_equal(compute_ai(r, rep(1e-13, n), pn_der, ai_lims), rep(0, n))

  ai_lims <- c(-10, 0)
  expect_equal(compute_ai(r, pn, pn_der, ai_lims), rep(-5, n))
  expect_equal(compute_ai(r, pn, rep(0, n), ai_lims), rep(0, n))
  expect_equal(compute_ai(r, rep(1e-13, n), pn_der, ai_lims), rep(-10, n))
})


test_that("compute_ai() fails when type of main arguments is illegal", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  expect_error(compute_ai(rep("a", n), pn, pn_der), class = "ai_compute_error_invalid_arg")
  expect_error(compute_ai(r, rep("a", n), pn_der), class = "ai_compute_error_invalid_arg")
  expect_error(compute_ai(r, pn, rep("a", n)), class = "ai_compute_error_invalid_arg")
})


test_that("compute_ai() fails when length of main arguments is different", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  expect_error(compute_ai(r[2:n], pn, pn_der), class = "ai_compute_error_arg_length_mismatch")
  expect_error(compute_ai(r, pn[2:n], pn_der), class = "ai_compute_error_arg_length_mismatch")
  expect_error(compute_ai(r, pn, pn_der[2:n]), class = "ai_compute_error_arg_length_mismatch")
})


test_that("compute_ai() fails when illegal ai_lims argument is provided", {
  n <- 100
  poisson_data <- gen_poisson_data(n)
  r <- poisson_data$r
  pn <- poisson_data$pn
  pn_der <- poisson_data$pn_der

  expect_error(compute_ai(r, pn, pn_der, c(-1, 1, 2)), class = "ai_compute_error_invalid_ai_lims")
  expect_error(compute_ai(r, pn, pn_der, c("a", "b")), class = "ai_compute_error_invalid_ai_lims")
  expect_error(compute_ai(r, pn, pn_der, c(1, -1)), class = "ai_compute_error_unordered_ai_lims")
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
  expect_error(single_mpi_derivative(model, data=list(x=x)), class = "ai_derivative_error_unsupported_model")

  model <- scam::scam(y ~ s(x))
  expect_error(single_mpi_derivative(model, data=list(x=x)), class = "ai_derivative_error_unsupported_model_formula")

  x1 <- sort(runif(n)*3-1)
  f1 <- exp(-1.3 * x1)
  y <- f + f1 + + rnorm(n) * 0.2
  model <- scam::scam(y ~ s(x, k=20, bs="mpi") + s(x1, k=15, bs="mpd"))
  expect_error(single_mpi_derivative(model, data=list(x=x, x1=x1)), class = "ai_derivative_error_unsupported_model_formula")

  model <- scam::scam(y ~ s(x, bs="mpi") + x1)
  expect_error(single_mpi_derivative(model, data=list(x=x, x1=x1)), class = "ai_derivative_error_unsupported_model_formula")
})


test_that("single_mpi_derivative() throws error when illegal data is provided", {
  n <- 100
  sigmoid_data <- gen_sigmoid_data(n)
  x <- sigmoid_data$r
  f <- sigmoid_data$pn
  y <- f + rnorm(n) * 0.1

  model <- scam::scam(y ~ s(x, bs="mpi"))
  x_out <- seq(-2, 0, length.out= 2 * n)
  expect_error(single_mpi_derivative(model, data=list(x=x_out)), class = "ai_derivative_error_invalid_argument")

  model <- scam::scam(y ~ s(x, bs="mpi"))
  expect_error(single_mpi_derivative(model, data=list(r=x)), class = "ai_derivative_error_missing_argument")
})


# AIest() tests
test_that("AIest() returns an aifv object with expected attribute values", {
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  ai <- AIest(rpp)
  expect_equal(ncol(ai), 3)
  expect_named(ai, c("r", "theo", "iso"))
  expect_s3_class(ai, "aifv")
  expect_s3_class(ai, "fv")
  expect_equal(attr(ai, "fname"), "AI")
  expect_equal(attr(ai, "fmla"), ".~r")
  expect_equal(attr(ai, "ylab"), quote(AI(r)))
  expect_equal(attr(ai, "yexp"), quote(AI(r)))
  expect_gte(min(ai[[3]]), -1)
  expect_lte(max(ai[[3]]), 1)
})


test_that("AIest() applies provided correction attribute", {
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  ai <- AIest(rpp, "border")
  expect_named(ai, c("r", "theo", "border"))

  ai <- AIest(rpp, "translate")
  expect_named(ai, c("r", "theo", "trans"))
})


test_that("AIest() applies provided rmax attribute", {
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  ai <- AIest(rpp, rmax=0.1)
  expect_gte(min(ai$r), 0)
  expect_lte(max(ai$r), 0.1)

  ai <- AIest(rpp, rmax=0.5)
  expect_gte(min(ai$r), 0)
  expect_lte(max(ai$r), 0.5)

  # Passed rmax is too large
  ai <- AIest(rpp, rmax=1)
  first_na_ind <- which(is.na(ai$iso))[1]
  expect_equal(sum(is.na(ai$iso)), nrow(ai) - first_na_ind + 1)
})


test_that("AIest() value for clustered pattern is greater than 0 at domain size", {
  parent_num <- 20
  clust_pn <- 25
  rad <- 0.05
  mat_clust_pp <- spatstat.random::rMatClust(parent_num,
                                             rad,
                                             clust_pn)

  ai <- AIest(mat_clust_pp)

  domain_size <- rad * 2
  domain_size_ind <- which.min(abs(ai$r - domain_size))
  expect_gt(ai[[3]][domain_size_ind], 0)
})


test_that("AIest() value for the Hardcore pattern equals to -1 until the inhibition distance", {
  id <- 0.05
  point_num <- 300
  hardcore_pp <- spatstat.random::rHardcore(point_num, R=id)

  ai <- AIest(hardcore_pp)
  id_ind <- which.min(abs(ai$r - id))
  expect_equal(max(ai[[3]][1:id_ind]), -1)
})


test_that("AIest() value reaches 1 and plateaues for maximum clustering pattern", {
  rad <- 0.05
  point_num <- 500

  seed <- 1511
  withr::with_seed(seed, {
    pp_mc <- gen_max_clustered_pattern(rad, point_num)
  })

  ai <- AIest(pp_mc)

  domain_size <- rad * 2
  domain_size_ind <- which.min(abs(ai$r - domain_size))
  expect_gt(min(ai[[3]][domain_size_ind:length(ai$r)]), 0.975)
  expect_equal(max(ai[[3]][domain_size_ind:length(ai$r)]), 1)
})


test_that("AIest() reaches expected values for the pattern with three clusters", {
  rad <- 0.03
  cluster_dist <- 0.15
  point_num <- 200

  seed <- 1510
  withr::with_seed(seed, {
    pp_three <- gen_three_clusters_pattern(rad, point_num, cluster_dist)
  })

  ai <- AIest(pp_three)

  # Find indices of the r vector that corresponds to distances related to
  # spatial arrangement of the point pattern
  domain_size <- rad * 2
  domain_size_ind <- which.min(abs(ai$r - domain_size))

  min_point_dist <- cluster_dist - domain_size
  min_point_dist_ind <- which.min(abs(ai$r - min_point_dist))

  cluster_dist_ind <- which.min(abs(ai$r - cluster_dist))

  max_point_dist <- cluster_dist + domain_size
  max_point_dist_ind <- which.min(abs(ai$r - max_point_dist))

  # TODO: perhaps check that function grows and decreases at the
  # corresponding intervals
  expect_gt(min(ai[[3]][domain_size_ind:min_point_dist_ind]), 0.95)
  expect_equal(max(ai[[3]][domain_size_ind:min_point_dist_ind]), 1)

  expect_lt(ai[[3]][cluster_dist_ind], 0)

  expect_gt(min(ai[[3]][max_point_dist_ind:length(ai$r)]), 0.95)
  expect_equal(max(ai[[3]][max_point_dist_ind:length(ai$r)]), 1)
})


test_that("AIest() throws error when the first argument is not of spatstat.geom::ppp type", {
  # First object has other type than spatstat.geom::ppp
  n <- 100
  df <- data.frame(x = runif(n), y = runif(n))
  expect_error(AIest(df))
})


test_that("AIest() throws error when illegal correctiion argument is provided", {
  # First object has other type than spatstat.geom::ppp
  point_num <- 500
  rpp <- spatstat.random::rpoispp(point_num)

  # Illegal correction parameter
  expect_error(AIest(rpp, 1), class = "ai_error_bad_correction")
  expect_error(AIest(rpp, c("Ripley", "border")), class = "ai_error_bad_correction")
  expect_error(AIest(rpp, c("all")), class = "ai_error_inefficient_correction")
  expect_error(AIest(rpp, c("qwerty")), "unrecognised correction")
})



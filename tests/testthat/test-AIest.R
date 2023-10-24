test_that("choose_basis_dim() returns sqrt of the sample size as dim", {
  expect_equal(choose_basis_dim(100), 10)
})

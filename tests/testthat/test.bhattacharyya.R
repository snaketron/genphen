context("Test Bhattacharyya")

test_that("getBhattacharyya", {
  expect_is(class = "numeric",
            object = getBhattacharyya(x = rnorm(10^4, -5, 1),
                                      y = rnorm(10^4, 5, 1))$bc)
  expect_equal(object = length(getBhattacharyya(x = rnorm(10^4, -5, 1),
                                                y = rnorm(10^4, 5, 1))$bc), 
               expected = 1)
  expect_equal(object = getBhattacharyya(x = rnorm(10^4, -5, 1),
                                         y = rnorm(10^4, 5, 1))$bc, 
               expected = 0, 0+0.01)
  expect_is(class = "numeric",
            object = getBhattacharyya(x = rnorm(10^4, 0, 2),
                                      y = rnorm(10^4, 0, 2))$bc)
  expect_equal(object = length(getBhattacharyya(x = rnorm(10^4, 0, 2),
                                                y = rnorm(10^4, 0, 2))$bc), 
               expected = 1)
  expect_equal(object = getBhattacharyya(x = rnorm(10^4, 0, 2),
                                         y = rnorm(10^4, 0, 2))$bc, 
               expected = 1, 1+0.01)
})
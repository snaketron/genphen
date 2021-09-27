context("Test HDI")


test_that("getHdi output is numeric vector of length 2", {
  expect_is(class = "numeric",
            object = getHdi(vec = runif(n = 10^3, min = 0, max = 1), 
                            hdi.level = 0.99))
  expect_equal(object = length(getHdi(vec = runif(n = 10^3, min = 0, max = 1), 
                                      hdi.level = 0.99)), 
               expected = 2)
  expect_equal(object = getHdi(vec = rnorm(n = 10^6, mean = 0, sd = 1), 
                               hdi.level = 0.95)[1], 
               expected = -1.96, tolerance = 0.1)
  expect_equal(object = getHdi(vec = rnorm(n = 10^6, mean = 0, sd = 1), 
                               hdi.level = 0.95)[2], 
               expected = 1.96, tolerance = 0.1)
})


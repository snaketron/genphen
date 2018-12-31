context("Test Kappa")




# Case 1:
# Real : 3A, 2B
# Pred: 4A, 1B
# CA_obs = 4/5 = 0.8
# CA_exp = (3*4+2*1)/5^2 = 0.56
# kappa = (CA_obs - CA_exp)/(1-CA_exp) = (0.8 - 0.56)/(1 - 0.56) = 0.5454545
# --------------
# Case 2:
# Real : 5A, 0B
# Pred: 5A, 0B
# CA_obs = 5/5 = 1
# CA_exp = (5*5+0*0)/5^2 = 1.0
# kappa = (CA_obs - CA_exp)/(1-CA_exp) = (1.0 - 1.0)/(1.0 - 1.0) = 0

test_that("getKappa", {
  expect_is(class = "numeric",
            object = getKappa(predicted = c("A", "A", "A", "A", "B"), 
                              real = c("A", "A", "A", "B", "B"), 
                              aas = c("A", "B")))
  expect_is(class = "numeric",
            object = getKappa(predicted = c("A", "A", "A", "A", "A"), 
                              real = c("A", "A", "A", "A", "A"), 
                              aas = c("A")))
  expect_equal(object = length(getKappa(predicted = c("A", "A", "A", "A", "B"), 
                                        real = c("A", "A", "A", "B", "B"), 
                                        aas = c("A", "B"))), 
               expected = 1)
  expect_equal(object = length(getKappa(predicted = c("A", "A", "A", "A", "A"), 
                                        real = c("A", "A", "A", "A", "A"), 
                                        aas = c("A"))), 
               expected = 1)
  expect_equal(object = getKappa(predicted = c("A", "A", "A", "A", "B"), 
                                 real = c("A", "A", "A", "B", "B"), 
                                 aas = c("A", "B")), 
               expected = 0.5454545, 0.5454545+10^-6)
  expect_equal(object = getKappa(predicted = c("A", "A", "A", "A", "A"), 
                                 real = c("A", "A", "A", "A", "A"), 
                                 aas = c("A")), 
               expected = 0, 0+10^-6)
})

context("Test model compilation")

# 1
testthat::expect_is(class = "stanmodel", 
                    object = compileModel(phenotype.type = "dichotomous",
                                          model.type = "univariate"))

# 2
testthat::expect_is(class = "stanmodel", 
                    object = compileModel(phenotype.type = "dichotomous",
                                          model.type = "hierarchical"))

# 3
testthat::expect_is(class = "stanmodel", 
                    object = compileModel(phenotype.type = "continuous",
                                          model.type = "univariate"))

# 4
testthat::expect_is(class = "stanmodel", 
                    object = compileModel(phenotype.type = "continuous",
                                          model.type = "hierarchical"))



p <- assignmentParameters('FIE')

test_that("technique can be returned", {
  expect_type(technique(p),'character')
})

test_that("limit can be returned", {
  expect_type(limit(p),'double')
})

test_that("limit can be set", {
  new_limit <- 0.002
  limit(p) <- new_limit
  
  expect_identical(limit(p),new_limit)
})

test_that("adducts can be returned", {
  expect_type(adducts(p),'list')
})

test_that("adducts can be set", {
  new_adducts <- list(n = c("[M-H]1-","[M+Cl]1-"),
                      p = c("[M+H]1+","[M+K]1+"))
  adducts(p) <- new_adducts
  
  expect_identical(adducts(p),new_adducts)
})

test_that("isotopes can be returned", {
  expect_type(isotopes(p),'character')
})

test_that("isotopes can be set", {
  new_isotopes <- '13C'
  isotopes(p) <- new_isotopes

  expect_identical(isotopes(p),new_isotopes)
})

test_that("transformations can be returned", {
  expect_type(transformations(p),'character')
})

test_that("transformations can be set", {
  new_transformations <- "M - [O] + [NH2]"
  transformations(p) <- new_transformations
  
  expect_identical(transformations(p),new_transformations)
})

test_that("adduct rules can be returned", {
  expect_s3_class(adductRules(p),'tbl_df')
})

test_that("adducts rules can be set", {
  new_adduct_rules <- mzAnnotation::adduct_rules()[1,]
  adductRules(p) <- new_adduct_rules
  
  expect_identical(adductRules(p),new_adduct_rules)
})

test_that("isotope rules can be returned", {
  expect_s3_class(isotopeRules(p),'tbl_df')
})

test_that("isotope rules can be set", {
  new_isotope_rules <- mzAnnotation::isotope_rules()[1,]
  isotopeRules(p) <- new_isotope_rules
  
  expect_identical(isotopeRules(p),new_isotope_rules)
})

test_that("transformation rules can be returned", {
  expect_s3_class(transformationRules(p),'tbl_df')
})

test_that("transformation rules can be set", {
  new_transformation_rules <- mzAnnotation::transformation_rules()[1,]
  transformationRules(p) <- new_transformation_rules
  
  expect_identical(transformationRules(p),new_transformation_rules)
})
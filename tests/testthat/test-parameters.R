
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
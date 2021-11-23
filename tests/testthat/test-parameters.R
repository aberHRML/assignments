
p <- assignmentParameters('FIE')

test_that("isotopes can be returned", {
  expect_type(isotopes(p),'character')
})

test_that("isotopes can be set", {
  new_isotopes <- '13C'
  isotopes(p) <- new_isotopes

  expect_identical(isotopes(p),new_isotopes)
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

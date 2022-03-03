test_that("Digits correctly set upon package load", {
  options('digits' = 7)
  
  MFassign:::.onLoad()
  
  expect_equal(getOption('digits'),10)
})

test_that('Digits not set upon package load if already above 10', {
  options('digits' = 11)
  
  MFassign:::.onLoad()
  
  expect_equal(getOption('digits'),11)
})

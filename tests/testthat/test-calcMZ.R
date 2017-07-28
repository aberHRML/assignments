
context('calcMZ')

test_that('calcMZ returns correctly',{
  res <- MFassign:::calcM(134.0215,'C13','[M-H]1-')
  expect_true(class(res) == 'numeric')
})

context('calcM')

test_that('calcM returns correctly',{
  res <- MFassign:::calcM(134.0176,'C13','[M-H]1-')
  expect_true(class(res) == 'numeric')
})
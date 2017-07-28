
context('MFscore')

test_that('MFscore returns correctly',{
  res <- MFassign:::MFscore('C4H6O5PS')
  expect_true(class(res) == 'numeric')
})

context('MFgen')

test_that('MFgen returns correctly',{
  res <- MFassign:::MFgen(117.07898,118.08626)
  
  expect_false(F %in% (class(res) == c('tbl_df','tbl','data.frame')))
  expect_false(F %in% (colnames(res) == c('MF','Theoretical M','PPM Error','Measured M','Measured m/z')))
  expect_true(class(res$MF) == 'character')
  expect_true(class(res$`Theoretical M`) == 'numeric')
  expect_true(class(res$`PPM Error`) == 'numeric')
  expect_true(class(res$`Measured M`) == 'numeric')
  expect_true(class(res$`Measured m/z`) == 'numeric')
})
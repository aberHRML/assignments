
context('assignMFs')

p <- assignmentParameters('FIE')

assignment <- assignMFs(peakData,p,verbose = TRUE)

test_that('assignMFs works',{
  expect_true(class(assignment) == "Assignment")
})
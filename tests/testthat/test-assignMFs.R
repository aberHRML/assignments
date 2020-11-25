
context('assignMFs')

p <- assignmentParameters('FIE')
p@nCores <- 2

assignment <- assignMFs(peakData,p,verbose = TRUE)

test_that('assignMFs works',{
  expect_true(class(assignment) == "Assignment")
})
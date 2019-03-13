
context('assignMFs')

p <- assignmentParameters('FIE')
p@nCores <- 2

cors <- correlations[correlations$r > 0.998 | correlations$r < -0.998,]

assignment <- assignMFs(cors,p,verbose = FALSE)

test_that('assignMFs works',{
  expect_true(class(assignment) == "Assignment")
})
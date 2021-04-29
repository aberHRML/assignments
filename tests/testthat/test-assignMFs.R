
context('assignMFs')

p <- assignmentParameters('FIE')

plan(future::multisession,workers = 2)

assignment <- assignMFs(peakData,
                        p,
                        verbose = TRUE)

test_that('assignMFs works',{
  expect_true(class(assignment) == "Assignment")
})
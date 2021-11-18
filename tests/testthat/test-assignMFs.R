
context('assignMFs')

p <- assignmentParameters('FIE')

plan(future::multisession,workers = 2)

assignment <- assignMFs(peakData,
                        p,
                        verbose = TRUE)

test_that('assignMFs works',{
  expect_s4_class(assignment,"Assignment")
})

test_that('feature solutions can be plotted',{
  pl <- plotFeatureSolutions(assignment,
                             'n191.01962',
                             maxComponents = 2)
  
  expect_s3_class(pl,'patchwork')
})
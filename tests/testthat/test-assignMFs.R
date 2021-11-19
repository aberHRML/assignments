
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

test_that('assignment network can be plotted',{
  pl <- plotNetwork(assignment)

  expect_s3_class(pl,'ggraph')
})

test_that('adduct distributions can be plotted',{
  pl <- plotAdductDist(assignment)
  
  expect_s3_class(pl,'patchwork')
})

test_that('assignment spectrum can be plotted',{
  pl <- plotSpectrum(assignment,'C6H8O7')
  
  expect_s3_class(pl,'ggplot')
})


assignment_parameters_FIE <- assignmentParameters()
assignment_parameters_LC <- assignmentParameters('RP-LC-HRMS')

test_adducts <- list(
  n = c("[M-H]1-",
        "[M+Cl]1-",
        "[M+K-2H]1-",
        "[2M-H]1-",
        "[M+Cl37]1-"),
  p = c("[M+H]1+"))
  
adducts(assignment_parameters_FIE) <- test_adducts
adducts(assignment_parameters_LC) <- test_adducts

LC_features <- new('Analysis')
metabolyseR::preTreated(LC_features) <- metabolyseR::analysisData(
  feature_data %>% 
    {magrittr::set_colnames(.,
                            paste0(colnames(.),
                                   '@1.00'))},
  info = tibble::tibble(
    ID = feature_data %>% 
      nrow() %>% 
      seq_len()))


assignment_FIE <- assignMFs(feature_data,
                            assignment_parameters_FIE,
                            verbose = TRUE)

assignment_LC <- assignMFs(LC_features,
                           assignment_parameters_LC,
                           verbose = TRUE)

test_that('assignment works for FIE technique',{
  expect_s4_class(assignment_FIE,"Assignment")
})

test_that('assignment works for LC techniques',{
  expect_s4_class(assignment_LC,"Assignment")
})

test_that('assignment class show method works',{
  expect_output(print(assignment_FIE),
                'Assignment:')
})

test_that('assignment data can be returned',{
  expect_s3_class(assignmentData(assignment_FIE),'tbl_df')
})

test_that('data with assigned feature names can be returned',{
  expect_s3_class(assignedData(assignment_FIE),'tbl_df')
})

test_that('a summary of assignments can be returned',{
  expect_s3_class(summariseAssignment(assignment_FIE),'tbl_df')
})
test_that('feature solutions can be plotted',{
  pl <- plotFeatureSolutions(assignment_FIE,
                             'n191.01962',
                             maxComponents = 2)
  
  expect_s3_class(pl,'patchwork')
})

test_that('feature solutions plotting throws an error if an incorrect feature is provided',{
  expect_error(plotFeatureSolutions(assignment_FIE,
                                    'test'))
})

test_that('assignment network can be plotted',{
  pl <- plotNetwork(assignment_FIE)
  
  expect_s3_class(pl,'ggraph')
})

test_that('adduct distributions can be plotted',{
  pl <- plotAdductDist(assignment_FIE)
  
  expect_s3_class(pl,'patchwork')
})

test_that('assignment spectrum can be plotted',{
  pl <- plotSpectrum(assignment_FIE,'C6H8O7')
  
  expect_s3_class(pl,'ggplot')
})

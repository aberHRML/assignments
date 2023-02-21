
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
transformations(assignment_parameters_LC) <- character()

FIE_features <-  new('Analysis')
metabolyseR::raw(FIE_features) <- metabolyseR::analysisData(
  feature_data,
  info = tibble::tibble(
    ID = feature_data %>% 
      nrow() %>% 
      seq_len()))

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


assignment_FIE <- assignMFs(FIE_features,
                            assignment_parameters_FIE,
                            type = 'raw',
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
  expect_s3_class(featureData(assignment_FIE),'tbl_df')
})

test_that('assignment correlations can be returned',{
  expect_s3_class(correlations(assignment_FIE),'tbl_df')
})

test_that('graph method throws an error if incorrect iteration specified',{
  expect_error(graph(assignment_FIE,'incorrect'))
})

test_that('component method throws an error if an incorrect component is specified',{
  expect_error(component(assignment_FIE,'incorrect','A&I1'))
})

test_that('data with assigned feature names can be returned',{
  expect_s3_class(assignedData(assignment_FIE),'tbl_df')
})

test_that('a summary of assignments can be returned',{
  expect_s3_class(summariseAssignments(assignment_FIE),'tbl_df')
})

test_that('feature components can be plotted',{
  pl <- plotFeatureComponents(assignment_FIE,
                             'n191.01962',
                             'A&I1')
  
  expect_s3_class(pl,'patchwork')
})

test_that('plotFeatureComponents throws an error if incorrect feature specified',{
  expect_error(plotFeatureComponents(assignment_FIE,
                              'incorrect',
                              'A&I1'))
})

test_that('plotFeatureComponents throws an error if no components are found',{
  expect_error(plotFeatureComponents(assignment_FIE,"n228.97636",'A&I1'))
})

test_that('a component can be plotted',{
  pl <- plotComponent(assignment_FIE,
                      1,
                      'A&I1')
  
  expect_s3_class(pl,'ggplot')
})

test_that('plotComponent throws an error if incorrect feature specified for highlighting',{
  expect_error(plotComponent(assignment_FIE,
                      1,
                      'A&I1',
                      highlight = 'incorrect'))
})

test_that('feature solutions plotting throws an error if an incorrect feature is provided',{
  expect_error(plotFeatureSolutions(assignment_FIE,
                                    'test'))
})

test_that('adduct distributions can be plotted',{
  pl <- plotAdductDist(assignment_FIE)
  
  expect_s3_class(pl,'patchwork')
})

test_that('assignment spectrum can be plotted',{
  pl <- plotSpectrum(assignment_FIE,'C6H8O7')
  
  expect_s3_class(pl,'ggplot')
})

test_that('Assignment class object can be created from a tibble',{
  expect_s4_class(assignment(feature_data,
                             assignment_parameters_FIE),
                  'Assignment')
})

test_that('Assignment class object can be created from an AnalysisData class object',{
  expect_s4_class(assignment(FIE_features %>% 
                               raw(),
                             assignment_parameters_FIE),
                  'Assignment')
})

test_that('Assignment class object can be created from an Analysis class object',{
  expect_s4_class(assignment(FIE_features,
                             assignment_parameters_FIE,
                             type = 'raw'),
                  'Assignment')
  
  expect_s4_class(assignment(LC_features,
                             assignment_parameters_FIE,
                             type = 'pre-treated'),
                  'Assignment')
})

test_that('assignment methods error correctly when slots are empty',{
  mf_assignments <- assignment(feature_data,
                               assignment_parameters_FIE)
  
  expect_error(calcRelationships(mf_assignments))
  expect_error(addIsoAssign(mf_assignments))
  expect_error(transformationAssign(mf_assignments))
})

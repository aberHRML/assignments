
setMethod('doAssignment',signature = 'Assignment',
          function(assignment){
            assignmentMethod <- assignMethods(assignment@parameters@technique)
            
            elements <- names(assignmentMethod())
            elements <- elements[!(elements %in% assignment@flags)]
            
            for(i in elements){
              method <- assignmentMethod(i)
              flag <- 'fail'
              try({
                assignment <- method(assignment)
                assignment@flags <- c(assignment@flags,i)
                flag <- 'success'
              })
              if (flag == 'fail') {
                cat('Failed at assignment step',i,'\n')
                return(assignment)
              }
            }
            return(assignment)
          }
)
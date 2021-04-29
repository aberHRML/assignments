
setMethod('doAssignment',signature = 'Assignment',
          function(assignment){
            assignmentMethod <- assignMethods(assignment@parameters@technique)
            
            elements <- names(assignmentMethod())
            elements <- elements[!(elements %in% assignment@flags)]
            
            for(i in elements){
              method <- assignmentMethod(i)
              assignment <- method(assignment)
              assignment@flags <- c(assignment@flags,i)
            }
            
            return(assignment)
          })
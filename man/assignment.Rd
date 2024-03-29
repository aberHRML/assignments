% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assignment.R
\name{assignment}
\alias{assignment}
\alias{assignment,tbl_df,AssignmentParameters-method}
\alias{assignment,AnalysisData,AssignmentParameters-method}
\alias{assignment,Analysis,AssignmentParameters-method}
\title{Create an Assignment S4 class object}
\usage{
assignment(feature_data, parameters, ...)

\S4method{assignment}{tbl_df,AssignmentParameters}(feature_data, parameters)

\S4method{assignment}{AnalysisData,AssignmentParameters}(feature_data, parameters)

\S4method{assignment}{Analysis,AssignmentParameters}(feature_data, parameters, type = c("pre-treated", "raw"))
}
\arguments{
\item{feature_data}{a tibble or an object of S4 class \code{AnalysisData} or \code{Analysis} containing the feature intensity matrix of m/z for which to assign molecular formulas. See details.}

\item{parameters}{an S4 object of class \code{AssignmentParamters} containing the parameters for molecular formula assignment}

\item{...}{arguments to pass to the relevant method}

\item{type}{type \code{pre-treated} or \code{raw} data on which to perform assignment when argument \code{feature_data} is of class \code{Analysis}}
}
\value{
An object of S4 class \code{Assignment}.
}
\description{
Constructor methods for creating an object of S4 class \code{Assignment}.
}
\examples{
mf_assignments <- assignment(feature_data,assignmentParameters('FIE-HRMS'))
mf_assignments
}

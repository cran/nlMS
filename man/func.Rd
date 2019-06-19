\name{auxiliary functions}
\alias{func}
\title{Auxiliary functions in Model Selection (MS)}
\usage{
func(mod, X, Y, name_mod, k=2, SST, N, file_name=NULL, form=NULL, make_plot=TRUE)
}
\description{
Auxiliary functions.
}
\arguments{
\item{form}{An object of class formula for the regression analysis.}
\item{mod}{Output of gnls.}
\item{X, Y}{Independent and dependent variables.}
\item{name_mod}{Name of the current regressionmodel preformed.}
\item{k}{Number of parameters used in the regression analysis.}
\item{SST}{Total Sum of Squares.}
\item{N}{X length.}
\item{file_name}{Name of the file where to save results.}
\item{make_plot}{Boolean value.}
}
\author{Carme Font <carme.font.moragon@gmail.com>}

\name{incidence}
\docType{data}
\alias{incidence}
\title{Age-specific incidence rates of disease}
\description{
A tab-separated text-file storing age-specific incidence rates for a defined calendar period.
}
\usage{incidence}
\format{A matrix.

\itemize{
\item{(a) The first row:}{file header consisting
of the keyword \samp{age} and the years for which incidence rates 
are available. \cr For example, it could look like
\tabular{rrrrrr}{
\samp{age} \tab \samp{1955} \tab \samp{1956} \tab \samp{1957} \tab ...  \tab \samp{2002} \cr 
}
} 
\item{(b) The first column:}{age groups (e.g. from 0 to 100)}
\item{(c) The remaining columns:}{age-specific incidence rates of disease}
} 
}

\source{CANCERMondial Statistical Information System [database online]: 
International Agency for Research on Cancer. Available at:
\emph{http://www-dep.iarc.fr/}
}
\keyword{datasets}

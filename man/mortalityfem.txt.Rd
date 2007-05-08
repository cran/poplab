\name{mortalityfem}
\docType{data}
\alias{mortalityfem}
\title{Age-specific mortality rates}
\description{
A tab-separated text-file storing age-specific mortality rates for a defined calendar period.
}
\usage{mortalityfem}
\format{A matrix.

\itemize{
\item{(a) The first row:}{file header consisting
of the keyword \samp{age} and the years for which mortality rates 
are available. \cr For example, it could look like
\tabular{rrrrrr}{
\samp{age} \tab \samp{1955} \tab \samp{1956} \tab \samp{1957} \tab ...  \tab \samp{2002} \cr 
}
} 
\item{(b) The first column:}{age groups (e.g. from 0 to 100)}
\item{(c) The remaining columns:}{age-specific mortality rates}
} 
}

\source{Statistics Sweden home-page:
\emph{http://www.scb.se}
}
\keyword{datasets}

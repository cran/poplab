\name{fertility}
\docType{data}
\alias{fertility}
\title{Age-specific female fertility rates}
\description{
A tab-separated text-file storing age-specific fertility rates for a defined calendar period.
}
\usage{fertility}
\format{A matrix.

\itemize{
\item{(a) The first row:}{file header consisting
of the keyword \samp{age} and the years for which fertility rates 
are available. \cr For example, it could look like
\tabular{rrrrrr}{
\samp{age} \tab \samp{1955} \tab \samp{1956} \tab \samp{1957} \tab ...  \tab \samp{2002} \cr 
}
} 
\item{(b) The first column:}{age groups (e.g. from 0 to 100)}
\item{(c) The remaining columns:}{age-specific fertility rates}
} 
}

\source{Johansson L. and Finnas F. \emph{Fertility of Swedish women born 1927-1960.}
Orebro: Statistics Sweden; 1983
}
\keyword{datasets}

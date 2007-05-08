\name{plot.poplab}
\alias{plot.poplab}
\title{Plot the age profile of a simulated population}

\description{
Plot a comparison between the simulated population and real population age profiles
for a specified calendar year, separately for females and males. 

The function looks in the specified folder for 2 tab-separated files containing
age-specific counts from the real female and male population. 
}

\usage{
\method{plot}{poplab} (x, option, population.fem, 
population.male, year, folder, \dots)
}
\arguments{
\item{x}{an object returned by the \code{\link[poplab]{create.baseline.complete}} or the
\code{\link[poplab]{simped}} functions.}
\item{option}{simulated population source. 2 possible values \samp{"base"} and \samp{"current"} (see Details). 
The default is \samp{"current"}.}
\item{population.fem}{name (quoted) of the tab-separated file storing age-specific counts from the real female
population (e.g. \samp{"populationfem.txt"}).}
\item{population.male}{name (quoted) of the tab-separated file storing age-specific counts from the real male
population (e.g. \samp{"populationmale.txt"})}
\item{year}{the year for which to plot the comparison}
\item{folder}{folder (quoted) where the function looks for the real population files. For example, 
it might be your R\_HOME/library/poplab/data/ directory, 
where R\_HOME is the location where R is installed on your computer (see Examples below).}
\item{\dots}{further arguments passed to or from other methods}
}

\details{
\code{option} specifies which simulated population to be displayed:\cr
\samp{"base"} requires a virtual baseline population (an object returned by the 
\code{\link[poplab]{create.baseline.complete}} function or the first list component 
returned by the \code{\link[poplab]{simped}} function). \cr
\samp{"current"} requires the population alive during the year specified by \code{year};
this is the second component of the list returned by the \code{\link[poplab]{simped}}
function.\cr
}

\examples{

\dontrun{
For the following examples it is necessary to load the simpop
object or to use the simped function and assign the result to
a variable called 'simpop'.
}

## load simulate the population for the calendar period 1955-2002:

load (file = file.path(system.file("data", package = "poplab"), 
"simpop.Rdata"))

## plot the simulated population versus the real population
## for the calendar year 1990:
plot.poplab(x = simpop, option = "current", population.fem =
"populationfem.txt", population.male = "populationmale.txt",
year = 1990, folder =
system.file("data", package = "poplab")) 

## plot the baseline population versus the real population
## for the calendar year 1955:
plot.poplab(x = simpop, option = "base", population.fem =
"populationfem.txt", population.male = "populationmale.txt",
year = 1955, folder =
system.file("data", package = "poplab")) 

}

\note{
(1) If \samp{yyyy} is the value of \code{year} and is smaller
than the \bold{last simulated year} (either the \samp{baseyear} when plotting a baseline
population or the \code{endyear} when plotting an evolved population), the
function will first exclude from the simulated population those individuals born after \samp{yyyy}. \cr

(2) If \samp{yyyy} is greater than the last simulated year, the function will plot
a comparison of the real population of the year \samp{yyyy} and the simulated population
that was still alive during the last simulated year, only older with the calendar
interval between the last simulated year and \samp{yyyy}. This is not a recommended use
of the function.
}


\seealso{
\code{\link[poplab]{create.baseline.complete}}, \code{\link[poplab]{simped}},
\code{\link[poplab]{print.poplab}}
}

\keyword{aplot}
\keyword{methods}

\name{simped}
\alias{simped}
\title{Simulate the evolution of a population of related individuals over a calendar period}

\description{
Simulate the evolution of a baseline population of related individuals over
a specified calendar period. Requires the period-specific fertility, mortality
and disease incidence rates, and the real population age profile for the
first simulated year.

The function looks in the folder specified in its call for 6 input data files: age-specific population
counts for males and females separately, female fertility, mortality rates for males and females separately,
and disease (e.g. breast cancer) incidence rates. 
}

\usage{
simped (baseyear, basehealth, basefamrisk, sex.a, basetotal,
warmuptime, simyears, endyear, healthy, famrisk, mortratio, 
d.mod, fam.rel, print.option, population.fem, population.male, 
mortality.fem, mortality.male, fertility, incidence, seed, 
folder, name.base)
}

\arguments{
\item{baseyear}{first calendar year to be simulated (the year for which to create
a baseline population - see \code{\link[poplab]{create.baseline.complete}}
help for more details).}
\item{basehealth}{0/1 values - disease incidence simulated/not simulated during the creation
of the baseline population (see \code{\link[poplab]{create.baseline.complete}} help). The 
default is 1 (i.e. disease incidence not simulated).}
\item{basefamrisk}{parameter for familial aggregation of disease applied in the creation
of the baseline population; used only if the \code{basehealth} is set to 0 (see
\code{\link[poplab]{create.baseline.complete}} help). The default is 1 (i.e. no familial
aggregation). See Details below.}
\item{sex.a}{1/2 values - gender specification (1 = males, 2 = females) for familial
aggregation of disease; used only if \code{healthy} or \code{basehealth} is set to 0.
The default is 2 (mother/daughter or sister/sister disease aggregation).}
\item{basetotal}{this factor multiplied by 100 is the sex-specific number of individuals in
the starting population (see \code{\link[poplab]{create.baseline.complete}} help). The
default is 500.}
\item{warmuptime}{number of years the \code{baseyear} real population age profile
is set back in time (see \code{\link[poplab]{create.baseline.complete}} help). The 
default is 100.}
\item{simyears}{number of years to be simulated.}
\item{endyear}{last calendar year to be simulated.}
\item{healthy}{0/1 values - disease incidence simulated/not simulated during the evolution
of the population. The default is 1.}
\item{famrisk}{parameter for familial aggregation of disease; used only if \code{healthy}
is set to 0. The default is 1 (i.e. no familial aggregation). See Details below.}
\item{mortratio}{mortality ratio for cases compared to disease free individuals. The
default is 1 (i.e. cases have the same risk of death as disease free individuals).}
\item{d.mod}{familial model of disease aggregation to choose from (see Details below); used 
only if \code{healthy} is set to 0.}
\item{fam.rel}{familial relationship of disease aggregation to choose from (see Details below);
used only if \code{healthy} is set to 0.}
\item{print.option}{TRUE/FALSE values - option for saving the simulated pedigree to file. The
default is FALSE.}
\item{population.fem}{name (quoted) of the tab-separated file storing age-specific counts
from the real female population (e.g. \samp{"populationfem.txt"}).}
\item{population.male}{name (quoted) of the tab-separated file storing age-specific counts
from the real male population (e.g. \samp{"populationmale.txt"}).}
\item{mortality.fem}{name (quoted) of the tab-separated file storing age-specific female
mortality rates (e.g. \samp{"mortalityfem.txt"}).}
\item{mortality.male}{name (quoted) of the tab-separated file storing age-specific male
mortality rates (e.g. \samp{"mortalitymale.txt"}).}
\item{fertility}{name (quoted) of the tab-separated file storing female fertility rates
(e.g. \samp{"fertility.txt"}).}
\item{incidence}{name (quoted) of the tab-separated file storing incidence rates (e.g. 
\samp{"incidence.txt"}).}
\item{seed}{seed used in the creation of the population. The default is \samp{NULL}. Useful when
the user wants to reproduce the exact same simulation results (population).}
\item{folder}{folder (quoted) where the function looks for input files such as population counts,
mortality, fertility and incidence rates. For example, it might be your R\_HOME/library/poplab/data/
directory, where R\_HOME is the location where R is installed on your computer (see Examples below).}
\item{name.base}{Ending string of the file name storing the baseline population. The default is "".
See Details below.}
}

\details{
\item{\code{baseyear}, \code{simyears} and \code{endyear}}{give the calendar interval
to be simulated; the user may specify any combination of 2 out of these 3 arguments and omit
the other argument.

The function will first identify if a virtual baseline population has been created and
saved to file, by searching into the specified \code{folder}
for a file named, for example, 'base\_pop\_1955.txt' if \code{baseyear} is 1955. 
\samp{base\_pop\_yyyy.txt}, where yyyy is a \code{baseyear}, is the name format that the 
\code{\link[poplab]{print.poplab}} function uses to save to file virtual baseline 
populations, after their creation with the \code{\link[poplab]{create.baseline.complete}}
function. This name format should not be used for other purposes. If the baseline
population is found, the \code{\link[poplab]{simped}} function reads in this
population and simulates further its evolution. If it is not found, a call to the 
\code{\link[poplab]{create.baseline.complete}} function is implicitly made by the simped
function, as follows:

    \preformatted{create.baseline.complete (baseyear, healthy = basehealth, 
    risk = basefamrisk, sex.a = sex.a, mortratio = mortratio,
    print.option = print.option, base.scale = basetotal, 
    runintime = warmuptime, population.fem = population.fem, 
    population.male = population.male, mortality.fem = mortality.fem,
    mortality.male = mortality.male, fertility = fertility,
    incidence = incidence, seed = seed, folder = folder)
}
Immediately after creating or reading in the baseline population, a plot of the simulated
baseline age profile versus the real population age profile for the baseyear year is displayed
and the function starts simulating the evolution of this population over the specified calendar
period.

}

\item{\code{d.mod}}{The options to choose from are: familial relative risk ("rr"),
familial odds-ratio ("or"), familial relative risk depending on the age at
incidence of the affected relative ("agesprr") and familial odds-ratio
depending on age at incidence of the affected relative ("agespor").
}

\item{\code{fam.rel}}{The options to choose from are: familial aggregation occurs
through a parent ("p") and familial aggregation occurs thorugh a sibling ("s"). 
}

\item{\code{basefamrisk}, \code{famrisk}}{If this is a single value, it is expected that the value for \code{d.mod}
is either "rr" or "or". For the familial models of disease where the risk and the odds,
respectively, depend on the age at incidence of the affected relative, a numerical 
array is specified. The first part of the array represents the upper limit of the age 
categories and the last part are the actual coefficients of increase. The age-intervals
are open to the left (see also Examples).}
\item{name.base}{When reading in the simulated baseline population the default of the file name
is base\_pop\_yyyy.txt (where yyyy is the baseline year). Other specifications (e.g., "fam\_risk2")
should be included in this parameter.}
}
\value{
A 2-component list. The first component (matrix) is the baseline population (either
read-in from the specified folder or just created). The second component (matrix) is
the evolved population. 
}

\examples{

\dontrun{
The syntax below will search into the specified folder for a
file called 'base_pop_1955.txt' (see Details). If the file
is found, it is read-in as the baseline population. If it is
not found, an implicit call to the create.baseline.complete
function is made as follows:

create.baseline.complete(baseyear = 1955, healthy = 0, risk =
2, sex.a = 2, base.scale = 50, runintime = 70, d.mod = 
"rr", fam.rel = "p", print.option = FALSE, population.fem = 
"populationfem.txt", population.male = "populationmale.txt", 
mortality.fem = "mortalityfem.txt", mortality.male = 
"mortalitymale.txt", fertility = "fertility.txt", incidence = 
"incidence.txt", seed = 52482337, folder =
system.file("data", package = "poplab"))

Next the evolution of this baseline population is simulated
over the calendar period 1955 - 1985 ([baseyear, endyear]).
Incident cancer is assigned (healthy) to the female population
(sex.a) via a model of familial relative risk (d.mod), where the 
risk ratio for women with an affected mother (fam.rel) is 2 (famrisk). 
The mortality ratio for cases is 2 (mortratio). The resulting 
pedigree is not saved to file (print.option). The input files 
are "populationfem.txt", "populationmale.txt", "mortalityfem.txt", 
"mortalitymale.txt", "fertility.txt" and "incidence.txt" from the "data"
folder of the poplab package.
}

simpop <- 
simped(baseyear = 1955, basehealth = 0, basefamrisk = 2,
sex.a = 2, basetotal = 50, warmuptime = 70, endyear = 1985,
healthy = 0, famrisk = 2, mortratio = 2, d.mod = "rr", 
fam.rel = "p", print.option = FALSE, population.fem = 
"populationfem.txt",  population.male = "populationmale.txt", 
mortality.fem = "mortalityfem.txt", mortality.male = 
"mortalitymale.txt", fertility = "fertility.txt", incidence = 
"incidence.txt", seed = 52482337,
folder = system.file("data", package = "poplab"))

## resulting baseline population
base <- simpop[[1]]

## resulting evolved population
evolved <- simpop[[2]]

\dontrun{
In the syntax below the familial relative risk (d.mod) depends
on maternal (fam.rel) age at incidence: the risk of disease (risk) 
is increased by factor of 4 for women whose mothers where strictly 
younger than 50 years of age at diagnosis compared to daughters 
of unaffected mothers, and by a factor of 2 for daughters of 
women diagnosed between the age of 50 to 100 years.

simpop <- 
simped(baseyear = 1955, basehealth = 0, basefamrisk = 
c(50, 101, 4, 2), sex.a = 2, basetotal = 50, warmuptime = 
70, endyear = 1985, healthy = 0, famrisk = c(50, 101, 4, 2), 
mortratio = 2, d.mod = "agesprr", fam.rel = "p", print.option = 
FALSE, population.fem = "populationfem.txt", population.male = 
"populationmale.txt", mortality.fem = "mortalityfem.txt", 
mortality.male = "mortalitymale.txt", fertility = 
"fertility.txt", incidence = "incidence.txt", 
seed = 52482337,
folder = system.file("data", package = "poplab"))
}

}

\note{
If a call to this function is made without assigning its value to a variable, the result will
be printed on the screen (most probably, a very large matrix)! Thus it is advisable to always assign
the object returned by this function to a variable.
}

\seealso{
\code{\link[poplab]{create.baseline.complete}}, \code{\link[poplab]{plot.poplab}},
\code{\link[poplab]{print.poplab}}
}

\keyword{datagen}

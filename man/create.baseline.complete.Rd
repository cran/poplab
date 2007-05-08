\name{create.baseline.complete}
\alias{create.baseline.complete}
\title{Create a baseline population of related individuals}

\description{
Create a baseline population for a given calendar year. Require fertility and mortality rates,
and age-profile for the specified year in order to create a population of related individuals of 
the correct age profile (the "trimmed" population), to which there are added the parents of these 
individuals who were excluded by the trimming process. Disease incidence is simulated when desired.
  
The function looks in the folder specified in its call for 6 input data files: age-specific population
counts for males and females separately, female fertility, mortality rates for males and females separately,
and disease (e.g. breast cancer) incidence rates.
}

\usage{
create.baseline.complete (baseyear, healthy, risk,
sex.a, mortratio, base.scale, runintime,
d.mod, fam.rel, print.option, population.fem, 
population.male, mortality.fem, mortality.male, 
fertility, incidence, seed, folder)
}

\arguments{
\item{baseyear}{year for which to create a baseline population.}
\item{healthy}{0/1 values - disease incidence simulated/not simulated during the creation
of the baseline population. The default is 1.}
\item{risk}{parameter for familial aggregation of disease; used only if \code{healthy}
is set to 0. The default is 1 (i.e. no familial aggregation). See Details below.}
\item{sex.a}{1/2 values - gender specification (1 = males, 2 = females) for familial
aggregation of disease; used only if \code{healthy} is set to 0. The default is 2 
(mother/daughter or sister/sister disease aggregation).}
\item{mortratio}{mortality ratio for cases compared to disease free individuals. The 
default is 1 (i.e. cases have the same risk of death as disease free individuals).}
\item{base.scale}{this factor multiplied by 100 is the sex-specific number of individuals
in the starting population. The default is 500 i.e. 50000 females and approximately
50000 males in the starting population (see Details below).}
\item{runintime}{number of years the \code{baseyear} real population age profile is set back
in time (see Details below). The default is 100.}
\item{d.mod}{familial model of disease aggregation to choose from (see Details below); used 
only if \code{healthy} is set to 0.}
\item{fam.rel}{familial relationship of disease aggregation to choose from (see Details below);
used only if \code{healthy} is set to 0.}
\item{print.option}{TRUE/FALSE values - option for saving to file the baseline population.
The default is FALSE.}
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
\item{incidence}{name (quoted) of the tab-separated file storing disease incidence rates
(e.g. \samp{"incidence.txt"}).}
\item{seed}{seed used in the creation of the population. The default is \samp{NULL}. Useful when
the user wants to reproduce the exact same simulation results (population).}
\item{folder}{folder (quoted) where the function looks for input files such as population counts,
mortality, fertility and incidence rates. For example, it might be your R\_HOME\/library\/poplab\/data\/
directory, where R\_HOME is the location where R is installed on your computer (see Examples below).}
}

\value{the matrix containing the baseline population}

\details{
\code{base.scale} The exact number of male individuals is calculated in such a way
as to maintain the real ratio between the number of male individuals and the
number of female individuals.

\code{runintime} To create a baseline population of related
individuals for the \code{baseyear}, a population of unrelated individuals
(the starting population) is constructed for the calendar year calculated
as \samp{(baseyear - runintime)}; this starting population has the age profile
of the \samp{baseyear} real population. For the number of years specified
by \code{runintime}, the \code{baseyear} fertility and mortality
rates are applied to this starting population and their descendants,
"arriving" in the baseyear with a complete population of related
individuals.

\code{d.mod} The options to choose from are: familial relative risk ("rr"),
familial odds-ratio ("or"), familial relative risk depending on the age at
incidence of the affected relative ("agesprr") and familial odds-ratio
depending on age at incidence of the affected relative ("agespor").

\code{fam.rel} The options to choose from are: familial aggregation occurs
through a parent ("p") and familial aggregation occurs thorugh a sibling ("s"). 

\code{risk} If this is a single value, it is expected that the value for \code{d.mod}
is either "rr" or "or". For the familial models of disease where the risk and the odds,
respectively, depend on the age at incidence of the affected relative, a numerical 
array is specified. The first part of the array represents the upper limit of the age 
categories and the last part are the actual coefficients of increase. The age-intervals
are open to the left (see also Examples).  
}
\examples{

\dontrun{
The syntax below constructs the baseline population for the
calendar year 1955 (baseyear). It applies the 1955 fertility
and mortality rates for 70 years (runintime) to the virtual
population of 1855 constructed from the 1955 age profile; the
1855 population consists of 5000 females and approximately
5000 males (base.scale), all unrelated individuals. Incident
cancer is assigned during the creation (healthy), but there is
no familial risk of cancer (risk). The mortality ratio for cases
is 2 (mortratio). The resulting population is not saved to file
(print.option). The input files are "populationfem.txt",
"populationmale.txt", "mortalityfem.txt", "mortalitymale.txt",
"fertility.txt" and "incidence.txt" from the "data" folder of the
poplab package.}

base <- 
create.baseline.complete(baseyear = 1955, healthy = 0, risk =
1, mortratio = 2, d.mod = "rr", fam.rel = "p", base.scale = 
50, runintime = 70, print.option = FALSE, population.fem = 
"populationfem.txt", population.male = "populationmale.txt", 
mortality.fem = "mortalityfem.txt", mortality.male = 
"mortalitymale.txt", fertility = "fertility.txt", incidence = 
"incidence.txt", seed = 52482337, folder =
system.file("data", package = "poplab"))

\dontrun{
In the syntax below the familial relative risk (d.mod) depends
on maternal (fam.rel) age at incidence: the risk of disease (risk)
is increased by factor of 4 for women whose mothers where strictly 
younger than 50 years of age at diagnosis compared to daughters 
of unaffected mothers, and by a factor of 2 for daughters of 
women diagnosed between the age of 50 to 100 years.

base <- 
create.baseline.complete(baseyear = 1955, healthy = 0, risk =
c(50, 101, 4, 2), mortratio = 2, d.mod = "agesprr", fam.rel =
"p", base.scale = 70, runintime = 50, print.option = FALSE, 
population.fem = "populationfem.txt", population.male = 
"populationmale.txt", mortality.fem = "mortalityfem.txt", 
mortality.male = "mortalitymale.txt", fertility = 
"fertility.txt", incidence = "incidence.txt",
seed = 52482337, folder =
system.file("data", package = "poplab"))
}

}

\note{
If a call to this function is made without assigning its value to a variable,
the result will be printed on the screen (most probably, a very large matrix!).
Thus it is advisable to always assign the object returned by this function to a
variable.
}

\seealso{
\code{\link[poplab]{simped}}, \code{\link[poplab]{plot.poplab}},
\code{\link[poplab]{print.poplab}}
}

\keyword{datagen}

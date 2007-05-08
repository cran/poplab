\name{assigncancer}
\alias{assigncancer}
\title{Assign incident disease (cancer)}

\description{
Not intended to be called directly.  
} 

\usage{
assigncancer(ped, inc, risk, yr, sex, type, fam.rel) 
}

\arguments{
\item{ped}{an object storing pedigree information}
\item{inc}{an object specifying population (age-specific) incidence rates for the simulated year}
\item{risk}{multiplication value of the population incidence rates}
\item{yr}{year for which to assign incident disease}
\item{sex}{1/2 values - gender specification (1 = males, 2 = females)}
\item{type}{familial model of disease aggregation to choose from ("rr", "or", "agesprr", "agespor")}
\item{fam.rel}{familial relationship of disease aggregation to choose from}
}

\keyword{internal}

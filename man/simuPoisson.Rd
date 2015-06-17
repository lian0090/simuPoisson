\name{simuPoisson}
\alias{simuPoisson}
\title{simulate genotypes from Poisson process}
\description{Simulate recombination from Poisson process}
\usage{
	simuPoisson(parentsGeno,chr,cM,N)
	}
\arguments{
\item{parentsGeno}{ genotype matrix for parents on all chromosomes: a data.frame or matrix with two rows, each row represent the marker gentoypes for one parent. Genotypes must be coded additively: the heterozygotes genotype must be coded as half the two homozygotes gentoypes }
\item{chr}{  a vector of size p for the chromosome numbers of all markers }
\item{cM}{a vector of size p for the centiMorgan map positions for each marker}
\item{N}{ total number of lines to simulate }
\item{mapchr}{map for a single chromosome}
}
\details{
simulation of genotypes is done for each individual chromosome.
minimum and maximum number of recombinations was set up based on the 95\% quantile of Poisson distribution
minrecmb=qpois(0.025,chrlength/100,lower.tail=T)
maxrecmb=qpois(0.975,chrlength/100,lower.tail=T)
	
}
\value{
\item{simuPoisson}{returns a matrix with N rows, each row is the simulated genotype for an individual} 
}
\author{Lian Lian}
\examples{
	data(parentsGeno)
	data(map)
dat=simuPoisson(parentsGeno,map$chr,map$cM,10)
}

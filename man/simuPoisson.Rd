\name{simuPoisson}
\alias{simuPoisson}
\alias{samplegamete}
\alias{samplegametechr}
\title{simulate genotypes from Poisson process}
\description{Simulate recombination form Poisson process}
\usage{
	simuPoisson(fpgenomat,map,gtype,nptotal,...)
	samplegamete(fpgenomat,map,...)
	samplegametechr(fpchrmat,mapchr,...)
	}
\arguments{
\item{fpgenomat}{ genotype matrix for parents on all chromosomes: a data.frame or matrix with 2*nmar, colnames should be mknames }
\item{map}{ a data.frame with nmar*(3-4), colnames should include mkname, chr, pos }
\item{gtype}{ the cross type of the project:F2,DH }
\item{nptotal}{ total number of lines to simulate }
\item{fpchrmat}{fpgenomat for a single chromosome}
\item{mapchr}{map for a single chromosome}
\item{\dots}{ savedir: the save directory of the output file, if NULL or not specified, the results will not be saved;
projectname:if savedir is true, the file will be saved as projectname.rds;
}
}
\details{
simulation of genotypes is done for each individual chromosome.
minimum and maximum number of recombinations cut off based on Morgan map distance is based on the 95\% quantile of Poisson distribution
minrecmb=qpois(0.025,chrlength/100,lower.tail=T)
maxrecmb=qpois(0.975,chrlength/100,lower.tail=T)
	
}
\value{
\item{samplegametechr}{returns a gamete on a single chromosome}
\item{samplegamete}{returns gamete on all chromosomes}
\item{simuPoisson}{returns a combinedgeno object first four rows are map_chr and map_pos,parent genotypes, colnames are "LINE" and marker names, other rows are marker genotypes for the simulated progenies} 
}
\author{Lian Lian}
\examples{
	data(fpgenomat)
	data(map)
dat=simuPoisson(fpgenomat,map,"F2",10)
dat[1:10,1:10]
samplegamete(fpgenomat[,1:300],map[1:300,])
samplegametechr(fpgenomat[,1:30],map[1:30,])
}
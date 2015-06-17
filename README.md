# simuPoisson
Simulate F2 populations based on Poisson distribution for the number of recombinants

This is currently only for F2 populations.

     (1) For each chrosome, number of cross-over events is sampled from a truncated Poisson distribution 
     (2) For each coss-over events, the position of this event on the chrosome is sampled from a uniform distribution along the chrosome intervals. The restrition is that two cross-over events need to be 20 cM distance from each other. 
 

## Install
If `devtools` is not installed, install devtools first. 

```R
install.packages("devtools")
library(devtools)
install_github("lian0090/simuPoisson")
```


## Usage
`simuPoisson(parentsGeno,chr,cM,N)`

- Arguments
    - `parentsGeno`: marker genotypes for two parents. Each row is an individual. Genotypes must be coded additively: the coded value for heterozygotes must be half the coded value of the two homozygotes. For example, `-1,0,1` or `0,1,2`. 
    - `chr`: a vector of size p for the chromosome numbers of all markers
    - `cM`: a vector of size p for the centiMorgan map positions for each marker
    - `N`: total number of lines to simulate

- Return Values
    - returns a matrix with N rows, each row is the simulated genotype for an individual


## Example
```R
library(simuPoisson)
data(parentsGeno)
data(map)
#simulate 10 individuals
progeny=simuPoisson(parentsGeno,map$chr,map$cM,10)
```

# simuPoisson
Simulate F2 populations based on Poisson distribution for the number of cross-over events

For each chromosome, a random haplotype was sampled as either from parent A or parent B, with crossing over occurring at random. The expected number of crossovers (L) in this haplotype was the length of the chromosome in Morgans, whereas the observed number of crossovers for the haplotype was sampled from a Poisson distribution with a mean of L (with restrictions that the number of cross-overs should be confined to the 0.025 and 0.975 quantile of the Poisson distribution). Crossover positions were randomly sampled along a chromosome according to a uniform distribution. To account for interference, two adjacent crossovers were arbitrarily assumed to be at least 10 cM from each other. 


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

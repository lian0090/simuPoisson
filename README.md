# simuPoisson
Simulate F2 populations based on Poisson distribution for the number of recombinants

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

library(simuPoisson)
data(parentsGeno)
data(map)

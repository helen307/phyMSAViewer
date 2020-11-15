# phyMSAViewer 

<!-- badges: start -->
<!-- badges: end -->

## Description
The goal of phyMSAViewer is to provide a general overview of a phylogenetic tree with its multiple sequence alignment on the side. The users can either provide a sequence alignment file in fasta or can just type in the Uniprot IDs manually to get the results. Another unique function is that the users are able to zoom on a particular range of the multiple sequence alignment from the R ShinyApp.

## Installation

You can install the released version of phyMSAViewer from GitHub with:

``` r
library(devtools)
devtools::install_github("helen307/phyMSAViewer", build_vignettes=TRUE)
library(phyMSAViewer)
```

## Example

This is a basic example which shows you how to get a phylogenetic graph with MSA from the Uniprot ID entered:

``` r
library(phyMSAViewer)
treeWithMSA <- uniprotToPhy("AC=Q9H9K5 OR AC=P04439 OR AC=P01889")
treeWithMSA
```

## Overview
ls("package:phyMSAViewer")
data(package="phyMSAViewer")
browseVignettes("phyMSAViewer")

There are mainly two functions in this package, including producing phylogenetic trees and multiple sequence alignment from Uniprot IDs and from fasta files. A sample dataset is provided as an example, and we will require the users to provide a file of the same format as input to `seqToPhy`. For `uniproToPhy`, please provide a correct character string, separated by " OR ".

An overview of the package is shown below:
![Overview of phyMSAViewer](/Users/helending/Desktop/phyMSAViewer/man/figures/pic.png){width=70%}

## Contributions
The author of the package is Yining Ding. The uniprotToPhy function makes use of msaplot from ggtree package to generate the phylogenetic tree with multiple sequence alignment. The msa R package is used for performing MUSCLE multiple sequence alignment. The R shiny package was used to produce the Shiny App.


## References


## Acknowledgements


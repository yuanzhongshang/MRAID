# MRAID

MRAID(Mendelian Randomization with Automated Instrument Determination),is an R package for efficient statistical inference of two-sample Mendelian Randomization. MRAID takes GWAS summary statistics as inputs to estimate causal effects of one trait on another.  MRAID is able to model an initial set of candidate SNP instruments that are in high LD with each other and perform automated instrument selection to identify suitable SNPs to serve as instrumental variables. MRAID simultaneously accounts for both uncorrelated and correlated horizontal pleiotropy, relies on a scalable sampling-based inference algorithm to perform numerical integration, circumventing the difficulty in likelihood function, leading to calibrated p-values that enable reasonably large-scale exposure screening.

This approach is described in,
> Yuan Z, Liu L, Guo P, Yan R, Xue F, Zhou X.  [Likelihood-based Mendelian randomization analysis with automated instrument selection and horizontal pleiotropic modeling](https://xzlab.org/papers/2022_Yuan_etal_SA.pdf). Science Advances 8, eabl5744.

# Installation
It is easy to install the development version of MRAID package using the 'devtools' package. 

```
# install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/MRAID")
```
# Usage
The main function in the package is MRAID, you can find the instructions by '?MRAID'.
```
library(MRAID)

?MRAID
```

# Example
One simple example to use the package can be found at https://github.com/yuanzhongshang/MRAID/tree/master/example

# Development
This R package is developed by Zhongshang Yuan and Lu Liu.

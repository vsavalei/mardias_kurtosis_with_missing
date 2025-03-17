
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Computing Mardia’s Kurtosis with Missing Data

The provided functions compute two versions of Mardia’s Kurtosis (MK):
one due to Yuan, Lambert, and Fouladi (2004) and one due to Savalei
(2010, unpublished, in prep for 2025). The original YLF equation is only
valid for MCAR data. The updated version is valid for MCAR and MAR data.
The YLF version requires sorting the data by missing data patterns and
looping through and adding computations across patterns. The updated
MAR-consistent version is almost trivial to compute because it depends
on `lavaan`’s internal matrices that are now accessible to the user
(Savalei and Rosseel, 2020). In the current code, both functions require
raw data as input, and they fit the saturated model to the data to
compute MK.

To load the functions, download the `Mardias_with_missing_functions.R`
file, set the working directory to where this file is stored, and source
the code:

``` r
source("Mardias_with_missing_functions.R")
#> This is lavaan 0.6-20.2276
#> lavaan is FREE software! Please report any bugs.
```

We illustrate the functions on a simulated normally distributed dataset
(also available for download) that contains 600 observations on 18
variables, where some have missing data (simulated to be MAR):

``` r
datai<-read.table("n600MAR15_rep1.dat",sep=",", na.strings = "NA",header=TRUE)
```

The following function computes the YLF MK value, its squared standard
error (sampling variance), and the associated z-test:

``` r
mylf <- mardia_ylf(datai)
mylf
#> $mYLF
#> [1] -1.625388
#> 
#> $mYLF_var
#> [1] 4.368889
#> 
#> $mYLF_ztest
#> [1] -0.7776273
```

<!-- old code gave: mYLF = -1.625496; old EQS code gave: -1.6241. current code gives: -1.625388 -->

The following function computes the updated MAR-consistent MK value, its
approximate sampling variance, and the associated z-test:

``` r
mmar <- mardia_mar(datai)
mmar
#> $mmar
#> [1] -0.2393642
#> 
#> $mmarvar
#> [1] 4.8
#> 
#> $ztest
#> [1] -0.1092543
```

<!-- old code: -.2393470; current code:  -0.2393642 -->

The raw value of MK is centered: It is expected to be zero when data are
normal.

Caution should be exercised with taking the standard error and z-test
for MK too literally, especially in small samples. For the
MAR-consistent version, the standard error is only approximate.
Bootstrapping can be used instead (functions to be added shortly). As
with all significance tests, in small samples nonnormality may not be
detectable due to low power, but it can still affect the results of SEM
analyses. In large samples trivial deviations from nonnormality may be
detected. It’s good advice to simply use `estimator='MLR'` when in doubt
or when any degree of nonnormality is present, at least as a robustness
check.

When data are complete, the two version of Mardia’s kurtosis are equal:

``` r
datac<-datai[complete.cases(datai),] #N = 272, down from 600
mylfc <- mardia_ylf(datac)
mylfc
#> $mYLF
#> [1] 2.570061
#> 
#> $mYLF_var
#> [1] 10.58824
#> 
#> $mYLF_ztest
#> [1] 0.7898262
mmarc <- mardia_mar(datac)
mmarc
#> $mmar
#> [1] 2.570061
#> 
#> $mmarvar
#> [1] 10.58824
#> 
#> $ztest
#> [1] 0.7898262
```

<!-- m = 2.570061, var=10.58824, z=.7898262 -->

They both reduce to the usual Mardia’s kurtosis for complete data. Below
is output from `mardiaKurtosis` function `semTools`:

``` r
library(semTools)
#> 
#> ###############################################################################
#> This is semTools 0.5-6
#> All users of R (or SEM) are invited to submit functions or ideas for functions.
#> ###############################################################################
mcomp <- mardiaKurtosis(datac) #oops
mcomp
#>          b2d            z            p 
#> 359.90900484  -0.02796447   0.97769049
```

The implementation in `semTools` computes the kurtosis value that is
uncentered; i.e., it is not expected to be zero when data are normal.
Its expected value under normality is $p(p+2)=18*20=360$, so the
centered value would be $359.9090048 - 360 = -0.0909952$. The remaining
small difference is because the function in `semTools` uses a sample
covariance matrix with $N-1$ rather than $N$ divisor. Below is
`semTools` source code with this change made:

``` r
#getting semTools to match values 
dat <- datac
centeredDat <- scale(dat, center = TRUE, scale = FALSE)
S <- cov(dat)*(nrow(dat)-1)/nrow(dat) #added this line 
invS <- solve(S)  #changed cov(dat) to S
    FUN <- function(vec, invS) {
        as.numeric(t(as.matrix(vec)) %*% invS %*% as.matrix(vec))
    }
indivTerm <- sapply(as.list(data.frame(t(centeredDat))), 
        FUN, invS = invS)
b2d <- sum(indivTerm^2)/nrow(dat)
d <- ncol(dat)
m <- d * (d + 2)
v <- 8 * d * (d + 2)/nrow(dat)
z <- (b2d - m)/sqrt(v)
p <- pnorm(-abs(z)) * 2

b2d-m #centered value
#> [1] 2.570061
```

In the future version of the MK functions, the function for
MAR-consistent MK will also work with just a `lavaan` object (any model
run on the data) as input.

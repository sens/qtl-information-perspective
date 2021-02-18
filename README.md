# Quantitative trait locus study design from an information perspective S Sen, JM Satagopan, and GA Churchill
-----------------------------------------------------------------------------------------------------------

[paper](https://www.genetics.org/content/170/1/447.long) ~ [symbolic](#symbolic) ~ [figures](#figures)


This page contains supplementary material for Sen S, Satagopan JM,
Churchill GA (2005) "Quantitative trait locus study design from an
information perspective," Genetics, 170:447-464. Please email Saunak
Sen if you have any problems or questions about the contents of this
webpage.


## Symbolic computation code <a name="symbolic"></a>

Some of the results in the paper were derived using symbolic
calculations in [Maxima](http://maxima.sourceforge.net). Start Maxima
in your system, and then you can cut and paste the contents of the
files below into the command window. The files are commented, so you
should be able to follow the steps.

-   Formula for missing information in backcross:
    [bc-missing.max](bc-missing.max)
-   Calculating the determinant and inverse of the information matrix
    for F~2~\'s:
    [f2det.max](f2det.max)
-   Formula for missing information in backcross in the presence of a
    second QTL, assuming that first QTL has small effect:
    [2qtl.max](2qtl.max)

## Figures <a name="figures"></a>

-   Figure 1:
    [genopat.m](genopat.m);
    this uses [Pseudomarker version
    0.9](http://www.jax.org/staff/churchill/labsite/software/pseudomarker/pseudomarker0_9/index.html)
    written in Matlab, and the [salt-induced hypertension
    data](hypertension1/index.html)
    from Sugiyama et.al. (2001)
-   Figure 2:
    [chr4.R](chr4.R);
    this uses the [R/qtl](http://www.rqtl.org)
    package
-   Figure 3:
    [numerical.R](numerical.R)
-   Figures 4 and 5:
    [optimal.R](optimal.R)
-   Figure 6:
    [opt-alpha.R](opt-alpha.R);
    this uses the R/qtlDesign package version 0.32 ([see
    below](#qtldesign)).
-   Figure 7:
    [replication.R](replication.R)
-   Figure 8:
    [2qtl.R](2qtl.R)

## R/qtlDesign <a id="qtldesign"><a/> 

This package performs power calculations and minimum effect size
determinations for backcross and F~2~ intercross populations. These
calculations take into account selective genotyping of the extreme
phenotypic individuals and marker spacing. It is an add-on package to
the [R](http://www.r-project.org) programming language. To install
version 0.32 of the package (in UNIX or OS X) download the file
[qtlDesign\_0.32.tar.gz](qtlDesign_0.32.tar.gz),
and give the type in a command window:

       R CMD INSTALL qtlDesign_0.32.tar.gz

For more recent versions of the package see the
[software](http://www.epibiostat.ucsf.edu/biostat/sen/software.html)
page.

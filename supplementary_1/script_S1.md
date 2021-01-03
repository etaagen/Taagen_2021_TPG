script\_S1: SynOpDH analysis
================

All packages, data, and statistical analysis for reproducing SynOpDH
population results reported in Taagen et al. 2021. Please see
`script_S1.Rmd` for full R script.

**Load packages**

``` r
library(tidyverse) # R/tidyverse version 1.3.0
library(lme4) #R/lme4 package version 1.1-21 
library(knitr) # R/knitr package version 1.28 
library(kableExtra) # R/kableExtra package version 1.1.0
library(qtl) #R/qtl package version 1.46-2
library(tibble) # R/tible package version 2.1.3 
library(ggpubr) # R/ggpubr package version 0.2.5 
library(scales) # R/scales package version 1.1.0
library(svglite) #R/svglite version 1.2.3.2
library(gt) #R/gt version 0.2.2
library(car) # R/car package version 3.0-10
library(emmeans) # R/emmeans package version 1.4.6
library(rstatix) # R/rstatix package version 0.6.0
```

### SynOpDH broad sense heritability and BLUP phenotypes

A 162-entry subset of 215 SynOpDH entries were grown in headrows in four
field-year combinations, with up to six replicates per entry, from
2016-2018 in Ithaca, NY. Univariate mixed linear models with random
environment and genotype effects were fitted with the R/lme4 package to
obtain best linear unbiased predictions (BLUPs) for TGW, GL, GW and HD
phenotypes across the four environments.

**Load data**

``` r
SynOpDH_Phenotypes <- read_csv("https://raw.githubusercontent.com/etaagen/Taagen_2021_TPG/main/supplementary_1/file_S1.1.csv", na = "")
SynOpDH_Phenotypes$SynOpDH_entry = as.factor(SynOpDH_Phenotypes$SynOpDH_entry)
SynOpDH_Phenotypes$Year = as.factor(SynOpDH_Phenotypes$Year)
SynOpDH_Phenotypes$Location = as.factor(SynOpDH_Phenotypes$Location)
SynOpDH_Phenotypes$Rep = as.factor(SynOpDH_Phenotypes$Rep)
SynOpDH_Phenotypes$Environment = as.factor(SynOpDH_Phenotypes$Environment)
```

**Falconer & Mackay `H^2` estimate function**

\[ 
H^2 = \frac{\sigma^2_G}{\sigma^2_G + \frac{\sigma^2_{GxE}}{l} + \frac{\sigma^2_r}{rl}}
\]

Where \(H^2\) is the broad-sense heritability estimate, \(\sigma^2_G\)
is the genetic variance, \(\sigma^2_{GxE}\) is the genotype x
environment variance, \(\sigma^2_r\) is the residual variance, \(l\) is
the number of environments, and \(r\) is the number of unique
observations.

**TGW**

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

Mixed model<sup>a</sup>

</th>

<th style="text-align:right;">

TGW H^2<sup>b</sup>

</th>

<th style="text-align:right;">

AIC<sup>c</sup>

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

TGW \~ (1|Entry)

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

4638.53

</td>

</tr>

<tr>

<td style="text-align:left;">

TGW \~ (1|Entry) + (1|Env)

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

3710.14

</td>

</tr>

<tr>

<td style="text-align:left;">

TGW \~ (1|Entry) + (1|Env) + Rep

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

3710.86

</td>

</tr>

<tr>

<td style="text-align:left;">

TGW \~ (1|Entry) + (1|Rep:Env)

</td>

<td style="text-align:right;">

0.69

</td>

<td style="text-align:right;">

3709.09

</td>

</tr>

<tr>

<td style="text-align:left;">

TGW \~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env)

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

3704.04

</td>

</tr>

<tr>

<td style="text-align:left;">

TGW \~ (1|Entry) + (1|Rep:Env) + HD

</td>

<td style="text-align:right;">

0.68

</td>

<td style="text-align:right;">

2993.01

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<sup>a</sup> Mixed model variables are defined as: thousand grain weight
(TGW), SynOpDH line entry (Entry), year environment unique to location
and year (Env), replicate (Rep), heading date (HD), <sup>b</sup> Broad
sense heritability calculated with Falconer & Mackay 1996 formula
<sup>c</sup> AIC was calculated using `AIC` in R/stats package

</td>

</tr>

</tfoot>

</table>

**GL**

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

Mixed model<sup>a</sup>

</th>

<th style="text-align:right;">

GL H^2<sup>b</sup>

</th>

<th style="text-align:right;">

AIC<sup>c</sup>

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GL \~ (1|Entry)

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

805.15

</td>

</tr>

<tr>

<td style="text-align:left;">

GL \~ (1|Entry) + (1|Env)

</td>

<td style="text-align:right;">

0.82

</td>

<td style="text-align:right;">

70.58

</td>

</tr>

<tr>

<td style="text-align:left;">

GL \~ (1|Entry) + (1|Env) + Rep

</td>

<td style="text-align:right;">

0.82

</td>

<td style="text-align:right;">

70.93

</td>

</tr>

<tr>

<td style="text-align:left;">

GL \~ (1|Entry) + (1|Rep:Env)

</td>

<td style="text-align:right;">

0.83

</td>

<td style="text-align:right;">

79.33

</td>

</tr>

<tr>

<td style="text-align:left;">

GL \~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env)

</td>

<td style="text-align:right;">

0.75

</td>

<td style="text-align:right;">

\-50.85

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<sup>a</sup> Mixed model variables are defined as: grain length (GL),
SynOpDH line entry (Entry), year environment unique to location and year
(Env), replicate (Rep), heading date (HD), <sup>b</sup> Broad sense
heritability calculated with Falconer & Mackay 1996 formula <sup>c</sup>
AIC was calculated using `AIC` in R/stats package

</td>

</tr>

</tfoot>

</table>

**GW**

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

Mixed model<sup>a</sup>

</th>

<th style="text-align:right;">

GW H^2<sup>b</sup>

</th>

<th style="text-align:right;">

AIC<sup>c</sup>

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry)

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

\-192.52

</td>

</tr>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry) + (1|Env)

</td>

<td style="text-align:right;">

0.72

</td>

<td style="text-align:right;">

\-1063.32

</td>

</tr>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry) + (1|Env) + Rep

</td>

<td style="text-align:right;">

0.72

</td>

<td style="text-align:right;">

\-1062.49

</td>

</tr>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry) + (1|Rep:Env)

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

\-1055.86

</td>

</tr>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env)

</td>

<td style="text-align:right;">

0.64

</td>

<td style="text-align:right;">

\-1130.19

</td>

</tr>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env) + HD

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:right;">

\-1136.80

</td>

</tr>

<tr>

<td style="text-align:left;">

GW \~ (1|Entry) + (1|Rep:Env) + HD

</td>

<td style="text-align:right;">

0.81

</td>

<td style="text-align:right;">

\-1137.19

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<sup>a</sup> Mixed model variables are defined as: grain width (GW) ,
SynOpDH line entry (Entry), year environment unique to location and year
(Env), replicate (Rep), heading date (HD), <sup>b</sup> Broad sense
heritability calculated with Falconer & Mackay 1996 formula <sup>c</sup>
AIC was calculated using `AIC` in R/stats package

</td>

</tr>

</tfoot>

</table>

**HD**

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

Mixed model<sup>a</sup>

</th>

<th style="text-align:right;">

HD ^2<sup>b</sup>

</th>

<th style="text-align:right;">

AIC<sup>c</sup>

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

HD \~ (1|Entry)

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

2951.50

</td>

</tr>

<tr>

<td style="text-align:left;">

HD \~ (1|Entry) + (1|Env)

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

2764.46

</td>

</tr>

<tr>

<td style="text-align:left;">

HD \~ (1|Entry) + (1|Env) + Rep

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

2764.58

</td>

</tr>

<tr>

<td style="text-align:left;">

HD \~ (1|Entry) + (1|Rep:Env)

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

2769.31

</td>

</tr>

<tr>

<td style="text-align:left;">

HD \~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env)

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:right;">

2756.02

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<sup>a</sup> Mixed model variables are defined as: heading date (HD) ,
SynOpDH line entry (Entry), year environment unique to location and year
(Env), replicate (Rep) <sup>b</sup> Broad sense heritability calculated
with Falconer & Mackay 1996 formula <sup>c</sup> AIC was calculated
using `AIC` in R/stats package

</td>

</tr>

</tfoot>

</table>

### Select models, extract BLUPs

The models were selected based on low AIC value and improved H^2
estimate:

  - `TGW ~ (1|Entry) + (1|Rep:Env) + HD`

  - `GL ~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env)`

  - `GW ~ (1|Entry) + (1|Rep:Env) + HD`

  - `HD ~ (1|Entry) + (1|Rep:Env) + (1|Entry:Env)`

**Extract random effects from univariate mixed linear models for entries
(obtain BLUPs):**

``` r
BLUP_TGW <- ranef(TGW_model5)$SynOpDH_entry 
BLUP_GL <- ranef(GL_model3.5)$SynOpDH_entry 
BLUP_GW <- ranef(GW_model5)$SynOpDH_entry 
BLUP_HD <- ranef(HD_model3.5)$SynOpDH_entry 
```

**Phenotype means, add to BLUP for scaled phenotype value:**

``` r
TGW_mean <- mean(SynOpDH_Phenotypes$TGW, na.rm = TRUE) 
BLUP_TGW <- BLUP_TGW + TGW_mean
names(BLUP_TGW)[1] <- "TGW"

GL_mean <- mean(SynOpDH_Phenotypes$GL, na.rm = TRUE)
BLUP_GL <- BLUP_GL + GL_mean
names(BLUP_GL)[1] <- "GL"

GW_mean <- mean(SynOpDH_Phenotypes$GW, na.rm = TRUE)
BLUP_GW <- BLUP_GW + GW_mean
names(BLUP_GW)[1] <- "GW"

HD_mean <- mean(SynOpDH_Phenotypes$HD, na.rm = TRUE)
BLUP_HD <- BLUP_HD + HD_mean
names(BLUP_HD)[1] <- "HD"
```

**Create csv files**  
*eval = FALSE, BLUP files in repo, `file_S1.8.csv`*

``` r
write.csv(BLUP_TGW, "SynOpDH_BLUP_TGW.csv")
write.csv(BLUP_GL, "SynOpDH_BLUP_GL.csv")
write.csv(BLUP_GW, "SynOpDH_BLUP_GW.csv")
write.csv(BLUP_HD, "SynOpDH_BLUP_HD.csv")
```

### SynOpDH Genetic Linkage Map

A genetic linkage map of the SynOpDH population was constructed from a
subset of 1,551 [GBS](https://doi.org/10.1371/journal.pone.0032253) and
[SSR](https://wheat.pw.usda.gov/cgi-bin/GG3/browse.cgi?class=marker)
markers using the maximum likelihood algorithm in `R/qtl`.

The `SynOpDH_Import_GeneticMap.csv` data set contains SynOpDH genotype
data and genetic positions of markers on all chromosomes, and an index
phenotype (TGW BLUP). The `SynOpDH_Import_GeneticMap.csv` marker order
has not been constructed using maximum likelihood, and markers are
spaced every cM. “A” genotype represent W7984 variants of genetic
markers, “B” genotype represent Opata variants of genetic markers.

**Load data**

``` r
SynOpDH <- read.cross("csvr",  file = "https://raw.githubusercontent.com/etaagen/Taagen_2021_TPG/main/supplementary_1/file_S1.2.csv", map.function = "kosambi", na.strings=c("-","NA"))
```

    ##  --Read the following data:
    ##   162  individuals
    ##   1551  markers
    ##   1  phenotypes
    ##  --Cross type: f2

``` r
SynOpDH<- convert2riself(SynOpDH) #changes cross type from f2 to inbred population

summary(SynOpDH)
```

    ##     RI strains via selfing
    ## 
    ##     No. individuals:    162 
    ## 
    ##     No. phenotypes:     1 
    ##     Percent phenotyped: 92.6 
    ## 
    ##     No. chromosomes:    21 
    ##         Autosomes:      1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 
    ##                         7B 7D 
    ## 
    ##     Total markers:      1551 
    ##     No. markers:        51 85 33 76 99 36 87 108 63 78 60 26 63 141 90 66 87 57 
    ##                         97 88 60 
    ##     Percent genotyped:  99.4 
    ##     Genotypes (%):      AA:49.6  BB:50.4

**Maximum likelihood genetic map**  
The import file contains genetic markers evenly spaced every cM. The
code below compares and selects the best marker order by maximum
likelihood analysis.

*Note: there is a small level of randomness when calculating marker
positions, even when receiving identical inputs and parameters*  
*Eval set to FALSE because computationally intensive, please set to TRUE
to validate results.*

``` r
est.map(SynOpDH, map.function = "kosambi") #restimate cM 
#1A
SynOpDH_ripple <- ripple(SynOpDH, chr = "1A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple) 
SynOpDH <- switch.order(SynOpDH, chr = "1A", SynOpDH_ripple[1,])
#1B
SynOpDH_ripple <- ripple(SynOpDH, chr = "1B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "1B", SynOpDH_ripple[1,])
#1D
SynOpDH_ripple <- ripple(SynOpDH, chr = "1D", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "1D", SynOpDH_ripple[1,])
#2A
SynOpDH_ripple <- ripple(SynOpDH, chr = "2A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "2A", SynOpDH_ripple[1,])
#2B
SynOpDH_ripple <- ripple(SynOpDH, chr = "2B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "2B", SynOpDH_ripple[1,])
#2D
SynOpDH_ripple <- ripple(SynOpDH, chr = "2D", window = 2, method = "likelihood", map.function = "kosambi", maxit=10000)
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "2D", SynOpDH_ripple[1,])
#3A
SynOpDH_ripple <- ripple(SynOpDH, chr = "3A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "3A", SynOpDH_ripple[1,])
#3B
SynOpDH_ripple <- ripple(SynOpDH, chr = "3B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "3B", SynOpDH_ripple[1,])
#3D
SynOpDH_ripple <- ripple(SynOpDH, chr = "3D", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "3D", SynOpDH_ripple[1,])
#4A
SynOpDH_ripple <- ripple(SynOpDH, chr = "4A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "4A", SynOpDH_ripple[1,])
#4B
SynOpDH_ripple <- ripple(SynOpDH, chr = "4B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "4B", SynOpDH_ripple[1,])
#4D
SynOpDH_ripple <- ripple(SynOpDH, chr = "4D", window = 2, method = "likelihood", map.function = "kosambi", maxit = 10000)
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "4D", SynOpDH_ripple[1,])
#5A
SynOpDH_ripple <- ripple(SynOpDH, chr = "5A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "5A", SynOpDH_ripple[1,])
#5B
SynOpDH_ripple <- ripple(SynOpDH, chr = "5B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "5B", SynOpDH_ripple[1,])
#5D
SynOpDH_ripple <- ripple(SynOpDH, chr = "5D", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "5D", SynOpDH_ripple[1,])
#6A
SynOpDH_ripple <- ripple(SynOpDH, chr = "6A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "6A", SynOpDH_ripple[1,])
#6B
SynOpDH_ripple <- ripple(SynOpDH, chr = "6B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "6B", SynOpDH_ripple[1,])
#6D
SynOpDH_ripple <- ripple(SynOpDH, chr = "6D", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "6D", SynOpDH_ripple[1,])
#7A
SynOpDH_ripple <- ripple(SynOpDH, chr = "7A", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "7A", SynOpDH_ripple[1,])
#7B
SynOpDH_ripple <- ripple(SynOpDH, chr = "7B", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "7B", SynOpDH_ripple[1,])
#7D
SynOpDH_ripple <- ripple(SynOpDH, chr = "7D", window = 2, method = "likelihood", map.function = "kosambi")
summary(SynOpDH_ripple)
SynOpDH <- switch.order(SynOpDH, chr = "7D", SynOpDH_ripple[1,])

plot.map(SynOpDH)
```

**Export map file:**  
*Eval set to False, `file_S1.9.csv` available in repo, and is the marker
order used for QTL mapping, import file `file_ S1.3.csv` and `file_
S1.4.csv`. Please set eval to TRUE to validate results*

``` r
SynOpDHMap<-pull.map(SynOpDH, as.table=T)
write.csv(SynOpDHMap, "file_S1.9.csv")
```

### QTL mapping

In order to associate SynOpDH genotype and phenotype data, QTL mapping
was performed with the R package “qtl”. The `file_ S1.3.csv` data set
contains all SynOpDH phenotype data, genotype data and genetic linkage
map positions (from `file_S1.9.csv`) of markers on all chromosomes. The
markers on chromosome 5A are in order of physical position, and a single
marker representing the chromosome 5AS structural variant has been
included as well (1,552 markers total). “A” genotype represent W7984
variants of genetic markers, “B” genotype represent Opata variants of
genetic markers.

**Load data**

``` r
#SynOpDH 5A physical positions
SynOpDH <- read.cross("csvr", file = "https://raw.githubusercontent.com/etaagen/Taagen_2021_TPG/main/supplementary_1/file_S1.3.csv", map.function = "kosambi", na.strings=c("-","NA"))
```

    ##  --Read the following data:
    ##   162  individuals
    ##   1552  markers
    ##   26  phenotypes
    ##  --Cross type: f2

``` r
#SynOpDH 5A cM
SynOpDH_cM <- read.cross("csvr", file = "https://raw.githubusercontent.com/etaagen/Taagen_2021_TPG/main/supplementary_1/file_S1.4.csv", map.function = "kosambi", na.strings=c("-","NA"))
```

    ##  --Read the following data:
    ##   162  individuals
    ##   1552  markers
    ##   26  phenotypes
    ##  --Cross type: f2

``` r
SynOpDH<-convert2riself(SynOpDH) #changes cross type from f2 to inbred selfing population
SynOpDH <- jittermap(SynOpDH) #some markers at same position 

SynOpDH_cM<-convert2riself(SynOpDH_cM) #changes cross type from f2 to inbred selfing population
SynOpDH_cM <- jittermap(SynOpDH_cM) #some markers at same position 

summary(SynOpDH)
```

    ##     RI strains via selfing
    ## 
    ##     No. individuals:    162 
    ## 
    ##     No. phenotypes:     26 
    ##     Percent phenotyped: 87 91.4 92 92 90.1 92.6 81.5 92 91.4 92 92 90.1 95.7 
    ##                         81.5 92 91.4 92 92 90.1 92.6 92.6 92.6 92.6 91.4 94.4 
    ##                         93.8 
    ## 
    ##     No. chromosomes:    21 
    ##         Autosomes:      1A 1B 1D 2A 2B 2D 3A 3B 3D 4A 4B 4D 5A 5B 5D 6A 6B 6D 7A 
    ##                         7B 7D 
    ## 
    ##     Total markers:      1552 
    ##     No. markers:        51 85 33 76 99 36 87 108 63 78 60 26 64 141 90 66 87 57 
    ##                         97 88 60 
    ##     Percent genotyped:  99.4 
    ##     Genotypes (%):      AA:49.6  BB:50.4

``` r
SynOpDH_map <- pull.map(SynOpDH, as.table=T)
```

**Scan for QTL**

``` r
#impute missing genotype
SynOpDH <- fill.geno(SynOpDH, method = "imp")
SynOpDH_cM <- fill.geno(SynOpDH_cM, method = "imp")

SynOpData <- calc.genoprob(SynOpDH, step = 0, map.function = "kosambi") 
SynOpData_cM <- calc.genoprob(SynOpDH_cM, step = 0, map.function = "kosambi") 

#warning, scanone will drop individuals with missing phenotypes
TGW_BLUP <- scanone(SynOpData, pheno.col = "BLUP_TGW", method = "hk")
GL_BLUP <- scanone(SynOpData, pheno.col = "BLUP_GL", method = "hk")
GW_BLUP <- scanone(SynOpData, pheno.col = "BLUP_GW", method = "hk") 
HD_BLUP <- scanone(SynOpData, pheno.col = "BLUP_HD", method = "hk") 

summary(TGW_perm <- scanone(SynOpData, pheno.col = "BLUP_TGW", method="hk", n.perm=1000))#5% 3.14
```

    ## Doing permutation in batch mode ...

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 5%  3.22
    ## 10% 2.89

``` r
summary(GL_perm <- scanone(SynOpData, pheno.col = "BLUP_GL", method="hk", n.perm=1000))#5% 3.19
```

    ## Doing permutation in batch mode ...

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 5%  3.13
    ## 10% 2.83

``` r
summary(GW_perm <- scanone(SynOpData, pheno.col = "BLUP_GW", method="hk", n.perm=1000))#5% 3.10
```

    ## Doing permutation in batch mode ...

    ## LOD thresholds (1000 permutations)
    ##     lod
    ## 5%  3.1
    ## 10% 2.8

``` r
summary(HD_perm <- scanone(SynOpData, pheno.col = "BLUP_HD", method="hk", n.perm=1000))#5% 3.07
```

    ## Doing permutation in batch mode ...

    ## LOD thresholds (1000 permutations)
    ##      lod
    ## 5%  3.21
    ## 10% 2.84

``` r
LOD_BLUP <- cbind(TGW_BLUP, GL_BLUP, GW_BLUP, HD_BLUP) 

#check cM 5A
TGW_BLUP_cM <- scanone(SynOpData_cM, pheno.col = "BLUP_TGW", method = "hk")
GL_BLUP_cM <- scanone(SynOpData_cM, pheno.col = "BLUP_GL", method = "hk")
GW_BLUP_cM <- scanone(SynOpData_cM, pheno.col = "BLUP_GW", method = "hk") 
HD_BLUP_cM <- scanone(SynOpData_cM, pheno.col = "BLUP_HD", method = "hk") 
```

### QTL Summary

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#bqtkdzqdzc .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 3px;
  border-bottom-color: white;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#bqtkdzqdzc .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bqtkdzqdzc .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#bqtkdzqdzc .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#bqtkdzqdzc .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bqtkdzqdzc .gt_col_headings {
  border-top-style: solid;
  border-top-width: 3px;
  border-top-color: white;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: black;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bqtkdzqdzc .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#bqtkdzqdzc .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#bqtkdzqdzc .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bqtkdzqdzc .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bqtkdzqdzc .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: black;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#bqtkdzqdzc .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#bqtkdzqdzc .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#bqtkdzqdzc .gt_from_md > :first-child {
  margin-top: 0;
}

#bqtkdzqdzc .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bqtkdzqdzc .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: white;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#bqtkdzqdzc .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#bqtkdzqdzc .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bqtkdzqdzc .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#bqtkdzqdzc .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bqtkdzqdzc .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bqtkdzqdzc .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bqtkdzqdzc .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bqtkdzqdzc .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bqtkdzqdzc .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#bqtkdzqdzc .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bqtkdzqdzc .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#bqtkdzqdzc .gt_left {
  text-align: left;
}

#bqtkdzqdzc .gt_center {
  text-align: center;
}

#bqtkdzqdzc .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bqtkdzqdzc .gt_font_normal {
  font-weight: normal;
}

#bqtkdzqdzc .gt_font_bold {
  font-weight: bold;
}

#bqtkdzqdzc .gt_font_italic {
  font-style: italic;
}

#bqtkdzqdzc .gt_super {
  font-size: 65%;
}

#bqtkdzqdzc .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="bqtkdzqdzc" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1">

</th>

<th class="gt_col_heading gt_center gt_columns_bottom_border" rowspan="2" colspan="1" style="font-weight: bold;">

position<sup class="gt_footnote_marks">1</sup>

</th>

<th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="4">

<span class="gt_column_spanner">SynOpDH LOD score</span>

</th>

</tr>

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

TGW

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

GL

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

GW

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

HD

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr class="gt_group_heading_row">

<td colspan="6" class="gt_group_heading">

2D

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS248

</td>

<td class="gt_row gt_right">

6.82

</td>

<td class="gt_row gt_right">

1.55

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.47

</td>

<td class="gt_row gt_right">

0.27

</td>

<td class="gt_row gt_right">

0.19

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS5

</td>

<td class="gt_row gt_right">

6.82

</td>

<td class="gt_row gt_right">

1.55

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.47

</td>

<td class="gt_row gt_right">

0.27

</td>

<td class="gt_row gt_right">

0.19

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS972

</td>

<td class="gt_row gt_right">

7.45

</td>

<td class="gt_row gt_right">

1.55

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.18

</td>

<td class="gt_row gt_right">

0.21

</td>

<td class="gt_row gt_right">

0.08

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS973

</td>

<td class="gt_row gt_right">

7.45

</td>

<td class="gt_row gt_right">

1.55

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.18

</td>

<td class="gt_row gt_right">

0.21

</td>

<td class="gt_row gt_right">

0.08

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS221

</td>

<td class="gt_row gt_right">

8.08

</td>

<td class="gt_row gt_right">

1.85

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.47

</td>

<td class="gt_row gt_right">

0.14

</td>

<td class="gt_row gt_right">

0.06

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS941

</td>

<td class="gt_row gt_right">

8.08

</td>

<td class="gt_row gt_right">

1.85

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.47

</td>

<td class="gt_row gt_right">

0.14

</td>

<td class="gt_row gt_right">

0.06

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS110

</td>

<td class="gt_row gt_right">

9.03

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.15

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1128

</td>

<td class="gt_row gt_right">

9.03

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.15

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1146

</td>

<td class="gt_row gt_right">

9.03

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.15

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS343

</td>

<td class="gt_row gt_right">

9.03

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.15

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS666

</td>

<td class="gt_row gt_right">

9.03

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.15

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS560

</td>

<td class="gt_row gt_right">

9.34

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.85

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS58

</td>

<td class="gt_row gt_right">

9.34

</td>

<td class="gt_row gt_right">

1.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.85

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.10

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS618

</td>

<td class="gt_row gt_right">

9.65

</td>

<td class="gt_row gt_right">

1.83

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.38

</td>

<td class="gt_row gt_right">

0.07

</td>

<td class="gt_row gt_right">

0.17

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS632

</td>

<td class="gt_row gt_right">

9.65

</td>

<td class="gt_row gt_right">

1.83

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.38

</td>

<td class="gt_row gt_right">

0.07

</td>

<td class="gt_row gt_right">

0.17

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS579

</td>

<td class="gt_row gt_right">

19.82

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.70

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.16

</td>

<td class="gt_row gt_right">

0.00

</td>

<td class="gt_row gt_right">

0.68

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1212

</td>

<td class="gt_row gt_right">

31.19

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.22

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.02

</td>

<td class="gt_row gt_right">

0.10

</td>

<td class="gt_row gt_right">

1.69

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1310

</td>

<td class="gt_row gt_right">

31.19

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.22

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.02

</td>

<td class="gt_row gt_right">

0.10

</td>

<td class="gt_row gt_right">

1.69

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS745

</td>

<td class="gt_row gt_right">

31.58

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.83

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.99

</td>

<td class="gt_row gt_right">

0.09

</td>

<td class="gt_row gt_right">

1.83

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS250

</td>

<td class="gt_row gt_right">

32.93

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.31

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.00

</td>

<td class="gt_row gt_right">

0.14

</td>

<td class="gt_row gt_right">

2.18

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS421

</td>

<td class="gt_row gt_right">

32.93

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.31

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.00

</td>

<td class="gt_row gt_right">

0.14

</td>

<td class="gt_row gt_right">

2.18

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS216

</td>

<td class="gt_row gt_right">

39.61

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.52

</td>

<td class="gt_row gt_right">

0.27

</td>

<td class="gt_row gt_right">

1.23

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1434

</td>

<td class="gt_row gt_right">

39.92

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.52

</td>

<td class="gt_row gt_right">

0.27

</td>

<td class="gt_row gt_right">

0.95

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1288

</td>

<td class="gt_row gt_right">

39.92

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.52

</td>

<td class="gt_row gt_right">

0.27

</td>

<td class="gt_row gt_right">

0.95

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1436

</td>

<td class="gt_row gt_right">

39.92

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.52

</td>

<td class="gt_row gt_right">

0.27

</td>

<td class="gt_row gt_right">

0.95

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS757

</td>

<td class="gt_row gt_right">

40.23

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.61

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.70

</td>

<td class="gt_row gt_right">

0.25

</td>

<td class="gt_row gt_right">

1.03

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS156

</td>

<td class="gt_row gt_right">

40.23

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.61

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.70

</td>

<td class="gt_row gt_right">

0.25

</td>

<td class="gt_row gt_right">

1.03

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS171

</td>

<td class="gt_row gt_right">

40.23

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.61

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.70

</td>

<td class="gt_row gt_right">

0.25

</td>

<td class="gt_row gt_right">

1.03

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS756

</td>

<td class="gt_row gt_right">

40.23

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.61

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.70

</td>

<td class="gt_row gt_right">

0.25

</td>

<td class="gt_row gt_right">

1.03

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1286

</td>

<td class="gt_row gt_right">

41.48

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.72

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.09

</td>

<td class="gt_row gt_right">

0.23

</td>

<td class="gt_row gt_right">

0.45

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS505

</td>

<td class="gt_row gt_right">

42.48

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.91

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.29

</td>

<td class="gt_row gt_right">

0.58

</td>

<td class="gt_row gt_right">

0.63

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS929

</td>

<td class="gt_row gt_right">

43.83

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.81

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.03

</td>

<td class="gt_row gt_right">

0.66

</td>

<td class="gt_row gt_right">

0.36

</td>

</tr>

<tr class="gt_group_heading_row">

<td colspan="6" class="gt_group_heading">

5A

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

chr5AS\_status

</td>

<td class="gt_row gt_right">

0.00

</td>

<td class="gt_row gt_right">

2.61

</td>

<td class="gt_row gt_right">

0.05

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

13.64

</td>

<td class="gt_row gt_right">

1.55

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

wmc705

</td>

<td class="gt_row gt_right">

290.04

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.61

</td>

<td class="gt_row gt_right">

0.22

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

15.92

</td>

<td class="gt_row gt_right">

1.88

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

wmc805aD

</td>

<td class="gt_row gt_right">

364.42

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.35

</td>

<td class="gt_row gt_right">

0.29

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

15.32

</td>

<td class="gt_row gt_right">

1.28

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS429

</td>

<td class="gt_row gt_right">

380.82

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.62

</td>

<td class="gt_row gt_right">

0.23

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

16.81

</td>

<td class="gt_row gt_right">

1.36

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS147

</td>

<td class="gt_row gt_right">

436.89

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.49

</td>

<td class="gt_row gt_right">

0.45

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

15.35

</td>

<td class="gt_row gt_right">

1.30

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS902

</td>

<td class="gt_row gt_right">

463.77

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS440

</td>

<td class="gt_row gt_right">

467.37

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS70

</td>

<td class="gt_row gt_right">

470.19

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS422

</td>

<td class="gt_row gt_right">

471.00

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1198

</td>

<td class="gt_row gt_right">

473.01

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

wPt-0958

</td>

<td class="gt_row gt_right">

473.24

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS380

</td>

<td class="gt_row gt_right">

474.15

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS57

</td>

<td class="gt_row gt_right">

476.90

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.82

</td>

<td class="gt_row gt_right">

0.64

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

12.91

</td>

<td class="gt_row gt_right">

1.04

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS687

</td>

<td class="gt_row gt_right">

480.51

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.49

</td>

<td class="gt_row gt_right">

0.96

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

10.68

</td>

<td class="gt_row gt_right">

1.15

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS284

</td>

<td class="gt_row gt_right">

487.76

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.49

</td>

<td class="gt_row gt_right">

0.96

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

10.68

</td>

<td class="gt_row gt_right">

1.15

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS473

</td>

<td class="gt_row gt_right">

496.49

</td>

<td class="gt_row gt_right">

2.35

</td>

<td class="gt_row gt_right">

0.70

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.32

</td>

<td class="gt_row gt_right">

1.32

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS822

</td>

<td class="gt_row gt_right">

496.49

</td>

<td class="gt_row gt_right">

2.35

</td>

<td class="gt_row gt_right">

0.70

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.32

</td>

<td class="gt_row gt_right">

1.32

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS823

</td>

<td class="gt_row gt_right">

496.49

</td>

<td class="gt_row gt_right">

2.35

</td>

<td class="gt_row gt_right">

0.70

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.32

</td>

<td class="gt_row gt_right">

1.32

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1057

</td>

<td class="gt_row gt_right">

512.60

</td>

<td class="gt_row gt_right">

2.55

</td>

<td class="gt_row gt_right">

0.69

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.88

</td>

<td class="gt_row gt_right">

1.51

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

wmc415aD

</td>

<td class="gt_row gt_right">

535.15

</td>

<td class="gt_row gt_right">

1.75

</td>

<td class="gt_row gt_right">

0.50

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

8.37

</td>

<td class="gt_row gt_right">

1.70

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS836

</td>

<td class="gt_row gt_right">

548.99

</td>

<td class="gt_row gt_right">

0.54

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.35

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.41

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS837

</td>

<td class="gt_row gt_right">

548.99

</td>

<td class="gt_row gt_right">

0.54

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.35

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.41

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1419

</td>

<td class="gt_row gt_right">

549.41

</td>

<td class="gt_row gt_right">

0.69

</td>

<td class="gt_row gt_right">

0.03

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.42

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS578

</td>

<td class="gt_row gt_right">

550.10

</td>

<td class="gt_row gt_right">

0.69

</td>

<td class="gt_row gt_right">

0.03

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.42

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS285

</td>

<td class="gt_row gt_right">

553.34

</td>

<td class="gt_row gt_right">

0.69

</td>

<td class="gt_row gt_right">

0.03

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.42

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS463

</td>

<td class="gt_row gt_right">

554.08

</td>

<td class="gt_row gt_right">

0.69

</td>

<td class="gt_row gt_right">

0.00

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.88

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.42

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1347

</td>

<td class="gt_row gt_right">

554.98

</td>

<td class="gt_row gt_right">

0.55

</td>

<td class="gt_row gt_right">

0.03

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.52

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.81

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS998

</td>

<td class="gt_row gt_right">

555.04

</td>

<td class="gt_row gt_right">

0.55

</td>

<td class="gt_row gt_right">

0.03

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.52

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.81

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1143

</td>

<td class="gt_row gt_right">

561.12

</td>

<td class="gt_row gt_right">

0.23

</td>

<td class="gt_row gt_right">

0.19

</td>

<td class="gt_row gt_right">

2.63

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.50

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS984

</td>

<td class="gt_row gt_right">

562.79

</td>

<td class="gt_row gt_right">

0.19

</td>

<td class="gt_row gt_right">

0.22

</td>

<td class="gt_row gt_right">

2.58

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.62

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS354

</td>

<td class="gt_row gt_right">

564.56

</td>

<td class="gt_row gt_right">

0.19

</td>

<td class="gt_row gt_right">

0.22

</td>

<td class="gt_row gt_right">

2.58

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.62

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1174

</td>

<td class="gt_row gt_right">

565.40

</td>

<td class="gt_row gt_right">

0.19

</td>

<td class="gt_row gt_right">

0.22

</td>

<td class="gt_row gt_right">

2.58

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.62

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1176

</td>

<td class="gt_row gt_right">

565.40

</td>

<td class="gt_row gt_right">

0.19

</td>

<td class="gt_row gt_right">

0.22

</td>

<td class="gt_row gt_right">

2.58

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.62

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS707

</td>

<td class="gt_row gt_right">

565.61

</td>

<td class="gt_row gt_right">

0.19

</td>

<td class="gt_row gt_right">

0.22

</td>

<td class="gt_row gt_right">

2.58

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.62

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1046

</td>

<td class="gt_row gt_right">

567.50

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_right">

0.56

</td>

<td class="gt_row gt_right">

2.48

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

4.68

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1008

</td>

<td class="gt_row gt_right">

576.11

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.50

</td>

<td class="gt_row gt_right">

1.12

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.67

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS294

</td>

<td class="gt_row gt_right">

577.80

</td>

<td class="gt_row gt_right">

0.04

</td>

<td class="gt_row gt_right">

0.50

</td>

<td class="gt_row gt_right">

1.12

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.67

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS746

</td>

<td class="gt_row gt_right">

582.09

</td>

<td class="gt_row gt_right">

0.03

</td>

<td class="gt_row gt_right">

0.52

</td>

<td class="gt_row gt_right">

1.01

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

8.40

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1405

</td>

<td class="gt_row gt_right">

585.20

</td>

<td class="gt_row gt_right">

0.01

</td>

<td class="gt_row gt_right">

0.61

</td>

<td class="gt_row gt_right">

0.89

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.91

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS539

</td>

<td class="gt_row gt_right">

585.22

</td>

<td class="gt_row gt_right">

0.01

</td>

<td class="gt_row gt_right">

0.61

</td>

<td class="gt_row gt_right">

0.89

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.91

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS540

</td>

<td class="gt_row gt_right">

585.22

</td>

<td class="gt_row gt_right">

0.01

</td>

<td class="gt_row gt_right">

0.61

</td>

<td class="gt_row gt_right">

0.89

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.91

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS541

</td>

<td class="gt_row gt_right">

585.22

</td>

<td class="gt_row gt_right">

0.01

</td>

<td class="gt_row gt_right">

0.61

</td>

<td class="gt_row gt_right">

0.89

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

9.91

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS293

</td>

<td class="gt_row gt_right">

605.25

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_right">

0.36

</td>

<td class="gt_row gt_right">

0.56

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.28

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

cfa2141aD

</td>

<td class="gt_row gt_right">

620.17

</td>

<td class="gt_row gt_right">

0.06

</td>

<td class="gt_row gt_right">

0.36

</td>

<td class="gt_row gt_right">

0.56

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

7.28

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

barc232b

</td>

<td class="gt_row gt_right">

623.48

</td>

<td class="gt_row gt_right">

0.11

</td>

<td class="gt_row gt_right">

0.40

</td>

<td class="gt_row gt_right">

0.45

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

6.94

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1249

</td>

<td class="gt_row gt_right">

624.21

</td>

<td class="gt_row gt_right">

0.08

</td>

<td class="gt_row gt_right">

0.52

</td>

<td class="gt_row gt_right">

0.44

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.86

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS532

</td>

<td class="gt_row gt_right">

626.19

</td>

<td class="gt_row gt_right">

0.08

</td>

<td class="gt_row gt_right">

0.52

</td>

<td class="gt_row gt_right">

0.44

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.86

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS499

</td>

<td class="gt_row gt_right">

631.68

</td>

<td class="gt_row gt_right">

0.08

</td>

<td class="gt_row gt_right">

0.52

</td>

<td class="gt_row gt_right">

0.44

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

5.86

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS826

</td>

<td class="gt_row gt_right">

640.68

</td>

<td class="gt_row gt_right">

0.01

</td>

<td class="gt_row gt_right">

1.25

</td>

<td class="gt_row gt_right">

0.24

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.39

</td>

</tr>

<tr class="gt_group_heading_row">

<td colspan="6" class="gt_group_heading">

6A

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS4

</td>

<td class="gt_row gt_right">

35.42

</td>

<td class="gt_row gt_right">

2.03

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.12

</td>

<td class="gt_row gt_right">

0.53

</td>

<td class="gt_row gt_right">

0.18

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS1114

</td>

<td class="gt_row gt_right">

41.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.18

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.14

</td>

<td class="gt_row gt_right">

0.99

</td>

<td class="gt_row gt_right">

0.12

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS42

</td>

<td class="gt_row gt_right">

41.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.18

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.14

</td>

<td class="gt_row gt_right">

0.99

</td>

<td class="gt_row gt_right">

0.12

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

synopGBS811

</td>

<td class="gt_row gt_right">

41.68

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.18

</td>

<td class="gt_row gt_right" style="background-color: rgba(0,0,255,0.7); color: white; font-weight: bold;">

3.14

</td>

<td class="gt_row gt_right">

0.99

</td>

<td class="gt_row gt_right">

0.12

</td>

</tr>

</tbody>

<tfoot>

<tr class="gt_footnotes">

<td colspan="6">

<p class="gt_footnote">

<sup class="gt_footnote_marks"> <em>1</em> </sup>

Mbp for chr 5A <br />

</p>

</td>

</tr>

</tfoot>

</table>

</div>

<!--/html_preserve-->

### Phenotype variation explained by QTL

`% var = 1 - 10^(-2 * LOD / n)` <!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#rkrhyifmra .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rkrhyifmra .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rkrhyifmra .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rkrhyifmra .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rkrhyifmra .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rkrhyifmra .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rkrhyifmra .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rkrhyifmra .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rkrhyifmra .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rkrhyifmra .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rkrhyifmra .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rkrhyifmra .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#rkrhyifmra .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rkrhyifmra .gt_from_md > :first-child {
  margin-top: 0;
}

#rkrhyifmra .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rkrhyifmra .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rkrhyifmra .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#rkrhyifmra .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rkrhyifmra .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#rkrhyifmra .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rkrhyifmra .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rkrhyifmra .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rkrhyifmra .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rkrhyifmra .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rkrhyifmra .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#rkrhyifmra .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rkrhyifmra .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#rkrhyifmra .gt_left {
  text-align: left;
}

#rkrhyifmra .gt_center {
  text-align: center;
}

#rkrhyifmra .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rkrhyifmra .gt_font_normal {
  font-weight: normal;
}

#rkrhyifmra .gt_font_bold {
  font-weight: bold;
}

#rkrhyifmra .gt_font_italic {
  font-style: italic;
}

#rkrhyifmra .gt_super {
  font-size: 65%;
}

#rkrhyifmra .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="rkrhyifmra" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

QTL

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1">

% variation

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

QTL2D\_TGW

</td>

<td class="gt_row gt_right">

13.790415

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTL2D\_GL

</td>

<td class="gt_row gt_right">

19.658808

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTL5A\_TGW

</td>

<td class="gt_row gt_right">

10.392228

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTL5A\_GW

</td>

<td class="gt_row gt_right">

37.988943

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTL6A\_TGW

</td>

<td class="gt_row gt_right">

8.643229

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTL6A\_GL

</td>

<td class="gt_row gt_right">

8.539290

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

### QTL plots

![](script_S1_files/figure-gfm/QTL%20plot-1.png)<!-- -->![](script_S1_files/figure-gfm/QTL%20plot-2.png)<!-- -->![](script_S1_files/figure-gfm/QTL%20plot-3.png)<!-- -->

### Figure 1B

![](script_S1_files/figure-gfm/figure%201B-1.png)<!-- -->

### QTL physical postions

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#fyngjabqex .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#fyngjabqex .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#fyngjabqex .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#fyngjabqex .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#fyngjabqex .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fyngjabqex .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#fyngjabqex .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#fyngjabqex .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#fyngjabqex .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#fyngjabqex .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#fyngjabqex .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#fyngjabqex .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#fyngjabqex .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#fyngjabqex .gt_from_md > :first-child {
  margin-top: 0;
}

#fyngjabqex .gt_from_md > :last-child {
  margin-bottom: 0;
}

#fyngjabqex .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#fyngjabqex .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#fyngjabqex .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#fyngjabqex .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#fyngjabqex .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#fyngjabqex .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#fyngjabqex .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#fyngjabqex .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fyngjabqex .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#fyngjabqex .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#fyngjabqex .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#fyngjabqex .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#fyngjabqex .gt_left {
  text-align: left;
}

#fyngjabqex .gt_center {
  text-align: center;
}

#fyngjabqex .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#fyngjabqex .gt_font_normal {
  font-weight: normal;
}

#fyngjabqex .gt_font_bold {
  font-weight: bold;
}

#fyngjabqex .gt_font_italic {
  font-style: italic;
}

#fyngjabqex .gt_super {
  font-size: 65%;
}

#fyngjabqex .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="fyngjabqex" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

QTL<sup class="gt_footnote_marks">1</sup>

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Flanking markers

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

cM

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

RefSeqV1.0 flanking

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Peak marker

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

RefSeqV1.0 peak

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

QTgw.cnl-2D

</td>

<td class="gt_row gt_left">

synopGBS929 - synopGBS579

</td>

<td class="gt_row gt_left">

19.82 - 43.83

</td>

<td class="gt_row gt_left">

10,718,874 - 56,072,170

</td>

<td class="gt_row gt_left">

synopGBS1212

</td>

<td class="gt_row gt_left">

35,003,633

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTgw.cnl-5A

</td>

<td class="gt_row gt_left">

wmc705 - synopGBS284

</td>

<td class="gt_row gt_left">

0.1 - 4.49

</td>

<td class="gt_row gt_left">

290,035,203 - 487,755,214

</td>

<td class="gt_row gt_left">

synopGBS429

</td>

<td class="gt_row gt_left">

380,823,821

</td>

</tr>

<tr>

<td class="gt_row gt_left">

QTgw.cnl-6A

</td>

<td class="gt_row gt_left">

synopGBS42 - synopGBS1114

</td>

<td class="gt_row gt_left">

41.68 - 41.68

</td>

<td class="gt_row gt_left">

314,315,653 - 440,255,621

</td>

<td class="gt_row gt_left">

synopGBS811

</td>

<td class="gt_row gt_left">

358,607,834

</td>

</tr>

</tbody>

<tfoot>

<tr class="gt_footnotes">

<td colspan="6">

<p class="gt_footnote">

<sup class="gt_footnote_marks"> <em>1</em> </sup>

chromosome 2D and chromosome 6A markers reordered to match physical map
order <br />

</p>

</td>

</tr>

</tfoot>

</table>

</div>

<!--/html_preserve-->

### Test for interaction between 2D, 5A and 6A TGW QTL

    ## Analysis of Deviance Table (Type II Wald chisquare tests)
    ## 
    ## Response: TGW
    ##                                       Chisq Df Pr(>Chisq)    
    ## QTL2D_synopGBS1212                  29.1425  1  6.725e-08 ***
    ## QTL5A_wmc705                        28.0831  1  1.162e-07 ***
    ## QTL6A_synopGBS42                    20.9094  1  4.815e-06 ***
    ## QTL2D_synopGBS1212:QTL5A_wmc705      0.1849  1    0.66722    
    ## QTL5A_wmc705:QTL6A_synopGBS42        3.2000  1    0.07364 .  
    ## QTL2D_synopGBS1212:QTL6A_synopGBS42  0.0159  1    0.89953    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### Phenotypes across years, Table 1

``` r
cald16 <- SynOpDH_Phenotypes %>% 
  filter(Location == "Caldwell" & Year == "2016") 
summary(aov(TGW ~ QTL5A_wmc705, data = cald16))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2  567.1  283.55   19.72 2.93e-08 ***
    ## Residuals    138 1984.3   14.38                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 21 observations deleted due to missingness

``` r
summary(aov(GL ~ QTL5A_wmc705, data = cald16))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## QTL5A_wmc705   1  0.172  0.1715   1.214  0.273
    ## Residuals    130 18.371  0.1413               
    ## 30 observations deleted due to missingness

``` r
summary(aov(GW ~ QTL5A_wmc705, data = cald16))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   1  1.154  1.1536   49.11 1.18e-10 ***
    ## Residuals    130  3.054  0.0235                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 30 observations deleted due to missingness

``` r
cald17 <- SynOpDH_Phenotypes %>% 
  filter(Location == "Caldwell" & Year == "2017") 
summary(aov(GL ~ QTL5A_wmc705, data = cald17))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## QTL5A_wmc705   2  0.159 0.07956   0.402   0.67
    ## Residuals    146 28.873 0.19776               
    ## 13 observations deleted due to missingness

``` r
summary(aov(GW ~ QTL5A_wmc705, data = cald17))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2  0.963  0.4817   10.22 7.03e-05 ***
    ## Residuals    146  6.883  0.0471                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 13 observations deleted due to missingness

``` r
cald18A <- SynOpDH_Phenotypes %>% 
  filter(Location == "Caldwell" & Year == "2018" & Rep == "A") 
summary(aov(TGW ~ QTL5A_wmc705, data = cald18A))
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## QTL5A_wmc705   2  176.8   88.42   5.287 0.00607 **
    ## Residuals    145 2424.7   16.72                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 14 observations deleted due to missingness

``` r
summary(aov(GL ~ QTL5A_wmc705, data = cald18A))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## QTL5A_wmc705   2  0.084 0.04191   0.292  0.747
    ## Residuals    145 20.842 0.14374               
    ## 14 observations deleted due to missingness

``` r
summary(aov(GW ~ QTL5A_wmc705, data = cald18A))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2  0.543 0.27152   24.05 9.53e-10 ***
    ## Residuals    145  1.637 0.01129                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 14 observations deleted due to missingness

``` r
summary(aov(HD ~ QTL5A_wmc705, data = cald18A))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)  
    ## QTL5A_wmc705   2  122.2   61.10   3.277 0.0405 *
    ## Residuals    147 2740.4   18.64                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 12 observations deleted due to missingness

``` r
cald18B <- SynOpDH_Phenotypes %>% 
  filter(Location == "Caldwell" & Year == "2018" & Rep == "B") 
summary(aov(TGW ~ QTL5A_wmc705, data = cald18B))
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## QTL5A_wmc705   2  226.5  113.25   6.417 0.00213 **
    ## Residuals    146 2576.7   17.65                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 13 observations deleted due to missingness

``` r
summary(aov(GL ~ QTL5A_wmc705, data = cald18B))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## QTL5A_wmc705   2   0.17 0.08502    0.59  0.555
    ## Residuals    146  21.02 0.14401               
    ## 13 observations deleted due to missingness

``` r
summary(aov(GW ~ QTL5A_wmc705, data = cald18B))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2 0.6795  0.3397   32.21 2.58e-12 ***
    ## Residuals    146 1.5398  0.0105                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 13 observations deleted due to missingness

``` r
summary(aov(HD ~ QTL5A_wmc705, data = cald18B))
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## QTL5A_wmc705   2  187.3   93.65   5.223 0.00644 **
    ## Residuals    147 2636.0   17.93                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 12 observations deleted due to missingness

``` r
helf18A <- SynOpDH_Phenotypes %>% 
  filter(Location == "Helfer" & Year == "2018" & Rep == "A") 
summary(aov(TGW ~ QTL5A_wmc705, data = helf18A))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2  315.5  157.74   10.17 7.33e-05 ***
    ## Residuals    146 2264.6   15.51                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 13 observations deleted due to missingness

``` r
summary(aov(GL ~ QTL5A_wmc705, data = helf18A))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## QTL5A_wmc705   2   0.63  0.3150   2.244   0.11
    ## Residuals    146  20.50  0.1404               
    ## 13 observations deleted due to missingness

``` r
summary(aov(GW ~ QTL5A_wmc705, data = helf18A))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2  1.128  0.5639   44.43 8.47e-16 ***
    ## Residuals    146  1.853  0.0127                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 13 observations deleted due to missingness

``` r
summary(aov(HD ~ QTL5A_wmc705, data = helf18A))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)  
    ## QTL5A_wmc705   2  159.1   79.54   4.221 0.0165 *
    ## Residuals    147 2770.2   18.85                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 12 observations deleted due to missingness

``` r
helf18B <- SynOpDH_Phenotypes %>% 
  filter(Location == "Helfer" & Year == "2018" & Rep == "B") 
summary(aov(TGW ~ QTL5A_wmc705, data = helf18B))
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## QTL5A_wmc705   2  167.7   83.86   6.214 0.00258 **
    ## Residuals    143 1929.9   13.50                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 16 observations deleted due to missingness

``` r
summary(aov(GL ~ QTL5A_wmc705, data = helf18B))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## QTL5A_wmc705   2  0.391  0.1957   1.547  0.216
    ## Residuals    143 18.089  0.1265               
    ## 16 observations deleted due to missingness

``` r
summary(aov(GW ~ QTL5A_wmc705, data = helf18B))
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## QTL5A_wmc705   2 0.7841  0.3921    37.7 7.07e-14 ***
    ## Residuals    143 1.4870  0.0104                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 16 observations deleted due to missingness

``` r
summary(aov(HD ~ QTL5A_wmc705, data = helf18B))
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)  
    ## QTL5A_wmc705   2  132.9   66.46   3.619 0.0293 *
    ## Residuals    145 2662.9   18.36                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 14 observations deleted due to missingness

``` r
#BLUP
summary(aov(BLUP_TGW ~ KASP_341510829, data = chr5Apheno))
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## KASP_341510829   1  181.4  181.40   17.77 4.33e-05 ***
    ## Residuals      148 1511.3   10.21                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 12 observations deleted due to missingness

``` r
summary(aov(BLUP_GL ~ KASP_341510829, data = chr5Apheno))
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)
    ## KASP_341510829   1   0.16 0.15978   1.616  0.206
    ## Residuals      153  15.13 0.09887               
    ## 7 observations deleted due to missingness

``` r
summary(aov(BLUP_GW ~ KASP_341510829, data = chr5Apheno))
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## KASP_341510829   1 0.5982  0.5982   85.56 2.33e-16 ***
    ## Residuals      148 1.0347  0.0070                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 12 observations deleted due to missingness

``` r
summary(aov(BLUP_HD ~ KASP_341510829, data = chr5Apheno))
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## KASP_341510829   1  194.5  194.45   13.46 0.000337 ***
    ## Residuals      151 2181.1   14.44                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 9 observations deleted due to missingness

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#bvdqqwtjqc .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: black;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 3px;
  border-bottom-color: white;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#bvdqqwtjqc .gt_heading {
  background-color: #FFFFFF;
  text-align: left;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bvdqqwtjqc .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#bvdqqwtjqc .gt_subtitle {
  color: #333333;
  font-size: 16px;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#bvdqqwtjqc .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bvdqqwtjqc .gt_col_headings {
  border-top-style: solid;
  border-top-width: 3px;
  border-top-color: white;
  border-bottom-style: solid;
  border-bottom-width: 3px;
  border-bottom-color: black;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bvdqqwtjqc .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#bvdqqwtjqc .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#bvdqqwtjqc .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bvdqqwtjqc .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bvdqqwtjqc .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 3px;
  border-bottom-color: black;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#bvdqqwtjqc .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#bvdqqwtjqc .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#bvdqqwtjqc .gt_from_md > :first-child {
  margin-top: 0;
}

#bvdqqwtjqc .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bvdqqwtjqc .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: white;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#bvdqqwtjqc .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#bvdqqwtjqc .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bvdqqwtjqc .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#bvdqqwtjqc .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bvdqqwtjqc .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bvdqqwtjqc .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bvdqqwtjqc .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bvdqqwtjqc .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bvdqqwtjqc .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#bvdqqwtjqc .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bvdqqwtjqc .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#bvdqqwtjqc .gt_left {
  text-align: left;
}

#bvdqqwtjqc .gt_center {
  text-align: center;
}

#bvdqqwtjqc .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bvdqqwtjqc .gt_font_normal {
  font-weight: normal;
}

#bvdqqwtjqc .gt_font_bold {
  font-weight: bold;
}

#bvdqqwtjqc .gt_font_italic {
  font-style: italic;
}

#bvdqqwtjqc .gt_super {
  font-size: 65%;
}

#bvdqqwtjqc .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="bvdqqwtjqc" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_header">

<tr>

<th colspan="7" class="gt_heading gt_title gt_font_normal" style>

</th>

</tr>

<tr>

<th colspan="7" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style>

<strong>Table 1</strong> Mean thousand grain weight (TGW), grain length
(GL), grain width (GW), & heading date (HD) of SynOpDH entries with
<em>QTgw.cnl-5A+</em> vs <em>QTgw.cnl-5A-</em> alleles

</th>

</tr>

</thead>

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

Year

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

Location

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

Genotype

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

TGW (g)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

GL (mm)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

GW (mm)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

HD (julian)

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

2016

</td>

<td class="gt_row gt_left">

Caldwell

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left">

39.19

</td>

<td class="gt_row gt_left">

7.16

</td>

<td class="gt_row gt_left">

3.26

</td>

<td class="gt_row gt_left">

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

36.49

</td>

<td class="gt_row gt_left">

7.09

</td>

<td class="gt_row gt_left">

3.08

</td>

<td class="gt_row gt_left">

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

7.4% \*\*\*

</td>

<td class="gt_row gt_left">

0.98%

</td>

<td class="gt_row gt_left">

5.7% \*\*\*

</td>

<td class="gt_row gt_left">

</td>

</tr>

<tr>

<td class="gt_row gt_left">

2017

</td>

<td class="gt_row gt_left">

Caldwell

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

6.58

</td>

<td class="gt_row gt_left">

2.84

</td>

<td class="gt_row gt_left">

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

6.53

</td>

<td class="gt_row gt_left">

2.68

</td>

<td class="gt_row gt_left">

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

0.81%

</td>

<td class="gt_row gt_left">

6.1% \*\*\*

</td>

<td class="gt_row gt_left">

</td>

</tr>

<tr>

<td class="gt_row gt_left">

2018

</td>

<td class="gt_row gt_left">

Caldwell (A)

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left">

40.36

</td>

<td class="gt_row gt_left">

7.25

</td>

<td class="gt_row gt_left">

3.26

</td>

<td class="gt_row gt_left">

184.71

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

38.37

</td>

<td class="gt_row gt_left">

7.28

</td>

<td class="gt_row gt_left">

3.14

</td>

<td class="gt_row gt_left">

182.64

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

5.2% \*\*

</td>

<td class="gt_row gt_left">

\-0.39%

</td>

<td class="gt_row gt_left">

3.7% \*\*\*

</td>

<td class="gt_row gt_left">

1.1% \*

</td>

</tr>

<tr>

<td class="gt_row gt_left">

2018

</td>

<td class="gt_row gt_left">

Caldwell (B)

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left">

39.97

</td>

<td class="gt_row gt_left">

7.24

</td>

<td class="gt_row gt_left">

3.26

</td>

<td class="gt_row gt_left">

185.2

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

37.8

</td>

<td class="gt_row gt_left">

7.25

</td>

<td class="gt_row gt_left">

3.13

</td>

<td class="gt_row gt_left">

182.89

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

5.7% \*\*

</td>

<td class="gt_row gt_left">

\-0.21%

</td>

<td class="gt_row gt_left">

4.04% \*\*\*

</td>

<td class="gt_row gt_left">

1.3% \*\*

</td>

</tr>

<tr>

<td class="gt_row gt_left">

2018

</td>

<td class="gt_row gt_left">

Helfer (A)

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left">

31.04

</td>

<td class="gt_row gt_left">

7.08

</td>

<td class="gt_row gt_left">

3.08

</td>

<td class="gt_row gt_left">

187.05

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

28.55

</td>

<td class="gt_row gt_left">

6.97

</td>

<td class="gt_row gt_left">

2.92

</td>

<td class="gt_row gt_left">

184.8

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

8.7% \*\*\*

</td>

<td class="gt_row gt_left">

1.70%

</td>

<td class="gt_row gt_left">

5.4% \*\*\*

</td>

<td class="gt_row gt_left">

1.2% \*

</td>

</tr>

<tr>

<td class="gt_row gt_left">

2018

</td>

<td class="gt_row gt_left">

Helfer (B)

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left">

31.61

</td>

<td class="gt_row gt_left">

7.03

</td>

<td class="gt_row gt_left">

3.09

</td>

<td class="gt_row gt_left">

186.91

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

29.77

</td>

<td class="gt_row gt_left">

6.97

</td>

<td class="gt_row gt_left">

2.96

</td>

<td class="gt_row gt_left">

184.79

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

6.2% \*\*

</td>

<td class="gt_row gt_left">

0.81%

</td>

<td class="gt_row gt_left">

4.5% \*\*\*

</td>

<td class="gt_row gt_left">

1.1% \*

</td>

</tr>

<tr>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

BLUP

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black; font-style: italic;">

QTgw.cnl-5A+

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

36.46

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

7.07

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

3.12

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

186.02

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

QTgw.cnl-5A-

</td>

<td class="gt_row gt_left">

34.26

</td>

<td class="gt_row gt_left">

7.004

</td>

<td class="gt_row gt_left">

2.996

</td>

<td class="gt_row gt_left">

183.77

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left" style="font-style: italic;">

</td>

<td class="gt_row gt_left">

6.4% \*\*\*

</td>

<td class="gt_row gt_left">

0.92%

</td>

<td class="gt_row gt_left">

4.2% \*\*\*

</td>

<td class="gt_row gt_left">

1.2% \*\*

</td>

</tr>

<tr>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

H^2

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black; font-style: italic;">

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

0.68

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

0.75

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

0.81

</td>

<td class="gt_row gt_left" style="border-top-width: 2px; border-top-style: solid; border-top-color: black;">

0.78

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

### Parent and check entry averages across all observations

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#msroixerxk .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#msroixerxk .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#msroixerxk .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#msroixerxk .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#msroixerxk .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#msroixerxk .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#msroixerxk .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#msroixerxk .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#msroixerxk .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#msroixerxk .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#msroixerxk .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#msroixerxk .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#msroixerxk .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#msroixerxk .gt_from_md > :first-child {
  margin-top: 0;
}

#msroixerxk .gt_from_md > :last-child {
  margin-bottom: 0;
}

#msroixerxk .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#msroixerxk .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#msroixerxk .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#msroixerxk .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#msroixerxk .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#msroixerxk .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#msroixerxk .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#msroixerxk .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#msroixerxk .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#msroixerxk .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#msroixerxk .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#msroixerxk .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#msroixerxk .gt_left {
  text-align: left;
}

#msroixerxk .gt_center {
  text-align: center;
}

#msroixerxk .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#msroixerxk .gt_font_normal {
  font-weight: normal;
}

#msroixerxk .gt_font_bold {
  font-weight: bold;
}

#msroixerxk .gt_font_italic {
  font-style: italic;
}

#msroixerxk .gt_super {
  font-size: 65%;
}

#msroixerxk .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="msroixerxk" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

Genotype

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

TGW (g)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

GL (mm)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

GW (mm)

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

HD (julian)

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

Opata

</td>

<td class="gt_row gt_left">

37.985

</td>

<td class="gt_row gt_left">

6.686

</td>

<td class="gt_row gt_left">

3.363

</td>

<td class="gt_row gt_left">

171

</td>

</tr>

<tr>

<td class="gt_row gt_left">

W7984

</td>

<td class="gt_row gt_left">

43.001

</td>

<td class="gt_row gt_left">

8.106

</td>

<td class="gt_row gt_left">

3.444

</td>

<td class="gt_row gt_left">

175.667

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

13.20%

</td>

<td class="gt_row gt_left">

21.20%

</td>

<td class="gt_row gt_left">

2.40%

</td>

<td class="gt_row gt_left">

2.80%

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Tom

</td>

<td class="gt_row gt_left">

42.882

</td>

<td class="gt_row gt_left">

6.827

</td>

<td class="gt_row gt_left">

3.532

</td>

<td class="gt_row gt_left">

169.833

</td>

</tr>

<tr>

<td class="gt_row gt_left">

Glenn

</td>

<td class="gt_row gt_left">

32.69

</td>

<td class="gt_row gt_left">

5.879

</td>

<td class="gt_row gt_left">

3.256

</td>

<td class="gt_row gt_left">

169.5

</td>

</tr>

<tr>

<td class="gt_row gt_left">

</td>

<td class="gt_row gt_left">

31.20%

</td>

<td class="gt_row gt_left">

16.10%

</td>

<td class="gt_row gt_left">

8.50%

</td>

<td class="gt_row gt_left">

0.19%

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

### Boxplots chromosome 5AS vs TGW

![](script_S1_files/figure-gfm/5AS%20plot-1.png)<!-- -->

### SynOpDH chr5AS and QTgw.cnl-5A linkage disequilibrium

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#whmracclwe .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#whmracclwe .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#whmracclwe .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#whmracclwe .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#whmracclwe .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#whmracclwe .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#whmracclwe .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#whmracclwe .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#whmracclwe .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#whmracclwe .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#whmracclwe .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#whmracclwe .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#whmracclwe .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#whmracclwe .gt_from_md > :first-child {
  margin-top: 0;
}

#whmracclwe .gt_from_md > :last-child {
  margin-bottom: 0;
}

#whmracclwe .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#whmracclwe .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#whmracclwe .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#whmracclwe .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#whmracclwe .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#whmracclwe .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#whmracclwe .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#whmracclwe .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#whmracclwe .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#whmracclwe .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#whmracclwe .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#whmracclwe .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#whmracclwe .gt_left {
  text-align: left;
}

#whmracclwe .gt_center {
  text-align: center;
}

#whmracclwe .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#whmracclwe .gt_font_normal {
  font-weight: normal;
}

#whmracclwe .gt_font_bold {
  font-weight: bold;
}

#whmracclwe .gt_font_italic {
  font-style: italic;
}

#whmracclwe .gt_super {
  font-size: 65%;
}

#whmracclwe .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="whmracclwe" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="3">

<span class="gt_column_spanner">Recombinants</span>

</th>

<th class="gt_col_heading gt_center gt_columns_bottom_border" rowspan="2" colspan="1" style="font-weight: bold;">

TGW

</th>

<th class="gt_col_heading gt_center gt_columns_bottom_border" rowspan="2" colspan="1" style="font-weight: bold;">

GL

</th>

<th class="gt_col_heading gt_center gt_columns_bottom_border" rowspan="2" colspan="1" style="font-weight: bold;">

GW

</th>

<th class="gt_col_heading gt_center gt_columns_bottom_border" rowspan="2" colspan="1" style="font-weight: bold;">

HD

</th>

</tr>

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

Entry

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

chr5AS

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" style="font-weight: bold;">

QTgw.cnl\_5A

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

HC05.B110

</td>

<td class="gt_row gt_center">

present

</td>

<td class="gt_row gt_center">

W7984

</td>

<td class="gt_row gt_right">

32.61

</td>

<td class="gt_row gt_right">

6.14

</td>

<td class="gt_row gt_right">

3.06

</td>

<td class="gt_row gt_right">

182.55

</td>

</tr>

<tr>

<td class="gt_row gt_left">

HC06.C514

</td>

<td class="gt_row gt_center">

present

</td>

<td class="gt_row gt_center">

W7984

</td>

<td class="gt_row gt_right">

35.55

</td>

<td class="gt_row gt_right">

7.27

</td>

<td class="gt_row gt_right">

3.11

</td>

<td class="gt_row gt_right">

178.04

</td>

</tr>

<tr>

<td class="gt_row gt_left">

HC06.H939

</td>

<td class="gt_row gt_center">

present

</td>

<td class="gt_row gt_center">

W7984

</td>

<td class="gt_row gt_right">

35.52

</td>

<td class="gt_row gt_right">

7.07

</td>

<td class="gt_row gt_right">

3.08

</td>

<td class="gt_row gt_right">

178.28

</td>

</tr>

<tr>

<td class="gt_row gt_left">

HC06.H775

</td>

<td class="gt_row gt_center">

absent

</td>

<td class="gt_row gt_center">

Opata

</td>

<td class="gt_row gt_right">

40.38

</td>

<td class="gt_row gt_right">

7.27

</td>

<td class="gt_row gt_right">

3.21

</td>

<td class="gt_row gt_right">

188.25

</td>

</tr>

<tr>

<td class="gt_row gt_left">

HC06.K593

</td>

<td class="gt_row gt_center">

absent

</td>

<td class="gt_row gt_center">

Opata

</td>

<td class="gt_row gt_right">

33.89

</td>

<td class="gt_row gt_right">

6.95

</td>

<td class="gt_row gt_right">

3.00

</td>

<td class="gt_row gt_right">

184.69

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#wmblqilhtj .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#wmblqilhtj .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#wmblqilhtj .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#wmblqilhtj .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#wmblqilhtj .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmblqilhtj .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#wmblqilhtj .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#wmblqilhtj .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#wmblqilhtj .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#wmblqilhtj .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#wmblqilhtj .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#wmblqilhtj .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#wmblqilhtj .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#wmblqilhtj .gt_from_md > :first-child {
  margin-top: 0;
}

#wmblqilhtj .gt_from_md > :last-child {
  margin-bottom: 0;
}

#wmblqilhtj .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#wmblqilhtj .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#wmblqilhtj .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmblqilhtj .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#wmblqilhtj .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#wmblqilhtj .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#wmblqilhtj .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#wmblqilhtj .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#wmblqilhtj .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#wmblqilhtj .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#wmblqilhtj .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#wmblqilhtj .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#wmblqilhtj .gt_left {
  text-align: left;
}

#wmblqilhtj .gt_center {
  text-align: center;
}

#wmblqilhtj .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#wmblqilhtj .gt_font_normal {
  font-weight: normal;
}

#wmblqilhtj .gt_font_bold {
  font-weight: bold;
}

#wmblqilhtj .gt_font_italic {
  font-style: italic;
}

#wmblqilhtj .gt_super {
  font-size: 65%;
}

#wmblqilhtj .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="wmblqilhtj" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="font-weight: bold;">

Allele frequency

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left gt_stub">

5AS, present

</td>

<td class="gt_row gt_right">

0.48

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

5AS, absent

</td>

<td class="gt_row gt_right">

0.46

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

QTgw.cnl-5A, Opata

</td>

<td class="gt_row gt_right">

0.49

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

QTgw.cnl-5A, W7984

</td>

<td class="gt_row gt_right">

0.49

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#fhovegbmlc .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#fhovegbmlc .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#fhovegbmlc .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#fhovegbmlc .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#fhovegbmlc .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fhovegbmlc .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#fhovegbmlc .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#fhovegbmlc .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#fhovegbmlc .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#fhovegbmlc .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#fhovegbmlc .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#fhovegbmlc .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#fhovegbmlc .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#fhovegbmlc .gt_from_md > :first-child {
  margin-top: 0;
}

#fhovegbmlc .gt_from_md > :last-child {
  margin-bottom: 0;
}

#fhovegbmlc .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#fhovegbmlc .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#fhovegbmlc .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#fhovegbmlc .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#fhovegbmlc .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#fhovegbmlc .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#fhovegbmlc .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#fhovegbmlc .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#fhovegbmlc .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#fhovegbmlc .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#fhovegbmlc .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#fhovegbmlc .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#fhovegbmlc .gt_left {
  text-align: left;
}

#fhovegbmlc .gt_center {
  text-align: center;
}

#fhovegbmlc .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#fhovegbmlc .gt_font_normal {
  font-weight: normal;
}

#fhovegbmlc .gt_font_bold {
  font-weight: bold;
}

#fhovegbmlc .gt_font_italic {
  font-style: italic;
}

#fhovegbmlc .gt_super {
  font-size: 65%;
}

#fhovegbmlc .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="fhovegbmlc" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1">

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="font-weight: bold;">

Haplotype frequency

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left gt_stub">

present:Opata

</td>

<td class="gt_row gt_right">

0.46

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

present:W7984

</td>

<td class="gt_row gt_right">

0.02

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

absent:Opata

</td>

<td class="gt_row gt_right">

0.01

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

absent:W7984

</td>

<td class="gt_row gt_right">

0.44

</td>

</tr>

<tr>

<td class="gt_row gt_left gt_stub">

coefficient of correlation

</td>

<td class="gt_row gt_right" style="font-weight: bold;">

0.95

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

### T-test: phenotype vs chromosome 5AS presence / absence OR QTL Opata / W7984 allele

<!--html_preserve-->

<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#adikfltnct .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: 0;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#adikfltnct .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#adikfltnct .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#adikfltnct .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#adikfltnct .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#adikfltnct .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#adikfltnct .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#adikfltnct .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#adikfltnct .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#adikfltnct .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#adikfltnct .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#adikfltnct .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#adikfltnct .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#adikfltnct .gt_from_md > :first-child {
  margin-top: 0;
}

#adikfltnct .gt_from_md > :last-child {
  margin-bottom: 0;
}

#adikfltnct .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#adikfltnct .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#adikfltnct .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#adikfltnct .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#adikfltnct .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#adikfltnct .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#adikfltnct .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#adikfltnct .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#adikfltnct .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#adikfltnct .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#adikfltnct .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#adikfltnct .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#adikfltnct .gt_left {
  text-align: left;
}

#adikfltnct .gt_center {
  text-align: center;
}

#adikfltnct .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#adikfltnct .gt_font_normal {
  font-weight: normal;
}

#adikfltnct .gt_font_bold {
  font-weight: bold;
}

#adikfltnct .gt_font_italic {
  font-style: italic;
}

#adikfltnct .gt_super {
  font-size: 65%;
}

#adikfltnct .gt_footnote_marks {
  font-style: italic;
  font-size: 65%;
}
</style>

<div id="adikfltnct" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">

<table class="gt_table">

<thead class="gt_col_headings">

<tr>

<th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" style="font-weight: bold;">

t.test

</th>

<th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" style="font-weight: bold;">

P value

</th>

</tr>

</thead>

<tbody class="gt_table_body">

<tr>

<td class="gt_row gt_left">

TGW \~ 5AS

</td>

<td class="gt_row gt_right">

5.910 × 10<sup class='gt_super'>−4</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

TGW \~ QTL

</td>

<td class="gt_row gt_right">

5.450 × 10<sup class='gt_super'>−5</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

GL \~ 5AS

</td>

<td class="gt_row gt_right">

4.420 × 10<sup class='gt_super'>−1</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

GL \~ QTL

</td>

<td class="gt_row gt_right">

2.050 × 10<sup class='gt_super'>−1</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

GW \~ 5AS

</td>

<td class="gt_row gt_right">

4.490 × 10<sup class='gt_super'>−15</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

GW \~ QTL

</td>

<td class="gt_row gt_right">

4.620 × 10<sup class='gt_super'>−16</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

HD \~ 5AS

</td>

<td class="gt_row gt_right">

9.000 × 10<sup class='gt_super'>−3</sup>

</td>

</tr>

<tr>

<td class="gt_row gt_left">

HD \~ QTL

</td>

<td class="gt_row gt_right">

3.390 × 10<sup class='gt_super'>−4</sup>

</td>

</tr>

</tbody>

</table>

</div>

<!--/html_preserve-->

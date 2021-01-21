script\_S3: DPA analysis
================

All packages, data, and statistical analysis for reproducing SynOp HIF
DPA results reported in Taagen et al. 2021. Please see `script_S3.Rmd`
for full R script.

## 2019 field environmnet

A subset of ten HIF haplotypes (five *QTgw.cnl-5A-* and five
*QTgw.cnl-5A+*) were selected for a days post anthesis (DPA) study that
tracked fresh and dry grain weight (`FW`, `DW`) and grain width and
grain length (`GW`, `GL`) development over 0, 4, 10, 16 and 22 DPA.

**Load packages:**

``` r
library(tidyverse) # R/tidyverse package version 1.3.0  
library(knitr) # R/knitr package version 1.28 
library(kableExtra) # R/kableExtra package version 1.1.0
library(tibble) # R/tible package version 2.1.3
library(lme4) # R/lme4 package version 1.1-26
library(emmeans) # R/emmeans package version 1.4.6
library(ggpubr) # R/ggpubr package version 0.4.0 
library(inauguration) #R/inauguration version 0.0.0.90
```

**Load data:**

``` r
#when public repo change  
DPA <- read_csv("https://raw.githubusercontent.com/etaagen/Taagen_2021_TPG/main/supplementary_3/file_S3.1.csv", na = "")
DPA$family =  as.factor(DPA$family)              
DPA$entry =  as.factor(DPA$entry)
DPA$genotype =  as.factor(DPA$genotype)    
DPA$plot =  as.factor(DPA$plot)  
DPA$tube =  as.factor(DPA$tube)
DPA$DPA =  as.factor(DPA$DPA)  
```

**10 HIF haplotype entries**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

Group

</th>

<th style="text-align:left;">

Entry

</th>

<th style="text-align:left;">

QTgw.cnl-5A

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="2">

1

</td>

<td style="text-align:left;">

7-956-2-19-1-44

</td>

<td style="text-align:left;">

W7984 (-)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-956-2-19-1-54

</td>

<td style="text-align:left;">

Opata (+)

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="2">

2

</td>

<td style="text-align:left;">

7-956-2-50-1-50

</td>

<td style="text-align:left;">

W7984 (-)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-956-2-50-1-01

</td>

<td style="text-align:left;">

Opata (+)

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="2">

3

</td>

<td style="text-align:left;">

7-1201-10-36-2-18

</td>

<td style="text-align:left;">

W7984 (-)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-1201-10-36-2-04

</td>

<td style="text-align:left;">

Opata (+)

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="2">

4

</td>

<td style="text-align:left;">

7-1201-10-36-2-10

</td>

<td style="text-align:left;">

W7984 (-)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-1201-10-36-2-25

</td>

<td style="text-align:left;">

Opata (+)

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="2">

5

</td>

<td style="text-align:left;">

7-1201-10-45-2-23

</td>

<td style="text-align:left;">

W7984 (-)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-1201-10-45-2-23

</td>

<td style="text-align:left;">

Opata (+)

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> QTgw.cnl-5A: QTL genotype (impact on trait)

</td>

</tr>

</tfoot>

</table>

### Statistical analysis

Single-trait mixed models with a fixed interaction between `genotype`
and `DPA`, and a random effect of `plot` replicate nested within
`genotype`, and post hoc comparison of least-squares means for genotype
and DPA performed within each model. To avoid repetition, the code is
shown only for the first group. Please see the .Rmd file for all code.

**Group 1 W7984 - Opata contrast over 22 DPA**

``` r
group1 <- DPA[DPA$family=="31", ]

#grain width
model1_GW <- lmer(GW ~ genotype*DPA + (1|plot:genotype), 
                  data = group1)
#summary(model1_GW) 
GW_1 <- emmeans(model1_GW, pairwise ~ genotype|DPA)
#bonferroni multiply p value by 5 (because 5 DPA multiple testing conditions, keeps alpha < 0.05)
#Or bonferroni for five DPA multiple testing condtions, 0.05/5 = 0.01
GW_1_df <- as.data.frame(summary(GW_1$contrasts))[,c('DPA', 'p.value')] 

#grain length
model1_GL <- lmer(GL ~ genotype*DPA + (1|plot:genotype), 
                  data = group1)
#summary(model1_GL) 
GL_1 <- emmeans(model1_GL, pairwise ~ genotype|DPA)
#bonferroni for five DPA multiple testing condtions, 0.05/5 = 0.01
GL_1_df <- as.data.frame(summary(GL_1$contrasts))[,c('DPA', 'p.value')] 

#fresh weight
model1_FW <- lmer(FW ~ genotype*DPA + (1|plot:genotype), 
                  data = group1)
#summary(model1_FW) 
FW_1 <- emmeans(model1_FW, pairwise ~ genotype|DPA)
#bonferroni for five DPA multiple testing condtions, 0.05/5 = 0.01
FW_1_df <- as.data.frame(summary(FW_1$contrasts))[,c('DPA', 'p.value')] 

#dry weight
model1_DW <- lmer(DW ~ genotype*DPA + (1|plot:genotype), 
                  data = group1)
#summary(model1_DW) 
DW_1 <- emmeans(model1_DW, pairwise ~ genotype|DPA)
#bonferroni for five DPA multiple testing condtions, 0.05/5 = 0.01
DW_1_df <- as.data.frame(summary(DW_1$contrasts))[,c('DPA', 'p.value')] 

summary <- GW_1_df %>% 
  left_join(GL_1_df, by = "DPA") %>% 
  mutate(p.value_GW = p.value.x, p.value.x = NULL, 
         p.value_GL = p.value.y, p.value.y = NULL) %>% 
  left_join(FW_1_df, by = "DPA") %>% 
  mutate(p.value_FW = p.value, p.value = NULL) %>% 
  left_join(DW_1_df, by = "DPA") %>% 
  mutate(p.value_DW = p.value, p.value = NULL)
  
kable(summary, 
        escape = F,
        digits = 8,
        format.args = list(scientific = FALSE)) %>% 
  kable_styling(bootstrap_options = "striped", full_width = FALSE,  position = "left") %>% 
  footnote(general = "Bonferroni correction for five DPA multiple testing condtions, 0.05/5, p-values less than 0.01 are significant") 
```

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.53447160

</td>

<td style="text-align:right;">

0.34836809

</td>

<td style="text-align:right;">

0.89025966

</td>

<td style="text-align:right;">

0.97253688

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.25828285

</td>

<td style="text-align:right;">

0.42067562

</td>

<td style="text-align:right;">

0.49900294

</td>

<td style="text-align:right;">

0.76260961

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.02377324

</td>

<td style="text-align:right;">

0.00629607

</td>

<td style="text-align:right;">

0.00999664

</td>

<td style="text-align:right;">

0.18486105

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00637367

</td>

<td style="text-align:right;">

0.00168532

</td>

<td style="text-align:right;">

0.00001194

</td>

<td style="text-align:right;">

0.00053771

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.00055109

</td>

<td style="text-align:right;">

0.00390687

</td>

<td style="text-align:right;">

0.00000032

</td>

<td style="text-align:right;">

0.00000136

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for five DPA multiple testing
condtions, 0.05/5, p-values less than 0.01 are significant

</td>

</tr>

</tfoot>

</table>

``` r
#set color palette
pal2 <- inauguration("inauguration_2021", 2)

p1 <- ggplot(data = group1, aes(x = DPA, y=GW, color = genotype, group=genotype))+
  geom_jitter(width = 0.1, alpha = 0.4, na.rm = TRUE)+ 
  geom_line(stat="summary", na.rm = TRUE)+
  scale_color_manual(values = c(pal2), 
                     labels = c("Opata (+)", "W7984 (-)"))+
  theme_classic2()+
  ylab("grain width (mm)")+
  labs(color = "QTgw.cnl-5A") +
  theme(text = element_text(size=12), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        legend.title = element_text(face="bold")) 

p2 <- ggplot(data = group1, aes(x = DPA, y=GL, color = genotype, group=genotype))+
  geom_jitter(width = 0.1, alpha = 0.4, na.rm = TRUE)+ 
  geom_line(stat="summary", na.rm = TRUE)+
  scale_color_manual(values = pal2, 
                     labels = c("Opata (+)", "W7984 (-)"))+
  theme_classic2()+
  ylab("grain length (mm)")+
  labs(color = "QTgw.cnl-5A") +
  theme(text = element_text(size=12), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        legend.title = element_text(face="bold")) 

p3 <- ggplot(data = group1, aes(x = DPA, y=FW, color = genotype, group=genotype))+
  geom_jitter(width = 0.1, alpha = 0.4, na.rm = TRUE)+ 
  geom_line(stat="summary", na.rm = TRUE)+
  scale_color_manual(values = pal2, 
                     labels = c("Opata (+)", "W7984 (-)"))+
  theme_classic2()+
  ylab("fresh weight (g)")+
  labs(color = "QTgw.cnl-5A") +
  theme(text = element_text(size=12), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        legend.title = element_text(face="bold"))

p4 <- ggplot(data = group1, aes(x = DPA, y=DW, color = genotype, group=genotype))+
  geom_jitter(width = 0.1, alpha = 0.4, na.rm = TRUE)+ 
  geom_line(stat="summary", na.rm = TRUE)+
  scale_color_manual(values = pal2, 
                     labels = c("Opata (+)", "W7984 (-)"))+
  theme_classic2()+
  ylab("dry weight (g)")+
  labs(color = "QTgw.cnl-5A") +
  theme(text = element_text(size=12), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        legend.title = element_text(face="bold"))

ggarrange(p1, p2, p3, p4, 
          ncol = 4,
          common.legend = TRUE) 
```

![](script_S3_files/figure-gfm/group%201%20figures-1.png)<!-- -->

**Group 2 W7984 - Opata contrast over 22 DPA**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.47302271

</td>

<td style="text-align:right;">

0.82262608

</td>

<td style="text-align:right;">

0.96582857

</td>

<td style="text-align:right;">

0.9400937

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.34649778

</td>

<td style="text-align:right;">

0.13427210

</td>

<td style="text-align:right;">

0.90523113

</td>

<td style="text-align:right;">

0.9188970

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.01127709

</td>

<td style="text-align:right;">

0.81599806

</td>

<td style="text-align:right;">

0.03384936

</td>

<td style="text-align:right;">

0.3639281

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00055474

</td>

<td style="text-align:right;">

0.00014767

</td>

<td style="text-align:right;">

0.00000018

</td>

<td style="text-align:right;">

0.0000000

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.00606014

</td>

<td style="text-align:right;">

0.72831687

</td>

<td style="text-align:right;">

0.00000408

</td>

<td style="text-align:right;">

0.0000000

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for five DPA multiple testing
condtions, 0.05/5, p-values less than 0.01 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/group%202%20figures-1.png)<!-- -->

**Group 3 W7984 - Opata contrast over 22 DPA**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.01068069

</td>

<td style="text-align:right;">

0.04053759

</td>

<td style="text-align:right;">

0.84024103

</td>

<td style="text-align:right;">

0.91457048

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.02380156

</td>

<td style="text-align:right;">

0.73871627

</td>

<td style="text-align:right;">

0.00467872

</td>

<td style="text-align:right;">

0.86031250

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.02837364

</td>

<td style="text-align:right;">

0.03122914

</td>

<td style="text-align:right;">

0.09494213

</td>

<td style="text-align:right;">

0.14477112

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00008333

</td>

<td style="text-align:right;">

0.01305512

</td>

<td style="text-align:right;">

0.00545166

</td>

<td style="text-align:right;">

0.00000003

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.01116197

</td>

<td style="text-align:right;">

0.79382199

</td>

<td style="text-align:right;">

0.03199825

</td>

<td style="text-align:right;">

0.00000002

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for five DPA multiple testing
condtions, 0.05/5, p-values less than 0.01 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/group%203%20figures-1.png)<!-- -->

**Group 4 W7984 - Opata contrast over 22 DPA**  

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.01068069

</td>

<td style="text-align:right;">

0.04053759

</td>

<td style="text-align:right;">

0.84024103

</td>

<td style="text-align:right;">

0.91457048

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.02380156

</td>

<td style="text-align:right;">

0.73871627

</td>

<td style="text-align:right;">

0.00467872

</td>

<td style="text-align:right;">

0.86031250

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.02837364

</td>

<td style="text-align:right;">

0.03122914

</td>

<td style="text-align:right;">

0.09494213

</td>

<td style="text-align:right;">

0.14477112

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00008333

</td>

<td style="text-align:right;">

0.01305512

</td>

<td style="text-align:right;">

0.00545166

</td>

<td style="text-align:right;">

0.00000003

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.01116197

</td>

<td style="text-align:right;">

0.79382199

</td>

<td style="text-align:right;">

0.03199825

</td>

<td style="text-align:right;">

0.00000002

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for five DPA multiple testing
condtions, 0.05/5, p-values less than 0.01 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/group%204%20figures-1.png)<!-- -->

**Group 5 W7984 - Opata contrast over 22 DPA**  

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.73878234

</td>

<td style="text-align:right;">

0.30132072

</td>

<td style="text-align:right;">

0.93329627

</td>

<td style="text-align:right;">

0.96829286

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.00300203

</td>

<td style="text-align:right;">

0.12317939

</td>

<td style="text-align:right;">

0.47883846

</td>

<td style="text-align:right;">

0.84430377

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.00000731

</td>

<td style="text-align:right;">

0.02071055

</td>

<td style="text-align:right;">

0.00029872

</td>

<td style="text-align:right;">

0.09193257

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00045298

</td>

<td style="text-align:right;">

0.01274794

</td>

<td style="text-align:right;">

0.00000001

</td>

<td style="text-align:right;">

0.00000000

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.00031106

</td>

<td style="text-align:right;">

0.63624089

</td>

<td style="text-align:right;">

0.00001702

</td>

<td style="text-align:right;">

0.00000005

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for five DPA multiple testing
condtions, 0.05/5, p-values less than 0.01 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/group%205%20figures-1.png)<!-- -->

**All groups W7984 - Opata contrast over 22 DPA**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.358681321728

</td>

<td style="text-align:right;">

0.931607558092

</td>

<td style="text-align:right;">

0.894344585479

</td>

<td style="text-align:right;">

0.953358712

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.003563114740

</td>

<td style="text-align:right;">

0.191845485565

</td>

<td style="text-align:right;">

0.059596821208

</td>

<td style="text-align:right;">

0.853929328

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.000000004181

</td>

<td style="text-align:right;">

0.000000854124

</td>

<td style="text-align:right;">

0.000000328478

</td>

<td style="text-align:right;">

0.002787326

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.000000000000

</td>

<td style="text-align:right;">

0.000000002517

</td>

<td style="text-align:right;">

0.000000000000

</td>

<td style="text-align:right;">

0.000000000

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.000000000000

</td>

<td style="text-align:right;">

0.000417597869

</td>

<td style="text-align:right;">

0.000000000000

</td>

<td style="text-align:right;">

0.000000000

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for five DPA multiple testing
condtions across 5 groups, 0.05/25, p-values less than 0.002 are
significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/all%20groups%20figures-1.png)<!-- -->

**All group mean trait values**

``` r
W7984_10 <- DPA %>% 
  filter(genotype=="W7984", DPA == 10) 
Opata_10 <- DPA %>% 
  filter(genotype=="Opata", DPA == 10) 
W7984_16 <- DPA %>% 
  filter(genotype=="W7984", DPA == 16) 
Opata_16 <- DPA %>% 
  filter(genotype=="Opata", DPA == 16) 

GL_opata <- mean(Opata_10$GL, na.rm = TRUE)
GL_W7984 <- mean(W7984_10$GL, na.rm = TRUE)
GW_opata <- mean(Opata_10$GW, na.rm = TRUE)
GW_W7984 <- mean(W7984_10$GW, na.rm = TRUE)
FW_opata <- mean(Opata_10$FW, na.rm = TRUE)
FW_W7984 <- mean(W7984_10$FW, na.rm = TRUE)
DW_opata <- mean(Opata_16$FW, na.rm = TRUE)
DW_W7984 <- mean(W7984_16$FW, na.rm = TRUE) 

GL_opata/GL_W7984
```

    ## [1] 1.045804

``` r
GW_opata/GW_W7984
```

    ## [1] 1.057114

``` r
FW_opata/FW_W7984
```

    ## [1] 1.189599

``` r
DW_opata/DW_W7984
```

    ## [1] 1.252364

## 2020 greenhouse environment

A subset of four *QTgw.cnl-5A* HIF haplotypes were selected to validate
a days post anthesis (DPA) study in a greenhouse environment. The four
*QTgw.cnl-5A* haplotypes were also used for the RNA-seq experiment, see
table below. For plants grown in the greenhouse we measured fresh and
dry grain weight (`FW`, `DW`) and grain width and grain length (`GW`,
`GL`) development over 0, 4, 10, 16, 22 and 28 DPA, and senescence.

**Load data:**

``` r
#when public repo change  
DPA <- read.csv("https://raw.githubusercontent.com/etaagen/Taagen_2021_TPG/main/supplementary_3/file_S3.2.csv")
DPA$family <- as.factor(DPA$family)
DPA$entry <- as.factor(DPA$entry)
DPA$haplotype <- as.factor(DPA$haplotype)
DPA$gh_stake <- as.factor(DPA$gh_stake)
DPA$tube <- as.factor(DPA$tube)
#DPA$DPA <- as.factor(DPA$DPA)
DPA$DPA <- ordered(DPA$DPA, levels = c("0", "4", "10", "16", "22", "28", "sen."))
```

**SynOpHIF *QTgw.cnl-5A* Haplotype DPA entries**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

Entry

</th>

<th style="text-align:left;">

QTgw.cnl-5A haplotype

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

7-956-2-19-1-31-03

</td>

<td style="text-align:left;">

Opata (+)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-956-2-19-1-44

</td>

<td style="text-align:left;">

W7984 (-)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-956-2-19-1-31-05

</td>

<td style="text-align:left;">

Crossover one, O:W (+)

</td>

</tr>

<tr>

<td style="text-align:left;">

7-956-2-12-1-69-07

</td>

<td style="text-align:left;">

Crossover two, W:O:W (+)

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> QTgw.cnl-5A haplotype (impact on trait)

</td>

</tr>

</tfoot>

</table>

**DPA trait means by haplotype**

    ## `summarise()` regrouping output by 'DPA' (override with `.groups` argument)

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:left;">

haplotype

</th>

<th style="text-align:right;">

meanGW

</th>

<th style="text-align:right;">

meanGL

</th>

<th style="text-align:right;">

meanFW

</th>

<th style="text-align:right;">

meanDW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

0

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

2.21646

</td>

<td style="text-align:right;">

2.93536

</td>

<td style="text-align:right;">

0.0232400

</td>

<td style="text-align:right;">

0.0055200

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

2.10956

</td>

<td style="text-align:right;">

2.75162

</td>

<td style="text-align:right;">

0.0235000

</td>

<td style="text-align:right;">

0.0049600

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

2.26960

</td>

<td style="text-align:right;">

2.78472

</td>

<td style="text-align:right;">

0.0262600

</td>

<td style="text-align:right;">

0.0059800

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

2.24291

</td>

<td style="text-align:right;">

3.01327

</td>

<td style="text-align:right;">

0.0237600

</td>

<td style="text-align:right;">

0.0053100

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

4

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

3.41913

</td>

<td style="text-align:right;">

4.89661

</td>

<td style="text-align:right;">

0.1203000

</td>

<td style="text-align:right;">

0.0258000

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

3.63647

</td>

<td style="text-align:right;">

5.43507

</td>

<td style="text-align:right;">

0.1675500

</td>

<td style="text-align:right;">

0.0357100

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

3.57069

</td>

<td style="text-align:right;">

5.29482

</td>

<td style="text-align:right;">

0.1440400

</td>

<td style="text-align:right;">

0.0312000

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

3.45613

</td>

<td style="text-align:right;">

5.26223

</td>

<td style="text-align:right;">

0.1283500

</td>

<td style="text-align:right;">

0.0284200

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

10

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

4.23637

</td>

<td style="text-align:right;">

8.12747

</td>

<td style="text-align:right;">

0.4825100

</td>

<td style="text-align:right;">

0.1164200

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

4.28629

</td>

<td style="text-align:right;">

8.36212

</td>

<td style="text-align:right;">

0.5088800

</td>

<td style="text-align:right;">

0.1326000

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

4.46837

</td>

<td style="text-align:right;">

8.44918

</td>

<td style="text-align:right;">

0.5459400

</td>

<td style="text-align:right;">

0.1297600

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

4.21068

</td>

<td style="text-align:right;">

8.18203

</td>

<td style="text-align:right;">

0.4780000

</td>

<td style="text-align:right;">

0.1201000

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

16

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

5.04687

</td>

<td style="text-align:right;">

8.48762

</td>

<td style="text-align:right;">

0.8523800

</td>

<td style="text-align:right;">

0.3184800

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

4.81401

</td>

<td style="text-align:right;">

8.64505

</td>

<td style="text-align:right;">

0.7523300

</td>

<td style="text-align:right;">

0.3136400

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

5.19493

</td>

<td style="text-align:right;">

8.70518

</td>

<td style="text-align:right;">

0.9005500

</td>

<td style="text-align:right;">

0.3304600

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

4.73633

</td>

<td style="text-align:right;">

8.30035

</td>

<td style="text-align:right;">

0.7149500

</td>

<td style="text-align:right;">

0.2857000

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

22

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

4.92366

</td>

<td style="text-align:right;">

8.46900

</td>

<td style="text-align:right;">

0.8429500

</td>

<td style="text-align:right;">

0.4523700

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

4.72085

</td>

<td style="text-align:right;">

8.44185

</td>

<td style="text-align:right;">

0.7615500

</td>

<td style="text-align:right;">

0.4136800

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

5.09685

</td>

<td style="text-align:right;">

8.64028

</td>

<td style="text-align:right;">

0.8885400

</td>

<td style="text-align:right;">

0.4583500

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

4.70967

</td>

<td style="text-align:right;">

8.39573

</td>

<td style="text-align:right;">

0.7378800

</td>

<td style="text-align:right;">

0.4022200

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

28

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

4.29229

</td>

<td style="text-align:right;">

8.17815

</td>

<td style="text-align:right;">

0.6420800

</td>

<td style="text-align:right;">

0.4758300

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

4.83285

</td>

<td style="text-align:right;">

8.49445

</td>

<td style="text-align:right;">

0.8057800

</td>

<td style="text-align:right;">

0.4919900

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

4.18611

</td>

<td style="text-align:right;">

8.05834

</td>

<td style="text-align:right;">

0.5529600

</td>

<td style="text-align:right;">

0.4701000

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

3.63358

</td>

<td style="text-align:right;">

7.57576

</td>

<td style="text-align:right;">

0.4262900

</td>

<td style="text-align:right;">

0.3892800

</td>

</tr>

<tr>

<td style="text-align:left;vertical-align: top !important;" rowspan="4">

sen.

</td>

<td style="text-align:left;">

CO1

</td>

<td style="text-align:right;">

3.63674

</td>

<td style="text-align:right;">

7.06497

</td>

<td style="text-align:right;">

0.4511643

</td>

<td style="text-align:right;">

0.4511643

</td>

</tr>

<tr>

<td style="text-align:left;">

CO2

</td>

<td style="text-align:right;">

3.63839

</td>

<td style="text-align:right;">

7.14985

</td>

<td style="text-align:right;">

0.4269971

</td>

<td style="text-align:right;">

0.4269971

</td>

</tr>

<tr>

<td style="text-align:left;">

Opata

</td>

<td style="text-align:right;">

3.71801

</td>

<td style="text-align:right;">

7.15965

</td>

<td style="text-align:right;">

0.4488600

</td>

<td style="text-align:right;">

0.4488600

</td>

</tr>

<tr>

<td style="text-align:left;">

W7984

</td>

<td style="text-align:right;">

3.14918

</td>

<td style="text-align:right;">

6.81561

</td>

<td style="text-align:right;">

0.3236329

</td>

<td style="text-align:right;">

0.3236329

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> mean of 10 samples of 10 grains per haplotype DPA, and sen.
FW and DW are the same observation

</td>

</tr>

</tfoot>

</table>

### GH statistical analysis

Single-trait mixed models with a fixed interaction between `haplotype`
and `DPA`, and a random effect of `tube` nested within `haplotype`, and
post hoc comparison of least-squares means for haplotype and DPA
performed within each model. For equations please see .Rmd file.

**Opata vs W7984 haplotype contrasts over 28 DPA, and senescence**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.61383786

</td>

<td style="text-align:right;">

0.01868546

</td>

<td style="text-align:right;">

0.85852014

</td>

<td style="text-align:right;">

0.9178007

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.03184258

</td>

<td style="text-align:right;">

0.73461337

</td>

<td style="text-align:right;">

0.26440210

</td>

<td style="text-align:right;">

0.6685969

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.00000320

</td>

<td style="text-align:right;">

0.00617947

</td>

<td style="text-align:right;">

0.00000351

</td>

<td style="text-align:right;">

0.1384661

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00004629

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.0000000

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.01198988

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.0000000

</td>

</tr>

<tr>

<td style="text-align:left;">

28

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00000164

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.0000000

</td>

</tr>

<tr>

<td style="text-align:left;">

sen.

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00047783

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.0000000

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for 7 DPA multiple testing condtions,
0.05/7, p-values less than 0.007 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/O%20vs%20W%20figures-1.png)<!-- -->

**CO1 vs W7984 haplotype contrasts over 28 DPA, and senescence**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.62818263

</td>

<td style="text-align:right;">

0.47202198

</td>

<td style="text-align:right;">

0.9735517

</td>

<td style="text-align:right;">

0.97747334

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.49831161

</td>

<td style="text-align:right;">

0.00094862

</td>

<td style="text-align:right;">

0.6079625

</td>

<td style="text-align:right;">

0.72468793

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.63808095

</td>

<td style="text-align:right;">

0.61432346

</td>

<td style="text-align:right;">

0.7737281

</td>

<td style="text-align:right;">

0.62090012

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00000008

</td>

<td style="text-align:right;">

0.08537789

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.00002137

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.00014064

</td>

<td style="text-align:right;">

0.49876053

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.00000000

</td>

</tr>

<tr>

<td style="text-align:left;">

28

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00000014

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.00000000

</td>

</tr>

<tr>

<td style="text-align:left;">

sen.

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.02258274

</td>

<td style="text-align:right;">

0.0000000

</td>

<td style="text-align:right;">

0.00000000

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for 7 DPA multiple testing condtions,
0.05/7, p-values less than 0.007 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/CO1%20vs%20W%20figures-1.png)<!-- -->

**CO2 vs W7984 haplotype contrasts over 28 DPA, and senescence**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.00937467

</td>

<td style="text-align:right;">

0.01021618

</td>

<td style="text-align:right;">

0.98602158

</td>

<td style="text-align:right;">

0.96369774

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.00051017

</td>

<td style="text-align:right;">

0.08741740

</td>

<td style="text-align:right;">

0.00916361

</td>

<td style="text-align:right;">

0.34399416

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.13706030

</td>

<td style="text-align:right;">

0.07507666

</td>

<td style="text-align:right;">

0.03908809

</td>

<td style="text-align:right;">

0.10586809

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.12672234

</td>

<td style="text-align:right;">

0.00080176

</td>

<td style="text-align:right;">

0.01284992

</td>

<td style="text-align:right;">

0.00039572

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.82523587

</td>

<td style="text-align:right;">

0.64656105

</td>

<td style="text-align:right;">

0.11250660

</td>

<td style="text-align:right;">

0.13787971

</td>

</tr>

<tr>

<td style="text-align:left;">

28

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00000000

</td>

</tr>

<tr>

<td style="text-align:left;">

sen.

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00113514

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00000000

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for 7 DPA multiple testing condtions,
0.05/7, p-values less than 0.007 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/CO2%20vs%20W7984%20figures-1.png)<!-- -->

**CO1 vs CO2 haplotype contrasts over 28 DPA, and senescence**

<table class="table table-striped" style="width: auto !important; border-bottom: 0;">

<thead>

<tr>

<th style="text-align:left;">

DPA

</th>

<th style="text-align:right;">

p.value\_GW

</th>

<th style="text-align:right;">

p.value\_GL

</th>

<th style="text-align:right;">

p.value\_FW

</th>

<th style="text-align:right;">

p.value\_DW

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

0

</td>

<td style="text-align:right;">

0.05009473

</td>

<td style="text-align:right;">

0.05580622

</td>

<td style="text-align:right;">

0.98815333

</td>

<td style="text-align:right;">

0.94405870

</td>

</tr>

<tr>

<td style="text-align:left;">

4

</td>

<td style="text-align:right;">

0.00009888

</td>

<td style="text-align:right;">

0.00000010

</td>

<td style="text-align:right;">

0.00780391

</td>

<td style="text-align:right;">

0.21572489

</td>

</tr>

<tr>

<td style="text-align:left;">

10

</td>

<td style="text-align:right;">

0.35738226

</td>

<td style="text-align:right;">

0.01503802

</td>

<td style="text-align:right;">

0.13381845

</td>

<td style="text-align:right;">

0.04431307

</td>

</tr>

<tr>

<td style="text-align:left;">

16

</td>

<td style="text-align:right;">

0.00003273

</td>

<td style="text-align:right;">

0.10062304

</td>

<td style="text-align:right;">

0.00000007

</td>

<td style="text-align:right;">

0.54449611

</td>

</tr>

<tr>

<td style="text-align:left;">

22

</td>

<td style="text-align:right;">

0.00026557

</td>

<td style="text-align:right;">

0.77592899

</td>

<td style="text-align:right;">

0.00000799

</td>

<td style="text-align:right;">

0.00000346

</td>

</tr>

<tr>

<td style="text-align:left;">

28

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.00116580

</td>

<td style="text-align:right;">

0.00000000

</td>

<td style="text-align:right;">

0.04457205

</td>

</tr>

<tr>

<td style="text-align:left;">

sen.

</td>

<td style="text-align:right;">

0.97569060

</td>

<td style="text-align:right;">

0.37422444

</td>

<td style="text-align:right;">

0.16914405

</td>

<td style="text-align:right;">

0.00292935

</td>

</tr>

</tbody>

<tfoot>

<tr>

<td style="padding: 0; " colspan="100%">

<span style="font-style: italic;">Note: </span>

</td>

</tr>

<tr>

<td style="padding: 0; " colspan="100%">

<sup></sup> Bonferroni correction for 7 DPA multiple testing condtions,
0.05/7, p-values less than 0.007 are significant

</td>

</tr>

</tfoot>

</table>

![](script_S3_files/figure-gfm/CO1%20vs%20CO2%20figures-1.png)<!-- -->

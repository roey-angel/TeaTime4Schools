---
title: "TeaTime4schools: Joint analysis - bacteria"
subtitle: "TBI S and k comparison"
author: "Anne Daebeler Roey Angel"
date: "2021-05-19"
bibliography: references.bib
link-citations: yes
csl: fems-microbiology-ecology.csl
output:
  rmarkdown::html_document:
    toc: true
    toc_float: true
    toc_depth: 5
    keep_md: true
    number_sections: false
    highlight: "pygments"
    theme: "flatly"
    dev: "png"
    df_print: "kable"
    fig_caption: true
    code_folding: "show"
---







[roey.angel@bc.cas.cz](mailto: roey.angel@bc.cas.cz)  

## TBI stabilisation factor (S) and decomposition rate (k) comparison according to Keuskamp et al. [-@keuskamp_tea_2013]


```r
read.delim("TBI.csv", sep = "\t") %>% 
  mutate(across(c("Soil",
                  "Season"), 
                ~factor(.))) %>% 
  mutate(across("Season", ~fct_relevel(., "Winter", "Spring", "Summer", "Autumn"))) %>% 
  unite("Soil_Season", Soil:Season, remove = FALSE) ->  
  TBI 
```

### Test the differences in the stabilisation factor (S)

```r
(mod_S <- TestAlphaV3(data2test = TBI,
                        response_name = "S",
                        factor_names = c("Soil", "Season"),
                        boxcox.trans = FALSE))
```

```
## Call:
##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
## 
## Terms:
##                       Soil     Season Soil:Season  Residuals
## Sum of Squares  0.01600006 0.22325539  0.06312455 0.08150000
## Deg. of Freedom          2          3           6         33
## 
## Residual standard error: 0.04969605
## Estimated effects may be unbalanced
## 6 observations deleted due to missingness
```

![](TBI_figures/test diffs in S-1.svg)<!-- -->

```
## [1] "Unequal group sizes - showing SS type III"
## Anova Table (Type III tests)
## 
## Response: S
##              Sum Sq Df  F value    Pr(>F)    
## (Intercept) 0.65610  1 265.6601 < 2.2e-16 ***
## Soil        0.01552  2   3.1414  0.056381 .  
## Season      0.18505  3  24.9763 1.276e-08 ***
## Soil:Season 0.06312  6   4.2599  0.002772 ** 
## Residuals   0.08150 33                       
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Tables of means
## Grand mean
##       
## 0.334 
## 
##  Soil 
##     Cambisol Fluvisol Luvisol
##        0.357    0.312   0.334
## rep   15.000   16.000  14.000
## 
##  Season 
##     Winter Spring Summer Autumn
##      0.355  0.379   0.38  0.204
## rep 12.000 12.000  11.00 10.000
## 
##  Soil:Season 
##           Season
## Soil       Winter Spring Summer Autumn
##   Cambisol 0.41   0.40   0.43   0.14  
##   rep      4.00   4.00   4.00   3.00  
##   Fluvisol 0.32   0.36   0.31   0.26  
##   rep      4.00   4.00   4.00   4.00  
##   Luvisol  0.34   0.38   0.40   0.20  
##   rep      4.00   4.00   3.00   3.00
```

```
## Call:
##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
## 
## Terms:
##                       Soil     Season Soil:Season  Residuals
## Sum of Squares  0.01600006 0.22325539  0.06312455 0.08150000
## Deg. of Freedom          2          3           6         33
## 
## Residual standard error: 0.04969605
## Estimated effects may be unbalanced
## 6 observations deleted due to missingness
```

```r
# Post-hoc test
marginal <- emmeans(mod_S,
                   ~ Soil : Season)
summary(marginal)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Soil </th>
   <th style="text-align:left;"> Season </th>
   <th style="text-align:right;"> emmean </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> lower.CL </th>
   <th style="text-align:right;"> upper.CL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.4050000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3544463 </td>
   <td style="text-align:right;"> 0.4555537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.3200000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2694463 </td>
   <td style="text-align:right;"> 0.3705537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.3425000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2919463 </td>
   <td style="text-align:right;"> 0.3930537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.4000000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3494463 </td>
   <td style="text-align:right;"> 0.4505537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.3600000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3094463 </td>
   <td style="text-align:right;"> 0.4105537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.3775000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3269463 </td>
   <td style="text-align:right;"> 0.4280537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.4325000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3819463 </td>
   <td style="text-align:right;"> 0.4830537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.3125000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2619463 </td>
   <td style="text-align:right;"> 0.3630537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.4033333 </td>
   <td style="text-align:right;"> 0.028692 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3449590 </td>
   <td style="text-align:right;"> 0.4617077 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.1366667 </td>
   <td style="text-align:right;"> 0.028692 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.0782923 </td>
   <td style="text-align:right;"> 0.1950410 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.2550000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2044463 </td>
   <td style="text-align:right;"> 0.3055537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.1966667 </td>
   <td style="text-align:right;"> 0.028692 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.1382923 </td>
   <td style="text-align:right;"> 0.2550410 </td>
  </tr>
</tbody>
</table>

</div>

```r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

```
##  contrast                          estimate     SE df t.ratio p.value
##  Cambisol Winter - Fluvisol Winter  0.08500 0.0351 33  2.419  0.4219 
##  Cambisol Winter - Luvisol Winter   0.06250 0.0351 33  1.779  0.8177 
##  Cambisol Winter - Cambisol Spring  0.00500 0.0351 33  0.142  1.0000 
##  Cambisol Winter - Fluvisol Spring  0.04500 0.0351 33  1.281  0.9761 
##  Cambisol Winter - Luvisol Spring   0.02750 0.0351 33  0.783  0.9996 
##  Cambisol Winter - Cambisol Summer -0.02750 0.0351 33 -0.783  0.9996 
##  Cambisol Winter - Fluvisol Summer  0.09250 0.0351 33  2.632  0.3030 
##  Cambisol Winter - Luvisol Summer   0.00167 0.0380 33  0.044  1.0000 
##  Cambisol Winter - Cambisol Autumn  0.26833 0.0380 33  7.070  <.0001 
##  Cambisol Winter - Fluvisol Autumn  0.15000 0.0351 33  4.269  0.0073 
##  Cambisol Winter - Luvisol Autumn   0.20833 0.0380 33  5.489  0.0002 
##  Fluvisol Winter - Luvisol Winter  -0.02250 0.0351 33 -0.640  0.9999 
##  Fluvisol Winter - Cambisol Spring -0.08000 0.0351 33 -2.277  0.5111 
##  Fluvisol Winter - Fluvisol Spring -0.04000 0.0351 33 -1.138  0.9902 
##  Fluvisol Winter - Luvisol Spring  -0.05750 0.0351 33 -1.636  0.8829 
##  Fluvisol Winter - Cambisol Summer -0.11250 0.0351 33 -3.201  0.1009 
##  Fluvisol Winter - Fluvisol Summer  0.00750 0.0351 33  0.213  1.0000 
##  Fluvisol Winter - Luvisol Summer  -0.08333 0.0380 33 -2.196  0.5638 
##  Fluvisol Winter - Cambisol Autumn  0.18333 0.0380 33  4.830  0.0016 
##  Fluvisol Winter - Fluvisol Autumn  0.06500 0.0351 33  1.850  0.7797 
##  Fluvisol Winter - Luvisol Autumn   0.12333 0.0380 33  3.249  0.0909 
##  Luvisol Winter - Cambisol Spring  -0.05750 0.0351 33 -1.636  0.8829 
##  Luvisol Winter - Fluvisol Spring  -0.01750 0.0351 33 -0.498  1.0000 
##  Luvisol Winter - Luvisol Spring   -0.03500 0.0351 33 -0.996  0.9968 
##  Luvisol Winter - Cambisol Summer  -0.09000 0.0351 33 -2.561  0.3402 
##  Luvisol Winter - Fluvisol Summer   0.03000 0.0351 33  0.854  0.9992 
##  Luvisol Winter - Luvisol Summer   -0.06083 0.0380 33 -1.603  0.8960 
##  Luvisol Winter - Cambisol Autumn   0.20583 0.0380 33  5.423  0.0003 
##  Luvisol Winter - Fluvisol Autumn   0.08750 0.0351 33  2.490  0.3800 
##  Luvisol Winter - Luvisol Autumn    0.14583 0.0380 33  3.842  0.0222 
##  Cambisol Spring - Fluvisol Spring  0.04000 0.0351 33  1.138  0.9902 
##  Cambisol Spring - Luvisol Spring   0.02250 0.0351 33  0.640  0.9999 
##  Cambisol Spring - Cambisol Summer -0.03250 0.0351 33 -0.925  0.9983 
##  Cambisol Spring - Fluvisol Summer  0.08750 0.0351 33  2.490  0.3800 
##  Cambisol Spring - Luvisol Summer  -0.00333 0.0380 33 -0.088  1.0000 
##  Cambisol Spring - Cambisol Autumn  0.26333 0.0380 33  6.938  <.0001 
##  Cambisol Spring - Fluvisol Autumn  0.14500 0.0351 33  4.126  0.0107 
##  Cambisol Spring - Luvisol Autumn   0.20333 0.0380 33  5.357  0.0003 
##  Fluvisol Spring - Luvisol Spring  -0.01750 0.0351 33 -0.498  1.0000 
##  Fluvisol Spring - Cambisol Summer -0.07250 0.0351 33 -2.063  0.6500 
##  Fluvisol Spring - Fluvisol Summer  0.04750 0.0351 33  1.352  0.9649 
##  Fluvisol Spring - Luvisol Summer  -0.04333 0.0380 33 -1.142  0.9900 
##  Fluvisol Spring - Cambisol Autumn  0.22333 0.0380 33  5.884  0.0001 
##  Fluvisol Spring - Fluvisol Autumn  0.10500 0.0351 33  2.988  0.1575 
##  Fluvisol Spring - Luvisol Autumn   0.16333 0.0380 33  4.303  0.0067 
##  Luvisol Spring - Cambisol Summer  -0.05500 0.0351 33 -1.565  0.9095 
##  Luvisol Spring - Fluvisol Summer   0.06500 0.0351 33  1.850  0.7797 
##  Luvisol Spring - Luvisol Summer   -0.02583 0.0380 33 -0.681  0.9999 
##  Luvisol Spring - Cambisol Autumn   0.24083 0.0380 33  6.345  <.0001 
##  Luvisol Spring - Fluvisol Autumn   0.12250 0.0351 33  3.486  0.0530 
##  Luvisol Spring - Luvisol Autumn    0.18083 0.0380 33  4.764  0.0019 
##  Cambisol Summer - Fluvisol Summer  0.12000 0.0351 33  3.415  0.0626 
##  Cambisol Summer - Luvisol Summer   0.02917 0.0380 33  0.768  0.9997 
##  Cambisol Summer - Cambisol Autumn  0.29583 0.0380 33  7.794  <.0001 
##  Cambisol Summer - Fluvisol Autumn  0.17750 0.0351 33  5.051  0.0008 
##  Cambisol Summer - Luvisol Autumn   0.23583 0.0380 33  6.213  <.0001 
##  Fluvisol Summer - Luvisol Summer  -0.09083 0.0380 33 -2.393  0.4376 
##  Fluvisol Summer - Cambisol Autumn  0.17583 0.0380 33  4.633  0.0027 
##  Fluvisol Summer - Fluvisol Autumn  0.05750 0.0351 33  1.636  0.8829 
##  Fluvisol Summer - Luvisol Autumn   0.11583 0.0380 33  3.052  0.1384 
##  Luvisol Summer - Cambisol Autumn   0.26667 0.0406 33  6.572  <.0001 
##  Luvisol Summer - Fluvisol Autumn   0.14833 0.0380 33  3.908  0.0188 
##  Luvisol Summer - Luvisol Autumn    0.20667 0.0406 33  5.093  0.0007 
##  Cambisol Autumn - Fluvisol Autumn -0.11833 0.0380 33 -3.118  0.1207 
##  Cambisol Autumn - Luvisol Autumn  -0.06000 0.0406 33 -1.479  0.9362 
##  Fluvisol Autumn - Luvisol Autumn   0.05833 0.0380 33  1.537  0.9189 
## 
## P value adjustment: tukey method for comparing a family of 12 estimates
```

```r
(S_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Soil </th>
   <th style="text-align:left;"> Season </th>
   <th style="text-align:right;"> emmean </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> lower.CL </th>
   <th style="text-align:right;"> upper.CL </th>
   <th style="text-align:left;"> .group </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.1366667 </td>
   <td style="text-align:right;"> 0.028692 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.0485897 </td>
   <td style="text-align:right;"> 0.2247437 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.1966667 </td>
   <td style="text-align:right;"> 0.028692 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.1085897 </td>
   <td style="text-align:right;"> 0.2847437 </td>
   <td style="text-align:left;"> ab </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 11 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.2550000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.1787231 </td>
   <td style="text-align:right;"> 0.3312769 </td>
   <td style="text-align:left;"> abc </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.3125000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2362231 </td>
   <td style="text-align:right;"> 0.3887769 </td>
   <td style="text-align:left;"> bcd </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.3200000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2437231 </td>
   <td style="text-align:right;"> 0.3962769 </td>
   <td style="text-align:left;"> bcd </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.3425000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2662231 </td>
   <td style="text-align:right;"> 0.4187769 </td>
   <td style="text-align:left;"> cd </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.3600000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.2837231 </td>
   <td style="text-align:right;"> 0.4362769 </td>
   <td style="text-align:left;"> cd </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.3775000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3012231 </td>
   <td style="text-align:right;"> 0.4537769 </td>
   <td style="text-align:left;"> cd </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.4000000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3237231 </td>
   <td style="text-align:right;"> 0.4762769 </td>
   <td style="text-align:left;"> d </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.4033333 </td>
   <td style="text-align:right;"> 0.028692 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3152563 </td>
   <td style="text-align:right;"> 0.4914103 </td>
   <td style="text-align:left;"> d </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.4050000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3287231 </td>
   <td style="text-align:right;"> 0.4812769 </td>
   <td style="text-align:left;"> d </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.4325000 </td>
   <td style="text-align:right;"> 0.024848 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.3562231 </td>
   <td style="text-align:right;"> 0.5087769 </td>
   <td style="text-align:left;"> d </td>
  </tr>
</tbody>
</table>

</div>

```r
(mod_S %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_S_ANOVA)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Sum Sq </th>
   <th style="text-align:right;"> Mean Sq </th>
   <th style="text-align:right;"> F value </th>
   <th style="text-align:right;"> Pr(&gt;F) </th>
   <th style="text-align:right;"> Part Eta Sq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Soil </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0160001 </td>
   <td style="text-align:right;"> 0.0080000 </td>
   <td style="text-align:right;"> 3.239276 </td>
   <td style="text-align:right;"> 0.0519422 </td>
   <td style="text-align:right;"> 0.0416798 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Season </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.2232554 </td>
   <td style="text-align:right;"> 0.0744185 </td>
   <td style="text-align:right;"> 30.132630 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.5815760 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Soil:Season </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.0631245 </td>
   <td style="text-align:right;"> 0.0105208 </td>
   <td style="text-align:right;"> 4.259939 </td>
   <td style="text-align:right;"> 0.0027724 </td>
   <td style="text-align:right;"> 0.1644382 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 0.0815000 </td>
   <td style="text-align:right;"> 0.0024697 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.2123059 </td>
  </tr>
</tbody>
</table>

</div>

```r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_S, Soil ~ Season)
```

![](TBI_figures/test diffs in S-2.svg)<!-- -->

### Test the differences in the decomposition rate (k)

```r
(mod_k <- TestAlphaV3(data2test = TBI,
                        response_name = "k",
                        factor_names = c("Soil", "Season"),
                        boxcox.trans = FALSE))
```

```
## Call:
##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
## 
## Terms:
##                         Soil       Season  Soil:Season    Residuals
## Sum of Squares  0.0000087045 0.0002068740 0.0005405618 0.0017833333
## Deg. of Freedom            2            3            6           26
## 
## Residual standard error: 0.008281893
## Estimated effects may be unbalanced
## 13 observations deleted due to missingness
```

![](TBI_figures/test diffs in k-1.svg)<!-- -->

```
## [1] "Unequal group sizes - showing SS type III"
## Anova Table (Type III tests)
## 
## Response: k
##                Sum Sq Df  F value    Pr(>F)    
## (Intercept) 0.0084262  1 122.8499 2.386e-11 ***
## Soil        0.0000014  2   0.0102    0.9898    
## Season      0.0002289  3   1.1126    0.3619    
## Soil:Season 0.0005406  6   1.3135    0.2864    
## Residuals   0.0017833 26                       
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## Tables of means
## Grand mean
##            
## 0.01447368 
## 
##  Soil 
##     Cambisol Fluvisol Luvisol
##        0.015   0.0138  0.0146
## rep   12.000  13.0000 13.0000
## 
##  Season 
##     Winter  Spring  Summer Autumn
##     0.0167  0.0109  0.0151 0.0162
## rep 9.0000 11.0000 10.0000 8.0000
## 
##  Soil:Season 
##           Season
## Soil       Winter Spring Summer Autumn
##   Cambisol 0.02   0.01   0.02   0.01  
##   rep      3.00   3.00   3.00   3.00  
##   Fluvisol 0.02   0.01   0.01   0.03  
##   rep      3.00   4.00   4.00   2.00  
##   Luvisol  0.01   0.01   0.02   0.02  
##   rep      3.00   4.00   3.00   3.00
```

```
## Call:
##    aov(formula = as.formula(paste(response_name, paste(factor_names[1], 
##     factor_names[2], sep = " * "), sep = " ~ ")), data = data2test)
## 
## Terms:
##                         Soil       Season  Soil:Season    Residuals
## Sum of Squares  0.0000087045 0.0002068740 0.0005405618 0.0017833333
## Deg. of Freedom            2            3            6           26
## 
## Residual standard error: 0.008281893
## Estimated effects may be unbalanced
## 13 observations deleted due to missingness
```

```r
# Post-hoc test
marginal <- emmeans(mod_k,
                   ~ Soil : Season)
summary(marginal)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Soil </th>
   <th style="text-align:left;"> Season </th>
   <th style="text-align:right;"> emmean </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> lower.CL </th>
   <th style="text-align:right;"> upper.CL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.0200000 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0101714 </td>
   <td style="text-align:right;"> 0.0298286 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.0166667 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0068380 </td>
   <td style="text-align:right;"> 0.0264953 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.0133333 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0035047 </td>
   <td style="text-align:right;"> 0.0231620 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.0133333 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0035047 </td>
   <td style="text-align:right;"> 0.0231620 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0041409 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0014882 </td>
   <td style="text-align:right;"> 0.0185118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0041409 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0014882 </td>
   <td style="text-align:right;"> 0.0185118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.0166667 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0068380 </td>
   <td style="text-align:right;"> 0.0264953 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0041409 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0014882 </td>
   <td style="text-align:right;"> 0.0185118 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.0200000 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0101714 </td>
   <td style="text-align:right;"> 0.0298286 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0001714 </td>
   <td style="text-align:right;"> 0.0198286 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.0250000 </td>
   <td style="text-align:right;"> 0.0058562 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0129624 </td>
   <td style="text-align:right;"> 0.0370376 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.0166667 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0068380 </td>
   <td style="text-align:right;"> 0.0264953 </td>
  </tr>
</tbody>
</table>

</div>

```r
contrast(marginal, 
         method = "pairwise", 
         adjust = "tukey")
```

```
##  contrast                          estimate      SE df t.ratio p.value
##  Cambisol Winter - Fluvisol Winter  0.00333 0.00676 26  0.493  1.0000 
##  Cambisol Winter - Luvisol Winter   0.00667 0.00676 26  0.986  0.9967 
##  Cambisol Winter - Cambisol Spring  0.00667 0.00676 26  0.986  0.9967 
##  Cambisol Winter - Fluvisol Spring  0.01000 0.00633 26  1.581  0.9010 
##  Cambisol Winter - Luvisol Spring   0.01000 0.00633 26  1.581  0.9010 
##  Cambisol Winter - Cambisol Summer  0.00333 0.00676 26  0.493  1.0000 
##  Cambisol Winter - Fluvisol Summer  0.01000 0.00633 26  1.581  0.9010 
##  Cambisol Winter - Luvisol Summer   0.00000 0.00676 26  0.000  1.0000 
##  Cambisol Winter - Cambisol Autumn  0.01000 0.00676 26  1.479  0.9336 
##  Cambisol Winter - Fluvisol Autumn -0.00500 0.00756 26 -0.661  0.9999 
##  Cambisol Winter - Luvisol Autumn   0.00333 0.00676 26  0.493  1.0000 
##  Fluvisol Winter - Luvisol Winter   0.00333 0.00676 26  0.493  1.0000 
##  Fluvisol Winter - Cambisol Spring  0.00333 0.00676 26  0.493  1.0000 
##  Fluvisol Winter - Fluvisol Spring  0.00667 0.00633 26  1.054  0.9943 
##  Fluvisol Winter - Luvisol Spring   0.00667 0.00633 26  1.054  0.9943 
##  Fluvisol Winter - Cambisol Summer  0.00000 0.00676 26  0.000  1.0000 
##  Fluvisol Winter - Fluvisol Summer  0.00667 0.00633 26  1.054  0.9943 
##  Fluvisol Winter - Luvisol Summer  -0.00333 0.00676 26 -0.493  1.0000 
##  Fluvisol Winter - Cambisol Autumn  0.00667 0.00676 26  0.986  0.9967 
##  Fluvisol Winter - Fluvisol Autumn -0.00833 0.00756 26 -1.102  0.9918 
##  Fluvisol Winter - Luvisol Autumn   0.00000 0.00676 26  0.000  1.0000 
##  Luvisol Winter - Cambisol Spring   0.00000 0.00676 26  0.000  1.0000 
##  Luvisol Winter - Fluvisol Spring   0.00333 0.00633 26  0.527  1.0000 
##  Luvisol Winter - Luvisol Spring    0.00333 0.00633 26  0.527  1.0000 
##  Luvisol Winter - Cambisol Summer  -0.00333 0.00676 26 -0.493  1.0000 
##  Luvisol Winter - Fluvisol Summer   0.00333 0.00633 26  0.527  1.0000 
##  Luvisol Winter - Luvisol Summer   -0.00667 0.00676 26 -0.986  0.9967 
##  Luvisol Winter - Cambisol Autumn   0.00333 0.00676 26  0.493  1.0000 
##  Luvisol Winter - Fluvisol Autumn  -0.01167 0.00756 26 -1.543  0.9140 
##  Luvisol Winter - Luvisol Autumn   -0.00333 0.00676 26 -0.493  1.0000 
##  Cambisol Spring - Fluvisol Spring  0.00333 0.00633 26  0.527  1.0000 
##  Cambisol Spring - Luvisol Spring   0.00333 0.00633 26  0.527  1.0000 
##  Cambisol Spring - Cambisol Summer -0.00333 0.00676 26 -0.493  1.0000 
##  Cambisol Spring - Fluvisol Summer  0.00333 0.00633 26  0.527  1.0000 
##  Cambisol Spring - Luvisol Summer  -0.00667 0.00676 26 -0.986  0.9967 
##  Cambisol Spring - Cambisol Autumn  0.00333 0.00676 26  0.493  1.0000 
##  Cambisol Spring - Fluvisol Autumn -0.01167 0.00756 26 -1.543  0.9140 
##  Cambisol Spring - Luvisol Autumn  -0.00333 0.00676 26 -0.493  1.0000 
##  Fluvisol Spring - Luvisol Spring   0.00000 0.00586 26  0.000  1.0000 
##  Fluvisol Spring - Cambisol Summer -0.00667 0.00633 26 -1.054  0.9943 
##  Fluvisol Spring - Fluvisol Summer  0.00000 0.00586 26  0.000  1.0000 
##  Fluvisol Spring - Luvisol Summer  -0.01000 0.00633 26 -1.581  0.9010 
##  Fluvisol Spring - Cambisol Autumn  0.00000 0.00633 26  0.000  1.0000 
##  Fluvisol Spring - Fluvisol Autumn -0.01500 0.00717 26 -2.091  0.6325 
##  Fluvisol Spring - Luvisol Autumn  -0.00667 0.00633 26 -1.054  0.9943 
##  Luvisol Spring - Cambisol Summer  -0.00667 0.00633 26 -1.054  0.9943 
##  Luvisol Spring - Fluvisol Summer   0.00000 0.00586 26  0.000  1.0000 
##  Luvisol Spring - Luvisol Summer   -0.01000 0.00633 26 -1.581  0.9010 
##  Luvisol Spring - Cambisol Autumn   0.00000 0.00633 26  0.000  1.0000 
##  Luvisol Spring - Fluvisol Autumn  -0.01500 0.00717 26 -2.091  0.6325 
##  Luvisol Spring - Luvisol Autumn   -0.00667 0.00633 26 -1.054  0.9943 
##  Cambisol Summer - Fluvisol Summer  0.00667 0.00633 26  1.054  0.9943 
##  Cambisol Summer - Luvisol Summer  -0.00333 0.00676 26 -0.493  1.0000 
##  Cambisol Summer - Cambisol Autumn  0.00667 0.00676 26  0.986  0.9967 
##  Cambisol Summer - Fluvisol Autumn -0.00833 0.00756 26 -1.102  0.9918 
##  Cambisol Summer - Luvisol Autumn   0.00000 0.00676 26  0.000  1.0000 
##  Fluvisol Summer - Luvisol Summer  -0.01000 0.00633 26 -1.581  0.9010 
##  Fluvisol Summer - Cambisol Autumn  0.00000 0.00633 26  0.000  1.0000 
##  Fluvisol Summer - Fluvisol Autumn -0.01500 0.00717 26 -2.091  0.6325 
##  Fluvisol Summer - Luvisol Autumn  -0.00667 0.00633 26 -1.054  0.9943 
##  Luvisol Summer - Cambisol Autumn   0.01000 0.00676 26  1.479  0.9336 
##  Luvisol Summer - Fluvisol Autumn  -0.00500 0.00756 26 -0.661  0.9999 
##  Luvisol Summer - Luvisol Autumn    0.00333 0.00676 26  0.493  1.0000 
##  Cambisol Autumn - Fluvisol Autumn -0.01500 0.00756 26 -1.984  0.6993 
##  Cambisol Autumn - Luvisol Autumn  -0.00667 0.00676 26 -0.986  0.9967 
##  Fluvisol Autumn - Luvisol Autumn   0.00833 0.00756 26  1.102  0.9918 
## 
## P value adjustment: tukey method for comparing a family of 12 estimates
```

```r
(k_pairwise <- cld(marginal,
                      alpha = 0.05,
                      Letters = letters,
                      adjust = "tukey")) # works with lm but not with two-factor ART
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Soil </th>
   <th style="text-align:left;"> Season </th>
   <th style="text-align:right;"> emmean </th>
   <th style="text-align:right;"> SE </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> lower.CL </th>
   <th style="text-align:right;"> upper.CL </th>
   <th style="text-align:left;"> .group </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0041409 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> -0.0029681 </td>
   <td style="text-align:right;"> 0.0229681 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0041409 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> -0.0029681 </td>
   <td style="text-align:right;"> 0.0229681 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0041409 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> -0.0029681 </td>
   <td style="text-align:right;"> 0.0229681 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.0100000 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> -0.0049742 </td>
   <td style="text-align:right;"> 0.0249742 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Spring </td>
   <td style="text-align:right;"> 0.0133333 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> -0.0016409 </td>
   <td style="text-align:right;"> 0.0283076 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.0133333 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> -0.0016409 </td>
   <td style="text-align:right;"> 0.0283076 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.0166667 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0016924 </td>
   <td style="text-align:right;"> 0.0316409 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.0166667 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0016924 </td>
   <td style="text-align:right;"> 0.0316409 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.0166667 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0016924 </td>
   <td style="text-align:right;"> 0.0316409 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> Luvisol </td>
   <td style="text-align:left;"> Summer </td>
   <td style="text-align:right;"> 0.0200000 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0050258 </td>
   <td style="text-align:right;"> 0.0349742 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Cambisol </td>
   <td style="text-align:left;"> Winter </td>
   <td style="text-align:right;"> 0.0200000 </td>
   <td style="text-align:right;"> 0.0047816 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0050258 </td>
   <td style="text-align:right;"> 0.0349742 </td>
   <td style="text-align:left;"> a </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 11 </td>
   <td style="text-align:left;"> Fluvisol </td>
   <td style="text-align:left;"> Autumn </td>
   <td style="text-align:right;"> 0.0250000 </td>
   <td style="text-align:right;"> 0.0058562 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0066604 </td>
   <td style="text-align:right;"> 0.0433396 </td>
   <td style="text-align:left;"> a </td>
  </tr>
</tbody>
</table>

</div>

```r
(mod_k %>% 
  anova() %>% 
  mutate(`Part Eta Sq`=`Sum Sq`/sum(`Sum Sq`) ) ->
  mod_k_ANOVA)
```

<div class="kable-table">

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Sum Sq </th>
   <th style="text-align:right;"> Mean Sq </th>
   <th style="text-align:right;"> F value </th>
   <th style="text-align:right;"> Pr(&gt;F) </th>
   <th style="text-align:right;"> Part Eta Sq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Soil </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0000087 </td>
   <td style="text-align:right;"> 4.40e-06 </td>
   <td style="text-align:right;"> 0.063453 </td>
   <td style="text-align:right;"> 0.9386631 </td>
   <td style="text-align:right;"> 0.0034277 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Season </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.0002069 </td>
   <td style="text-align:right;"> 6.90e-05 </td>
   <td style="text-align:right;"> 1.005369 </td>
   <td style="text-align:right;"> 0.4061663 </td>
   <td style="text-align:right;"> 0.0814634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Soil:Season </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.0005406 </td>
   <td style="text-align:right;"> 9.01e-05 </td>
   <td style="text-align:right;"> 1.313515 </td>
   <td style="text-align:right;"> 0.2863652 </td>
   <td style="text-align:right;"> 0.2128637 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 0.0017833 </td>
   <td style="text-align:right;"> 6.86e-05 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.7022453 </td>
  </tr>
</tbody>
</table>

</div>

```r
# pwpp(marginal) # Pairwise P-value plot. Fails for unbalanced design
emmip(mod_k, Soil ~ Season)
```

![](TBI_figures/test diffs in k-2.svg)<!-- -->

### Plot TBI factors

```r
TBI %>% 
  pivot_longer(cols = c("S", "k"), 
               names_to = c("Factor"), 
               values_to = c("Value")) %>% 
  mutate(across(Factor, ~dplyr::recode(., S = "Stabilisation factor (S)", k = "Decomposition rate (k)"))) ->
  TBI2plot
  
p_alpha <- ggplot() +
  geom_violin(data = TBI2plot,
             aes(
               x = Soil,
               y = Value,
               # ymin = lerr,
               # ymax = herr
             ), colour = "grey",
              fill = "grey",
              alpha = 1 / 3) +
  geom_jitter(data = TBI2plot,
               aes(x = Soil,
               y = Value,
               # ymin = lerr,
               # ymax = herr,
               colour = Soil), size = 3, width = 0.2, alpha = 2/3) +
  # scale_colour_manual(values = Gradient.colours[c(5, 6, 11)], name = "") +
  scale_colour_locuszoom(name = "") +
  scale_fill_locuszoom(name = "") +
  theme(legend.position = "none") +
  # geom_errorbar(alpha = 1 / 2, width = 0.3) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(
    angle = 45,
    vjust = 0.9,
    hjust = 1
  )) +
  facet_grid(Factor ~ Season, scale = "free") +
  theme(strip.text = element_text(size = f_size - 4)) +
  background_grid(major = "y",
                  minor = "none",
                  size.major = 0.8) 

bind_rows(`Stabilisation factor (S)` = S_pairwise, 
          `Decomposition rate (k)` = k_pairwise, 
          .id = "Factor") %>% 
  bind_cols(., y = rep(c(0.58, 0.048), each = 12)) %>% 
  mutate(across(Soil, ~fct_inorder(.x))) ->
  dat_text

p_alpha <- p_alpha + geom_text(
  data = dat_text,
  mapping = aes(x = Soil, y = y, label = .group),
  nudge_x = 0,
  nudge_y = 0
)
print(p_alpha)
```

![](TBI_figures/plot alpha-1.svg)<!-- -->



```r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
 )
```

<details open>
<summary> <span title='Click to Expand'> Current session info </span> </summary>

```r

─ Session info ─────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.3 (2020-10-10)
 os       Ubuntu 18.04.5 LTS          
 system   x86_64, linux-gnu           
 ui       X11                         
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Prague               
 date     2021-05-19                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package      * version date       lib source        
 abind          1.4-5   2016-07-21 [1] CRAN (R 4.0.2)
 assertthat     0.2.1   2019-03-21 [1] CRAN (R 4.0.2)
 backports      1.2.1   2020-12-09 [1] CRAN (R 4.0.2)
 broom          0.7.6   2021-04-05 [1] CRAN (R 4.0.3)
 bslib          0.2.5   2021-05-12 [1] CRAN (R 4.0.3)
 car          * 3.0-10  2020-09-29 [1] CRAN (R 4.0.2)
 carData      * 3.0-4   2020-05-22 [1] CRAN (R 4.0.2)
 cellranger     1.1.0   2016-07-27 [1] CRAN (R 4.0.2)
 cli            2.5.0   2021-04-26 [1] CRAN (R 4.0.3)
 clipr          0.7.1   2020-10-08 [1] CRAN (R 4.0.2)
 coda           0.19-4  2020-09-30 [1] CRAN (R 4.0.2)
 codetools      0.2-18  2020-11-04 [1] CRAN (R 4.0.2)
 colorspace     2.0-1   2021-05-04 [1] CRAN (R 4.0.3)
 cowplot      * 1.1.1   2020-12-30 [1] CRAN (R 4.0.2)
 crayon         1.4.1   2021-02-08 [1] CRAN (R 4.0.3)
 curl           4.3.1   2021-04-30 [1] CRAN (R 4.0.3)
 data.table     1.14.0  2021-02-21 [1] CRAN (R 4.0.3)
 DBI            1.1.1   2021-01-15 [1] CRAN (R 4.0.3)
 dbplyr         2.1.1   2021-04-06 [1] CRAN (R 4.0.3)
 desc           1.3.0   2021-03-05 [1] CRAN (R 4.0.3)
 details        0.2.1   2020-01-12 [1] CRAN (R 4.0.2)
 digest         0.6.27  2020-10-24 [1] CRAN (R 4.0.2)
 dplyr        * 1.0.6   2021-05-05 [1] CRAN (R 4.0.3)
 ellipsis       0.3.2   2021-04-29 [1] CRAN (R 4.0.3)
 emmeans      * 1.6.0   2021-04-24 [1] CRAN (R 4.0.3)
 estimability   1.3     2018-02-11 [1] CRAN (R 4.0.2)
 evaluate       0.14    2019-05-28 [1] CRAN (R 4.0.2)
 extrafont    * 0.17    2014-12-08 [1] CRAN (R 4.0.2)
 extrafontdb    1.0     2012-06-11 [1] CRAN (R 4.0.2)
 fansi          0.4.2   2021-01-15 [1] CRAN (R 4.0.3)
 farver         2.1.0   2021-02-28 [1] CRAN (R 4.0.3)
 forcats      * 0.5.1   2021-01-27 [1] CRAN (R 4.0.3)
 foreign        0.8-81  2020-12-22 [4] CRAN (R 4.0.3)
 fs             1.5.0   2020-07-31 [1] CRAN (R 4.0.2)
 generics       0.1.0   2020-10-31 [1] CRAN (R 4.0.2)
 ggplot2      * 3.3.3   2020-12-30 [1] CRAN (R 4.0.2)
 ggsci        * 2.9     2018-05-14 [1] CRAN (R 4.0.2)
 glue           1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
 gtable         0.3.0   2019-03-25 [1] CRAN (R 4.0.2)
 haven          2.4.1   2021-04-23 [1] CRAN (R 4.0.3)
 highr          0.9     2021-04-16 [1] CRAN (R 4.0.3)
 hms            1.1.0   2021-05-17 [1] CRAN (R 4.0.3)
 htmltools      0.5.1.1 2021-01-22 [1] CRAN (R 4.0.3)
 httr           1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
 jquerylib      0.1.4   2021-04-26 [1] CRAN (R 4.0.3)
 jsonlite       1.7.2   2020-12-09 [1] CRAN (R 4.0.2)
 knitr        * 1.33    2021-04-24 [1] CRAN (R 4.0.3)
 labeling       0.4.2   2020-10-20 [1] CRAN (R 4.0.2)
 lattice        0.20-44 2021-05-02 [1] CRAN (R 4.0.3)
 lifecycle      1.0.0   2021-02-15 [1] CRAN (R 4.0.3)
 lubridate      1.7.10  2021-02-26 [1] CRAN (R 4.0.3)
 magrittr       2.0.1   2020-11-17 [1] CRAN (R 4.0.2)
 MASS         * 7.3-54  2021-05-03 [1] CRAN (R 4.0.3)
 Matrix         1.3-3   2021-05-04 [1] CRAN (R 4.0.3)
 modelr         0.1.8   2020-05-19 [1] CRAN (R 4.0.2)
 multcomp     * 1.4-17  2021-04-29 [1] CRAN (R 4.0.3)
 munsell        0.5.0   2018-06-12 [1] CRAN (R 4.0.2)
 mvtnorm      * 1.1-1   2020-06-09 [1] CRAN (R 4.0.2)
 openxlsx       4.2.3   2020-10-27 [1] CRAN (R 4.0.2)
 pillar         1.6.1   2021-05-16 [1] CRAN (R 4.0.3)
 pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.0.2)
 png            0.1-7   2013-12-03 [1] CRAN (R 4.0.2)
 purrr        * 0.3.4   2020-04-17 [1] CRAN (R 4.0.2)
 R6             2.5.0   2020-10-28 [1] CRAN (R 4.0.2)
 Rcpp           1.0.6   2021-01-15 [1] CRAN (R 4.0.3)
 readr        * 1.4.0   2020-10-05 [1] CRAN (R 4.0.2)
 readxl         1.3.1   2019-03-13 [1] CRAN (R 4.0.2)
 reprex         2.0.0   2021-04-02 [1] CRAN (R 4.0.3)
 rio            0.5.26  2021-03-01 [1] CRAN (R 4.0.3)
 rlang          0.4.11  2021-04-30 [1] CRAN (R 4.0.3)
 rmarkdown      2.8     2021-05-07 [1] CRAN (R 4.0.3)
 rprojroot      2.0.2   2020-11-15 [1] CRAN (R 4.0.2)
 rstudioapi     0.13    2020-11-12 [1] CRAN (R 4.0.2)
 Rttf2pt1       1.3.8   2020-01-10 [1] CRAN (R 4.0.2)
 rvest          1.0.0   2021-03-09 [1] CRAN (R 4.0.3)
 sandwich       3.0-0   2020-10-02 [1] CRAN (R 4.0.2)
 sass           0.4.0   2021-05-12 [1] CRAN (R 4.0.3)
 scales         1.1.1   2020-05-11 [1] CRAN (R 4.0.2)
 sessioninfo    1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
 stringi        1.6.2   2021-05-17 [1] CRAN (R 4.0.3)
 stringr      * 1.4.0   2019-02-10 [1] CRAN (R 4.0.2)
 survival     * 3.2-11  2021-04-26 [1] CRAN (R 4.0.3)
 svglite        2.0.0   2021-02-20 [1] CRAN (R 4.0.3)
 systemfonts    1.0.2   2021-05-11 [1] CRAN (R 4.0.3)
 TH.data      * 1.0-10  2019-01-21 [1] CRAN (R 4.0.2)
 tibble       * 3.1.2   2021-05-16 [1] CRAN (R 4.0.3)
 tidyr        * 1.1.3   2021-03-03 [1] CRAN (R 4.0.3)
 tidyselect     1.1.1   2021-04-30 [1] CRAN (R 4.0.3)
 tidyverse    * 1.3.1   2021-04-15 [1] CRAN (R 4.0.3)
 utf8           1.2.1   2021-03-12 [1] CRAN (R 4.0.3)
 vctrs          0.3.8   2021-04-29 [1] CRAN (R 4.0.3)
 withr          2.4.2   2021-04-18 [1] CRAN (R 4.0.3)
 xfun           0.23    2021-05-15 [1] CRAN (R 4.0.3)
 xml2           1.3.2   2020-04-23 [1] CRAN (R 4.0.2)
 xtable         1.8-4   2019-04-21 [1] CRAN (R 4.0.2)
 yaml           2.2.1   2020-02-01 [1] CRAN (R 4.0.2)
 zip            2.1.1   2020-08-27 [1] CRAN (R 4.0.2)
 zoo            1.8-9   2021-03-09 [1] CRAN (R 4.0.3)

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library

```

</details>
<br>

## References


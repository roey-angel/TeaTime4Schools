TeaTime4schools: Joint analysis - fungi
================
Roey Angel
2021-05-05

-   [Beta diversity analysis](#beta-diversity-analysis)
    -   [Load phyloseq object](#load-phyloseq-object)
    -   [Standardize abundances to the median sequencing depth (and
        convert to
        proportion)](#standardize-abundances-to-the-median-sequencing-depth-and-convert-to-proportion)
    -   [Variance partitioning models and
        ordinations](#variance-partitioning-models-and-ordinations)
        -   [Partitioning the data using discrete
            distance](#partitioning-the-data-using-discrete-distance)
            -   [Calculate ordinations](#calculate-ordinations)
    -   [Variance partitioning models and ordinations - soils
        excluded](#variance-partitioning-models-and-ordinations---soils-excluded)
        -   [Partitioning the data using discrete
            distance](#partitioning-the-data-using-discrete-distance-1)
            -   [Calculate ordinations](#calculate-ordinations-1)
        -   [Test differences between samples on the phylum and order
            levels](#test-differences-between-samples-on-the-phylum-and-order-levels)
    -   [Differential abundance models](#differential-abundance-models)
-   [References](#references)

[roey.angel@bc.cas.cz](mailto:%20roey.angel@bc.cas.cz)

## Beta diversity analysis

This analysis explores the alpha-diversity ditribution patters in the
different samples, based on the DADA2-produced sequences. \#\#\# Setting
general parameters:

``` r
set.seed(1000)
min_lib_size <- 5000
data_path <- "./DADA2_pseudo/"
Ps_file <- "TeaTime4Schools_ITS_filt3.RDS"
Proj_name <- "TeaTime4Schools"
```

### Load phyloseq object

This phyloseq object was created in
[05\_Taxonomical\_analysis.html](05_Taxonomical_analysis.html). The
Ps\_obj\_filt object excludes contaminants or unknown but still includes
taxa with low prevalence

``` r
readRDS(file = paste0(data_path, Ps_file)) %>%
  subset_samples(., Field != "Unburied") %>% # drop unburied samples
  prune_samples(sample_sums(.) > min_lib_size, .) %>% # remove samples  < min_lib_size
  filter_taxa(., function(x) sum(x) > 0, TRUE) -> # drop taxa with 0 abundance
  Ps_obj_filt
```

``` r
qplot(rowSums(otu_table(Ps_obj_filt)), geom = "histogram") + 
  xlab("Library size")
```

![](07_Beta_figures/unnamed-chunk-1-1.svg)<!-- -->

``` r
qplot(log10(rowSums(otu_table(Ps_obj_filt)))) +
  xlab("Logged library size")
```

![](07_Beta_figures/unnamed-chunk-1-2.svg)<!-- -->

### Standardize abundances to the median sequencing depth (and convert to proportion)

``` r
adonis(
  otu_table(Ps_obj_filt) ~ Lib.size,
  data = get_variable(Ps_obj_filt),
  method = "horn",
  permutations = 9999
)
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt) ~ Lib.size, data = get_variable(Ps_obj_filt),      permutations = 9999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Lib.size    1      3.57  3.5697  10.438 0.08321  1e-04 ***
    ## Residuals 115     39.33  0.3420         0.91679           
    ## Total     116     42.90                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
Ps_obj_filt %>%
  otu_table(.) %>%
  as(., "matrix") %>%
  rowSums() %>% 
  median() ->
  total
standf = function(x, t = total) round(t * (x / sum(x)))
Ps_obj_filt_median <- transform_sample_counts(Ps_obj_filt, standf)  # Standardize abundances to median sequencing depth

Ps_obj_filt_median_rel <- transform_sample_counts(Ps_obj_filt_median, function(x) x / sum(x)) # convert to relative abundance (just incase it's explicitly needed)

sample_data(Ps_obj_filt_median)$Lib.size <- sample_sums(Ps_obj_filt_median)

qplot(rowSums(otu_table(Ps_obj_filt_median)), geom = "histogram") + 
  xlab("Library size")
```

![](07_Beta_figures/median-1.svg)<!-- -->

``` r
PlotLibDist(Ps_obj_filt_median)
```

![](07_Beta_figures/median-2.svg)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_filt_median), "matrix"))
```

![](07_Beta_figures/median%20diag%20plots-1.svg)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_filt_median))) > 0)
meanSdPlot(as(log2(t(otu_table(Ps_obj_filt_median))[notAllZero, ] + 1), "matrix"))
```

![](07_Beta_figures/median%20diag%20plots-2.svg)<!-- -->

### Variance partitioning models and ordinations

#### Partitioning the data using discrete distance

``` r
adonis(
  otu_table(Ps_obj_filt_median) ~ Lib.size,
  data = get_variable(Ps_obj_filt_median),
  method = "horn",
  permutations = 9999
)
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_median) ~ Lib.size, data = get_variable(Ps_obj_filt_median),      permutations = 9999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## Lib.size    1     0.122 0.12239 0.32903 0.00285 0.9683
    ## Residuals 115    42.778 0.37198         0.99715       
    ## Total     116    42.900                 1.00000

``` r
adonis(
  otu_table(Ps_obj_filt_median_rel) ~ Field + Sample.type * Season,
  data = get_variable(Ps_obj_filt_median_rel),
  method = "horn",
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_median_rel) ~ Field +      Sample.type * Season, data = get_variable(Ps_obj_filt_median_rel),      permutations = 999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field                2     2.327  1.1636   7.969 0.05425  0.001 ***
    ## Sample.type          2     9.741  4.8706  33.356 0.22707  0.001 ***
    ## Season               3     9.463  3.1543  21.602 0.22058  0.001 ***
    ## Sample.type:Season   6     6.329  1.0548   7.224 0.14752  0.001 ***
    ## Residuals          103    15.040  0.1460         0.35058           
    ## Total              116    42.900                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# sample_data(Ps_obj_filt_median)$Type.season <- paste(sample_data(Ps_obj_filt)$Sample.type, sample_data(Ps_obj_filt)$Season) # too many comparisons

# Ps_obj_filt_median_s <- prune_taxa(names(sort(taxa_sums(Ps_obj_filt_median), TRUE)[1:100]), Ps_obj_filt_median) # for testing

# Compare sample-type pairs
mod_pairwise1 <- PairwiseAdonis(
  x = otu_table(Ps_obj_filt_median),
  factors = get_variable(Ps_obj_filt_median, "Sample.type"),
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod_pairwise1)
```

    ##                  pairs total.DF   F.Model         R2 p.value p.adjusted sig
    ## 1 Green tea vs Rooibos       86  7.952102 0.08555053   0.001      0.001  **
    ## 2    Green tea vs Soil       74 29.422941 0.28726905   0.001      0.001  **
    ## 3      Rooibos vs Soil       71 18.848044 0.21213797   0.001      0.001  **

``` r
(sig_pairs1 <- as.character(mod_pairwise1$pairs[mod_pairwise1$p.adjusted < 0.05]))
```

    ## [1] "Green tea vs Rooibos" "Green tea vs Soil"    "Rooibos vs Soil"

``` r
# (meaningful_sig_pairs1 <- c("Green tea vs Rooibos", "Green tea vs Soil", "Rooibos vs Soil"))

# compare seasons
mod_pairwise2 <- PairwiseAdonis(
  otu_table(Ps_obj_filt_median),
  sample_data(Ps_obj_filt_median)$Season,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod_pairwise2)
```

    ##              pairs total.DF   F.Model         R2 p.value p.adjusted sig
    ## 1 Winter vs Spring       55  5.013270 0.08495157   0.001     0.0012   *
    ## 2 Winter vs Summer       58 22.846552 0.28613073   0.001     0.0012   *
    ## 3 Winter vs Autumn       61 17.165072 0.22244613   0.001     0.0012   *
    ## 4 Spring vs Summer       54 14.048243 0.20952440   0.001     0.0012   *
    ## 5 Spring vs Autumn       57  9.216003 0.14131505   0.001     0.0012   *
    ## 6 Summer vs Autumn       60  2.448077 0.03983977   0.017     0.0170   .

``` r
(sig_pairs2 <- as.character(mod_pairwise2$pairs[mod_pairwise2$p.adjusted < 0.05]))
```

    ## [1] "Winter vs Spring" "Winter vs Summer" "Winter vs Autumn" "Spring vs Summer"
    ## [5] "Spring vs Autumn" "Summer vs Autumn"

``` r
# (meaningful_sig_pairs2 <- c("Green tea vs Rooibos", "Green tea vs Soil", "Rooibos vs Soil"))
```

##### Calculate ordinations

``` r
Ps_obj_ord1 <- ordinate(Ps_obj_filt_median, "CAP", "horn", formula = Ps_obj_filt_median ~  Field + Sample.type * Season)
evals <- eigenvals(Ps_obj_ord1) # /sum(eigenvals(Ps_obj_ord)) * 100

Ps_obj_filt_median %>% 
  plot_ordination(., Ps_obj_ord1, type = "samples", shape = "Field", color = "Sample.type", justDF = TRUE) %>% 
  mutate_at(., "Season", ~fct_relevel(., "Winter", "Spring", "Summer", "Autumn")) %>% 
  mutate_at(., "Sample.type", ~fct_relevel(., "Soil", "Green tea", "Rooibos")) %>% 
  dplyr::rename(., `Sample type` = Sample.type) ->
  ord_df

# p_ord.file <- "PCoA_bray"
# svglite(paste0(p_ord.file, ".svg"),
#         width = 10, height = 7.2)

p_ord <- ggplot(ord_df,
             aes(
               x = CAP1,
               y = CAP2,
               shape = Field,
               color = `Sample type`
             )) +
  geom_point(size = 4, alpha = 2 / 3) +
  theme_bw(base_size = 14) +
  scale_colour_manual(values = pal_locuszoom("default")(3)[c(2, 3, 1)]) +
  stat_ellipse(
    aes(x = CAP1, y = CAP2, group = `Sample type`),
    alpha = 0.5,
    type = "norm",
    level = 0.95,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  labs(x = sprintf("CAP1 [%s%%]", round(evals[1], 1)), 
       y = sprintf("CAP2 [%s%%]", round(evals[2], 1))) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  facet_wrap(~ Season)

# # Now add the environmental variables as arrows
# arrowmat <- scores(Ps_obj_ord1, display = "bp") # bipplot arrows
# # Add labels, make a data.frame
# arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# 
# arrowdf %<>%
#   mutate(labels = fct_recode(labels,
#     "Braunerde" = "FieldBraunerde",
#     "Kolluvisol" = "FieldKolluvisol",
#     "Rooibos" = "Sample.typeRooibos",
#     "Soil" = "Sample.typeSoil"
#   ))
# arrowdf %<>% dplyr::slice(., c(1:4))
# 
# # Define the arrow aesthetic mapping
# arrow_map = aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, 
#     label = labels)
# label_map = aes(x = 1.2 * CAP1, y = 1.2 * CAP2, shape = NULL, color = NULL, 
#     label = labels)
# # Make a new graphic
# arrowhead = arrow(length = unit(0.05, "npc"))
# p_ord <- p_ord +
#   geom_segment(
#     arrow_map,
#     size = 0.5,
#     data = arrowdf,
#     color = "gray",
#     arrow = arrowhead
#   ) + 
#   geom_text(label_map, size = 2, data = arrowdf)
print(p_ord)
```

<img src="07_Beta_figures/ordinate all-1.svg" height="40%" />

``` r
# dev.off()

# ggsave(
#   paste0(p_ord.file, ".png"),
#   p_ord,
#   device = "png",
#   width = 10,
#   height = 6
# )
# gz(paste0(p_ord.file, ".svg"), paste0(p_ord.file, ".svgz"))
```

### Variance partitioning models and ordinations - soils excluded

#### Partitioning the data using discrete distance

``` r
Ps_obj_filt_median_tea <- subset_samples(Ps_obj_filt_median, Sample.type != "Soil")
adonis(
  otu_table(Ps_obj_filt_median_tea) ~ Lib.size,
  data = get_variable(Ps_obj_filt_median_tea),
  method = "horn",
  permutations = 9999
)
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_median_tea) ~ Lib.size,      data = get_variable(Ps_obj_filt_median_tea), permutations = 9999,      method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## Lib.size   1     0.292 0.29248 0.78766 0.00918  0.575
    ## Residuals 85    31.563 0.37133         0.99082       
    ## Total     86    31.855                 1.00000

``` r
adonis(
  otu_table(Ps_obj_filt_median_tea) ~ Field + Sample.type * Season,
  data = get_variable(Ps_obj_filt_median_tea),
  method = "horn",
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = otu_table(Ps_obj_filt_median_tea) ~ Field +      Sample.type * Season, data = get_variable(Ps_obj_filt_median_tea),      permutations = 999, method = "horn") 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field               2     2.254  1.1270  7.4101 0.07076  0.001 ***
    ## Sample.type         1     2.608  2.6084 17.1497 0.08188  0.001 ***
    ## Season              3    11.925  3.9751 26.1360 0.37436  0.001 ***
    ## Sample.type:Season  3     3.356  1.1188  7.3558 0.10536  0.001 ***
    ## Residuals          77    11.711  0.1521         0.36764           
    ## Total              86    31.855                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# sample_data(Ps_obj_filt_median)$Type.season <- paste(sample_data(Ps_obj_filt)$Sample.type, sample_data(Ps_obj_filt)$Season) # too many comparisons

# Ps_obj_filt_median_s <- prune_taxa(names(sort(taxa_sums(Ps_obj_filt_median), TRUE)[1:100]), Ps_obj_filt_median) # for testing

# Compare sample-type pairs
mod_pairwise3 <- PairwiseAdonis(
  otu_table(Ps_obj_filt_median_tea),
  sample_data(Ps_obj_filt_median_tea)$Sample.type,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod_pairwise3)
```

    ##                  pairs total.DF  F.Model         R2 p.value p.adjusted sig
    ## 1 Green tea vs Rooibos       86 7.952102 0.08555053   0.001      0.001  **

``` r
(sig_pairs3 <- as.character(mod_pairwise3$pairs[mod_pairwise3$p.adjusted < 0.05]))
```

    ## [1] "Green tea vs Rooibos"

``` r
# (meaningful_sig_pairs3 <- c("Green tea vs Rooibos", "Green tea vs Soil", "Rooibos vs Soil"))

# compare seasons
mod_pairwise4 <- PairwiseAdonis(
  otu_table(Ps_obj_filt_median_tea),
  sample_data(Ps_obj_filt_median_tea)$Season,
  sim.function = "vegdist",
  sim.method = "horn",
  p.adjust.m = "BH"
)
print(mod_pairwise4)
```

    ##              pairs total.DF   F.Model        R2 p.value p.adjusted sig
    ## 1 Winter vs Spring       43  8.533529 0.1688687   0.001     0.0012   *
    ## 2 Winter vs Summer       43 33.801961 0.4459246   0.001     0.0012   *
    ## 3 Winter vs Autumn       46 24.599278 0.3534416   0.001     0.0012   *
    ## 4 Spring vs Summer       39 20.786718 0.3535955   0.001     0.0012   *
    ## 5 Spring vs Autumn       42 13.464088 0.2472104   0.001     0.0012   *
    ## 6 Summer vs Autumn       42  3.578983 0.0802841   0.004     0.0040   *

``` r
(sig_pairs4 <- as.character(mod_pairwise4$pairs[mod_pairwise4$p.adjusted < 0.05]))
```

    ## [1] "Winter vs Spring" "Winter vs Summer" "Winter vs Autumn" "Spring vs Summer"
    ## [5] "Spring vs Autumn" "Summer vs Autumn"

``` r
# (meaningful_sig_pairs2 <- c("Green tea vs Rooibos", "Green tea vs Soil", "Rooibos vs Soil"))
```

##### Calculate ordinations

``` r
Ps_obj_ord3 <- ordinate(Ps_obj_filt_median_tea, "CAP", "horn", formula = Ps_obj_filt_median_tea ~  Field + Sample.type * Season)
evals <- eigenvals(Ps_obj_ord3) # /sum(eigenvals(Ps_obj_ord)) * 100

Ps_obj_filt_median_tea %>% 
  plot_ordination(., Ps_obj_ord3, type = "samples", shape = "Field", color = "Sample.type", justDF = TRUE) %>% 
  mutate_at(., "Season", ~fct_relevel(., "Winter", "Spring", "Summer", "Autumn")) %>% 
  mutate_at(., "Sample.type", ~fct_relevel(., "Green tea", "Rooibos")) %>% 
  dplyr::rename(., `Sample type` = Sample.type) ->
  ord_df3

# p_ord.file <- "PCoA_bray"
# svglite(paste0(p_ord.file, ".svg"),
#         width = 10, height = 7.2)

p_ord3 <- ggplot(ord_df3,
             aes(
               x = CAP1,
               y = CAP2,
               shape = Field,
               color = `Sample type`
             )) +
  geom_point(size = 4, alpha = 2 / 3) +
  theme_bw(base_size = 14) +
  scale_colour_manual(values = pal_locuszoom("default")(3)[c(3, 1)]) +
  stat_ellipse(
    aes(x = CAP1, y = CAP2, group = `Sample type`),
    alpha = 0.5,
    type = "norm",
    level = 0.95,
    linetype = 2,
    inherit.aes = FALSE
  ) +
  labs(x = sprintf("Axis1 [%s%%]", round(evals[1], 1)), 
       y = sprintf("Axis2 [%s%%]", round(evals[2], 1))) +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  facet_wrap(~ Season)

print(p_ord3)
```

<img src="07_Beta_figures/ordinate tea-1.svg" height="40%" />

``` r
# dev.off()

# ggsave(
#   paste0(p_ord.file, ".png"),
#   p_ord,
#   device = "png",
#   width = 10,
#   height = 6
# )
# gz(paste0(p_ord.file, ".svg"), paste0(p_ord.file, ".svgz"))
```

#### Test differences between samples on the phylum and order levels

``` r
Taxa_tests_phylum1 <- STAMPR(Ps_obj_filt_median, 
                             vars2test = "Sample.type", 
                             sig_pairs = sig_pairs1, 
                             outputfile = paste0(Proj_name, "_Sample.type"))

plotSTAMPR(Taxa_tests_phylum1, pair = "Green tea vs Rooibos")
```

![](07_Beta_figures/STAMPR-1.svg)<!-- -->

``` r
plotSTAMPR(Taxa_tests_phylum1, pair = "Green tea vs Soil")
```

![](07_Beta_figures/STAMPR-2.svg)<!-- -->

``` r
plotSTAMPR(Taxa_tests_phylum1, pair = "Rooibos vs Soil")
```

![](07_Beta_figures/STAMPR-3.svg)<!-- -->

``` r
Taxa_tests_order1 <- STAMPR(Ps_obj_filt_median, vars2test = "Sample.type", rank = "Order", sig_pairs = sig_pairs1, outputfile = paste0(Proj_name, "_Sample.type"))
plotSTAMPR(Taxa_tests_order1, pair = "Green tea vs Rooibos", tax_level = "Order")
```

![](07_Beta_figures/STAMPR-4.svg)<!-- -->

``` r
plotSTAMPR(Taxa_tests_order1, pair = "Green tea vs Soil", tax_level = "Order")
```

![](07_Beta_figures/STAMPR-5.svg)<!-- -->

``` r
plotSTAMPR(Taxa_tests_order1, pair = "Rooibos vs Soil", tax_level = "Order")
```

![](07_Beta_figures/STAMPR-6.svg)<!-- -->

``` r
Taxa_tests_phylum2 <- STAMPR(Ps_obj_filt_median, vars2test = "Season", sig_pairs = sig_pairs2, outputfile = paste0(Proj_name, "_Season"))
plotSTAMPR(Taxa_tests_phylum2, pair = "Winter vs Summer")
```

![](07_Beta_figures/STAMPR-7.svg)<!-- -->

``` r
Taxa_tests_order2 <- STAMPR(Ps_obj_filt_median, vars2test = "Season", rank = "Order", sig_pairs = sig_pairs2, outputfile = paste0(Proj_name, "_Season"))
plotSTAMPR(Taxa_tests_order2, pair = "Winter vs Summer", tax_level = "Order")
```

![](07_Beta_figures/STAMPR-8.svg)<!-- -->

``` r
Taxa_tests_phylum4 <- STAMPR(Ps_obj_filt_median_tea, vars2test = "Season", sig_pairs = sig_pairs4, outputfile = paste0(Proj_name, "Teabags_Season"))
plotSTAMPR(Taxa_tests_phylum4, pair = "Winter vs Summer")
```

![](07_Beta_figures/STAMPR-9.svg)<!-- -->

``` r
Taxa_tests_order4 <- STAMPR(Ps_obj_filt_median_tea, vars2test = "Season", rank = "Order", sig_pairs = sig_pairs4, outputfile = paste0(Proj_name, "Teabags_Season"))
plotSTAMPR(Taxa_tests_order4, pair = "Winter vs Summer", tax_level = "Order")
```

![](07_Beta_figures/STAMPR-10.svg)<!-- -->

### Differential abundance models

Tag rare phyla (for plotting purposes only)

``` r
Ps_obj_filt_media_glom <- tax_glom(Ps_obj_filt_median, 
                             "Order", 
                             NArm = TRUE) # glomerate to the Order level
Ps_obj_filt_media_glom_rel <- transform_sample_counts(Ps_obj_filt_media_glom, function(x) x / sum(x)) # transform to rel. ab.
Ps_obj_filt_media_glom_rel_DF <- psmelt(Ps_obj_filt_media_glom_rel) # generate a df
Ps_obj_filt_media_glom_rel_DF$Order %<>% as.character() # factor to char

# group dataframe by Order, calculate median rel. abundance
Ps_obj_filt_media_glom_rel_DF %>%
  group_by(Order) %>%
  summarise(median = median(Abundance)) ->
  medians

# find Phyla whose median rel. abund. is less than 0.5%
Rare_phyla <- medians[medians$median <= 0.005, ]$Order

# change their name to "Rare"
Ps_obj_filt_media_glom_rel_DF[Ps_obj_filt_media_glom_rel_DF$Order %in% Rare_phyla, ]$Order <- 'Rare'

# re-group
Ps_obj_filt_media_glom_rel_DF %>%
  group_by(Order) %>%
  summarise(Abundance = sum(Abundance)) %>% 
  arrange(desc(Abundance)) -> Taxa_rank
```

Detect differentially abundant OTUs using ALDEx2 ([Fernandes *et al.*
2013](#ref-fernandes_anova-like_2013))

**Sample type differences**

``` r
significance = 0.05

# run full model 
data2test <- t(otu_table(Ps_obj_filt_median))
# comparison <- as.character(get_variable(Ps_obj_filt_median, "Sample.type"))
#ALDEx_full <- aldex.clr(data2test, comparison, mc.samples = 128, denom = "iqlr", verbose = TRUE, useMC = TRUE) # iqlr for slight assymetry in composition
#ALDEx_full_glm <- aldex.glm(ALDEx_full, comparison, useMC = TRUE) # for more than two conditions
#sig_taxa <- rownames(ALDEx_full_glm)[ALDEx_full_glm$glm.eBH < 0.05] # save names of taxa that are significant under the full model
#write.csv(sig_taxa, file = "Aldex_full_significant_taxa.csv")

# Pairwise comparisons
# 
ALDEx_comparisons <- list()
ALDEx_comparisons$Comparisons <- sig_pairs1

# Ps_obj_filt_median_subset <- prune_taxa(names(sort(taxa_sums(Ps_obj_filt_median), TRUE)[1:100]), Ps_obj_filt_median)

for (j in seq(1, length(ALDEx_comparisons$Comparisons))) {
  # print(j)
  ALDEx_comparisons$Comparisons[j] %>% 
    str_split(., " vs ", simplify = FALSE) %>% 
    unlist() ->
    comparison_string
  Ps_obj_filt_median %>%
    subset_samples(Sample.type == comparison_string[1] | Sample.type == comparison_string[2]) ->
    # subset_samples(Treatment == comparison_string[2] | Treatment == comparison_string[4]) %>%
    # subset_samples(Year == "2016" | Year == "2017" | Year == "2018") ->
    Ps_obj_filt_median_pairwise

#  Remove species with prevalence < 10%
  Ps_obj_filt_median_pairwise_s <- DropRareSpecies(Ps_obj = Ps_obj_filt_median_pairwise, prevalence = 0.1)
    
  sample_data(Ps_obj_filt_median_pairwise_s)$Sample.type %<>% fct_relevel(., str_split(ALDEx_comparisons$Comparisons[j], " vs ", simplify = TRUE)[1])
  
  # make Joint.sample.name for matching OTUs between compared samples (for GGPlotTopOTUs)
  suppressWarnings(
  sample_data(Ps_obj_filt_median_pairwise_s) %<>%
    as(., "data.frame") %>% 
    rownames_to_column() %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(Joint.sample.name  = paste0(.$Sample.type, "_", .$Field, "_", .$Season, "_", .$Replicate)) %>% 
    mutate_at("Joint.sample.name", ~str_replace_all(., paste0(unique(comparison_string), collapse = "|"), "")) %>% # remove the levels participating in the comparison from the names
    mutate_at("Joint.sample.name", ~str_replace_all(., "^_", "")) %>% # remove first "_"
    column_to_rownames()
  )
    
  ALDEx2plot_pairwise <- CalcALDEx(
    physeq_obj = Ps_obj_filt_median_pairwise_s,
    vars2test = "Sample.type",
    rare_phyla = Rare_phyla,
    sig_level = significance,
    LFC = 0
    )
  ALDEx_comparisons$Results[[j]] <- ALDEx2plot_pairwise # store results
  ALDEx_comparisons$Results[[j]]$Var1 <- str_split(ALDEx_comparisons$Comparisons[j], " vs ", simplify = TRUE)[1]
  ALDEx_comparisons$Results[[j]]$Var2 <- str_split(ALDEx_comparisons$Comparisons[j], " vs ", simplify = TRUE)[2]
  ALDEX_summary <- tibble(Label = c(paste0("⬆", sum( ALDEx_comparisons$Results[[j]]$effect > 0 &  ALDEx_comparisons$Results[[j]]$Significance == "Pass"), " ⬇", sum( ALDEx_comparisons$Results[[j]]$effect < 0 &  ALDEx_comparisons$Results[[j]]$Significance == "Pass"), " (", nrow( ALDEx_comparisons$Results[[j]]), ")")))
  
  ALDEx2plot_pairwise %>%
    filter(Significance == "Pass") %>%
    dplyr::select(OTU, baseMean, effect, Order, Class, Order, Family, Genus) %>%
    arrange(desc(abs(effect))) ->
    ALDEx2plot_pairwise_results
    
  # print(ALDEx2plot_pairwise_results %>% 
  #   kable(., digits = c(2), caption = "Significantly different taxa:") %>%
  #   kable_styling(
  #     bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  #     full_width = F
  #   ))
  
  write.csv(ALDEx2plot_pairwise_results, file = paste0("Aldex", "_", paste0(comparison_string, collapse = "_"), ".csv"))
  
  # Plot ALDEX plot
  p1 <- GGPlotALDExTax(ALDEx2plot_pairwise, OTU_labels = FALSE, Taxa = "Order", sig_level = significance) +
    # ggtitle(ALDEx_comparisons$Comparisons[j]) +
    geom_text(
    data    = ALDEX_summary,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  ) 
  p1 <- p1 + labs(title = ALDEx_comparisons$Comparisons[j])
  print(p1)
  
  # Plot OTU plots
  # GGPlotOTU(Ps_obj_filt_subset_pairwise_s, vars2test = "Spill.Treatment", "Seq_8")
  p2 <- GGPlotTopOTUs(
    physeq_obj = Ps_obj_filt_median_pairwise_s,
    vars2test = "Sample.type",
    ALDEx_obj = ALDEx2plot_pairwise,
    tax_level = "Order",
    rank_by = "effect",
    Ntop = 12
  )
  print(p2)
}
```

phyloseq-class experiment-level object otu\_table() OTU Table: \[ 402
taxa and 87 samples \] sample\_data() Sample Data: \[ 87 samples by 25
sample variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6
taxonomic ranks \] phyloseq-class experiment-level object otu\_table()
OTU Table: \[ 266 taxa and 87 samples \] sample\_data() Sample Data: \[
87 samples by 25 sample variables \] tax\_table() Taxonomy Table: \[ 266
taxa by 6 taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20sample.type-1.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20sample.type-2.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 75
samples \] sample\_data() Sample Data: \[ 75 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 275 taxa and 75 samples \] sample\_data() Sample Data: \[ 75 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 275 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20sample.type-3.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20sample.type-4.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 72
samples \] sample\_data() Sample Data: \[ 72 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 297 taxa and 72 samples \] sample\_data() Sample Data: \[ 72 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 297 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20sample.type-5.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20sample.type-6.svg)<!-- -->

``` r
filter(ALDEx_comparisons$Results[[1]], Significance == "Pass")$OTU[common_taxa <- filter(ALDEx_comparisons$Results[[1]], Significance == "Pass")$OTU %in% filter(ALDEx_comparisons$Results[[2]], Significance == "Pass")$OTU]
```

\[1\] “Seq\_41” “Seq\_42” “Seq\_78” “Seq\_204”

\#\#\#\#Seasonal differences Comparing the seasonsal differences in
teabag samples only

``` r
significance = 0.05

data2test <- t(otu_table(Ps_obj_filt_median_tea))

# Pairwise comparisons
ALDEx_comparisons <- list()
ALDEx_comparisons$Comparisons <- sig_pairs2

# Ps_obj_filt_median_tea_subset <- prune_taxa(names(sort(taxa_sums(Ps_obj_filt_median_tea), TRUE)[1:100]), Ps_obj_filt_median_tea)

for (j in seq(1, length(ALDEx_comparisons$Comparisons))) {
  # print(j)
  ALDEx_comparisons$Comparisons[j] %>% 
    str_split(., " vs ", simplify = FALSE) %>% 
    unlist() ->
    comparison_string
  Ps_obj_filt_median_tea %>%
    subset_samples(Season == comparison_string[1] | Season == comparison_string[2]) ->
    # subset_samples(Treatment == comparison_string[2] | Treatment == comparison_string[4]) %>%
    # subset_samples(Year == "2016" | Year == "2017" | Year == "2018") ->
    Ps_obj_filt_median_tea_pairwise

#  Remove species with prevalence < 10%
  Ps_obj_filt_median_tea_pairwise_s <- DropRareSpecies(Ps_obj = Ps_obj_filt_median_tea_pairwise, prevalence = 0.1)
    
  sample_data(Ps_obj_filt_median_tea_pairwise_s)$Season %<>% fct_relevel(., str_split(ALDEx_comparisons$Comparisons[j], " vs ", simplify = TRUE)[1])
  
  # make Joint.sample.name for matching OTUs between compared samples (for GGPlotTopOTUs)
  suppressWarnings(
  sample_data(Ps_obj_filt_median_tea_pairwise_s) %<>% 
    as(., "data.frame") %>% 
    rownames_to_column() %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(Joint.sample.name  = paste0(.$Sample.type, "_", .$Field, "_", .$Season, "_", .$Replicate)) %>% 
    mutate_at("Joint.sample.name", ~str_replace_all(., paste0(unique(comparison_string), collapse = "|"), "")) %>% # remove the levels participating in the comparison from the names
    mutate_at("Joint.sample.name", ~str_replace_all(., "^_", "")) %>% # remove first "_"
    mutate_at("Joint.sample.name", ~str_replace_all(., "__", "_")) %>% # remove double "_"
    column_to_rownames()
  )
    
  ALDEx2plot_pairwise <- CalcALDEx(
    physeq_obj = Ps_obj_filt_median_tea_pairwise_s,
    vars2test = "Season",
    rare_phyla = Rare_phyla,
    sig_level = significance,
    LFC = 0
    )
  ALDEx_comparisons$Results[[j]] <- ALDEx2plot_pairwise # store results
  ALDEx_comparisons$Results[[j]]$Var1 <- str_split(ALDEx_comparisons$Comparisons[j], " vs ", simplify = TRUE)[1]
  ALDEx_comparisons$Results[[j]]$Var2 <- str_split(ALDEx_comparisons$Comparisons[j], " vs ", simplify = TRUE)[2]
  ALDEX_summary <- tibble(Label = c(paste0("⬆", sum( ALDEx_comparisons$Results[[j]]$effect > 0 &  ALDEx_comparisons$Results[[j]]$Significance == "Pass"), " ⬇", sum( ALDEx_comparisons$Results[[j]]$effect < 0 &  ALDEx_comparisons$Results[[j]]$Significance == "Pass"), " (", nrow( ALDEx_comparisons$Results[[j]]), ")")))
  
  ALDEx2plot_pairwise %>%
    filter(Significance == "Pass") %>%
    dplyr::select(OTU, baseMean, effect, Order, Class, Order, Family, Genus) %>%
    arrange(desc(abs(effect))) ->
    ALDEx2plot_pairwise_results
    
  # print(ALDEx2plot_pairwise_results %>% 
  #   kable(., digits = c(2), caption = "Significantly different taxa:") %>%
  #   kable_styling(
  #     bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  #     full_width = F
  #   ))
  
  write.csv(ALDEx2plot_pairwise_results, file = paste0("Aldex", "_", paste0(comparison_string, collapse = "_"), ".csv"))
  
  # Plot ALDEX plot
  p1 <- GGPlotALDExTax(ALDEx2plot_pairwise, OTU_labels = FALSE, sig_level = significance) +
    # ggtitle(ALDEx_comparisons$Comparisons[j]) +
    geom_text(
    data    = ALDEX_summary,
    mapping = aes(x = Inf, y = Inf, label = Label),
    hjust   = 1.1,
    vjust   = 1.6
  ) 
  p1 <- p1 + labs(title = ALDEx_comparisons$Comparisons[j])
  print(p1)
  
  # Plot OTU plots
  p2 <- GGPlotTopOTUs(
    physeq_obj = Ps_obj_filt_median_tea_pairwise_s,
    vars2test = "Season",
    ALDEx_obj = ALDEx2plot_pairwise,
    tax_level = "Order",
    rank_by = "effect",
    Ntop = 12
  )
  print(p2)
}
```

phyloseq-class experiment-level object otu\_table() OTU Table: \[ 402
taxa and 44 samples \] sample\_data() Sample Data: \[ 44 samples by 25
sample variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6
taxonomic ranks \] phyloseq-class experiment-level object otu\_table()
OTU Table: \[ 262 taxa and 44 samples \] sample\_data() Sample Data: \[
44 samples by 25 sample variables \] tax\_table() Taxonomy Table: \[ 262
taxa by 6 taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20season-1.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20season-2.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 44
samples \] sample\_data() Sample Data: \[ 44 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 293 taxa and 44 samples \] sample\_data() Sample Data: \[ 44 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 293 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20season-3.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20season-4.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 47
samples \] sample\_data() Sample Data: \[ 47 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 316 taxa and 47 samples \] sample\_data() Sample Data: \[ 47 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 316 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20season-5.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20season-6.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 40
samples \] sample\_data() Sample Data: \[ 40 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 196 taxa and 40 samples \] sample\_data() Sample Data: \[ 40 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 196 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20season-7.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20season-8.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 43
samples \] sample\_data() Sample Data: \[ 43 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 198 taxa and 43 samples \] sample\_data() Sample Data: \[ 43 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 198 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20season-9.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20season-10.svg)<!-- -->phyloseq-class
experiment-level object otu\_table() OTU Table: \[ 402 taxa and 43
samples \] sample\_data() Sample Data: \[ 43 samples by 25 sample
variables \] tax\_table() Taxonomy Table: \[ 402 taxa by 6 taxonomic
ranks \] phyloseq-class experiment-level object otu\_table() OTU Table:
\[ 186 taxa and 43 samples \] sample\_data() Sample Data: \[ 43 samples
by 25 sample variables \] tax\_table() Taxonomy Table: \[ 186 taxa by 6
taxonomic ranks \]
![](07_Beta_figures/ALDEx2%20-%20season-11.svg)<!-- -->![](07_Beta_figures/ALDEx2%20-%20season-12.svg)<!-- -->

``` r
filter(ALDEx_comparisons$Results[[1]], Significance == "Pass")$OTU[common_taxa <- filter(ALDEx_comparisons$Results[[1]], Significance == "Pass")$OTU %in% filter(ALDEx_comparisons$Results[[2]], Significance == "Pass")$OTU]
```

\[1\] “Seq\_3” “Seq\_17” “Seq\_20” “Seq\_27” “Seq\_31” “Seq\_33”
“Seq\_36” “Seq\_40” \[9\] “Seq\_43” “Seq\_44” “Seq\_50” “Seq\_55”
“Seq\_57” “Seq\_64” “Seq\_67” “Seq\_70” \[17\] “Seq\_72” “Seq\_74”
“Seq\_80” “Seq\_85” “Seq\_95” “Seq\_107” “Seq\_128” “Seq\_203” \[25\]
“Seq\_240”

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-fernandes_anova-like_2013" class="csl-entry">

Fernandes AD, Macklaim JM, Linn TG *et al.* ANOVA-Like Differential
Expression (ALDEx) Analysis for Mixed Population RNA-Seq. *PLOS ONE*
2013;**8**:e67019.

</div>

</div>

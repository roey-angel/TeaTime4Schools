TeaTime4Schools joint analysis - bacteria
================
Roey Angel
2021-05-21

-   [Preliminary analysis of MiSeq amplicon data for Evrona
    project](#preliminary-analysis-of-miseq-amplicon-data-for-evrona-project)
    -   [Setting general parameters:](#setting-general-parameters)
    -   [Reading in raw data and generate phyloseq
        object](#reading-in-raw-data-and-generate-phyloseq-object)
    -   [Inspect Library Sizes](#inspect-library-sizes)
    -   [Identify contaminants -
        Frequency](#identify-contaminants---frequency)
    -   [Identify contaminants -
        Prevalence](#identify-contaminants---prevalence)
    -   [Save contaminant sequence
        names](#save-contaminant-sequence-names)
-   [References](#references)

[roey.angel@bc.cas.cz](mailto:%20roey.angel@bc.cas.cz)

## Preliminary analysis of MiSeq amplicon data for Evrona project

Decontamination of sequence library based on [Introduction to
decontam](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html)
and Davis and colleagues ([2017](#ref-davis_simple_2017)).
Decontamination is based on correlating sequence abundance frequencies
to initial DNA concentrations used for PCR and also on examining
sequence abundance prevalence patterns in the negative control samples.

### Setting general parameters:

``` r
set.seed(1000)
metadata_path <- "./"
data_path <- "./DADA2_pseudo/"
Metadata_table <- "TeaTime_joint_Bacteria_metadata.csv"
Seq_table <- "DADA2.seqtab_nochim.tsv"
# Metadata_table <- ""
```

### Reading in raw data and generate phyloseq object

``` r
# read OTU mat from data file
Raw_data <- read_tsv(paste0(data_path, Seq_table), 
                        trim_ws = TRUE)

Raw_data[, 2:ncol(Raw_data)] %>% 
  t() %>% 
  as.data.frame() -> abundance_mat # convert to abundance matrix
colnames(abundance_mat) <- pull(Raw_data, "#OTU") # add sequence names

# get short names of samples
# abundance_mat %>% 
#   rownames() %>% 
#   str_remove("^Roey[0-9]{3,}-?") %>% 
#   str_split("_", simplify = T) %>% 
#   .[, 1] ->
#   short_names

# Read metadata file
Metadata <- read_csv(paste0(metadata_path, Metadata_table),
                        trim_ws = TRUE)
Metadata$Control %<>% fct_relevel("Biological", after = 1)

Metadata$Read1_file <- str_replace(Metadata$Read1_file, "(.*)_L001_R1_001.fastq.gz|\\.1\\.fastq.gz", "\\1")
Metadata <- Metadata[Metadata$Read1_file %in% rownames(abundance_mat), ] # remove metadata rows if the samples did not go through qual processing
Metadata %<>% mutate(`Conc. (ng/ul)` = replace(`Conc. (ng/ul)`, which(`Conc. (ng/ul)` == 0 | `Conc. (ng/ul)` < 0), 0.001)) # replace 0 conc. with 0.001

# Order abundance_mat samples according to the metadata
sample_order <- match(rownames(abundance_mat), Metadata$Read1_file)
abundance_mat %<>% arrange(., sample_order)
rownames(abundance_mat) <- Metadata$Read1_file # needed for pyhloseq

Metadata$LibrarySize <- rowSums(abundance_mat)
Metadata <- data.frame(row.names = Metadata$Read1_file, Metadata)

# generate phyloseq object
Ps_obj <- phyloseq(otu_table(abundance_mat, taxa_are_rows = FALSE),
                        sample_data(Metadata)
                        )
Ps_obj <- filter_taxa(Ps_obj, function(x) sum(x) > 0, TRUE) # remove 0 abundance taxa
Ps_obj <- subset_samples(Ps_obj, sample_sums(Ps_obj) > 0) # remove 0 abundance samples
```

### Inspect Library Sizes

``` r
Ps_obj_df <- get_variable(Ps_obj) # Put sample_data into a ggplot-friendly data.frame
Ps_obj_df <- Ps_obj_df[order(Ps_obj_df$LibrarySize), ]
Ps_obj_df$Index <- seq(nrow(Ps_obj_df))
ggplot(data = Ps_obj_df, 
       aes(x = Index, y = LibrarySize, color = Control)) + 
  geom_point() +
  scale_y_log10(breaks = c(
    min(Ps_obj_df$Lib.size),
    10,
    100,
    1000,
    5000,
    10000,
    ceiling(max(Ps_obj_df$Lib.size) / 10000) * 10000
    ))
```

![](02_Decontamination_figures/Library%20Sizes-1.png)<!-- -->

``` r
summary(sample_sums(Ps_obj))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     182   27498   34260   33388   39484   71171

``` r
summary(taxa_sums(Ps_obj))
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       1       4      15     336      88  155005

No sample has 0 counts and only two controls have more sequences than
any true sample.

This is how many reads remained per control sample after DADA2 pipeline

``` r
Ps_obj %>% 
  subset_samples(., Control == "TRUE") %>% 
  sample_sums() %>% 
  sort() %>% 
  as_tibble(rownames = "Sample") %>% 
  rename(Reads = value) %>% 
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Sample
</th>
<th style="text-align:right;">
Reads
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Angel022-blank-Extr\_S22
</td>
<td style="text-align:right;">
182
</td>
</tr>
<tr>
<td style="text-align:left;">
Angel034-ExtrCTRL\_S279
</td>
<td style="text-align:right;">
316
</td>
</tr>
<tr>
<td style="text-align:left;">
Eva6\_blank\_extr\_S215
</td>
<td style="text-align:right;">
345
</td>
</tr>
<tr>
<td style="text-align:left;">
Angel023-NTC\_S23
</td>
<td style="text-align:right;">
371
</td>
</tr>
<tr>
<td style="text-align:left;">
Angel035-CTRL\_S280
</td>
<td style="text-align:right;">
941
</td>
</tr>
<tr>
<td style="text-align:left;">
Eva12\_NTC\_S221
</td>
<td style="text-align:right;">
1432
</td>
</tr>
<tr>
<td style="text-align:left;">
Eva24-NTC\_S206
</td>
<td style="text-align:right;">
9476
</td>
</tr>
<tr>
<td style="text-align:left;">
Re-Eva24-NTC\_S24
</td>
<td style="text-align:right;">
25687
</td>
</tr>
</tbody>
</table>

Summary of all nono-control samples

``` r
# Ps_obj <-
  # subset_samples(Ps_obj, Sample.name != names(which(sample_sums(Ps_obj) == 0)))
Ps_obj %>% 
  subset_samples(., Control == "FALSE") %>% 
  sample_sums(.) %>% 
  summary()
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    1683   30351   35312   35617   40831   71171

### Identify contaminants - Frequency

Use the distribution of the frequency of each sequence feature as a
function of the input DNA concentration to identify contaminants.

``` r
contamdf.freq <-
  isContaminant(Ps_obj, method = "frequency", conc = "Conc...ng.ul.")
# print(contamdf.freq)
# How many contaminants are found?
table(contamdf.freq$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ## 15983   117

``` r
# Which ones?
contamdf.freq %>% 
  rownames_to_column(., var = "ESV") %>% 
  filter(contaminant == TRUE) %>% 
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
ESV
</th>
<th style="text-align:right;">
freq
</th>
<th style="text-align:right;">
prev
</th>
<th style="text-align:right;">
p.freq
</th>
<th style="text-align:left;">
p.prev
</th>
<th style="text-align:right;">
p
</th>
<th style="text-align:left;">
contaminant
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Seq\_445
</td>
<td style="text-align:right;">
0.0002763
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0387245
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0387245
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_2408
</td>
<td style="text-align:right;">
0.0000407
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0905435
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0905435
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_2409
</td>
<td style="text-align:right;">
0.0000427
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0538470
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0538470
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_2736
</td>
<td style="text-align:right;">
0.0000349
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.0925595
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0925595
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_2805
</td>
<td style="text-align:right;">
0.0000317
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.0796852
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0796852
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3126
</td>
<td style="text-align:right;">
0.0000341
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0521557
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0521557
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4011
</td>
<td style="text-align:right;">
0.0000192
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0261061
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0261061
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4152
</td>
<td style="text-align:right;">
0.0000151
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0592902
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0592902
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5080
</td>
<td style="text-align:right;">
0.0000105
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0991985
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0991985
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5172
</td>
<td style="text-align:right;">
0.0000096
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0808375
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0808375
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5268
</td>
<td style="text-align:right;">
0.0000088
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0608814
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0608814
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5317
</td>
<td style="text-align:right;">
0.0000095
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0.0490381
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0490381
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5457
</td>
<td style="text-align:right;">
0.0000159
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0370295
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0370295
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5523
</td>
<td style="text-align:right;">
0.0000116
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0519441
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0519441
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5671
</td>
<td style="text-align:right;">
0.0000077
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0555776
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0555776
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5804
</td>
<td style="text-align:right;">
0.0000074
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0234556
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0234556
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6066
</td>
<td style="text-align:right;">
0.0000052
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0413684
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0413684
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6227
</td>
<td style="text-align:right;">
0.0000057
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0838105
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0838105
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6394
</td>
<td style="text-align:right;">
0.0000057
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0147065
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0147065
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6510
</td>
<td style="text-align:right;">
0.0000073
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0423869
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0423869
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6559
</td>
<td style="text-align:right;">
0.0000058
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0306906
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0306906
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6602
</td>
<td style="text-align:right;">
0.0000075
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0748997
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0748997
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6617
</td>
<td style="text-align:right;">
0.0000058
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0561624
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0561624
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6823
</td>
<td style="text-align:right;">
0.0000041
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0895070
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0895070
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6839
</td>
<td style="text-align:right;">
0.0000044
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0120885
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0120885
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6876
</td>
<td style="text-align:right;">
0.0000030
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0098784
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0098784
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7013
</td>
<td style="text-align:right;">
0.0000038
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0195169
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0195169
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7189
</td>
<td style="text-align:right;">
0.0000041
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0874963
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0874963
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7333
</td>
<td style="text-align:right;">
0.0000039
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0886636
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0886636
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7348
</td>
<td style="text-align:right;">
0.0000041
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0052122
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0052122
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7455
</td>
<td style="text-align:right;">
0.0000028
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0077646
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0077646
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7496
</td>
<td style="text-align:right;">
0.0000039
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0297187
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0297187
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7505
</td>
<td style="text-align:right;">
0.0000038
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0107067
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0107067
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7507
</td>
<td style="text-align:right;">
0.0000039
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0968903
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0968903
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7529
</td>
<td style="text-align:right;">
0.0000025
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0750584
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0750584
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7535
</td>
<td style="text-align:right;">
0.0000040
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0100901
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0100901
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7561
</td>
<td style="text-align:right;">
0.0000034
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0300724
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0300724
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7609
</td>
<td style="text-align:right;">
0.0000028
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0750064
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0750064
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7632
</td>
<td style="text-align:right;">
0.0000031
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0249306
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0249306
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7648
</td>
<td style="text-align:right;">
0.0000028
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0781969
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0781969
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7804
</td>
<td style="text-align:right;">
0.0000041
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0510789
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0510789
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7864
</td>
<td style="text-align:right;">
0.0000032
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0555134
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0555134
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7865
</td>
<td style="text-align:right;">
0.0000027
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0881105
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0881105
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7881
</td>
<td style="text-align:right;">
0.0000032
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0565157
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0565157
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7896
</td>
<td style="text-align:right;">
0.0000029
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0788497
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0788497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7927
</td>
<td style="text-align:right;">
0.0000033
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0914994
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0914994
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7976
</td>
<td style="text-align:right;">
0.0000026
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0086166
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0086166
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8036
</td>
<td style="text-align:right;">
0.0000029
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0735807
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0735807
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8049
</td>
<td style="text-align:right;">
0.0000027
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0327067
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0327067
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8078
</td>
<td style="text-align:right;">
0.0000027
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0833882
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0833882
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8166
</td>
<td style="text-align:right;">
0.0000025
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0567015
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0567015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8241
</td>
<td style="text-align:right;">
0.0000023
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0930813
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0930813
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8254
</td>
<td style="text-align:right;">
0.0000028
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0876634
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0876634
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8509
</td>
<td style="text-align:right;">
0.0000024
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0296828
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0296828
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8578
</td>
<td style="text-align:right;">
0.0000027
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0163416
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0163416
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8781
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0577997
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0577997
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8782
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0577997
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0577997
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8833
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0329594
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0329594
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8844
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0783735
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0783735
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8944
</td>
<td style="text-align:right;">
0.0000018
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0937240
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0937240
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8956
</td>
<td style="text-align:right;">
0.0000017
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0060297
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0060297
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9007
</td>
<td style="text-align:right;">
0.0000019
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0371456
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0371456
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9300
</td>
<td style="text-align:right;">
0.0000015
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0482373
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0482373
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9324
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0232842
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0232842
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9442
</td>
<td style="text-align:right;">
0.0000018
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0242152
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0242152
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9443
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0994433
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0994433
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9791
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0728597
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0728597
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9857
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0215645
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0215645
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9975
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0351238
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0351238
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10037
</td>
<td style="text-align:right;">
0.0000017
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0978638
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0978638
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10085
</td>
<td style="text-align:right;">
0.0000026
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0587871
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0587871
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10109
</td>
<td style="text-align:right;">
0.0000018
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0000151
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0000151
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10132
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0780088
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0780088
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10171
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0278458
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0278458
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10179
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0850435
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0850435
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10211
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0040031
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0040031
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10370
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0985753
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0985753
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10570
</td>
<td style="text-align:right;">
0.0000024
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0996305
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0996305
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10600
</td>
<td style="text-align:right;">
0.0000015
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0932461
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0932461
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10608
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0702671
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0702671
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10674
</td>
<td style="text-align:right;">
0.0000018
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0404042
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0404042
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10715
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0702211
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0702211
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10716
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0108627
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0108627
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10718
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0182145
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0182145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10860
</td>
<td style="text-align:right;">
0.0000013
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0527437
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0527437
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11160
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0901958
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0901958
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11162
</td>
<td style="text-align:right;">
0.0000013
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0899744
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0899744
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11195
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0758771
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0758771
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11339
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0087656
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0087656
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11432
</td>
<td style="text-align:right;">
0.0000009
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0975414
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0975414
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11678
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0360607
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0360607
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11797
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0467104
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0467104
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11802
</td>
<td style="text-align:right;">
0.0000009
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0831606
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0831606
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11827
</td>
<td style="text-align:right;">
0.0000013
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0430464
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0430464
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11889
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0647571
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0647571
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11927
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0370623
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0370623
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12275
</td>
<td style="text-align:right;">
0.0000007
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0410935
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0410935
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12413
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0585262
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0585262
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12439
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0910046
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0910046
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12657
</td>
<td style="text-align:right;">
0.0000009
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0727111
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0727111
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12658
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0886334
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0886334
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12842
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0829197
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0829197
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12876
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0237253
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0237253
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12878
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0591997
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0591997
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12972
</td>
<td style="text-align:right;">
0.0000004
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0161293
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0161293
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13116
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0032925
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0032925
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13118
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0034260
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0034260
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13235
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0280197
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0280197
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13345
</td>
<td style="text-align:right;">
0.0000007
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0560180
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0560180
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13538
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0537444
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0537444
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13684
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0669854
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0669854
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13704
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0170558
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0170558
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13816
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0239575
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0239575
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13890
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0794512
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0794512
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14052
</td>
<td style="text-align:right;">
0.0000004
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0893503
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0893503
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14163
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0654242
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0654242
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14204
</td>
<td style="text-align:right;">
0.0000004
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0921764
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0921764
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
</tbody>
</table>

Plot the frequency of the first 20 non-contaminant ESVs against the DNA
concentration, as an example.

``` r
plot_frequency(Ps_obj, taxa_names(Ps_obj)[which(!contamdf.freq$contaminant)[1:20]], conc = "Conc...ng.ul.")
```

![](02_Decontamination_figures/plot%20frequency%201-1.png)<!-- -->

Plot the frequency of the the first 20 contaminant sequences against the
DNA concentration.

``` r
plot_frequency(Ps_obj, taxa_names(Ps_obj)[which(contamdf.freq$contaminant)[1:20]], conc = "Conc...ng.ul.")
```

![](02_Decontamination_figures/plot%20frequency%202-1.png)<!-- -->

The frequency analysis detected 117 sequences as contaminants.

### Identify contaminants - Prevalence

Use the prevalence of sequences found in the control samples
(no-template controls) to identify contaminants.

``` r
Ps_obj_noBioControl <- subset_samples(Ps_obj, Control != "Biological") # remove "Biological control" samples
sample_data(Ps_obj_noBioControl)$Control <- sample_data(Ps_obj_noBioControl)$Control == "TRUE" # convert to logical
contamdf.prev <- isContaminant(Ps_obj_noBioControl, method = "prevalence", neg = "Control")
# How many contaminants are found?
table(contamdf.prev$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ## 15931   169

``` r
# Which ones?
contamdf.prev %>% 
  rownames_to_column(., var = "ESV") %>% 
  filter(contaminant == TRUE) %>% 
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
ESV
</th>
<th style="text-align:right;">
freq
</th>
<th style="text-align:right;">
prev
</th>
<th style="text-align:left;">
p.freq
</th>
<th style="text-align:right;">
p.prev
</th>
<th style="text-align:right;">
p
</th>
<th style="text-align:left;">
contaminant
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Seq\_1426
</td>
<td style="text-align:right;">
0.0004878
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0002935
</td>
<td style="text-align:right;">
0.0002935
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_2604
</td>
<td style="text-align:right;">
0.0000802
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0280257
</td>
<td style="text-align:right;">
0.0280257
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_2850
</td>
<td style="text-align:right;">
0.0000370
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3024
</td>
<td style="text-align:right;">
0.0000472
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0195314
</td>
<td style="text-align:right;">
0.0195314
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3207
</td>
<td style="text-align:right;">
0.0000307
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3242
</td>
<td style="text-align:right;">
0.0000105
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3251
</td>
<td style="text-align:right;">
0.0000425
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0320554
</td>
<td style="text-align:right;">
0.0320554
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3470
</td>
<td style="text-align:right;">
0.0000265
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3671
</td>
<td style="text-align:right;">
0.0000203
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3725
</td>
<td style="text-align:right;">
0.0000213
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3758
</td>
<td style="text-align:right;">
0.0000217
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3894
</td>
<td style="text-align:right;">
0.0000188
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3972
</td>
<td style="text-align:right;">
0.0000216
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_3974
</td>
<td style="text-align:right;">
0.0000214
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4136
</td>
<td style="text-align:right;">
0.0000158
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4265
</td>
<td style="text-align:right;">
0.0000213
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4320
</td>
<td style="text-align:right;">
0.0000148
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4380
</td>
<td style="text-align:right;">
0.0000166
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4673
</td>
<td style="text-align:right;">
0.0000166
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4704
</td>
<td style="text-align:right;">
0.0000137
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4705
</td>
<td style="text-align:right;">
0.0000133
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4707
</td>
<td style="text-align:right;">
0.0000135
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4848
</td>
<td style="text-align:right;">
0.0000130
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4921
</td>
<td style="text-align:right;">
0.0000118
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4971
</td>
<td style="text-align:right;">
0.0000124
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4980
</td>
<td style="text-align:right;">
0.0000087
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4981
</td>
<td style="text-align:right;">
0.0000138
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4983
</td>
<td style="text-align:right;">
0.0000121
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_4984
</td>
<td style="text-align:right;">
0.0000112
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5024
</td>
<td style="text-align:right;">
0.0000130
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5048
</td>
<td style="text-align:right;">
0.0000120
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5064
</td>
<td style="text-align:right;">
0.0000119
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5077
</td>
<td style="text-align:right;">
0.0000122
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5145
</td>
<td style="text-align:right;">
0.0000095
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5152
</td>
<td style="text-align:right;">
0.0000111
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5267
</td>
<td style="text-align:right;">
0.0000104
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5344
</td>
<td style="text-align:right;">
0.0000108
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5349
</td>
<td style="text-align:right;">
0.0000107
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5362
</td>
<td style="text-align:right;">
0.0000107
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5394
</td>
<td style="text-align:right;">
0.0000113
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5402
</td>
<td style="text-align:right;">
0.0000090
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5468
</td>
<td style="text-align:right;">
0.0000133
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5469
</td>
<td style="text-align:right;">
0.0000107
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5472
</td>
<td style="text-align:right;">
0.0000102
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5602
</td>
<td style="text-align:right;">
0.0000122
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5617
</td>
<td style="text-align:right;">
0.0000107
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5620
</td>
<td style="text-align:right;">
0.0000094
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5623
</td>
<td style="text-align:right;">
0.0000093
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5628
</td>
<td style="text-align:right;">
0.0000094
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5669
</td>
<td style="text-align:right;">
0.0000096
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5670
</td>
<td style="text-align:right;">
0.0000096
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5724
</td>
<td style="text-align:right;">
0.0000086
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5755
</td>
<td style="text-align:right;">
0.0000076
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5774
</td>
<td style="text-align:right;">
0.0000091
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5780
</td>
<td style="text-align:right;">
0.0000090
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5821
</td>
<td style="text-align:right;">
0.0000086
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5825
</td>
<td style="text-align:right;">
0.0000086
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5882
</td>
<td style="text-align:right;">
0.0000087
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5897
</td>
<td style="text-align:right;">
0.0000098
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_5964
</td>
<td style="text-align:right;">
0.0000077
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6023
</td>
<td style="text-align:right;">
0.0000082
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6050
</td>
<td style="text-align:right;">
0.0000083
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6067
</td>
<td style="text-align:right;">
0.0000080
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6090
</td>
<td style="text-align:right;">
0.0000068
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6172
</td>
<td style="text-align:right;">
0.0000076
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6211
</td>
<td style="text-align:right;">
0.0000077
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:right;">
0.0989676
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6285
</td>
<td style="text-align:right;">
0.0000081
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6299
</td>
<td style="text-align:right;">
0.0000068
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6357
</td>
<td style="text-align:right;">
0.0000067
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:right;">
0.0291781
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6383
</td>
<td style="text-align:right;">
0.0000070
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6389
</td>
<td style="text-align:right;">
0.0001302
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0000204
</td>
<td style="text-align:right;">
0.0000204
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6408
</td>
<td style="text-align:right;">
0.0000066
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6433
</td>
<td style="text-align:right;">
0.0000069
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6451
</td>
<td style="text-align:right;">
0.0000069
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6517
</td>
<td style="text-align:right;">
0.0000086
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6605
</td>
<td style="text-align:right;">
0.0000070
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6614
</td>
<td style="text-align:right;">
0.0000056
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6624
</td>
<td style="text-align:right;">
0.0000077
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6723
</td>
<td style="text-align:right;">
0.0000071
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6865
</td>
<td style="text-align:right;">
0.0000053
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6907
</td>
<td style="text-align:right;">
0.0000056
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6938
</td>
<td style="text-align:right;">
0.0000051
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6952
</td>
<td style="text-align:right;">
0.0000052
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6956
</td>
<td style="text-align:right;">
0.0000049
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_6968
</td>
<td style="text-align:right;">
0.0000054
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7016
</td>
<td style="text-align:right;">
0.0000052
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7022
</td>
<td style="text-align:right;">
0.0000056
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7042
</td>
<td style="text-align:right;">
0.0000044
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7072
</td>
<td style="text-align:right;">
0.0000054
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:right;">
0.0870497
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7104
</td>
<td style="text-align:right;">
0.0000054
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7151
</td>
<td style="text-align:right;">
0.0000054
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7167
</td>
<td style="text-align:right;">
0.0000040
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7183
</td>
<td style="text-align:right;">
0.0000049
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7303
</td>
<td style="text-align:right;">
0.0000054
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7477
</td>
<td style="text-align:right;">
0.0000063
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7540
</td>
<td style="text-align:right;">
0.0000042
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:right;">
0.0757145
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7564
</td>
<td style="text-align:right;">
0.0000050
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:right;">
0.0650015
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7633
</td>
<td style="text-align:right;">
0.0000066
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0160635
</td>
<td style="text-align:right;">
0.0160635
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7653
</td>
<td style="text-align:right;">
0.0000040
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7670
</td>
<td style="text-align:right;">
0.0000046
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7673
</td>
<td style="text-align:right;">
0.0000049
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:right;">
0.0549505
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7846
</td>
<td style="text-align:right;">
0.0000039
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:right;">
0.0456021
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7873
</td>
<td style="text-align:right;">
0.0000038
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0065989
</td>
<td style="text-align:right;">
0.0065989
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_7953
</td>
<td style="text-align:right;">
0.0000042
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8084
</td>
<td style="text-align:right;">
0.0000051
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8230
</td>
<td style="text-align:right;">
0.0000027
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8493
</td>
<td style="text-align:right;">
0.0000042
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8513
</td>
<td style="text-align:right;">
0.0000029
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8661
</td>
<td style="text-align:right;">
0.0000028
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:right;">
0.0369974
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8764
</td>
<td style="text-align:right;">
0.0000047
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_8971
</td>
<td style="text-align:right;">
0.0000021
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9055
</td>
<td style="text-align:right;">
0.0000019
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9386
</td>
<td style="text-align:right;">
0.0000023
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9387
</td>
<td style="text-align:right;">
0.0000024
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9468
</td>
<td style="text-align:right;">
0.0000030
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9482
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:right;">
0.0221860
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9644
</td>
<td style="text-align:right;">
0.0000029
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9694
</td>
<td style="text-align:right;">
0.0000034
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_9858
</td>
<td style="text-align:right;">
0.0000017
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10313
</td>
<td style="text-align:right;">
0.0000015
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10478
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10524
</td>
<td style="text-align:right;">
0.0000020
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:right;">
0.0108535
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_10873
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11070
</td>
<td style="text-align:right;">
0.0000018
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11175
</td>
<td style="text-align:right;">
0.0000019
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11185
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11274
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11347
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11363
</td>
<td style="text-align:right;">
0.0000013
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11366
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11671
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11726
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11780
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11876
</td>
<td style="text-align:right;">
0.0000020
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11942
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_11966
</td>
<td style="text-align:right;">
0.0000014
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12047
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12085
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12105
</td>
<td style="text-align:right;">
0.0000007
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12172
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12266
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0065989
</td>
<td style="text-align:right;">
0.0065989
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12473
</td>
<td style="text-align:right;">
0.0000009
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12483
</td>
<td style="text-align:right;">
0.0000007
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12532
</td>
<td style="text-align:right;">
0.0000022
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12612
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12614
</td>
<td style="text-align:right;">
0.0000012
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12719
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12728
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12767
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12808
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_12811
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13127
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13129
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13152
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13193
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13237
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13253
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13346
</td>
<td style="text-align:right;">
0.0000013
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13347
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13356
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13377
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13390
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13414
</td>
<td style="text-align:right;">
0.0000011
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0033429
</td>
<td style="text-align:right;">
0.0033429
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_13435
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:right;">
0.0759060
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14051
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14361
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14362
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14374
</td>
<td style="text-align:right;">
0.0000004
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
<tr>
<td style="text-align:left;">
Seq\_14426
</td>
<td style="text-align:right;">
0.0000005
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:right;">
0.0506329
</td>
<td style="text-align:left;">
TRUE
</td>
</tr>
</tbody>
</table>

``` r
# And using a more aggrssive threshold
contamdf.prev05 <- isContaminant(Ps_obj_noBioControl, method = "prevalence", neg = "Control", threshold = 0.5)
table(contamdf.prev05$contaminant)
```

    ## 
    ## FALSE  TRUE 
    ## 14509  1591

``` r
# Make phyloseq object of presence-absence in negative controls
Ps_obj.pa <-
  transform_sample_counts(Ps_obj, function(abund)
    1 * (abund > 0))
Ps_obj.pa.neg <-
  prune_samples(sample_data(Ps_obj.pa)$Control == "TRUE", Ps_obj.pa)
Ps_obj.pa.pos <-
  prune_samples(sample_data(Ps_obj.pa)$Control == "FALSE", Ps_obj.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <-
  data.frame(
    pa.pos = taxa_sums(Ps_obj.pa.pos),
    pa.neg = taxa_sums(Ps_obj.pa.neg),
    contaminant = contamdf.prev$contaminant
  )
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
```

![](02_Decontamination_figures/prevalence-1.png)<!-- -->

The prevalence analysis detected 169 sequences as contaminants. In total
286 were detected as contaminants and will be removed.

### Save contaminant sequence names

``` r
c(taxa_names(Ps_obj)[which(contamdf.freq$contaminant)],
  taxa_names(Ps_obj)[which(contamdf.prev$contaminant)]) %>%
  data.frame() %>%
  write_csv(., 
            paste0(data_path, "decontam_contaminants.csv"), 
            col_names = FALSE)
```

``` r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
 )
```

<details open>
<summary>
<span title="Click to Expand"> Current session info </span>
</summary>

``` r
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
 date     2021-05-21                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package      * version date       lib source        
 ade4           1.7-16  2020-10-28 [1] CRAN (R 4.0.2)
 ape            5.5     2021-04-25 [1] CRAN (R 4.0.3)
 assertthat     0.2.1   2019-03-21 [1] CRAN (R 4.0.2)
 backports      1.2.1   2020-12-09 [1] CRAN (R 4.0.2)
 Biobase        2.50.0  2020-10-27 [1] Bioconductor  
 BiocGenerics   0.36.1  2021-04-16 [1] Bioconductor  
 biomformat     1.18.0  2020-10-27 [1] Bioconductor  
 Biostrings     2.58.0  2020-10-27 [1] Bioconductor  
 broom          0.7.6   2021-04-05 [1] CRAN (R 4.0.3)
 cellranger     1.1.0   2016-07-27 [1] CRAN (R 4.0.2)
 cli            2.5.0   2021-04-26 [1] CRAN (R 4.0.3)
 clipr          0.7.1   2020-10-08 [1] CRAN (R 4.0.2)
 cluster        2.1.2   2021-04-17 [1] CRAN (R 4.0.3)
 codetools      0.2-18  2020-11-04 [1] CRAN (R 4.0.2)
 colorspace     2.0-1   2021-05-04 [1] CRAN (R 4.0.3)
 cowplot      * 1.1.1   2020-12-30 [1] CRAN (R 4.0.2)
 crayon         1.4.1   2021-02-08 [1] CRAN (R 4.0.3)
 data.table     1.14.0  2021-02-21 [1] CRAN (R 4.0.3)
 DBI            1.1.1   2021-01-15 [1] CRAN (R 4.0.3)
 dbplyr         2.1.1   2021-04-06 [1] CRAN (R 4.0.3)
 decontam     * 1.10.0  2020-10-27 [1] Bioconductor  
 desc           1.3.0   2021-03-05 [1] CRAN (R 4.0.3)
 details        0.2.1   2020-01-12 [1] CRAN (R 4.0.2)
 digest         0.6.27  2020-10-24 [1] CRAN (R 4.0.2)
 dplyr        * 1.0.6   2021-05-05 [1] CRAN (R 4.0.3)
 ellipsis       0.3.2   2021-04-29 [1] CRAN (R 4.0.3)
 evaluate       0.14    2019-05-28 [1] CRAN (R 4.0.2)
 extrafont    * 0.17    2014-12-08 [1] CRAN (R 4.0.2)
 extrafontdb    1.0     2012-06-11 [1] CRAN (R 4.0.2)
 fansi          0.4.2   2021-01-15 [1] CRAN (R 4.0.3)
 farver         2.1.0   2021-02-28 [1] CRAN (R 4.0.3)
 forcats      * 0.5.1   2021-01-27 [1] CRAN (R 4.0.3)
 foreach        1.5.1   2020-10-15 [1] CRAN (R 4.0.2)
 fs             1.5.0   2020-07-31 [1] CRAN (R 4.0.2)
 generics       0.1.0   2020-10-31 [1] CRAN (R 4.0.2)
 ggplot2      * 3.3.3   2020-12-30 [1] CRAN (R 4.0.2)
 glue           1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
 gtable         0.3.0   2019-03-25 [1] CRAN (R 4.0.2)
 haven          2.4.1   2021-04-23 [1] CRAN (R 4.0.3)
 highr          0.9     2021-04-16 [1] CRAN (R 4.0.3)
 hms            1.1.0   2021-05-17 [1] CRAN (R 4.0.3)
 htmltools      0.5.1.1 2021-01-22 [1] CRAN (R 4.0.3)
 httr           1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
 igraph         1.2.6   2020-10-06 [1] CRAN (R 4.0.2)
 IRanges        2.24.1  2020-12-12 [1] Bioconductor  
 iterators      1.0.13  2020-10-15 [1] CRAN (R 4.0.2)
 jsonlite       1.7.2   2020-12-09 [1] CRAN (R 4.0.2)
 kableExtra   * 1.3.4   2021-02-20 [1] CRAN (R 4.0.3)
 knitr        * 1.33    2021-04-24 [1] CRAN (R 4.0.3)
 labeling       0.4.2   2020-10-20 [1] CRAN (R 4.0.2)
 lattice        0.20-44 2021-05-02 [1] CRAN (R 4.0.3)
 lifecycle      1.0.0   2021-02-15 [1] CRAN (R 4.0.3)
 lubridate      1.7.10  2021-02-26 [1] CRAN (R 4.0.3)
 magrittr     * 2.0.1   2020-11-17 [1] CRAN (R 4.0.2)
 MASS           7.3-54  2021-05-03 [1] CRAN (R 4.0.3)
 Matrix         1.3-3   2021-05-04 [1] CRAN (R 4.0.3)
 mgcv           1.8-35  2021-04-18 [1] CRAN (R 4.0.3)
 modelr         0.1.8   2020-05-19 [1] CRAN (R 4.0.2)
 multtest       2.46.0  2020-10-27 [1] Bioconductor  
 munsell        0.5.0   2018-06-12 [1] CRAN (R 4.0.2)
 my.tools     * 0.5     2020-09-30 [1] local         
 nlme           3.1-152 2021-02-04 [1] CRAN (R 4.0.3)
 permute        0.9-5   2019-03-12 [1] CRAN (R 4.0.2)
 phyloseq     * 1.34.0  2020-10-27 [1] Bioconductor  
 pillar         1.6.1   2021-05-16 [1] CRAN (R 4.0.3)
 pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.0.2)
 plyr           1.8.6   2020-03-03 [1] CRAN (R 4.0.2)
 png            0.1-7   2013-12-03 [1] CRAN (R 4.0.2)
 prettyunits    1.1.1   2020-01-24 [1] CRAN (R 4.0.2)
 progress       1.2.2   2019-05-16 [1] CRAN (R 4.0.2)
 purrr        * 0.3.4   2020-04-17 [1] CRAN (R 4.0.2)
 R6             2.5.0   2020-10-28 [1] CRAN (R 4.0.2)
 ragg         * 1.1.2   2021-03-17 [1] CRAN (R 4.0.3)
 Rcpp           1.0.6   2021-01-15 [1] CRAN (R 4.0.3)
 readr        * 1.4.0   2020-10-05 [1] CRAN (R 4.0.2)
 readxl         1.3.1   2019-03-13 [1] CRAN (R 4.0.2)
 reprex         2.0.0   2021-04-02 [1] CRAN (R 4.0.3)
 reshape2       1.4.4   2020-04-09 [1] CRAN (R 4.0.2)
 rhdf5          2.34.0  2020-10-27 [1] Bioconductor  
 rhdf5filters   1.2.1   2021-05-03 [1] Bioconductor  
 Rhdf5lib       1.12.1  2021-01-26 [1] Bioconductor  
 rlang          0.4.11  2021-04-30 [1] CRAN (R 4.0.3)
 rmarkdown    * 2.8     2021-05-07 [1] CRAN (R 4.0.3)
 rprojroot      2.0.2   2020-11-15 [1] CRAN (R 4.0.2)
 rstudioapi     0.13    2020-11-12 [1] CRAN (R 4.0.2)
 Rttf2pt1       1.3.8   2020-01-10 [1] CRAN (R 4.0.2)
 rvest          1.0.0   2021-03-09 [1] CRAN (R 4.0.3)
 S4Vectors      0.28.1  2020-12-09 [1] Bioconductor  
 scales       * 1.1.1   2020-05-11 [1] CRAN (R 4.0.2)
 sessioninfo    1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
 stringi        1.6.2   2021-05-17 [1] CRAN (R 4.0.3)
 stringr      * 1.4.0   2019-02-10 [1] CRAN (R 4.0.2)
 survival       3.2-11  2021-04-26 [1] CRAN (R 4.0.3)
 svglite      * 2.0.0   2021-02-20 [1] CRAN (R 4.0.3)
 systemfonts    1.0.2   2021-05-11 [1] CRAN (R 4.0.3)
 textshaping    0.3.4   2021-05-11 [1] CRAN (R 4.0.3)
 tibble       * 3.1.2   2021-05-16 [1] CRAN (R 4.0.3)
 tidyr        * 1.1.3   2021-03-03 [1] CRAN (R 4.0.3)
 tidyselect     1.1.1   2021-04-30 [1] CRAN (R 4.0.3)
 tidyverse    * 1.3.1   2021-04-15 [1] CRAN (R 4.0.3)
 utf8           1.2.1   2021-03-12 [1] CRAN (R 4.0.3)
 vctrs          0.3.8   2021-04-29 [1] CRAN (R 4.0.3)
 vegan          2.5-7   2020-11-28 [1] CRAN (R 4.0.3)
 viridisLite    0.4.0   2021-04-13 [1] CRAN (R 4.0.3)
 webshot        0.5.2   2019-11-22 [1] CRAN (R 4.0.2)
 withr          2.4.2   2021-04-18 [1] CRAN (R 4.0.3)
 xfun           0.23    2021-05-15 [1] CRAN (R 4.0.3)
 xml2           1.3.2   2020-04-23 [1] CRAN (R 4.0.2)
 XVector        0.30.0  2020-10-27 [1] Bioconductor  
 yaml           2.2.1   2020-02-01 [1] CRAN (R 4.0.2)
 zlibbioc       1.36.0  2020-10-27 [1] Bioconductor  

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```

</details>

<br>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-davis_simple_2017" class="csl-entry">

Davis NM, Proctor D, Holmes SP *et al.* <span class="nocase">Simple
statistical identification and removal of contaminant sequences in
marker-gene and metagenomics data</span>. *bioRxiv* 2017:221499.

</div>

</div>

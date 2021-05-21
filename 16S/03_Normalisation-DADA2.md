TeaTime4schools: joint analysis - bacteria
================
Roey Angel
2021-05-21

-   [Read-depth distribution and
    normalisation](#read-depth-distribution-and-normalisation)
    -   [Setting general parameters:](#setting-general-parameters)
    -   [Reading in raw data and inspecting read-depth
        variations](#reading-in-raw-data-and-inspecting-read-depth-variations)
    -   [Exploring Ps\_obj dataset
        features](#exploring-ps_obj-dataset-features)
    -   [Try various normalisation
        methods](#try-various-normalisation-methods)
        -   [Rarefaction](#rarefaction)
        -   [GMPR (Chen and Chen 2017)](#gmpr-chen_gmpr_2017)
        -   [Cumulative sum scaling normalization (Paulson *et al.*
            2013)](#cumulative-sum-scaling-normalization-paulson_differential_2013)
        -   [Standardize abundances to the median sequencing depth (and
            convert to
            proportion)](#standardize-abundances-to-the-median-sequencing-depth-and-convert-to-proportion)
        -   [Standardize abundances using log transformation for
            variance
            stabilisation](#standardize-abundances-using-log-transformation-for-variance-stabilisation)
    -   [Standardize abundances using Centered-Log-Ratio transformation
        for variance stabilisation (Fernandes *et al.*
        2013)](#standardize-abundances-using-centered-log-ratio-transformation-for-variance-stabilisation-fernandes_anova-like_2013)
-   [References](#references)

[roey.angel@bc.cas.cz](mailto:%20roey.angel@bc.cas.cz)

## Read-depth distribution and normalisation

This analysis tests the effect of library read-depth distribution on the
community composition. It then performs various read-depth normalisation
methods on the DADA2 zOTU-based dataset, for determining the optimal
strategy to handle the bias of uneven read distribution.

### Setting general parameters:

``` r
set.seed(1000)
min_lib_size <- 5000
metadata_path <- "./"
data_path <- "./DADA2_pseudo/"
Metadata_table <- "TeaTime_joint_Bacteria_metadata.csv"
Seq_table <- "DADA2.seqtab_nochim.tsv"
```

### Reading in raw data and inspecting read-depth variations

Read abundance table, taxonomic classification and metadata into a
phyloseq object. Also remove sequences detected as contaminants in
[03\_Decontamination.html](03_Decontamination.html).

``` r
# read OTU mat from data file
Raw_data <- read_tsv(paste0(data_path, Seq_table), 
                        trim_ws = TRUE)
contaminant_seqs <- read_csv(paste0(data_path, "decontam_contaminants.csv"), 
                        trim_ws = TRUE,
                        col_names = FALSE)

Raw_data %<>% # remove contaminant OTUs. 
  # .[, -grep("CTRL", colnames(.))] %>% # remove ext. cont. 
  .[!(Raw_data$`#OTU` %in% contaminant_seqs$X1), ] 

Raw_data[, 2:ncol(Raw_data)] %>% 
  t() %>% 
  as.data.frame() -> abundance_mat # convert to abundance matrix
colnames(abundance_mat) <- pull(Raw_data, "#OTU") # add sequence names

# Read metadata file
read_csv(paste0(metadata_path, Metadata_table),
         trim_ws = TRUE) %>%
  mutate_at(
    c(
      "Workshop",
      "Season",
      "Run",
      "Type",
      "Sample type",
      "Field",
      "Replicate",
      "Control",
      "Gene"
    ),
    funs(factor(.))
  ) %>% 
  mutate_at(c("Extr. Date", "PCR products_16S_send for seq"), ~as.Date(., "%d.%m.%Y")) ->
  Metadata
Metadata$Season %<>% fct_relevel("Winter", "Spring", "Summer", "Autumn")
Metadata$Read1_file <- str_replace(Metadata$Read1_file, "(.*)_L001_R1_001.fastq.gz|\\.1\\.fastq.gz", "\\1")
Metadata <- Metadata[Metadata$Read1_file %in% rownames(abundance_mat), ] # remove metadata rows if the samples did not go through qual processing

# Order abundance_mat samples according to the metadata
sample_order <- match(rownames(abundance_mat), Metadata$Read1_file)
abundance_mat %<>% 
  rownames_to_column('sample_name') %>% 
  arrange(., sample_order) %>% 
  column_to_rownames('sample_name') # needed for phyloseq

Metadata$Library.size <- rowSums(abundance_mat)
Metadata <- data.frame(row.names = Metadata$Read1_file, Metadata)

# generate phyloseq object
Ps_obj <- phyloseq(otu_table(abundance_mat, taxa_are_rows = FALSE),
                        sample_data(Metadata)
                        )
Ps_obj <- filter_taxa(Ps_obj, function(x) sum(x) > 0, TRUE) # remove 0 abundance taxa
Ps_obj <- subset_samples(Ps_obj, sample_sums(Ps_obj) > 0) # remove 0 abundance samples
# Remove mock and control samples
Ps_obj <- subset_samples(Ps_obj, Type != "Control" & Type != "Mock")

# Create a grouping variable for merging
sample_data(Ps_obj) %<>% 
  as(., "data.frame") %>% 
  # get_variable(., c("Sample.type", "Field", "Season", "Replicate")) %>% 
  unite(., "Description", c("Sample.type", "Field", "Season", "Replicate"), remove = FALSE)
sample_data(Ps_obj)$Description %<>% as_factor(.)

# merged_Ps_obj <- merge_samples(Ps_obj, "Description")
# merged_SD <- merge_samples(sample_data(Ps_obj), "Description")
Ps_obj_merged <- MergeSamples(Ps_obj, grouping_name = "Description")
```

### Exploring Ps\_obj dataset features

First let’s look at the count data distribution

``` r
PlotLibDist(Ps_obj)
```

![](03_Normalisation_DADA2_figures/plot%20abundance-1.png)<!-- -->

``` r
get_variable(Ps_obj) %>% 
  remove_rownames %>% 
  dplyr::select(Sample.name, Library.size) %>% 
  as(., "data.frame") %>% 
  format_engr() %>% 
  kable(., col.names = c("Sample name", "Library size")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Sample name
</th>
<th style="text-align:left;">
Library size
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GK1
</td>
<td style="text-align:left;">
28.94 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK1
</td>
<td style="text-align:left;">
50.44 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK2
</td>
<td style="text-align:left;">
26.81 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK2
</td>
<td style="text-align:left;">
71.08 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK3
</td>
<td style="text-align:left;">
32.51 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK3
</td>
<td style="text-align:left;">
32.52 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK4
</td>
<td style="text-align:left;">
32.20 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK4
</td>
<td style="text-align:left;">
35.47 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK1
</td>
<td style="text-align:left;">
32.99 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK1
</td>
<td style="text-align:left;">
25.69 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK2
</td>
<td style="text-align:left;">
32.45 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK2
</td>
<td style="text-align:left;">
45.07 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK3
</td>
<td style="text-align:left;">
33.51 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK3
</td>
<td style="text-align:left;">
16.63 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK4
</td>
<td style="text-align:left;">
34.32 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK4
</td>
<td style="text-align:left;">
32.02 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB1
</td>
<td style="text-align:left;">
35.71 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB1
</td>
<td style="text-align:left;">
49.37 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB2
</td>
<td style="text-align:left;">
38.45 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB2
</td>
<td style="text-align:left;">
22.98 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB3
</td>
<td style="text-align:left;">
36.92 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB3
</td>
<td style="text-align:left;">
30.32 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB4
</td>
<td style="text-align:left;">
34.06 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB4
</td>
<td style="text-align:left;">
25.20 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB1
</td>
<td style="text-align:left;">
30.92 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB1
</td>
<td style="text-align:left;">
36.23 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB2
</td>
<td style="text-align:left;">
36.95 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB2
</td>
<td style="text-align:left;">
23.29 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB3
</td>
<td style="text-align:left;">
36.15 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB3
</td>
<td style="text-align:left;">
29.30 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB4
</td>
<td style="text-align:left;">
39.52 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB4
</td>
<td style="text-align:left;">
19.58 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA1
</td>
<td style="text-align:left;">
30.44 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA1
</td>
<td style="text-align:left;">
23.58 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA2
</td>
<td style="text-align:left;">
39.63 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA2
</td>
<td style="text-align:left;">
19.50 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA3
</td>
<td style="text-align:left;">
36.05 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA3
</td>
<td style="text-align:left;">
12.77 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA4
</td>
<td style="text-align:left;">
34.87 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA4
</td>
<td style="text-align:left;">
19.85 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA1
</td>
<td style="text-align:left;">
31.64 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA1
</td>
<td style="text-align:left;">
33.09 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA2
</td>
<td style="text-align:left;">
39.14 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA2
</td>
<td style="text-align:left;">
14.18 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA3
</td>
<td style="text-align:left;">
37.75 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA3
</td>
<td style="text-align:left;">
19.79 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA4
</td>
<td style="text-align:left;">
38.75 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA4
</td>
<td style="text-align:left;">
21.64 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1
</td>
<td style="text-align:left;">
19.92 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1
</td>
<td style="text-align:left;">
24.69 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1
</td>
<td style="text-align:left;">
15.68 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1
</td>
<td style="text-align:left;">
20.69 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:left;">
34.63 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1
</td>
<td style="text-align:left;">
28.78 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1-1-2
</td>
<td style="text-align:left;">
34.10 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1-1-2
</td>
<td style="text-align:left;">
24.93 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1-1-2
</td>
<td style="text-align:left;">
27.59 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RF
</td>
<td style="text-align:left;">
2.080 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RF
</td>
<td style="text-align:left;">
8.196 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GF
</td>
<td style="text-align:left;">
2.980 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GF
</td>
<td style="text-align:left;">
14.27 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK5
</td>
<td style="text-align:left;">
45.22 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK7
</td>
<td style="text-align:left;">
34.08 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK8
</td>
<td style="text-align:left;">
32.18 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK5
</td>
<td style="text-align:left;">
35.68 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK6
</td>
<td style="text-align:left;">
46.20 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK7
</td>
<td style="text-align:left;">
33.14 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK8
</td>
<td style="text-align:left;">
31.62 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB5
</td>
<td style="text-align:left;">
30.27 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB6
</td>
<td style="text-align:left;">
33.33 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB7
</td>
<td style="text-align:left;">
32.68 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB8
</td>
<td style="text-align:left;">
37.00 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB5
</td>
<td style="text-align:left;">
27.74 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB6
</td>
<td style="text-align:left;">
26.78 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB7
</td>
<td style="text-align:left;">
32.31 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB8
</td>
<td style="text-align:left;">
38.73 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA5
</td>
<td style="text-align:left;">
33.86 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA6
</td>
<td style="text-align:left;">
37.63 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA7
</td>
<td style="text-align:left;">
35.71 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA8
</td>
<td style="text-align:left;">
41.33 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA5
</td>
<td style="text-align:left;">
32.83 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA6
</td>
<td style="text-align:left;">
31.90 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA7
</td>
<td style="text-align:left;">
37.63 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA8
</td>
<td style="text-align:left;">
37.33 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1-2-1
</td>
<td style="text-align:left;">
27.38 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1-2-1
</td>
<td style="text-align:left;">
34.18 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1-2-1
</td>
<td style="text-align:left;">
35.75 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1-2-2
</td>
<td style="text-align:left;">
32.31 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1-2-2
</td>
<td style="text-align:left;">
33.40 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1-2-2
</td>
<td style="text-align:left;">
32.97 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK10
</td>
<td style="text-align:left;">
41.93 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK11
</td>
<td style="text-align:left;">
44.16 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK9
</td>
<td style="text-align:left;">
48.02 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK10
</td>
<td style="text-align:left;">
43.50 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK12
</td>
<td style="text-align:left;">
40.76 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB9
</td>
<td style="text-align:left;">
38.22 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB10
</td>
<td style="text-align:left;">
38.48 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB11
</td>
<td style="text-align:left;">
33.14 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB12
</td>
<td style="text-align:left;">
31.39 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB9
</td>
<td style="text-align:left;">
46.35 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB10
</td>
<td style="text-align:left;">
44.48 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB11
</td>
<td style="text-align:left;">
43.94 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB12
</td>
<td style="text-align:left;">
6.680 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA9
</td>
<td style="text-align:left;">
36.23 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA10
</td>
<td style="text-align:left;">
35.75 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA11
</td>
<td style="text-align:left;">
38.70 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA12
</td>
<td style="text-align:left;">
36.52 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA9
</td>
<td style="text-align:left;">
36.78 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA10
</td>
<td style="text-align:left;">
36.15 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA11
</td>
<td style="text-align:left;">
37.26 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA12
</td>
<td style="text-align:left;">
53.28 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1-3-1
</td>
<td style="text-align:left;">
31.98 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1-3-2
</td>
<td style="text-align:left;">
34.79 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
K1-3-3
</td>
<td style="text-align:left;">
47.26 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1-3-1
</td>
<td style="text-align:left;">
32.99 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1-3-2
</td>
<td style="text-align:left;">
37.53 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
B1-3-3
</td>
<td style="text-align:left;">
36.56 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1-3-1
</td>
<td style="text-align:left;">
43.96 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1-3-2
</td>
<td style="text-align:left;">
48.38 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
A1-3-3
</td>
<td style="text-align:left;">
27.25 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK13
</td>
<td style="text-align:left;">
47.76 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK14
</td>
<td style="text-align:left;">
47.93 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK15
</td>
<td style="text-align:left;">
32.01 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RK16
</td>
<td style="text-align:left;">
51.30 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK13
</td>
<td style="text-align:left;">
58.42 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK14
</td>
<td style="text-align:left;">
50.94 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK15
</td>
<td style="text-align:left;">
59.82 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GK16
</td>
<td style="text-align:left;">
52.24 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB13
</td>
<td style="text-align:left;">
53.48 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB14
</td>
<td style="text-align:left;">
50.98 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB15
</td>
<td style="text-align:left;">
34.69 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RB16
</td>
<td style="text-align:left;">
56.18 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB13
</td>
<td style="text-align:left;">
42.41 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB14
</td>
<td style="text-align:left;">
35.62 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB15
</td>
<td style="text-align:left;">
41.01 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GB16
</td>
<td style="text-align:left;">
58.61 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA13
</td>
<td style="text-align:left;">
49.54 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA14
</td>
<td style="text-align:left;">
60.50 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA15
</td>
<td style="text-align:left;">
50.59 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
RA16
</td>
<td style="text-align:left;">
56.54 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA13
</td>
<td style="text-align:left;">
39.18 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA14
</td>
<td style="text-align:left;">
37.95 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA15
</td>
<td style="text-align:left;">
46.26 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
GA16
</td>
<td style="text-align:left;">
40.85 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 6-1
</td>
<td style="text-align:left;">
36.19 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 6-2
</td>
<td style="text-align:left;">
21.35 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 6-3
</td>
<td style="text-align:left;">
24.06 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 5-1
</td>
<td style="text-align:left;">
23.86 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 5-2
</td>
<td style="text-align:left;">
27.95 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 5-3
</td>
<td style="text-align:left;">
28.37 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 4-1
</td>
<td style="text-align:left;">
1.677 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 4-2
</td>
<td style="text-align:left;">
19.94 × 10<sup>3</sup>
</td>
</tr>
<tr>
<td style="text-align:left;">
soil 4-3
</td>
<td style="text-align:left;">
25.92 × 10<sup>3</sup>
</td>
</tr>
</tbody>
</table>

The figure and table indicate only a small deviation in the number of
reads per samples.

``` r
(mod1 <- adonis(vegdist(otu_table(Ps_obj), method = "bray") ~ Library.size,
  data = as(sample_data(Ps_obj), "data.frame"),
  permutations = 999
))
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj), method = "bray") ~      Library.size, data = as(sample_data(Ps_obj), "data.frame"),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Library.size   1     2.795 2.79458  8.4201 0.05282  0.001 ***
    ## Residuals    151    50.116 0.33189         0.94718           
    ## Total        152    52.911                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotReadHist(as(otu_table(Ps_obj), "matrix"))
```

![](03_Normalisation_DADA2_figures/mod%20abundance-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj))) > 0)
meanSdPlot(as.matrix(log2(t(otu_table(Ps_obj))[notAllZero, ] + 1)))
```

![](03_Normalisation_DADA2_figures/mod%20abundance-2.png)<!-- -->
Nevertheless, modelling library size shows a significant effect of read
depth on the community structure (explaining only 5%). The reads
histogram shows as expected a highly sparse and skewed sequence matrix.
The mean vs SD also shows as expected large dependency of SD on the mean
reads of a sequence across all samples.

### Try various normalisation methods

``` r
# subsample libraries from 1000 to max(sample_sums(Ps_obj)) and test
for (i in seq(1000, max(sample_sums(Ps_obj)), 1000)) {
  min_seqs <<- i
  Ps_obj_pruned <- subset_samples(Ps_obj, sample_sums(Ps_obj) > i)
  mod <-
    adonis2(vegdist(otu_table(Ps_obj_pruned), method = "bray") ~ sample_sums(Ps_obj_pruned), # I use adonis2 because it gives a data.frame
      data = get_variable(Ps_obj_pruned),
      permutations = 999
    )
  Pval <- tidy(mod)$p.value[1]
  if (Pval > 0.05)
    break()
}
```

Only by subsetting the samples to a minimum library size of
3.9 × 10<sup>4</sup> sequences do we get independence from library size
but this forces us to drop 111 out of 153 samples.

Let’s see the effect of this

``` r
Ps_obj_pruned_harsh <- subset_samples(Ps_obj_merged, sample_sums(Ps_obj) > i)
adonis(
  vegdist(otu_table(Ps_obj_pruned_harsh), method = "bray") ~ Library.size,
  data =
    get_variable(Ps_obj_pruned_harsh),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_harsh), method = "bray") ~      Library.size, data = get_variable(Ps_obj_pruned_harsh), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##              Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Library.size  1    1.2659 1.26593  3.8539 0.14351  0.001 ***
    ## Residuals    23    7.5551 0.32848         0.85649           
    ## Total        24    8.8210                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(
  vegdist(otu_table(Ps_obj_pruned_harsh), method = "bray") ~ Field + Sample.type * Season,
  data =
    get_variable(Ps_obj_pruned_harsh),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_harsh), method = "bray") ~      Field + Sample.type * Season, data = get_variable(Ps_obj_pruned_harsh),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                    Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field               3    1.1018 0.36726  1.8846 0.12490  0.005 ** 
    ## Sample.type         2    3.1843 1.59213  8.1701 0.36099  0.001 ***
    ## Season              3    1.3995 0.46651  2.3939 0.15866  0.001 ***
    ## Sample.type:Season  2    0.4072 0.20359  1.0447 0.04616  0.417    
    ## Residuals          14    2.7282 0.19487         0.30929           
    ## Total              24    8.8210                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_harsh), "matrix"))
```

![](03_Normalisation_DADA2_figures/effect%20s%20pruning-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_harsh))) > 0)
meanSdPlot(as.matrix(log2(t(otu_table(Ps_obj_pruned_harsh))[notAllZero, ] + 1)))
```

![](03_Normalisation_DADA2_figures/effect%20s%20pruning-2.png)<!-- -->

``` r
Pruned_ord <- ordinate(Ps_obj_pruned_harsh, "CAP", "bray", formula = Ps_obj_pruned_harsh ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_harsh, Pruned_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field") +
  facet_wrap(~ Season) #+
```

![](03_Normalisation_DADA2_figures/effect%20s%20pruning%20-%20ordinate-1.png)<!-- -->

``` r
  # geom_point(size = 5) +
  # geom_text(size = 5) 
```

Instead let’s drop all samples below 5000 reads and do try some
correction methods for library depths

``` r
Ps_obj_pruned_min <- subset_samples(Ps_obj_merged, sample_sums(Ps_obj_merged) > min_lib_size)
Ps_obj_pruned_min <- filter_taxa(Ps_obj_pruned_min, function(x) sum(x) > 0, TRUE)
```

#### Rarefaction

``` r
Ps_obj_pruned_rared <-
  rarefy_even_depth(
  Ps_obj_pruned_min,
  sample.size = min(sample_sums(Ps_obj_pruned_min)),
  rngseed = FALSE,
  replace = FALSE
  )
sample_data(Ps_obj_pruned_rared)$Library.size <- sample_sums(Ps_obj_pruned_rared)
# Ps_obj_pruned_rared <- Ps_obj_pruned_min
# Ps_obj_pruned_min %>%
#   otu_table() %>%
#   rowSums() %>%
#   min() %>%
#   rrarefy(otu_table(Ps_obj_pruned_min), .) ->
#   otu_table(Ps_obj_pruned_rared)

adonis(
  vegdist(otu_table(Ps_obj_pruned_rared), method = "bray") ~ Field + Sample.type * Season,
  data =
    get_variable(Ps_obj_pruned_rared),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_rared), method = "bray") ~      Field + Sample.type * Season, data = get_variable(Ps_obj_pruned_rared),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field                3     1.944  0.6481   4.208 0.04457  0.001 ***
    ## Sample.type          2    12.310  6.1552  39.963 0.28219  0.001 ***
    ## Season               3     7.903  2.6345  17.104 0.18117  0.001 ***
    ## Sample.type:Season   6     4.833  0.8054   5.229 0.11077  0.001 ***
    ## Residuals          108    16.635  0.1540         0.38130           
    ## Total              122    43.625                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotLibDist(Ps_obj_pruned_rared)
```

![](03_Normalisation_DADA2_figures/test%20rarefaction-1.png)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_rared), "matrix"))
```

![](03_Normalisation_DADA2_figures/test%20rarefaction%20diag%20plots-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_rared))) > 0)
meanSdPlot(as.matrix(log2(t(otu_table(Ps_obj_pruned_rared))[notAllZero, ] + 1)))
```

![](03_Normalisation_DADA2_figures/test%20rarefaction%20diag%20plots-2.png)<!-- -->

``` r
Rared_ord <- ordinate(Ps_obj_pruned_rared, "CAP", "bray", formula = Ps_obj_pruned_rared ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_rared, Rared_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field")  +
  facet_wrap(~ Season) #+
```

![](03_Normalisation_DADA2_figures/rarefaction%20-%20ordinate-1.png)<!-- -->

``` r
  # geom_point(size = 5) 
```

#### GMPR ([Chen and Chen 2017](#ref-chen_gmpr:_2017))

``` r
Ps_obj_pruned_GMPR <- Ps_obj_pruned_min
Ps_obj_pruned_min %>%
  otu_table(.) %>%
  t() %>%
  as(., "matrix") %>%
  GMPR() ->
  GMPR_factors
```

    ## Begin GMPR size factor calculation ...
    ## 50 
    ## 100 
    ## Completed!
    ## Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers!

``` r
Ps_obj_pruned_min %>%
  otu_table(.) %>%
  t() %*% diag(1 / GMPR_factors$gmpr) %>%
  t() %>%
  as.data.frame(., row.names = sample_names(Ps_obj_pruned_min)) %>%
  otu_table(., taxa_are_rows = FALSE) ->
  otu_table(Ps_obj_pruned_GMPR)
sample_data(Ps_obj_pruned_GMPR)$Library.size <- sample_sums(Ps_obj_pruned_GMPR)

adonis(
  vegdist(otu_table(Ps_obj_pruned_GMPR), method = "bray") ~ Library.size,
  data =
    get_variable(Ps_obj_pruned_GMPR),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_GMPR), method = "bray") ~      Library.size, data = get_variable(Ps_obj_pruned_GMPR), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Library.size   1     2.092  2.0921  6.1696 0.04852  0.001 ***
    ## Residuals    121    41.031  0.3391         0.95148           
    ## Total        122    43.123                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(
  vegdist(otu_table(Ps_obj_pruned_GMPR), method = "bray") ~ Field + Sample.type * Season,
  data =
    get_variable(Ps_obj_pruned_GMPR),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_GMPR), method = "bray") ~      Field + Sample.type * Season, data = get_variable(Ps_obj_pruned_GMPR),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field                3     1.973  0.6575   4.532 0.04574  0.001 ***
    ## Sample.type          2    12.528  6.2642  43.172 0.29052  0.001 ***
    ## Season               3     8.126  2.7086  18.667 0.18843  0.001 ***
    ## Sample.type:Season   6     4.826  0.8043   5.543 0.11190  0.001 ***
    ## Residuals          108    15.671  0.1451         0.36340           
    ## Total              122    43.123                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotLibDist(Ps_obj_pruned_GMPR)
```

![](03_Normalisation_DADA2_figures/GMPR-1.png)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_GMPR), "matrix"))
```

![](03_Normalisation_DADA2_figures/GMPR%20diag%20plots-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_GMPR))) > 0)
meanSdPlot(as.matrix(log2(t(otu_table(Ps_obj_pruned_GMPR))[notAllZero, ] + 1)))
```

![](03_Normalisation_DADA2_figures/GMPR%20diag%20plots-2.png)<!-- -->

``` r
GMPR_ord <- ordinate(Ps_obj_pruned_GMPR, "CAP", "bray", formula = Ps_obj_pruned_GMPR ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_GMPR, GMPR_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field") +
  facet_wrap(~ Season)
```

![](03_Normalisation_DADA2_figures/GMPR%20-%20ordinate-1.png)<!-- -->

#### Cumulative sum scaling normalization ([Paulson *et al.* 2013](#ref-paulson_differential_2013))

``` r
# Cumulative sum scaling normalization
Ps_obj_pruned_CSS <- Ps_obj_pruned_min

Ps_obj_pruned_CSS %>%
  otu_table(.) %>%
  t() %>%
  as(., "matrix") %>%
  newMRexperiment(.) ->
  mr_obj
p <- cumNormStatFast(mr_obj)
cumNormMat(mr_obj, p = p) %>% 
  otu_table(., taxa_are_rows = TRUE) ->
  otu_table(Ps_obj_pruned_CSS)

# Did any OTU produce Na?
Ps_obj_pruned_CSS %>% 
  otu_table() %>% 
  as(., "matrix") %>% 
  t(.) %>% 
  apply(., 2, function(x) !any(is.na(x))) %>% 
  unlist() ->
  OTUs2keep
Ps_obj_pruned_CSS <- prune_taxa(OTUs2keep, Ps_obj_pruned_CSS)
sample_data(Ps_obj_pruned_CSS)$Library.size <- sample_sums(Ps_obj_pruned_CSS)

# adonis(
#   vegdist(otu_table(Ps_obj_pruned_CS), method = "bray") ~ Library.size,
#   data =
#     get_variable(Ps_obj_pruned_CS), "data.frame"),
#   permutations = 999
# )
# adonis(
#   vegdist(otu_table(Ps_obj_pruned_CS), method = "bray") ~ Field + Sample.type * Season,
#   data =
#     get_variable(Ps_obj_pruned_CS), "data.frame"),
#   permutations = 999
# )

PlotLibDist(Ps_obj_pruned_CSS)
```

![](03_Normalisation_DADA2_figures/cumsum-1.png)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_CSS), "matrix"))
```

![](03_Normalisation_DADA2_figures/cumsum%20diag%20plots-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_CSS))) > 0)
meanSdPlot(as(log2(t(otu_table(Ps_obj_pruned_CSS))[notAllZero, ] + 1), "matrix"))
```

![](03_Normalisation_DADA2_figures/cumsum%20diag%20plots-2.png)<!-- -->

``` r
CSS_ord <- ordinate(Ps_obj_pruned_CSS, "CAP", "bray", formula = Ps_obj_pruned_CS ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_CSS, CSS_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field") +
  facet_wrap(~ Season)
```

![](03_Normalisation_DADA2_figures/CS%20-%20ordinate-1.png)<!-- -->

#### Standardize abundances to the median sequencing depth (and convert to proportion)

``` r
Ps_obj_pruned_min %>%
  otu_table(.) %>%
  as(., "matrix") %>%
  rowSums() %>% 
  median() ->
  total
standf = function(x, t=total) round(t * (x / sum(x)))
Ps_obj_pruned_median <- transform_sample_counts(Ps_obj_pruned_min, standf) # Standardize abundances to median sequencing depth
sample_data(Ps_obj_pruned_median)$Library.size <- sample_sums(Ps_obj_pruned_median)

adonis(
  vegdist(otu_table(Ps_obj_pruned_median), method = "bray") ~ Library.size,
  data =
    get_variable(Ps_obj_pruned_median),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_median), method = "bray") ~      Library.size, data = get_variable(Ps_obj_pruned_median),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
    ## Library.size   1     0.956 0.95566  2.7651 0.02234  0.005 **
    ## Residuals    121    41.820 0.34562         0.97766          
    ## Total        122    42.775                 1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(
  vegdist(otu_table(Ps_obj_pruned_median), method = "bray") ~ Field + Sample.type * Season,
  data =
    get_variable(Ps_obj_pruned_median),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_median), method = "bray") ~      Field + Sample.type * Season, data = get_variable(Ps_obj_pruned_median),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field                3     1.927  0.6423   4.584 0.04505  0.001 ***
    ## Sample.type          2    12.861  6.4306  45.889 0.30067  0.001 ***
    ## Season               3     8.066  2.6885  19.185 0.18856  0.001 ***
    ## Sample.type:Season   6     4.787  0.7978   5.693 0.11191  0.001 ***
    ## Residuals          108    15.135  0.1401         0.35382           
    ## Total              122    42.775                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotLibDist(Ps_obj_pruned_median)
```

![](03_Normalisation_DADA2_figures/median-1.png)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_median), "matrix"))
```

![](03_Normalisation_DADA2_figures/median%20diag%20plots-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_median))) > 0)
meanSdPlot(as(log2(t(otu_table(Ps_obj_pruned_median))[notAllZero, ] + 1), "matrix"))
```

![](03_Normalisation_DADA2_figures/median%20diag%20plots-2.png)<!-- -->

``` r
Median_ord <- ordinate(Ps_obj_pruned_median, "CAP", "bray", formula = Ps_obj_pruned_median ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_median, Median_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field") +
  facet_wrap(~ Season)
```

![](03_Normalisation_DADA2_figures/median%20-%20ordinate-1.png)<!-- -->

#### Standardize abundances using log transformation for variance stabilisation

``` r
Ps_obj_pruned_log <- transform_sample_counts(Ps_obj_pruned_min, function(x) log(1 + x))

# Ps_obj_pruned_rlog <- Ps_obj_pruned_min
# Ps_obj_pruned_min %>%
#   transform_sample_counts(., function(x) (1 + x)) %>% # add pseudocount
#   phyloseq_to_deseq2(., ~ Spill) %>%
#   rlog(., blind = TRUE , fitType = "parametric") %>%
#   assay() %>%
#   otu_table(, taxa_are_rows = TRUE) ->
#   otu_table(Ps_obj_pruned_rlog)
sample_data(Ps_obj_pruned_log)$Library.size <- sample_sums(Ps_obj_pruned_log)

adonis(
  vegdist(otu_table(Ps_obj_pruned_log), method = "bray") ~ Library.size,
  data =
    get_variable(Ps_obj_pruned_log),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_log), method = "bray") ~      Library.size, data = get_variable(Ps_obj_pruned_log), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Library.size   1     8.263  8.2626  37.894 0.23849  0.001 ***
    ## Residuals    121    26.383  0.2180         0.76151           
    ## Total        122    34.646                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adonis(
  vegdist(otu_table(Ps_obj_pruned_log), method = "bray") ~ Field + Sample.type * Season,
  data =
    get_variable(Ps_obj_pruned_log),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_log), method = "bray") ~      Field + Sample.type * Season, data = get_variable(Ps_obj_pruned_log),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field                3     1.640  0.5465   4.822 0.04733  0.001 ***
    ## Sample.type          2    10.161  5.0804  44.823 0.29328  0.001 ***
    ## Season               3     7.247  2.4156  21.312 0.20917  0.001 ***
    ## Sample.type:Season   6     3.357  0.5595   4.937 0.09690  0.001 ***
    ## Residuals          108    12.241  0.1133         0.35333           
    ## Total              122    34.646                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotLibDist(Ps_obj_pruned_log)
```

![](03_Normalisation_DADA2_figures/log-1.png)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_log), "matrix"), b.width = 0.1)
```

![](03_Normalisation_DADA2_figures/log%20diag%20plots-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_log))) > 0)
meanSdPlot(as(log2(t(otu_table(Ps_obj_pruned_log))[notAllZero, ] + 1), "matrix"))
```

![](03_Normalisation_DADA2_figures/log%20diag%20plots-2.png)<!-- -->

``` r
Rlog_ord <- ordinate(Ps_obj_pruned_log, "CAP", "bray", formula = Ps_obj_pruned_log ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_log, Rlog_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field") +
  facet_wrap(~ Season)
```

![](03_Normalisation_DADA2_figures/log%20-%20ordinate-1.png)<!-- -->

### Standardize abundances using Centered-Log-Ratio transformation for variance stabilisation ([Fernandes *et al.* 2013](#ref-fernandes_anova-like_2013))

``` r
Ps_obj_pruned_CLR <- phyloseq_CLR(Ps_obj_pruned_min)
```

    ## No. corrected values:  41013

``` r
sample_data(Ps_obj_pruned_CLR)$Library.size <- sample_sums(Ps_obj_pruned_CLR)

qplot(rowSums(otu_table(Ps_obj_pruned_CLR)), geom = "histogram") + 
  xlab("Library size")
```

![](03_Normalisation_DADA2_figures/clr-1.png)<!-- -->

``` r
adonis(
  vegdist(otu_table(Ps_obj_pruned_CLR), method = "euclidean") ~ Library.size,
  data =
    get_variable(Ps_obj_pruned_CLR),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_CLR), method = "euclidean") ~      Library.size, data = get_variable(Ps_obj_pruned_CLR), permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
    ## Library.size   1     25053   25053  1.2649 0.01035  0.197
    ## Residuals    121   2396557   19806         0.98965       
    ## Total        122   2421610                 1.00000

``` r
adonis(
  vegdist(otu_table(Ps_obj_pruned_CLR), method = "euclidean") ~ Field + Sample.type * Season,
  data =
    get_variable(Ps_obj_pruned_CLR),
  permutations = 999
)
```

    ## 
    ## Call:
    ## adonis(formula = vegdist(otu_table(Ps_obj_pruned_CLR), method = "euclidean") ~      Field + Sample.type * Season, data = get_variable(Ps_obj_pruned_CLR),      permutations = 999) 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## Field                3     94783   31594   3.236 0.03914  0.001 ***
    ## Sample.type          2    660984  330492  33.849 0.27295  0.001 ***
    ## Season               3    396633  132211  13.541 0.16379  0.001 ***
    ## Sample.type:Season   6    214726   35788   3.665 0.08867  0.001 ***
    ## Residuals          108   1054484    9764         0.43545           
    ## Total              122   2421610                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PlotLibDist(Ps_obj_pruned_CLR)
```

![](03_Normalisation_DADA2_figures/clr-2.png)<!-- -->

``` r
PlotReadHist(as(otu_table(Ps_obj_pruned_CLR), "matrix"), b.width = 0.1)
```

![](03_Normalisation_DADA2_figures/clr%20diag%20plots-1.png)<!-- -->

``` r
notAllZero <- (rowSums(t(otu_table(Ps_obj_pruned_CLR))) > 0)
meanSdPlot(as(log2(t(otu_table(Ps_obj_pruned_CLR))[notAllZero, ] + 1), "matrix"))
```

![](03_Normalisation_DADA2_figures/clr%20diag%20plots-2.png)<!-- -->

``` r
CLR_ord <- ordinate(Ps_obj_pruned_CLR, "RDA", formula = Ps_obj_pruned_CLR ~ Field + Sample.type * Season)
plot_ordination(Ps_obj_pruned_CLR, CLR_ord, type = "samples", color = "Sample.type", label = "Sample.name", shape = "Field") +
  facet_wrap(~ Season)
```

![](03_Normalisation_DADA2_figures/clr%20-%20ordinate-1.png)<!-- -->

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
 package        * version date       lib source        
 ade4             1.7-16  2020-10-28 [1] CRAN (R 4.0.2)
 affy             1.68.0  2020-10-27 [1] Bioconductor  
 affyio           1.60.0  2020-10-27 [1] Bioconductor  
 ape            * 5.5     2021-04-25 [1] CRAN (R 4.0.3)
 assertthat       0.2.1   2019-03-21 [1] CRAN (R 4.0.2)
 backports        1.2.1   2020-12-09 [1] CRAN (R 4.0.2)
 Biobase        * 2.50.0  2020-10-27 [1] Bioconductor  
 BiocGenerics   * 0.36.1  2021-04-16 [1] Bioconductor  
 BiocManager      1.30.15 2021-05-11 [1] CRAN (R 4.0.3)
 biomformat       1.18.0  2020-10-27 [1] Bioconductor  
 Biostrings     * 2.58.0  2020-10-27 [1] Bioconductor  
 bit              4.0.4   2020-08-04 [1] CRAN (R 4.0.2)
 bit64            4.0.5   2020-08-30 [1] CRAN (R 4.0.2)
 bitops           1.0-7   2021-04-24 [1] CRAN (R 4.0.3)
 blob             1.2.1   2020-01-20 [1] CRAN (R 4.0.2)
 broom          * 0.7.6   2021-04-05 [1] CRAN (R 4.0.3)
 cachem           1.0.5   2021-05-15 [1] CRAN (R 4.0.3)
 caTools          1.18.2  2021-03-28 [1] CRAN (R 4.0.3)
 cellranger       1.1.0   2016-07-27 [1] CRAN (R 4.0.2)
 cli              2.5.0   2021-04-26 [1] CRAN (R 4.0.3)
 clipr            0.7.1   2020-10-08 [1] CRAN (R 4.0.2)
 cluster          2.1.2   2021-04-17 [1] CRAN (R 4.0.3)
 codetools        0.2-18  2020-11-04 [1] CRAN (R 4.0.2)
 colorspace       2.0-1   2021-05-04 [1] CRAN (R 4.0.3)
 cowplot        * 1.1.1   2020-12-30 [1] CRAN (R 4.0.2)
 crayon           1.4.1   2021-02-08 [1] CRAN (R 4.0.3)
 data.table       1.14.0  2021-02-21 [1] CRAN (R 4.0.3)
 DBI              1.1.1   2021-01-15 [1] CRAN (R 4.0.3)
 dbplyr           2.1.1   2021-04-06 [1] CRAN (R 4.0.3)
 DECIPHER       * 2.18.1  2020-10-29 [1] Bioconductor  
 desc             1.3.0   2021-03-05 [1] CRAN (R 4.0.3)
 details          0.2.1   2020-01-12 [1] CRAN (R 4.0.2)
 digest           0.6.27  2020-10-24 [1] CRAN (R 4.0.2)
 docxtools      * 0.2.2   2020-06-03 [1] CRAN (R 4.0.3)
 doParallel     * 1.0.16  2020-10-16 [1] CRAN (R 4.0.2)
 dplyr          * 1.0.6   2021-05-05 [1] CRAN (R 4.0.3)
 ellipsis         0.3.2   2021-04-29 [1] CRAN (R 4.0.3)
 evaluate         0.14    2019-05-28 [1] CRAN (R 4.0.2)
 extrafont      * 0.17    2014-12-08 [1] CRAN (R 4.0.2)
 extrafontdb      1.0     2012-06-11 [1] CRAN (R 4.0.2)
 fansi            0.4.2   2021-01-15 [1] CRAN (R 4.0.3)
 farver           2.1.0   2021-02-28 [1] CRAN (R 4.0.3)
 fastmap          1.1.0   2021-01-25 [1] CRAN (R 4.0.3)
 fastmatch        1.1-0   2017-01-28 [1] CRAN (R 4.0.2)
 forcats        * 0.5.1   2021-01-27 [1] CRAN (R 4.0.3)
 foreach        * 1.5.1   2020-10-15 [1] CRAN (R 4.0.2)
 fs               1.5.0   2020-07-31 [1] CRAN (R 4.0.2)
 generics         0.1.0   2020-10-31 [1] CRAN (R 4.0.2)
 ggplot2        * 3.3.3   2020-12-30 [1] CRAN (R 4.0.2)
 ggsci          * 2.9     2018-05-14 [1] CRAN (R 4.0.2)
 glmnet         * 4.1-1   2021-02-21 [1] CRAN (R 4.0.3)
 glue             1.4.2   2020-08-27 [1] CRAN (R 4.0.2)
 gplots           3.1.1   2020-11-28 [1] CRAN (R 4.0.2)
 gridExtra      * 2.3     2017-09-09 [1] CRAN (R 4.0.2)
 gtable           0.3.0   2019-03-25 [1] CRAN (R 4.0.2)
 gtools           3.8.2   2020-03-31 [1] CRAN (R 4.0.2)
 haven            2.4.1   2021-04-23 [1] CRAN (R 4.0.3)
 hexbin           1.28.2  2021-01-08 [1] CRAN (R 4.0.2)
 highr            0.9     2021-04-16 [1] CRAN (R 4.0.3)
 hms              1.1.0   2021-05-17 [1] CRAN (R 4.0.3)
 htmltools        0.5.1.1 2021-01-22 [1] CRAN (R 4.0.3)
 httr             1.4.2   2020-07-20 [1] CRAN (R 4.0.2)
 igraph           1.2.6   2020-10-06 [1] CRAN (R 4.0.2)
 IRanges        * 2.24.1  2020-12-12 [1] Bioconductor  
 iterators      * 1.0.13  2020-10-15 [1] CRAN (R 4.0.2)
 jsonlite         1.7.2   2020-12-09 [1] CRAN (R 4.0.2)
 kableExtra     * 1.3.4   2021-02-20 [1] CRAN (R 4.0.3)
 KernSmooth       2.23-20 2021-05-03 [1] CRAN (R 4.0.3)
 knitr          * 1.33    2021-04-24 [1] CRAN (R 4.0.3)
 labeling         0.4.2   2020-10-20 [1] CRAN (R 4.0.2)
 lattice        * 0.20-44 2021-05-02 [1] CRAN (R 4.0.3)
 lifecycle        1.0.0   2021-02-15 [1] CRAN (R 4.0.3)
 limma          * 3.46.0  2020-10-27 [1] Bioconductor  
 locfit           1.5-9.4 2020-03-25 [1] CRAN (R 4.0.2)
 lubridate        1.7.10  2021-02-26 [1] CRAN (R 4.0.3)
 magrittr       * 2.0.1   2020-11-17 [1] CRAN (R 4.0.2)
 MASS           * 7.3-54  2021-05-03 [1] CRAN (R 4.0.3)
 Matrix         * 1.3-3   2021-05-04 [1] CRAN (R 4.0.3)
 matrixStats    * 0.58.0  2021-01-29 [1] CRAN (R 4.0.3)
 memoise          2.0.0   2021-01-26 [1] CRAN (R 4.0.3)
 metagenomeSeq  * 1.32.0  2020-10-27 [1] Bioconductor  
 mgcv             1.8-35  2021-04-18 [1] CRAN (R 4.0.3)
 modelr           0.1.8   2020-05-19 [1] CRAN (R 4.0.2)
 multtest       * 2.46.0  2020-10-27 [1] Bioconductor  
 munsell          0.5.0   2018-06-12 [1] CRAN (R 4.0.2)
 my.tools       * 0.5     2020-09-30 [1] local         
 NADA           * 1.6-1.1 2020-03-22 [1] CRAN (R 4.0.2)
 nlme             3.1-152 2021-02-04 [1] CRAN (R 4.0.3)
 permute        * 0.9-5   2019-03-12 [1] CRAN (R 4.0.2)
 phangorn       * 2.7.0   2021-05-03 [1] CRAN (R 4.0.3)
 phyloseq       * 1.34.0  2020-10-27 [1] Bioconductor  
 pillar           1.6.1   2021-05-16 [1] CRAN (R 4.0.3)
 pkgconfig        2.0.3   2019-09-22 [1] CRAN (R 4.0.2)
 plyr             1.8.6   2020-03-03 [1] CRAN (R 4.0.2)
 png              0.1-7   2013-12-03 [1] CRAN (R 4.0.2)
 preprocessCore   1.52.1  2021-01-08 [1] Bioconductor  
 prettyunits      1.1.1   2020-01-24 [1] CRAN (R 4.0.2)
 progress         1.2.2   2019-05-16 [1] CRAN (R 4.0.2)
 purrr          * 0.3.4   2020-04-17 [1] CRAN (R 4.0.2)
 quadprog         1.5-8   2019-11-20 [1] CRAN (R 4.0.2)
 R6               2.5.0   2020-10-28 [1] CRAN (R 4.0.2)
 ragg           * 1.1.2   2021-03-17 [1] CRAN (R 4.0.3)
 RColorBrewer   * 1.1-2   2014-12-07 [1] CRAN (R 4.0.2)
 Rcpp             1.0.6   2021-01-15 [1] CRAN (R 4.0.3)
 readr          * 1.4.0   2020-10-05 [1] CRAN (R 4.0.2)
 readxl           1.3.1   2019-03-13 [1] CRAN (R 4.0.2)
 reprex           2.0.0   2021-04-02 [1] CRAN (R 4.0.3)
 reshape2         1.4.4   2020-04-09 [1] CRAN (R 4.0.2)
 rhdf5            2.34.0  2020-10-27 [1] Bioconductor  
 rhdf5filters     1.2.1   2021-05-03 [1] Bioconductor  
 Rhdf5lib         1.12.1  2021-01-26 [1] Bioconductor  
 rlang            0.4.11  2021-04-30 [1] CRAN (R 4.0.3)
 rmarkdown      * 2.8     2021-05-07 [1] CRAN (R 4.0.3)
 rprojroot        2.0.2   2020-11-15 [1] CRAN (R 4.0.2)
 RSQLite        * 2.2.7   2021-04-22 [1] CRAN (R 4.0.3)
 rstudioapi       0.13    2020-11-12 [1] CRAN (R 4.0.2)
 Rttf2pt1         1.3.8   2020-01-10 [1] CRAN (R 4.0.2)
 rvest            1.0.0   2021-03-09 [1] CRAN (R 4.0.3)
 S4Vectors      * 0.28.1  2020-12-09 [1] Bioconductor  
 scales         * 1.1.1   2020-05-11 [1] CRAN (R 4.0.2)
 sessioninfo      1.1.1   2018-11-05 [1] CRAN (R 4.0.2)
 shape            1.4.6   2021-05-19 [1] CRAN (R 4.0.3)
 stringi          1.6.2   2021-05-17 [1] CRAN (R 4.0.3)
 stringr        * 1.4.0   2019-02-10 [1] CRAN (R 4.0.2)
 survival       * 3.2-11  2021-04-26 [1] CRAN (R 4.0.3)
 svglite        * 2.0.0   2021-02-20 [1] CRAN (R 4.0.3)
 systemfonts      1.0.2   2021-05-11 [1] CRAN (R 4.0.3)
 textshaping      0.3.4   2021-05-11 [1] CRAN (R 4.0.3)
 tibble         * 3.1.2   2021-05-16 [1] CRAN (R 4.0.3)
 tidyr          * 1.1.3   2021-03-03 [1] CRAN (R 4.0.3)
 tidyselect       1.1.1   2021-04-30 [1] CRAN (R 4.0.3)
 tidyverse      * 1.3.1   2021-04-15 [1] CRAN (R 4.0.3)
 truncnorm      * 1.0-8   2018-02-27 [1] CRAN (R 4.0.2)
 utf8             1.2.1   2021-03-12 [1] CRAN (R 4.0.3)
 vctrs            0.3.8   2021-04-29 [1] CRAN (R 4.0.3)
 vegan          * 2.5-7   2020-11-28 [1] CRAN (R 4.0.3)
 viridisLite      0.4.0   2021-04-13 [1] CRAN (R 4.0.3)
 vsn            * 3.58.0  2020-10-27 [1] Bioconductor  
 webshot          0.5.2   2019-11-22 [1] CRAN (R 4.0.2)
 withr            2.4.2   2021-04-18 [1] CRAN (R 4.0.3)
 Wrench           1.8.0   2020-10-27 [1] Bioconductor  
 xfun             0.23    2021-05-15 [1] CRAN (R 4.0.3)
 xml2             1.3.2   2020-04-23 [1] CRAN (R 4.0.2)
 XVector        * 0.30.0  2020-10-27 [1] Bioconductor  
 yaml             2.2.1   2020-02-01 [1] CRAN (R 4.0.2)
 zCompositions  * 1.3.4   2020-03-04 [1] CRAN (R 4.0.2)
 zlibbioc         1.36.0  2020-10-27 [1] Bioconductor  

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```

</details>

<br>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chen_gmpr:_2017" class="csl-entry">

Chen J, Chen L. <span class="nocase">GMPR: A novel normalization method
for microbiome sequencing data</span>. *bioRxiv* 2017:112565.

</div>

<div id="ref-fernandes_anova-like_2013" class="csl-entry">

Fernandes AD, Macklaim JM, Linn TG *et al.* ANOVA-Like Differential
Expression (ALDEx) Analysis for Mixed Population RNA-Seq. *PLOS ONE*
2013;**8**:e67019.

</div>

<div id="ref-paulson_differential_2013" class="csl-entry">

Paulson JN, Stine OC, Bravo HC *et al.* <span
class="nocase">Differential abundance analysis for microbial marker-gene
surveys</span>. *Nat Meth* 2013;**10**:1200–2.

</div>

</div>


# Finemapping pipeline

This repo provides a pipeline for performing finemapping of genome-wide
summary level data. The pipeline is primarily designed for proteins,
specifically many protein GWAS. Finemapping is performed in 4 different
ways:

1.  p-value and LD thresholds (not really “finemapping”)
    - identifies sentinel SNPs based on P-value thresholds and LD
      structure
2.  approximate bayes factor (single causal variant assumption)
    - uses P-values and summary statistics to compute SNP-level
      posterior probabilities
3.  sum of single effects (SuSiE)
    - iteratively models multiple causal variants, adjusting iterations
      until convergence
4.  product inverse-moment prior (Finimom)
    - fits a Bayesian hierarchical model to estimate posterior inclusion
      probabilities (PIPs) and credible sets

The idea of doing it 4 different ways is that you can run the pipeline
once and then decide which method best fits your data and assumptions.

## Limitations

1.  out-of-sample reference panels do not perform as well as in-sample -
    this is both an issue of alignment and coverage and is explained
    nicely
    [here](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html)
2.  

## Requirements

1.  GWAS data should have been downloaded and standardised according to:
    <https://github.com/mattlee821/000_GWAS_formatting>
2.  You should have downloaded and processed your reference panel of
    choice. I use the 1kG (phase 3) and code to do this is
    [here](https://github.com/mattlee821/000_GWAS_formatting/blob/main/references/1000genomes-phase3.sh)
3.  Gene region/window should have been created for each GWAS and saved
    as a single file per GWAS

## Gene region/window

To extract a window for a protein you first need to identify the gene
that codes for the protein and the chromosome and start/end position of
the gene. I have made files with mapping data for
[Somalogic](https://github.com/mattlee821/functions/blob/main/data/mapping_GRCh38_p14_somalogic.rda)
and
[Olink](https://github.com/mattlee821/functions/blob/main/data/mapping_GRCh38_p14_olink.rda)
using the following
[Somalogic](https://github.com/mattlee821/functions/blob/main/inst/scripts/mapping-somalogic.R)
and
[Olink](https://github.com/mattlee821/functions/blob/main/inst/scripts/mapping-olink.R)
specific code.

You then need to make a single file with the following structure:
`filepath ID CHR POS-start POS-end`. For example, we can load the
Somalogic mapping file as `data_map` and then list all GWAS files in a
directory and combine these:

``` r
# load mapping file
data_map <- functions::mapping_GRCh38_p14_somalogic
# create list of GWAS files
data <- data.frame(file = list.files("/data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS/", 
                                     pattern = "gz", all.files = TRUE, full.names = TRUE)) |>
# extract the mapping variable from the GWAS file name
  dplyr::mutate(SeqId = stringr::str_extract(file, "[0-9]+_[0-9]+") |>
                  stringr::str_replace("_", "-")) |> 
# join the mapping data to the GWAS file path
  dplyr::left_join(data_map, by = "SeqId") |>
# keep relevant columns
  dplyr::select(file, UNIPROT, CHR, START_hg38, END_hg38)
# save
write.table(x = data, file = "data/files.txt", sep = "\t", row.names = F, col.names = T, quote = F)
head(data)
```

    ##                                                                                         file
    ## 1 /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//.10903_50_STX8_STX8.txt.gz.I1HxeG
    ## 2      /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//10000_28_CRYBB2_CRBB2.txt.gz
    ## 3         /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//10001_7_RAF1_c_Raf.txt.gz
    ## 4       /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//10003_15_ZNF41_ZNF41.txt.gz
    ## 5         /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//10006_25_ELK1_ELK1.txt.gz
    ## 6      /data/GWAS_data/work/ferkingstad_2021_PMID34857953/GWAS//10008_43_GUCA1A_GUC1A.txt.gz
    ##   UNIPROT CHR START_hg38 END_hg38
    ## 1  Q9UNK0  17    9250471  9576591
    ## 2  P43320  22   25212085 25231870
    ## 3  P04049   3   12582101 12664201
    ## 4  P51814   X   47445178 47483222
    ## 5  P19419   X   47635521 47650604
    ## 6  P43080   6   42173364 42181338

This file can then be used to extract the specified window for each GWAS
using an array job [`window_array.sh`]() or as an individual file
[`window_individual.sh`](). Once you have a gene region/window for a
GWAS you can then perform the finemapping.

## Finemapping pipeline

The pipeline runs in R and can be run as an array or as an individual
job.

### Set-up

The first step is to assign variables based on your data and directory
structure. For example, i want to process some somalogic data from the
Ferkingstad 2021 study, which is a EUR population, and the location of
my reference data (bfile).

``` r
# VARS ====
VAR_data <- "somalogic"
VAR_study <- "ferkingstad_2021_PMID34857953"
VAR_population <- "EUR"
VAR_bfile <- paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/", VAR_population)
FILE <- data.table::fread(paste0("/data/IARC_Biostat/work/finemapping/src/", VAR_study, "/files.txt"), header = FALSE)[[1]][1612]
```

I then load my data and do some minor formatting of the POS column and
p-value column. I also run a check to see if any SNPs have P\<1E-5.

``` r
# data ====
data <- data.table::fread(FILE, header = TRUE, sep = "\t")
data <- data |>
  dplyr::mutate(
    POS = POS38,
    P = dplyr::if_else(
      P == 0,
      {
        min_nonzero <- min(P[P != 0], na.rm = TRUE)
        replacement_value <- ifelse(min_nonzero < 5e-8, 5e-8, min_nonzero)
        replacement_value
      },
      P
    )
  )
## Check if any value in the P column is less than X
if (!any(data$P < 1E-5, na.rm = TRUE)) {
  stop("### STOP: no value in column 'P' <1E-5")
}
```

I then combine the mapping data for the somalogic/olink file so that i
can extract an ID and UniProt for each GWAS. The ID is GWAS specific and
is taken directly from the GWAS file name and as such the
`str_extract()` and `str_replace()` commands below will vary based on
the study

``` r
# VAR_map ====
VAR_map <- if (VAR_data == "somalogic") {
  functions::mapping_GRCh38_p14_somalogic
} else if (VAR_data == "olink") {
  functions::mapping_GRCh38_p14_olink
} else {
  stop("Unsupported VAR_data: must be 'somalogic' or 'olink'")
}

# VAR_id ====
VAR_id <- VAR_map |>
  dplyr::filter(
    stringr::str_detect(SeqId, 
                        gsub("_", "-", sub("(^[^_]+_[^_]+).*", "\\1", unique(data$phenotype))))
  )

if (VAR_data == "somalogic") {
  VAR_id <- VAR_id |>
    dplyr::mutate(id = SeqId)
} else if (VAR_data == "olink") {
  VAR_id <- VAR_id |>
    dplyr::mutate(id = olinkID)
}
head(VAR_id)
```

    ##     SeqId   SomaId UNIPROT Target         TargetFullName EntrezGeneID
    ## 1 2192-63 SL003168  Q9Y4X3  CTACK C-C motif chemokine 27        10850
    ##   EntrezGeneSymbol uniprot_gn_id uniprot_gn_symbol entrezgene_id
    ## 1            CCL27        Q9Y4X3             CCL27         10850
    ##   entrezgene_accession    hgnc_id hgnc_symbol ensembl_gene_id
    ## 1                CCL27 HGNC:10626       CCL27 ENSG00000213927
    ##   external_gene_name   gene_biotype CHR START_hg19 END_hg19 strand_hg19
    ## 1              CCL27 protein_coding   9   34661893 34662689           -
    ##   START_hg38 END_hg38 strand_hg38      id
    ## 1   34661880 34664048           - 2192-63

Finally we load the reference data

``` r
reference <- data.table::fread(paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/stats.stats"))
```

### Alignment check

With all the data loaded and variables created we first start by
checking the alignment of the GWAS SNPs with those of the reference
panel. I checks if the alleles in the GWAS are aligned and flip them
when necessary to ensure consistency with the reference. The function
also computes an LD matrix, performs a kriging-based procedure to adjust
the z-scores, and generates a series of plots to visualize the
alignment. For information `?functions::alignment_check()`

``` r
# alignment check and allele flipping ====
list_data <- tryCatch({
  functions::alignment_check(
    df = data, 
    reference = reference, 
    bfile = VAR_bfile
  )
}, error = function(e) {
  stop("Error with alignment check and allele flipping: ", e$message)
})
```

    ## # START 
    ## ## DF made 
    ## ## LD matrix made 
    ## ## list_df made 
    ## ## plot_alignment made

    ## ## S: 0.2107239 
    ## ### larger S indicates inconsistency between Z and LD matrix

    ## ## kriging_rss made 
    ## ## plot_observed_expected made 
    ## ## Alleles to flip:  15 
    ## ## Starting iteration:  1

    ## ## Alleles to flip:  15 
    ## ## alleles flipped 
    ## ## LD matrix made 
    ## ## list_df made 
    ## ## plot_alignment made

    ## ## S: 0.1397395 
    ## ### larger S indicates inconsistency between Z and LD matrix

    ## ## kriging_rss made 
    ## ## plot_observed_expected made 
    ## ## Starting iteration:  2

    ## ## Alleles to flip:  1 
    ## ## alleles flipped 
    ## ## LD matrix made 
    ## ## list_df made 
    ## ## plot_alignment made

    ## ## S: 0.1352297 
    ## ### larger S indicates inconsistency between Z and LD matrix

    ## ## kriging_rss made 
    ## ## plot_observed_expected made 
    ## ## Starting iteration:  3

    ## No alleles to flip, exiting loop.

``` r
## Check if any value in the P column is less than X
if (!any(list_data$list_df$P < 1E-5, na.rm = TRUE)) {
  stop("### STOP: no value in column 'P' <1E-5")
}
```

This returns a list containing the data (`list_df`), alignemnt plots,
and the lambda estimated for the adjusted z-scores (i.e., the degree of
alignment). Again, we perform a check to see if there are any SNPs with
p\<1E-5. You can see from the below plots that the alignment check was
run 3 times and that the alignment improved after allele flipping but
its stil not great. But it is what it is….

``` r
# plot: alignment ====
tryCatch({
  plot_alignment <- plot_grid(
    plotlist = c(
      list_data$plots[seq(1, length(list_data$plots), by = 2)],  # Odd numbered plots
      list_data$plots[seq(2, length(list_data$plots), by = 2)]   # Even numbered plots
    ),
    nrow = 2
  )
}, error = function(e) {
  message("Error creating alignment plot: ", e$message)
})
plot_alignment
```

<img src="README_files/figure-gfm/plot-alignment-1.png" style="display: block; margin: auto;" />

If all of the above is satisfactory we then create a directory to save
the output of the finemapping steps in:

``` r
# directories ====
VAR_OUT <- paste0("data/finemapping/", VAR_study, "/", VAR_id$id, "_", VAR_id$UNIPROT, "/")
if (!dir.exists(VAR_OUT)) {
  dir.create(VAR_OUT, recursive = TRUE)
}
```

### Finemapping: p-value and LD threshold

This is a simple function that first extracts SNPs at a p-value
threshold and then performs clumping at multiple R2 values. Clumping is
performed locally.

``` r
# finemap: p-value and LD ====
tryCatch({
  # Attempt to run the fine-mapping function
  finemap_p_ld <- list_data$list_df[names(list_data$list_df) != "LD"] |>
    as.data.frame() |>
    rename(rsid = SNP, pval = P) |>
    functions::finemap_pval_LD(
      bfile = VAR_bfile, 
      clump_r2 = c(0.1, 0.01, 0.001),
      pval_threshold = 5E-8
    )
}, error = function(e) {
  # Handle error and extract SNPs with P < 5E-8
  message("Error with finemap p-value and LD: ", e$message)
  
  # Filter SNPs with P < 5E-8
  significant_snps <- list_data$list_df |>
    as.data.frame() |>
    filter(P < 5E-8) |>
    mutate(LD = NA, test = "p_ld") |>
    select(SNP, LD, test)
  
  # Assign fallback table
  if (nrow(significant_snps) > 0) {
    message("Extracting SNPs with P < 5E-8.")
    finemap_p_ld <- significant_snps
  } else {
    message("No SNPs with P < 5E-8.")
    finemap_p_ld <- data.frame(SNP = character(), LD = logical(), test = numeric())
  }
})
finemap_p_ld
```

    ## # A tibble: 3 × 3
    ##   SNP       test     LD
    ##   <chr>     <chr> <dbl>
    ## 1 rs2070074 p_ld  0.1  
    ## 2 rs2070074 p_ld  0.01 
    ## 3 rs2070074 p_ld  0.001

### Finemapping: approximate bayes factors

We perform finemapping under a single causal variant assumption using
approximate bayes factors, though only approximate if only the p-value
is used.

``` r
# finemap: approximate bayes factor ====
tryCatch({
  finemap_colocabf <- list_data$list_df |>
    purrr::set_names(~ case_when(
      . == "SNP" ~ "snp",
      . == "P" ~ "p", 
      . == "BETA" ~ "beta", 
      TRUE ~ .
    )) |>
    coloc::finemap.abf() |>
    tidyr::drop_na() |>
    dplyr::filter(SNP.PP == max(SNP.PP)) |>
    dplyr::select(snp, SNP.PP) |>
    dplyr::mutate(test = "coloc.abf") |>
    dplyr::rename(SNP = snp,
                  PIP = SNP.PP) 
}, error = function(e) {
  message("Error with finemap approximate bayes factor: ", e$message)
})
finemap_colocabf
```

    ##         SNP       PIP      test
    ## 1 rs2070074 0.7742873 coloc.abf

### Finemapping: SuSiE

We perform finemapping for multiple causal variants using the sum of
single effects model.

``` r
# finemap: susie ====
tryCatch({
  # Initial parameters
  max_iter <- 3000
  max_retries <- 2
  retry_count <- 0
  converged <- FALSE
  
  while (!converged && retry_count < max_retries) {
    retry_count <- retry_count + 1
    
    # Try running susie_rss
    finemap_susie <- tryCatch({
      susieR::susie_rss(
        bhat = list_data$list_df$BETA,
        shat = list_data$list_df$SE,
        R = list_data$list_df$LD,
        n = list_data$list_df$N,
        max_iter = max_iter,
        check_prior = TRUE
      )
    }, warning = function(w) {
      # Check if the warning indicates non-convergence
      if (grepl("IBSS algorithm did not converge", w$message)) {
        message("Warning: ", w$message)
        return(NULL)
      } else {
        warning(w) # Propagate other warnings
      }
    }, error = function(e) {
      # Stop if there is an error unrelated to convergence
      stop("Error in susie_rss: ", e$message)
    })
    
    # Check if the model converged
    if (!is.null(finemap_susie)) {
      converged <- TRUE
    } else {
      max_iter <- max_iter + 1000
      message("Retrying with max_iter = ", max_iter)
    }
  }
  
  # If convergence failed after max retries, raise an error
  if (!converged) {
    stop("SusieR failed to converge after ", retry_count, " retries with max_iter = ", max_iter)
  }
  
  # Generate finemap tables if susie_rss succeeded
  finemap_table_susie <- functions::susieR_cs_table(
    susieR_model = finemap_susie, 
    df = data.frame(
      SNP = list_data$list_df$SNP,
      POS = list_data$list_df$POS,
      P = list_data$list_df$P
    )
  )
  
  finemap_table_susie_cs <- finemap_table_susie |>
    filter(!is.na(cs)) |>
    select(-P,-POS)
  
}, error = function(e) {
  message("Error with susie finemap: ", e$message)
})
finemap_table_susie_cs
```

    ## # A tibble: 16 × 6
    ##    SNP            PIP cs    cs_snps                                  label test 
    ##    <chr>        <dbl> <chr> <chr>                                    <chr> <chr>
    ##  1 rs72729401  0.978  L4    <NA>                                     rs72… SuSiE
    ##  2 rs112600419 1      L3    <NA>                                     rs11… SuSiE
    ##  3 rs10814131  0.152  L5    rs2772546, rs2492365, rs2480560, rs2094… <NA>  SuSiE
    ##  4 rs2772546   0.0735 L5    rs10814131, rs2492365, rs2480560, rs209… <NA>  SuSiE
    ##  5 rs2492365   0.0732 L5    rs10814131, rs2772546, rs2480560, rs209… <NA>  SuSiE
    ##  6 rs144781957 0.199  L1    rs17355465, rs61746778, rs72737148, rs1… rs14… SuSiE
    ##  7 rs17355465  0.188  L1    rs144781957, rs61746778, rs72737148, rs… <NA>  SuSiE
    ##  8 rs2480560   0.0732 L5    rs10814131, rs2772546, rs2492365, rs209… <NA>  SuSiE
    ##  9 rs2094722   0.0732 L5    rs10814131, rs2772546, rs2492365, rs248… <NA>  SuSiE
    ## 10 rs1977494   0.377  L5    rs10814131, rs2772546, rs2492365, rs248… rs19… SuSiE
    ## 11 rs58087541  1      L2    <NA>                                     rs58… SuSiE
    ## 12 rs61746778  0.139  L1    rs144781957, rs17355465, rs72737148, rs… <NA>  SuSiE
    ## 13 rs72737148  0.142  L1    rs144781957, rs17355465, rs61746778, rs… <NA>  SuSiE
    ## 14 rs113129079 0.178  L5    rs10814131, rs2772546, rs2492365, rs248… <NA>  SuSiE
    ## 15 rs115291674 0.155  L1    rs144781957, rs17355465, rs61746778, rs… <NA>  SuSiE
    ## 16 rs72737151  0.177  L1    rs144781957, rs17355465, rs61746778, rs… <NA>  SuSiE

### Finemapping: finimom

We perform finemapping for multiple causal variants using product
inverse-moment prior for the effect sizes, and a beta-binomial prior for
the model dimension.

``` r
# finemap: finimom ====
tryCatch({
  # Run the finimom model
  finemap_finimom <- finimom::finimom(
    beta = list_data$list_df$BETA,
    se = list_data$list_df$SE,
    eaf = list_data$list_df$MAF,
    n = list_data$list_df$N,
    R = list_data$list_df$LD,
    insampleLD = FALSE,
    check_ld = TRUE,
    verbose = TRUE, 
    ala = TRUE, 
    cs = TRUE, 
    pip = TRUE,
    clump = TRUE,
    clump_r2 = 0.99^2,
    u = 2.25,
    purity = 0.5
  )
  
  # Create the credible set table
  finemap_table_finimom <- functions::finimom_cs_table(
    finimom_model = finemap_finimom, 
    df = data.frame(
      SNP = list_data$list_df$SNP,
      POS = list_data$list_df$POS,
      P = list_data$list_df$P
    )
  )
  
  # Filter the credible set table
  finemap_table_finimom_cs <- finemap_table_finimom |>
    filter(!is.na(cs)) |>
    select(-POS, -P)
  
}, error = function(e) {
  # Log an error message if anything goes wrong
  message("Error during finimom finemapping: ", e$message)
})
```

    ## Calculating tau based on the sample size.
    ## Clumping variants at r2=0.98
    ## Sampling from the posterior...
    ## 
    ## 12500 iterations done in 4.37 seconds

``` r
finemap_table_finimom_cs
```

    ## # A tibble: 13 × 6
    ##    SNP            PIP cs    cs_snps                                  label test 
    ##    <chr>        <dbl> <chr> <chr>                                    <chr> <chr>
    ##  1 rs112600419 0.783  L3    ""                                       rs11… Fini…
    ##  2 rs2070074   0.431  L1    ""                                       rs20… Fini…
    ##  3 rs144781957 0.0800 L2    "rs17355465, rs61746778, rs72737148, rs… <NA>  Fini…
    ##  4 rs17355465  0.0800 L2    "rs144781957, rs61746778, rs72737148, r… <NA>  Fini…
    ##  5 rs58087541  0.800  L4    ""                                       rs58… Fini…
    ##  6 rs61746778  0.0800 L2    "rs144781957, rs17355465, rs72737148, r… <NA>  Fini…
    ##  7 rs72737148  0.0800 L2    "rs144781957, rs17355465, rs61746778, r… <NA>  Fini…
    ##  8 rs115291674 0.0800 L2    "rs144781957, rs17355465, rs61746778, r… <NA>  Fini…
    ##  9 rs72737151  0.0800 L2    "rs144781957, rs17355465, rs61746778, r… <NA>  Fini…
    ## 10 rs3176823   0.0800 L2    "rs144781957, rs17355465, rs61746778, r… rs31… Fini…
    ## 11 rs72737155  0.0800 L2    "rs144781957, rs17355465, rs61746778, r… <NA>  Fini…
    ## 12 rs3176815   0.0800 L2    "rs144781957, rs17355465, rs61746778, r… <NA>  Fini…
    ## 13 rs3176814   0.0800 L2    "rs144781957, rs17355465, rs61746778, r… <NA>  Fini…

### Finemapping table

We combine all results into a single table

``` r
# finemap: table ====
tryCatch({
  # List all potential data frames
  finemap_candidates <- list(
    finemap_colocabf = if (exists("finemap_colocabf") && !is.null(finemap_colocabf)) finemap_colocabf else NULL,
    finemap_p_ld = if (exists("finemap_p_ld") && !is.null(finemap_p_ld)) finemap_p_ld else NULL,
    finemap_table_susie_cs = if (exists("finemap_table_susie_cs") && !is.null(finemap_table_susie_cs)) finemap_table_susie_cs else NULL,
    finemap_table_finimom_cs = if (exists("finemap_table_finimom_cs") && !is.null(finemap_table_finimom_cs)) finemap_table_finimom_cs else NULL
  )
  
  # Filter out missing or NULL data frames
  finemap_list <- purrr::compact(finemap_candidates)
  
  # Bind rows only for existing data frames
  if (length(finemap_list) > 0) {
    finemap <- dplyr::bind_rows(finemap_list) |>
      dplyr::mutate(id = VAR_id$id, UNIPROT = VAR_id$UNIPROT) |>
      dplyr::left_join(data, by = "SNP") |>
      dplyr::filter(!is.na(SNP))
    finemap_table <- finemap |>
      select(SNP,test,PIP,LD,cs,cs_snps,
             id,UNIPROT,CHR,POS19,POS38,EA,OA,EAF,BETA,SE,P,N,phenotype)
    # Write the final finemap table
    write.table(
      x = finemap_table,
      file = paste0(VAR_OUT, "finemap.txt"),
      row.names = FALSE, col.names = TRUE, sep = "\t"
    )
  } else {
    message("No data frames available to create the finemap table.")
  }
}, error = function(e) {
  message("Error creating finemap table: ", e$message)
})
finemap_table
```

    ##            SNP      test        PIP    LD   cs
    ## 1    rs2070074 coloc.abf 0.77428735    NA <NA>
    ## 2    rs2070074      p_ld         NA 0.100 <NA>
    ## 3    rs2070074      p_ld         NA 0.010 <NA>
    ## 4    rs2070074      p_ld         NA 0.001 <NA>
    ## 5   rs72729401     SuSiE 0.97770229    NA   L4
    ## 6  rs112600419     SuSiE 1.00000000    NA   L3
    ## 7   rs10814131     SuSiE 0.15220641    NA   L5
    ## 8    rs2772546     SuSiE 0.07354534    NA   L5
    ## 9    rs2492365     SuSiE 0.07319726    NA   L5
    ## 10 rs144781957     SuSiE 0.19878891    NA   L1
    ## 11  rs17355465     SuSiE 0.18755998    NA   L1
    ## 12   rs2480560     SuSiE 0.07319726    NA   L5
    ## 13   rs2094722     SuSiE 0.07319726    NA   L5
    ## 14   rs1977494     SuSiE 0.37654642    NA   L5
    ## 15  rs58087541     SuSiE 1.00000000    NA   L2
    ## 16  rs61746778     SuSiE 0.13949287    NA   L1
    ## 17  rs72737148     SuSiE 0.14221531    NA   L1
    ## 18 rs113129079     SuSiE 0.17813249    NA   L5
    ## 19 rs115291674     SuSiE 0.15514887    NA   L1
    ## 20  rs72737151     SuSiE 0.17697351    NA   L1
    ## 21 rs112600419   Finimom 0.78270000    NA   L3
    ## 22   rs2070074   Finimom 0.43090000    NA   L1
    ## 23 rs144781957   Finimom 0.07998000    NA   L2
    ## 24  rs17355465   Finimom 0.07998000    NA   L2
    ## 25  rs58087541   Finimom 0.79980000    NA   L4
    ## 26  rs61746778   Finimom 0.07998000    NA   L2
    ## 27  rs72737148   Finimom 0.07998000    NA   L2
    ## 28 rs115291674   Finimom 0.07998000    NA   L2
    ## 29  rs72737151   Finimom 0.07998000    NA   L2
    ## 30   rs3176823   Finimom 0.07998000    NA   L2
    ## 31  rs72737155   Finimom 0.07998000    NA   L2
    ## 32   rs3176815   Finimom 0.07998000    NA   L2
    ## 33   rs3176814   Finimom 0.07998000    NA   L2
    ##                                                                                                       cs_snps
    ## 1                                                                                                        <NA>
    ## 2                                                                                                        <NA>
    ## 3                                                                                                        <NA>
    ## 4                                                                                                        <NA>
    ## 5                                                                                                        <NA>
    ## 6                                                                                                        <NA>
    ## 7                                          rs2772546, rs2492365, rs2480560, rs2094722, rs1977494, rs113129079
    ## 8                                         rs10814131, rs2492365, rs2480560, rs2094722, rs1977494, rs113129079
    ## 9                                         rs10814131, rs2772546, rs2480560, rs2094722, rs1977494, rs113129079
    ## 10                                                rs17355465, rs61746778, rs72737148, rs115291674, rs72737151
    ## 11                                               rs144781957, rs61746778, rs72737148, rs115291674, rs72737151
    ## 12                                        rs10814131, rs2772546, rs2492365, rs2094722, rs1977494, rs113129079
    ## 13                                        rs10814131, rs2772546, rs2492365, rs2480560, rs1977494, rs113129079
    ## 14                                        rs10814131, rs2772546, rs2492365, rs2480560, rs2094722, rs113129079
    ## 15                                                                                                       <NA>
    ## 16                                               rs144781957, rs17355465, rs72737148, rs115291674, rs72737151
    ## 17                                               rs144781957, rs17355465, rs61746778, rs115291674, rs72737151
    ## 18                                          rs10814131, rs2772546, rs2492365, rs2480560, rs2094722, rs1977494
    ## 19                                                rs144781957, rs17355465, rs61746778, rs72737148, rs72737151
    ## 20                                               rs144781957, rs17355465, rs61746778, rs72737148, rs115291674
    ## 21                                                                                                           
    ## 22                                                                                                           
    ## 23   rs17355465, rs61746778, rs72737148, rs115291674, rs72737151, rs3176823, rs72737155, rs3176815, rs3176814
    ## 24  rs144781957, rs61746778, rs72737148, rs115291674, rs72737151, rs3176823, rs72737155, rs3176815, rs3176814
    ## 25                                                                                                           
    ## 26  rs144781957, rs17355465, rs72737148, rs115291674, rs72737151, rs3176823, rs72737155, rs3176815, rs3176814
    ## 27  rs144781957, rs17355465, rs61746778, rs115291674, rs72737151, rs3176823, rs72737155, rs3176815, rs3176814
    ## 28   rs144781957, rs17355465, rs61746778, rs72737148, rs72737151, rs3176823, rs72737155, rs3176815, rs3176814
    ## 29  rs144781957, rs17355465, rs61746778, rs72737148, rs115291674, rs3176823, rs72737155, rs3176815, rs3176814
    ## 30 rs144781957, rs17355465, rs61746778, rs72737148, rs115291674, rs72737151, rs72737155, rs3176815, rs3176814
    ## 31  rs144781957, rs17355465, rs61746778, rs72737148, rs115291674, rs72737151, rs3176823, rs3176815, rs3176814
    ## 32 rs144781957, rs17355465, rs61746778, rs72737148, rs115291674, rs72737151, rs3176823, rs72737155, rs3176814
    ## 33 rs144781957, rs17355465, rs61746778, rs72737148, rs115291674, rs72737151, rs3176823, rs72737155, rs3176815
    ##         id UNIPROT CHR    POS19    POS38                         EA OA     EAF
    ## 1  2192-63  Q9Y4X3   9 34649442 34649445                          G  A 0.09268
    ## 2  2192-63  Q9Y4X3   9 34649442 34649445                          G  A 0.09268
    ## 3  2192-63  Q9Y4X3   9 34649442 34649445                          G  A 0.09268
    ## 4  2192-63  Q9Y4X3   9 34649442 34649445                          G  A 0.09268
    ## 5  2192-63  Q9Y4X3   9 34094007 34094009                          G  A 0.01353
    ## 6  2192-63  Q9Y4X3   9 34641348 34641351                          A  G 0.02859
    ## 7  2192-63  Q9Y4X3   9 34666176 34666179                          A  G 0.98186
    ## 8  2192-63  Q9Y4X3   9 34668546 34668549                          G  A 0.98187
    ## 9  2192-63  Q9Y4X3   9 34672119 34672122                          C  T 0.98187
    ## 10 2192-63  Q9Y4X3   9 34673555 34673558                          T  C 0.06597
    ## 11 2192-63  Q9Y4X3   9 34674827 34674830                          A  G 0.06597
    ## 12 2192-63  Q9Y4X3   9 34678271 34678274                          T  G 0.98186
    ## 13 2192-63  Q9Y4X3   9 34678489 34678492                          G  A 0.98186
    ## 14 2192-63  Q9Y4X3   9 34680460 34680463                          C  T 0.98171
    ## 15 2192-63  Q9Y4X3   9 34680692 34680695 AAAATAAATAAATAAATAAAAAAATA  ! 0.02956
    ## 16 2192-63  Q9Y4X3   9 34680971 34680974                          T  C 0.06597
    ## 17 2192-63  Q9Y4X3   9 34681243 34681246                          G  A 0.06597
    ## 18 2192-63  Q9Y4X3   9 34681511 34681514                    ATTAAAC  ! 0.01795
    ## 19 2192-63  Q9Y4X3   9 34682243 34682246                          T  C 0.06596
    ## 20 2192-63  Q9Y4X3   9 34683912 34683915                          A  G 0.06596
    ## 21 2192-63  Q9Y4X3   9 34641348 34641351                          A  G 0.02859
    ## 22 2192-63  Q9Y4X3   9 34649442 34649445                          G  A 0.09268
    ## 23 2192-63  Q9Y4X3   9 34673555 34673558                          T  C 0.06597
    ## 24 2192-63  Q9Y4X3   9 34674827 34674830                          A  G 0.06597
    ## 25 2192-63  Q9Y4X3   9 34680692 34680695 AAAATAAATAAATAAATAAAAAAATA  ! 0.02956
    ## 26 2192-63  Q9Y4X3   9 34680971 34680974                          T  C 0.06597
    ## 27 2192-63  Q9Y4X3   9 34681243 34681246                          G  A 0.06597
    ## 28 2192-63  Q9Y4X3   9 34682243 34682246                          T  C 0.06596
    ## 29 2192-63  Q9Y4X3   9 34683912 34683915                          A  G 0.06596
    ## 30 2192-63  Q9Y4X3   9 34685765 34685768                          T  C 0.06596
    ## 31 2192-63  Q9Y4X3   9 34688093 34688096                          A  G 0.06599
    ## 32 2192-63  Q9Y4X3   9 34690636 34690639                          T  C 0.06599
    ## 33 2192-63  Q9Y4X3   9 34690642 34690645                          A  G 0.06599
    ##       BETA       SE           P     N           phenotype
    ## 1  -0.0768 0.013835 2.83700e-08 35350 2192_63_CCL27_CTACK
    ## 2  -0.0768 0.013835 2.83700e-08 35350 2192_63_CCL27_CTACK
    ## 3  -0.0768 0.013835 2.83700e-08 35350 2192_63_CCL27_CTACK
    ## 4  -0.0768 0.013835 2.83700e-08 35350 2192_63_CCL27_CTACK
    ## 5  -0.1430 0.035286 5.06500e-05 35349 2192_63_CCL27_CTACK
    ## 6   0.0069 0.023691 7.70858e-01 35350 2192_63_CCL27_CTACK
    ## 7   0.0904 0.029844 2.45300e-03 35350 2192_63_CCL27_CTACK
    ## 8   0.0910 0.029845 2.29500e-03 35350 2192_63_CCL27_CTACK
    ## 9   0.0910 0.029851 2.30000e-03 35350 2192_63_CCL27_CTACK
    ## 10 -0.0629 0.016083 9.19400e-05 35349 2192_63_CCL27_CTACK
    ## 11 -0.0629 0.016089 9.25300e-05 35349 2192_63_CCL27_CTACK
    ## 12  0.0910 0.029851 2.30000e-03 35349 2192_63_CCL27_CTACK
    ## 13  0.0910 0.029851 2.30000e-03 35349 2192_63_CCL27_CTACK
    ## 14  0.0969 0.029749 1.12500e-03 35349 2192_63_CCL27_CTACK
    ## 15  0.0037 0.023590 8.75366e-01 35349 2192_63_CCL27_CTACK
    ## 16 -0.0628 0.016094 9.53600e-05 35349 2192_63_CCL27_CTACK
    ## 17 -0.0628 0.016092 9.51500e-05 35349 2192_63_CCL27_CTACK
    ## 18 -0.0921 0.029911 2.07600e-03 35349 2192_63_CCL27_CTACK
    ## 19 -0.0628 0.016083 9.43600e-05 35349 2192_63_CCL27_CTACK
    ## 20 -0.0629 0.016095 9.30800e-05 35349 2192_63_CCL27_CTACK
    ## 21  0.0069 0.023691 7.70858e-01 35350 2192_63_CCL27_CTACK
    ## 22 -0.0768 0.013835 2.83700e-08 35350 2192_63_CCL27_CTACK
    ## 23 -0.0629 0.016083 9.19400e-05 35349 2192_63_CCL27_CTACK
    ## 24 -0.0629 0.016089 9.25300e-05 35349 2192_63_CCL27_CTACK
    ## 25  0.0037 0.023590 8.75366e-01 35349 2192_63_CCL27_CTACK
    ## 26 -0.0628 0.016094 9.53600e-05 35349 2192_63_CCL27_CTACK
    ## 27 -0.0628 0.016092 9.51500e-05 35349 2192_63_CCL27_CTACK
    ## 28 -0.0628 0.016083 9.43600e-05 35349 2192_63_CCL27_CTACK
    ## 29 -0.0629 0.016095 9.30800e-05 35349 2192_63_CCL27_CTACK
    ## 30 -0.0639 0.016096 7.18600e-05 35349 2192_63_CCL27_CTACK
    ## 31 -0.0633 0.016092 8.36700e-05 35349 2192_63_CCL27_CTACK
    ## 32 -0.0631 0.016089 8.78400e-05 35349 2192_63_CCL27_CTACK
    ## 33 -0.0631 0.016089 8.78500e-05 35349 2192_63_CCL27_CTACK

### Plots

``` r
# plot: manhattan ====
tryCatch({
  plot_manhattan <- data.frame(
    P = list_data$list_df$P, 
    log10p = -log10(list_data$list_df$P),
    POS = list_data$list_df$POS) |>
    ggplot(aes(x = POS / 1e6, y = log10p)) +
    geom_point() +
    geom_hline(yintercept = -log10(5E-8)) +
    labs(x = "Position (Mb)", y = "-log10(P)") +  
    cowplot::theme_cowplot() +
    scale_x_continuous(labels = scales::comma)
    tiff(filename = paste0(VAR_OUT, "manhattan.tiff"), 
       width = 600, height = 400, units = "px")
  print(plot_manhattan)
  dev.off()
}, error = function(e) {
  message("Error creating manhattan plot: ", e$message)
})
```

    ## png 
    ##   2

``` r
plot_manhattan
```

<img src="README_files/figure-gfm/plot-manhattan-1.png" style="display: block; margin: auto;" />

``` r
# plot_finimom ====
if (exists("finemap_table_finimom") && !is.null(finemap_table_finimom_cs) &&
    nrow(finemap_table_finimom_cs) > 0) {
  
  # Separate data
  df_no_cs <- finemap_table_finimom |> filter(is.na(cs))
  df_cs    <- finemap_table_finimom |> filter(!is.na(cs))
  
  # Generate the plot
  plot_finimom <- ggplot() +
    # SNPs not in a CS
    geom_point(data = df_no_cs,
               aes(x = POS / 1e6, y = -log10(P)),
               colour = "darkgrey", alpha = 0.8, size = 1) +
    # SNPs in a CS
    geom_point(data = df_cs,
               aes(x = POS / 1e6, y = -log10(P), colour = cs),
               alpha = 1, size = 2) +
    # SNP labels for CS
    ggrepel::geom_text_repel(data = df_cs,
                             aes(x = POS / 1e6, y = -log10(P), label = label),
                             min.segment.length = 0,
                             max.overlaps = 10,
                             box.padding = 0.5,
                             colour = "black") +
    labs(x = "Position (Mb)", y = "-log10(P)", title = "Finimom") +
    guides(colour = guide_legend(title = "Credible set"), size = NULL) +
    cowplot::theme_cowplot() +
    scale_x_continuous(labels = scales::comma)
  # Save to TIFF
  tiff(filename = paste0(VAR_OUT, "finemap-finimom.tiff"),
       width = 800, height = 400, units = "px")
  print(plot_finimom)
  dev.off()
} else {
  message("finemap_table_finimom is either missing or empty. Skipping plot generation.")
}
```

    ## png 
    ##   2

``` r
plot_finimom
```

<img src="README_files/figure-gfm/plot-finimom-1.png" style="display: block; margin: auto;" />

``` r
# plot_susie ====
if (exists("finemap_table_susie") && !is.null(finemap_table_susie_cs) &&
    nrow(finemap_table_susie_cs) > 0) {
  
  # Separate data
  df_no_cs <- finemap_table_susie |> filter(is.na(cs))
  df_cs    <- finemap_table_susie |> filter(!is.na(cs))
  
  # Generate the plot
  plot_susie <- ggplot() +
    # SNPs not in a CS
    geom_point(data = df_no_cs,
               aes(x = POS / 1e6, y = -log10(P)),
               colour = "darkgrey", alpha = 0.8, size = 1) +
    # SNPs in a CS
    geom_point(data = df_cs,
               aes(x = POS / 1e6, y = -log10(P), colour = cs),
               alpha = 1, size = 2) +
    # Add SNP labels for CS
    ggrepel::geom_text_repel(data = df_cs,
                             aes(x = POS / 1e6, y = -log10(P), label = label),
                             min.segment.length = 0,
                             max.overlaps = 100,
                             box.padding = 1,
                             colour = "black") +
    labs(x = "Position (Mb)", y = "-log10(P)", title = "SuSiE") +
    guides(colour = guide_legend(title = "Credible set"), size = NULL) +
    cowplot::theme_cowplot() +
    scale_x_continuous(labels = scales::comma)
  # Save to TIFF
  tiff(filename = paste0(VAR_OUT, "finemap-susie.tiff"),
       width = 800, height = 400, units = "px")
  print(plot_susie)
  dev.off()
} else {
  message("finemap_table_susie is either missing or empty. Skipping plot generation.")
}
```

    ## png 
    ##   2

``` r
plot_susie
```

<img src="README_files/figure-gfm/plot-susie-1.png" style="display: block; margin: auto;" />

``` r
# plot: regional association ====
## Check if 'finemap' exists and is not empty
if (exists("finemap") && !is.null(finemap) && nrow(finemap) > 0) {
  # Define label from unique SNPs in finemap
  label <- c(unique(finemap$SNP))
  label_p <- unique(finemap$SNP[finemap$P < 5E-8])
  plot_subtitle <- paste0(" Finemapped SNPs: ", length(label), "\n Finemapped SNPs & P<5E-8: ", length(label_p))
} else {
  # Set label and subtitle to NULL if 'finemap' is unavailable
  label <- NULL
  plot_subtitle <- "N finemapped SNPs: 0 (No finemap data available)"
}

if (exists("finemap") && !is.null(finemap) && nrow(finemap) > 0) {
  # Prepare data for the plot
  df <- list_data$list_df[names(list_data$list_df) != "LD"] |>
    as.data.frame() 
  
  # Create the plot
  plot_region <- functions::gg_regionplot(
    df = df, 
    rsid = "SNP", 
    chrom = "CHR", 
    pos = "POS", 
    p_value = "P", 
    trait = "phenotype", 
    plot_pvalue_threshold = 0.1, 
    genome_build = "GRCh38", 
    population = VAR_population, 
    plot_title = paste0("Region plot for: id = ", VAR_id$id, "; UNIPROT = ", VAR_id$UNIPROT, "; HGNC = ", VAR_id$hgnc_symbol),
    plot_subtitle = plot_subtitle, 
    label = label,
    labels = FALSE
  )
    tiff(filename = paste0(VAR_OUT, "region.tiff"), 
       width = 800, height = 800, units = "px")
  print(plot_region)
  dev.off()
}
```

    ## png 
    ##   2

``` r
plot_region
```

<img src="README_files/figure-gfm/plot-region-1.png" style="display: block; margin: auto;" />

``` r
# plot: locus ====
## identify snps to plot - cap at 25
if (exists("finemap") && !is.null(finemap) && nrow(finemap) > 0) {
  snps <- finemap |>
    mutate(
      # 1. If PIP is < 0.8, set label = NA
      label = ifelse(PIP < 0.8, NA, label),
      
      # 2. Conditions when 'test' is 'p_ld'
      label = case_when(
        test == "p_ld" & LD == 0.001 ~ SNP,                             # 2.1 LD = 0.001
        test == "p_ld" & LD == 0.01 & !any(LD == 0.001) ~ SNP,          # 2.2 LD = 0.01 if 0.001 is not present
        test == "p_ld" & LD == 0.1 & !any(LD %in% c(0.001, 0.01)) ~ SNP, # 2.3 LD = 0.1 if 0.001 and 0.01 are not present
        test == "coloc.abf" & PIP > 0.8 ~ SNP,                          # 3. test is 'coloc.abf' and PIP > 0.8
        TRUE ~ label                                                    # Otherwise, keep the existing label
      )
    ) |>
    filter(!is.na(label)) |>        
    distinct(SNP, .keep_all = TRUE) |>      
    arrange(P) |>      
    slice_head(n = 25) |>      
    pull(SNP) 
  
  ## make df
  df <- list_data$list_df[names(list_data$list_df) != "LD"] |>
    as.data.frame() |>
    dplyr::rename(rsid = SNP,
                  pos = POS) |>
    dplyr::mutate(logp = -log10(P))
  
  # make plot
  plot_locus <- lapply(snps, function(snp) {
    # Attempt to create the plot with plot_genes = TRUE
    plot <- tryCatch({
      functions::gg_locusplot(
        df = df,
        lead_snp = snp,
        rsid = rsid,
        chrom = CHR,
        pos = pos,
        ref = EA,
        alt = OA,
        effect = BETA,
        std_err = SE,
        trait = phenotype,
        p_value = P,
        plot_distance = ceiling((max(df$pos) - min(df$pos)) / 100000) * 100000,
        population = "ALL", 
        genome_build = "GRCh38",
        plot_genes = TRUE,  # Initially set plot_genes to TRUE
        plot_recombination = TRUE, 
        plot_subtitle = paste0("Sentinel SNP: ", snp),
        label = c(unique(finemap$SNP))
      )
    }, warning = function(w) {
      # Handle warnings related to genes not found
      if (grepl("No genes found in the specified region", w$message)) {
        message("No genes found in the specified region; rerunning with plot_genes = FALSE.")
        functions::gg_locusplot(
          df = df,
          lead_snp = snp,
          rsid = rsid,
          chrom = CHR,
          pos = pos,
          ref = EA,
          alt = OA,
          effect = BETA,
          std_err = SE,
          trait = phenotype,
          p_value = P,
          plot_distance = ceiling((max(df$pos) - min(df$pos)) / 100000) * 100000,
          population = "ALL", 
          genome_build = "GRCh38",
          plot_genes = FALSE,  # Set plot_genes to FALSE if no genes found
          plot_recombination = TRUE, 
          plot_subtitle = paste0("Sentinel SNP: ", snp),
          label = c(unique(finemap$SNP))
        )
      } else {
        warning(w)  # Rethrow any other warnings
      }
    })
    
    # Add title to each plot and return it
    plot <- patchwork::wrap_plots(plot)
    return(plot)
  })
  
  # sizes based on list size
  width = max(1000, length(plot_locus)*150)
  height = max(800, length(plot_locus)*150)
  
  # Combine all locus plots into a single plot
  plot_locus <- cowplot::plot_grid(plotlist = plot_locus)
  
  # Create a title for the plot
  plot_title <- ggdraw() + 
    draw_label(
      paste0("LD for each (capped at 25) finemapped SNP for: id = ", VAR_id$id, "; UNIPROT = ", VAR_id$UNIPROT, "; HGNC = ", VAR_id$hgnc_symbol),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  # Combine title and locus plots
  plot_locus <- plot_grid(
    plot_title, plot_locus,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  # Save the plot as a .tiff file
  tiff(filename = paste0(VAR_OUT, "locus.tiff"), 
       width = width, height = height+10, units = "px")
  print(plot_locus)
  dev.off()
} else {
  message("Finemap data is unavailable or empty; skipping locus plot generation.")
}
```

    ## png 
    ##   2

``` r
plot_locus
```

<img src="README_files/figure-gfm/plot-locus-1.png" style="display: block; margin: auto;" />

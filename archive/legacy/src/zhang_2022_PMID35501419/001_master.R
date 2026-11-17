rm(list = ls())
set.seed(821)

# environment ====
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# CRAN packages
cran_pkgs <- c(
  "dplyr",
  "tidyr",
  "cowplot",
  "ggplot2",
  "viridis",
  "data.table",
  "tidyverse",
  "curl",
  "Rcpp",
  "snpsettest"
)

# GitHub packages (mapping: pkg -> repo)
github_pkgs <- list(
  "functions"       = "mattlee821/functions",
  "Rfast"           = "RfastOfficial/Rfast",
  "susieR"          = "stephenslab/susieR",
  "ieugwasr"        = "MRCIEU/ieugwasr",
  "genetics.binaRies" = "MRCIEU/genetics.binaRies",
  "TwoSampleMR"     = "MRCIEU/TwoSampleMR",
  "plinkbinr"       = "explodecomputer/plinkbinr",
  "coloc"           = "chr1swallace/coloc@main",
  "locuscomparer"   = "boxiangliu/locuscomparer",
  "RcppParallel"    = "RcppCore/RcppParallel",
  "locusplotr"      = "mglev1n/locusplotr"
)

# Install + load CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Install + load GitHub packages
for (pkg in names(github_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing GitHub package: ", pkg, " from ", github_pkgs[[pkg]])
    remotes::install_github(github_pkgs[[pkg]], force = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# VARS ====
VAR_data <- "somalogic"
VAR_study <- "zhang_2022_PMID35501419"
VAR_population <- "EUR"
VAR_bfile <- paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/", VAR_population)

# VAR array ====
FILE <- data.table::fread(paste0("/data/IARC_Biostat/work/finemapping/src/", VAR_study, "/files.txt"), header = FALSE)
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
FILE <- FILE[i, ][[1]]

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

# VAR_map ====
VAR_map <- if (VAR_data == "somalogic") {
  functions::mapping_GRCh38_p14_somalogic
} else if (VAR_data == "olink") {
  functions::mapping_GRCh38_p14_olink
} else {
  stop("Unsupported VAR_data: must be 'somalogic' or 'olink'")
}

# VAR_id ====
VAR_id <- data$phenotype |>
  str_extract("(?<=SeqId_).*") |>  # Extract the part after "SeqId_"
  str_replace_all("_", "-")  |>
  unique()
VAR_id <- VAR_map |>
  dplyr::filter(
    SeqId %in% VAR_id
  )

if (VAR_data == "somalogic") {
  VAR_id <- VAR_id |>
    dplyr::mutate(id = SeqId)
} else if (VAR_data == "olink") {
  VAR_id <- VAR_id |>
    dplyr::mutate(id = olinkID)
}

# reference data ====
reference <- data.table::fread(paste0("/data/GWAS_data/work/references/1000genomes/phase3/", VAR_population, "/stats.stats"))

# alignment check and allele flipping ====
tryCatch({
  list_data <- functions::alignment_check(
    df = data, 
    reference = reference, 
    bfile = VAR_bfile
  )
}, error = function(e) {
  message("Error with alignment check and allele flipping: ", e$message)
})
## Check if any value in the P column is less than X
if (!any(list_data$list_df$P < 1E-5, na.rm = TRUE)) {
  stop("### STOP: no value in column 'P' <1E-5")
}

# directories ====
VAR_OUT <- paste0("/data/IARC_Biostat/work/finemapping/analysis/", VAR_study, "/", VAR_id$id, "_", VAR_id$UNIPROT, "/")
if (!dir.exists(VAR_OUT)) {
  dir.create(VAR_OUT, recursive = TRUE)
}

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

# finemap: approximate bayes factor (approximate only if p used) ====
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

# plot: alignment ====
tryCatch({
  plot_alignment <- plot_grid(
    plotlist = c(
      list_data$plots[seq(1, length(list_data$plots), by = 2)],  # Odd numbered plots
      list_data$plots[seq(2, length(list_data$plots), by = 2)]   # Even numbered plots
    ),
    nrow = 2
  )
  tiff(filename = paste0(VAR_OUT, "alignment.tiff"), 
       width = length(list_data$plots)*200, height = 600, units = "px")
  print(plot_alignment)
  dev.off()
}, error = function(e) {
  message("Error creating alignment plot: ", e$message)
})

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

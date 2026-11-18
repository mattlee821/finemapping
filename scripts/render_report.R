#!/usr/bin/env Rscript
# Render a self-contained report.html for one fine-mapped GWAS purely from
# {metadata,provenance,qc}.json and methods/{p-ld,susie-rss,susie-rrss,carma,polyfun}.{rda,rds}.

parse_options <- function(args) {
  if (length(args) %% 2L != 0L) {
    stop("Options must be supplied as --name value pairs.", call. = FALSE)
  }
  opts <- list()
  for (i in seq.int(1L, length(args), by = 2L)) {
    key <- gsub("-", "_", sub("^--", "", args[[i]]), fixed = TRUE)
    opts[[key]] <- args[[i + 1L]]
  }
  opts
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

html_escape <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

fmt_int <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return("&mdash;")
  format(as.numeric(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_num <- function(x, digits = 3L) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !is.finite(as.numeric(x))) {
    return("&mdash;")
  }
  formatC(as.numeric(x), digits = digits, format = "fg", flag = "#")
}

badge <- function(label, tone = "neutral") {
  sprintf('<span class="badge badge-%s">%s</span>', tone, html_escape(label))
}

load_rda_object <- function(path, obj_name) {
  if (!file.exists(path)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (exists(obj_name, envir = env, inherits = FALSE)) {
    get(obj_name, envir = env, inherits = FALSE)
  } else {
    objs <- ls(env)
    if (length(objs) == 1L) get(objs[[1L]], envir = env) else NULL
  }
}

load_method_object <- function(methods_dir, stem, obj_name) {
  rds_p <- file.path(methods_dir, paste0(stem, ".rds"))
  if (file.exists(rds_p)) return(readRDS(rds_p))
  rda_p <- file.path(methods_dir, paste0(stem, ".rda"))
  if (file.exists(rda_p)) return(load_rda_object(rda_p, obj_name))
  NULL
}

thin_indices <- function(n, keep_mask, max_bg = 2500L) {
  keep_idx <- which(keep_mask)
  bg_idx <- which(!keep_mask)
  if (length(bg_idx) > max_bg) {
    step <- ceiling(length(bg_idx) / max_bg)
    bg_idx <- bg_idx[seq.int(1L, length(bg_idx), by = step)]
  }
  sort(unique(c(bg_idx, keep_idx)))
}

r2_palette <- function(r2_vec) {
  cols <- rep("#64748b", length(r2_vec))
  valid <- !is.na(r2_vec)
  cols[valid & r2_vec < 0.2] <- "#1e3a8a"
  cols[valid & r2_vec >= 0.2 & r2_vec < 0.4] <- "#0284c7"
  cols[valid & r2_vec >= 0.4 & r2_vec < 0.6] <- "#16a34a"
  cols[valid & r2_vec >= 0.6 & r2_vec < 0.8] <- "#f59e0b"
  cols[valid & r2_vec >= 0.8] <- "#dc2626"
  cols
}

render_locus_png_base64 <- function(reg_id, rss_reg, rrss_reg, sig_lp = 7.30103) {
  if (is.null(rss_reg) || is.null(rss_reg$result) ||
      is.null(rss_reg$result$variants)) {
    return(NULL)
  }
  df_rss <- rss_reg$result$variants
  df_rrss <- if (!is.null(rrss_reg) && !is.null(rrss_reg$result)) {
    rrss_reg$result$variants
  } else {
    NULL
  }
  if (nrow(df_rss) == 0L) return(NULL)

  pos_mb <- as.numeric(
    if ("GWAS_POS" %in% names(df_rss)) df_rss$GWAS_POS else df_rss$POS
  ) / 1e6
  lp <- as.numeric(df_rss$LP)
  pip_rss <- as.numeric(df_rss$PIP)
  pip_rrss <- if (!is.null(df_rrss) && nrow(df_rrss) == nrow(df_rss)) {
    as.numeric(df_rrss$PIP)
  } else {
    rep(NA_real_, nrow(df_rss))
  }
  r2 <- if ("LD_R2_TO_LEAD" %in% names(df_rss)) {
    as.numeric(df_rss$LD_R2_TO_LEAD)
  } else {
    rep(NA_real_, nrow(df_rss))
  }
  cs_rss <- if ("credible_set_id" %in% names(df_rss)) {
    as.character(df_rss$credible_set_id)
  } else {
    rep(NA_character_, nrow(df_rss))
  }
  cs_rrss <- if (!is.null(df_rrss) && "credible_set_id" %in% names(df_rrss)) {
    as.character(df_rrss$credible_set_id)
  } else {
    rep(NA_character_, nrow(df_rss))
  }

  important <- (!is.na(lp) & lp >= sig_lp) |
    (!is.na(pip_rss) & pip_rss >= 0.02) |
    (!is.na(pip_rrss) & pip_rrss >= 0.02) |
    !is.na(cs_rss) | !is.na(cs_rrss)
  idx <- thin_indices(nrow(df_rss), important, max_bg = 2500L)

  tmp_png <- tempfile(pattern = "gwasfm_locus_", fileext = ".png")
  on.exit(unlink(tmp_png, force = TRUE), add = TRUE)

  btype <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")
  grDevices::png(
    filename = tmp_png,
    width = 1280,
    height = 760,
    res = 135,
    type = btype
  )
  graphics::layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE),
                   widths = c(1.55, 1), heights = c(1, 1))
  graphics::par(mar = c(3.8, 4.4, 2.5, 1.2), family = "sans", cex.axis = 0.9,
                cex.lab = 0.95, cex.main = 1.0)

  # Panel 1: Regional -log10(P) colored by Pan-UKB LD r^2 to lead
  cols_r2 <- r2_palette(r2[idx])
  ylim_lp <- c(0, max(c(lp[idx], sig_lp + 1), na.rm = TRUE) * 1.08)
  graphics::plot(
    pos_mb[idx], lp[idx],
    pch = 16, cex = 0.75, col = cols_r2,
    xlab = "GRCh38 position (Mb)",
    ylab = expression(-log[10](italic(P))),
    ylim = ylim_lp,
    main = paste0(reg_id, " — Marginal association & Pan-UKB LD r² to lead"),
    bty = "l"
  )
  graphics::abline(h = sig_lp, lty = 2, col = "#64748b")
  lead_i <- which.max(ifelse(is.na(lp), -Inf, lp))
  if (length(lead_i) == 1L) {
    graphics::points(pos_mb[lead_i], lp[lead_i], pch = 18, cex = 1.6, col = "#7c2d12")
  }
  graphics::legend(
    "topright",
    legend = c("r² >= 0.8", "0.6-0.8", "0.4-0.6", "0.2-0.4", "< 0.2"),
    col = c("#dc2626", "#f59e0b", "#16a34a", "#0284c7", "#1e3a8a"),
    pch = 16, pt.cex = 0.9, bty = "n", cex = 0.8, horiz = TRUE
  )

  # Panel 2: SuSiE-RSS vs SuSiE-rRSS PIPs & 95% Credible Sets
  cs_names <- unique( stats::na.omit(c(cs_rss, cs_rrss)) )
  cs_cols <- c("#dc2626", "#2563eb", "#16a34a", "#9333ea", "#ea580c", "#0891b2")
  col_pip <- rep("#94a3b8", length(idx))
  if (length(cs_names) > 0L) {
    for (k in seq_along(cs_names)) {
      in_k <- (!is.na(cs_rss[idx]) & cs_rss[idx] == cs_names[[k]]) |
        (!is.na(cs_rrss[idx]) & cs_rrss[idx] == cs_names[[k]])
      col_pip[in_k] <- cs_cols[((k - 1L) %% length(cs_cols)) + 1L]
    }
  }
  graphics::plot(
    pos_mb[idx], pip_rss[idx],
    pch = 16, cex = 0.85, col = col_pip,
    ylim = c(0, 1.05),
    xlab = "GRCh38 position (Mb)",
    ylab = "Posterior inclusion probability (PIP)",
    main = "SuSiE-RSS (filled) vs SuSiE-rRSS (open triangles) PIPs",
    bty = "l"
  )
  if (any(!is.na(pip_rrss[idx]))) {
    graphics::points(pos_mb[idx], pip_rrss[idx], pch = 2, cex = 0.85, col = col_pip)
  }
  graphics::abline(h = 0.5, lty = 3, col = "#cbd5e1")
  if (length(cs_names) > 0L) {
    graphics::legend(
      "topright",
      legend = c(paste0("95% CS ", cs_names), "RSS", "rRSS"),
      col = c(cs_cols[((seq_along(cs_names) - 1L) %% length(cs_cols)) + 1L],
              "#334155", "#334155"),
      pch = c(rep(15, length(cs_names)), 16, 2),
      bty = "n", cex = 0.8
    )
  }

  # Panel 3: Observed Z vs LD-expected Z diagnostic
  z_obs <- as.numeric(df_rss$EZ)
  z_exp <- if (!is.null(df_rrss) && "EXPECTED_Z" %in% names(df_rrss) &&
               any(!is.na(df_rrss$EXPECTED_Z))) {
    as.numeric(df_rrss$EXPECTED_Z)
  } else if ("EXPECTED_Z" %in% names(df_rss)) {
    as.numeric(df_rss$EXPECTED_Z)
  } else {
    rep(NA_real_, nrow(df_rss))
  }
  if (any(is.finite(z_obs[idx])) && any(is.finite(z_exp[idx]))) {
    lim <- range(c(z_obs[idx], z_exp[idx]), finite = TRUE)
    graphics::plot(
      z_exp[idx], z_obs[idx],
      pch = 16, cex = 0.75, col = col_pip,
      xlim = lim, ylim = lim,
      xlab = "SuSiE LD-expected Z",
      ylab = "Observed GWAS Z",
      main = "Observed vs LD-expected Z",
      bty = "l"
    )
    graphics::abline(0, 1, lty = 2, col = "#475569")
  } else {
    graphics::plot.new()
    graphics::title(main = "Observed vs LD-expected Z (unavailable)")
  }
  grDevices::dev.off()

  raw_bytes <- readBin(tmp_png, what = "raw", n = file.info(tmp_png)$size)
  paste0("data:image/png;base64,", gsub("\\s+", "", jsonlite::base64_enc(raw_bytes)))
}

options <- parse_options(commandArgs(trailingOnly = TRUE))
gwas_dir <- options$gwas_dir
if (is.null(gwas_dir) || !dir.exists(gwas_dir)) {
  stop("Supply an existing --gwas-dir.", call. = FALSE)
}

meta <- jsonlite::fromJSON(file.path(gwas_dir, "metadata.json"),
                           simplifyVector = FALSE)
prov <- jsonlite::fromJSON(file.path(gwas_dir, "provenance.json"),
                           simplifyVector = FALSE)
qc <- jsonlite::fromJSON(file.path(gwas_dir, "qc.json"),
                         simplifyVector = FALSE)

methods_dir <- file.path(gwas_dir, "methods")
pld_obj <- load_method_object(methods_dir, "p-ld", "p_ld")
susie_rss_obj <- load_method_object(methods_dir, "susie-rss", "susie_rss")
susie_rrss_obj <- load_method_object(methods_dir, "susie-rrss", "susie_rrss")
carma_obj <- load_method_object(methods_dir, "carma", "carma")
polyfun_obj <- load_method_object(methods_dir, "polyfun", "polyfun")

status <- as.character(qc$overall_status %||% "UNKNOWN")
status_tone <- switch(
  status,
  PASS = "ok",
  WARN_DIAGNOSTIC_FLAGS = "warn",
  NO_SIGNIFICANT_LOCI = "neutral",
  "bad"
)
anc_status <- as.character(qc$ancestry_ld_status %||% prov$ancestry_ld_status %||% "")
header_badges <- badge(status, status_tone)
if (anc_status == "multi_ancestry_gwas_single_panel_ld") {
  header_badges <- paste(
    header_badges,
    badge(paste0("ALL→", prov$ld_population, " LD"), "warn")
  )
}

loci_list <- qc$loci %||% list()
sig_lp <- as.numeric(prov$significance_threshold_lp %||% 7.30103)

# Build Locus Summary Table rows
locus_rows_html <- if (length(loci_list) == 0L) {
  '<tr><td colspan="13" class="muted">No genome-wide or cis regions passed the significance gate.</td></tr>'
} else {
  paste(vapply(loci_list, function(loc) {
    st <- as.character(loc$status %||% "")
    m_rss <- loc$methods$susie_rss %||% list()
    m_rrss <- loc$methods$susie_rrss %||% list()
    m_carma <- loc$methods$carma %||% list()
    m_pf <- loc$methods$polyfun %||% list()
    m_pld <- loc$methods$p_ld$selected_by_r2 %||% list()
    pld_str <- if (length(m_pld) > 0L) {
      paste(sprintf("%s:%s", names(m_pld), unlist(m_pld)), collapse = " / ")
    } else {
      "&mdash;"
    }
    flags <- character()
    if (isTRUE(loc$mhc_region)) flags <- c(flags, badge("MHC", "warn"))
    if (isTRUE(loc$mhc_core_grc)) flags <- c(flags, badge("GRC-HLA", "warn"))
    if (isTRUE(loc$cis_overlap)) flags <- c(flags, badge("cis", "info"))
    if (isTRUE(loc$split_at_panukb_radius)) flags <- c(flags, badge("9.5Mb-split", "warn"))
    if (isTRUE(loc$sample_size_heterogeneity_warning)) flags <- c(flags, badge("SS-het", "warn"))
    if (identical(m_rss$converged, FALSE)) flags <- c(flags, badge("RSS-nonconv", "bad"))
    if (identical(m_rrss$converged, FALSE)) flags <- c(flags, badge("rRSS-nonconv", "bad"))
    if (identical(m_carma$converged, FALSE)) flags <- c(flags, badge("CARMA-nonconv", "bad"))
    if (identical(m_pf$converged, FALSE)) flags <- c(flags, badge("PolyFun-nonconv", "bad"))
    if (isTRUE(m_rrss$SER_fallback) || isTRUE(m_rrss$R_reliability_flag)) {
      flags <- c(flags, badge("LD-mismatch", "warn"))
    }
    if (!is.null(m_carma$outlier_count) && as.integer(m_carma$outlier_count) > 0L) {
      flags <- c(flags, badge(paste0("CARMA-out:", as.integer(m_carma$outlier_count)), "warn"))
    }
    if (length(flags) == 0L) flags <- badge("clean", "ok")

    sprintf(
      paste0(
        "<tr>",
        "<td class='mono'>%s</td>",
        "<td>%s</td>",
        "<td class='mono'>%s</td>",
        "<td class='num'>%s</td>",
        "<td class='num'>%s / %s</td>",
        "<td class='mono'>%s (%s)</td>",
        "<td class='mono'>%s</td>",
        "<td class='num'>%s</td>",
        "<td class='num'>%s</td>",
        "<td class='num'>%s</td>",
        "<td class='num'>%s</td>",
        "<td>%s</td>",
        "<td class='num'>%ss</td>",
        "</tr>"
      ),
      html_escape(loc$region_id),
      html_escape(loc$region_type),
      html_escape(loc$source_region_grch38),
      fmt_num((as.numeric(loc$width_bp_grch38 %||% NA)) / 1e6, 2L),
      fmt_int(loc$variants_ld_matched),
      fmt_int(loc$variants_in_source_region),
      html_escape(loc$lead_variant_id %||% "-"),
      fmt_num(loc$lead_lp, 2L),
      html_escape(pld_str),
      fmt_int(m_rss$credible_sets),
      fmt_int(m_rrss$credible_sets),
      fmt_int(m_carma$credible_sets),
      fmt_int(m_pf$credible_sets),
      paste(flags, collapse = " "),
      fmt_num(loc$wall_seconds, 2L)
    )
  }, character(1)), collapse = "\n")
}

# Build Credible Set Top Variants Table & Per-Locus Embedded Figures
cs_rows <- character()
figure_blocks <- character()

bayesian_objs <- list(
  list(label = "SuSiE-RSS", obj = susie_rss_obj),
  list(label = "SuSiE-rRSS", obj = susie_rrss_obj),
  list(label = "CARMA", obj = carma_obj),
  list(label = "PolyFun+SuSiE", obj = polyfun_obj)
)
active_bayesian <- Filter(function(x) !is.null(x$obj) && length(x$obj$regions) > 0L, bayesian_objs)
all_reg_ids <- unique(unlist(lapply(active_bayesian, function(x) names(x$obj$regions))))

if (length(all_reg_ids) > 0L) {
  for (reg_id in all_reg_ids) {
    primary_reg <- if (!is.null(susie_rss_obj) && !is.null(susie_rss_obj$regions[[reg_id]])) {
      susie_rss_obj$regions[[reg_id]]
    } else {
      active_bayesian[[1L]]$obj$regions[[reg_id]]
    }
    secondary_reg <- if (!is.null(susie_rrss_obj) && !is.null(susie_rrss_obj$regions[[reg_id]])) {
      susie_rrss_obj$regions[[reg_id]]
    } else if (length(active_bayesian) >= 2L) {
      active_bayesian[[2L]]$obj$regions[[reg_id]]
    } else {
      NULL
    }
    img_uri <- render_locus_png_base64(reg_id, primary_reg, secondary_reg, sig_lp = sig_lp)
    if (!is.null(img_uri)) {
      figure_blocks <- c(
        figure_blocks,
        sprintf(
          '<div class="locus-card"><h3>%s <span class="muted">(%s)</span></h3><img src="%s" alt="%s diagnostic plot"/></div>',
          html_escape(reg_id),
          html_escape(primary_reg$region_type),
          img_uri,
          html_escape(reg_id)
        )
      )
    }

    for (pair in active_bayesian) {
      m_obj <- pair$obj$regions[[reg_id]]
      if (is.null(m_obj) || is.null(m_obj$result) ||
          is.null(m_obj$result$credible_sets) ||
          nrow(m_obj$result$credible_sets) == 0L) {
        next
      }
      cs_df <- m_obj$result$credible_sets
      v_df <- m_obj$result$variants
      for (k in seq_len(nrow(cs_df))) {
        cs_id <- cs_df$credible_set_id[[k]]
        members <- v_df[!is.na(v_df$credible_set_id) & v_df$credible_set_id == cs_id, , drop = FALSE]
        if (nrow(members) == 0L) next
        top_i <- which.max(members$PIP)
        top_v <- members[top_i, , drop = FALSE]
        cs_rows <- c(
          cs_rows,
          sprintf(
            paste0(
              "<tr>",
              "<td class='mono'>%s</td>",
              "<td>%s</td>",
              "<td class='mono'>%s</td>",
              "<td class='num'>%s</td>",
              "<td class='num'>%s</td>",
              "<td class='mono'>%s (%s)</td>",
              "<td class='num'>%s</td>",
              "<td class='num'>%s (%s)</td>",
              "<td class='num'>%s</td>",
              "<td class='num'>%s</td>",
              "</tr>"
            ),
            html_escape(reg_id),
            html_escape(pair$label),
            html_escape(cs_id),
            fmt_int(nrow(members)),
            fmt_num(cs_df$min_abs_corr[[k]], 3L),
            html_escape(top_v$GWAS_VARIANT_ID[[1L]]),
            html_escape(top_v$ID[[1L]] %||% "."),
            fmt_num(top_v$PIP[[1L]], 3L),
            fmt_num(top_v$ES[[1L]], 3L),
            fmt_num(top_v$SE[[1L]], 3L),
            fmt_num(top_v$LP[[1L]], 2L),
            fmt_num(top_v$AF[[1L]] %||% NA, 3L)
          )
        )
      }
    }
  }
}

cs_table_html <- if (length(cs_rows) == 0L) {
  '<tr><td colspan="10" class="muted">No 95% credible sets passed purity filtering.</td></tr>'
} else {
  paste(cs_rows, collapse = "\n")
}

figures_html <- if (length(figure_blocks) == 0L) {
  '<p class="muted">No analyzed loci to plot.</p>'
} else {
  paste(figure_blocks, collapse = "\n")
}

cis_info <- qc$cis_gene %||% list(status = "not_applicable")
mhc_info <- qc$mhc %||% list(status = "no_significant_variants")

html_doc <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>GWASfinemapping Report — %s / %s</title>
<style>
:root {
  --bg: #f8fafc; --card: #ffffff; --fg: #0f172a; --muted: #64748b;
  --border: #e2e8f0; --accent: #0f766e;
  --ok-bg: #dcfce7; --ok-fg: #166534;
  --warn-bg: #fef3c7; --warn-fg: #92400e;
  --bad-bg: #fee2e2; --bad-fg: #991b1b;
  --info-bg: #dbeafe; --info-fg: #1e40af;
  --neu-bg: #e2e8f0; --neu-fg: #334155;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #0f172a; --card: #1e293b; --fg: #f1f5f9; --muted: #94a3b8;
    --border: #334155; --accent: #2dd4bf;
    --ok-bg: #14532d; --ok-fg: #bbf7d0;
    --warn-bg: #78350f; --warn-fg: #fde68a;
    --bad-bg: #7f1d1d; --bad-fg: #fecaca;
    --info-bg: #1e3a8a; --info-fg: #bfdbfe;
    --neu-bg: #334155; --neu-fg: #e2e8f0;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 24px; background: var(--bg); color: var(--fg);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  line-height: 1.45; font-size: 14px;
}
.container { max-width: 1360px; margin: 0 auto; }
header {
  display: flex; justify-content: space-between; align-items: flex-start;
  gap: 16px; margin-bottom: 20px; padding-bottom: 14px; border-bottom: 1px solid var(--border);
}
h1 { margin: 0 0 4px 0; font-size: 22px; }
h2 { margin: 24px 0 10px 0; font-size: 16px; }
h3 { margin: 0 0 10px 0; font-size: 14px; }
.muted { color: var(--muted); }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12.5px; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
.grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
  gap: 12px; margin-bottom: 20px;
}
.stat-card, .panel, .locus-card {
  background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 14px 16px;
}
.stat-label { font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); }
.stat-val { font-size: 20px; font-weight: 600; margin-top: 4px; }
.badge {
  display: inline-block; padding: 2px 8px; border-radius: 999px;
  font-size: 11.5px; font-weight: 600; white-space: nowrap;
}
.badge-ok { background: var(--ok-bg); color: var(--ok-fg); }
.badge-warn { background: var(--warn-bg); color: var(--warn-fg); }
.badge-bad { background: var(--bad-bg); color: var(--bad-fg); }
.badge-info { background: var(--info-bg); color: var(--info-fg); }
.badge-neutral { background: var(--neu-bg); color: var(--neu-fg); }
.table-wrap { overflow-x: auto; }
table { width: 100%%; border-collapse: collapse; font-size: 13px; }
th, td { padding: 8px 10px; border-bottom: 1px solid var(--border); text-align: left; }
th { font-weight: 600; color: var(--muted); font-size: 12px; background: rgba(148, 163, 184, 0.08); }
.locus-card { margin-bottom: 16px; }
.locus-card img { width: 100%%; height: auto; border-radius: 4px; background: #ffffff; }
</style>
</head>
<body>
<div class="container">
  <header>
    <div>
      <h1>%s / %s</h1>
      <div class="muted">Trait: <strong>%s</strong> &middot; Build: GRCh38 &rarr; Pan-UKB %s (N = %s) &middot; Window: &plusmn;%s Mb merged &middot; MHC: %s</div>
    </div>
    <div>%s</div>
  </header>

  <div class="grid">
    <div class="stat-card">
      <div class="stat-label">Significant Variants</div>
      <div class="stat-val">%s</div>
      <div class="muted">Gate: LP &ge; %s</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Analyzed Regions</div>
      <div class="stat-val">%s / %s</div>
      <div class="muted">cis status: %s</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">95%% Credible Sets</div>
      <div class="stat-val">%s / %s / %s / %s</div>
      <div class="muted">RSS / rRSS / CARMA / PolyFun (L = %s)</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">MHC (chr6:25.5–34.0 Mb)</div>
      <div class="stat-val">%s</div>
      <div class="muted">Hits in MHC: %s</div>
    </div>
    <div class="stat-card">
      <div class="stat-label">Runtime</div>
      <div class="stat-val">%ss</div>
      <div class="muted">Self-contained .rda + .json</div>
    </div>
  </div>

  <h2>Region, Liftover &amp; Method Summary</h2>
  <div class="panel table-wrap">
    <table>
      <thead>
        <tr>
          <th>Region ID</th><th>Type</th><th>GRCh38 Span</th><th class="num">Width (Mb)</th>
          <th class="num">Matched / VCF</th><th>Lead (LP)</th><th>p-LD (0.1/0.01/0.001)</th>
          <th class="num">RSS CS</th><th class="num">rRSS CS</th><th class="num">CARMA CS</th><th class="num">PolyFun CS</th><th>Flags</th><th class="num">Time</th>
        </tr>
      </thead>
      <tbody>%s</tbody>
    </table>
  </div>

  <h2>95%% Credible Set Lead Variants</h2>
  <div class="panel table-wrap">
    <table>
      <thead>
        <tr>
          <th>Region ID</th><th>Method</th><th>CS</th><th class="num">CS Size</th>
          <th class="num">Min |r|</th><th>Top Variant (rsID)</th><th class="num">Top PIP</th>
          <th class="num">ES (SE)</th><th class="num">LP</th><th class="num">AF</th>
        </tr>
      </thead>
      <tbody>%s</tbody>
    </table>
  </div>

  <h2>Regional Association, PIP &amp; LD-Mismatch Diagnostics</h2>
  %s
</div>
</body>
</html>
',
  html_escape(prov$study_id), html_escape(prov$gwas_id),
  html_escape(prov$study_id), html_escape(prov$gwas_id),
  html_escape(prov$trait %||% prov$gwas_id),
  html_escape(prov$ld_population),
  fmt_int(prov$ld_reference_sample_size),
  fmt_num((as.numeric(prov$region_flank_bp %||% 1500000)) / 1e6, 2L),
  html_escape(prov$mhc_policy %||% "included_single_region"),
  header_badges,
  fmt_int(qc$genomewide_significant_variants),
  fmt_num(sig_lp, 2L),
  fmt_int(qc$analyzed_regions_count), fmt_int(qc$candidate_regions_total),
  html_escape(cis_info$status %||% "n/a"),
  fmt_int(qc$total_credible_sets_susie_rss),
  fmt_int(qc$total_credible_sets_susie_rrss),
  fmt_int(qc$total_credible_sets_carma),
  fmt_int(qc$total_credible_sets_polyfun),
  fmt_int(prov$max_causal_L %||% 10L),
  html_escape(mhc_info$status %||% "none"),
  fmt_int(mhc_info$significant_variant_count %||% 0L),
  fmt_num(prov$runtime_seconds, 1L),
  locus_rows_html,
  cs_table_html,
  figures_html
)

out_html <- file.path(gwas_dir, "report.html")
writeLines(html_doc, con = out_html, useBytes = TRUE)

read_npy_float64_matrix <- function(npy_path, ids = NULL, prevalidated = FALSE,
                                    scale = 1) {
  if (length(scale) != 1L || !is.numeric(scale) || is.na(scale) ||
      !is.finite(scale) || scale <= 0) {
    stop("`scale` must be one positive finite number.", call. = FALSE)
  }
  scale <- as.numeric(scale)

  con <- file(npy_path, "rb")
  on.exit(close(con), add = TRUE)

  magic <- readBin(con, what = "raw", n = 6L)
  expected_magic <- as.raw(c(0x93, charToRaw("NUMPY")))
  if (length(magic) != 6L || !identical(magic, expected_magic)) {
    stop("Invalid NumPy .npy file header: ", npy_path, call. = FALSE)
  }

  ver <- as.integer(readBin(con, what = "raw", n = 2L))
  header_len <- if (ver[[1L]] == 1L) {
    readBin(con, what = "integer", size = 2L, signed = FALSE, endian = "little")
  } else if (ver[[1L]] %in% c(2L, 3L)) {
    readBin(con, what = "integer", size = 4L, signed = TRUE, endian = "little")
  } else {
    stop("Unsupported .npy format version: ", paste(ver, collapse = "."),
         call. = FALSE)
  }

  header_raw <- readBin(con, what = "raw", n = header_len)
  header <- rawToChar(header_raw)

  descr_match <- regmatches(
    header,
    regexpr("'descr'\\s*:\\s*'[^']+'", header, perl = TRUE)
  )
  if (length(descr_match) != 1L || !grepl("(<f8|\\|f8|=f8)", descr_match)) {
    stop("Expected little-endian float64 ('<f8') .npy matrix: ", npy_path,
         call. = FALSE)
  }

  fortran_match <- regmatches(
    header,
    regexpr("'fortran_order'\\s*:\\s*(True|False)", header, perl = TRUE)
  )
  fortran_order <- length(fortran_match) == 1L && grepl("True", fortran_match)

  shape_match <- regmatches(
    header,
    regexpr("'shape'\\s*:\\s*\\([^\\)]*\\)", header, perl = TRUE)
  )
  if (length(shape_match) != 1L) {
    stop("Could not parse .npy matrix shape: ", npy_path, call. = FALSE)
  }
  dims_str <- gsub(".*\\(([^\\)]*)\\).*", "\\1", shape_match)
  dims <- as.integer(Filter(nzchar, trimws(strsplit(dims_str, ",")[[1L]])))
  if (length(dims) != 2L || anyNA(dims) || any(dims < 0L)) {
    stop("Expected a 2D .npy matrix: ", npy_path, call. = FALSE)
  }

  if (!is.null(ids) && (dims[[1L]] != length(ids) || dims[[2L]] != length(ids))) {
    stop(".npy LD matrix dimensions (", dims[[1L]], "x", dims[[2L]],
         ") do not match companion ID count (", length(ids), ").",
         call. = FALSE)
  }

  n_total <- as.numeric(dims[[1L]]) * as.numeric(dims[[2L]])
  has_c_fast <- is.loaded("gwasfm_frob_and_scale_inplace",
                          PACKAGE = "gwasfinemapping")
  # Base R readBin() coerces `n` to 32-bit integer (max 2^31 - 1 = 2,147,483,647,
  # exceeded at M >= 46,341 variants such as dense MHC). Chunk above 2^30 doubles.
  max_single_read <- 1073741824

  if (scale != 1 && !has_c_fast) {
    values <- numeric(n_total)
    frob_sq <- 0
    chunk_n <- 16777216L
    pos <- 1
    while (pos <= n_total) {
      take <- as.integer(min(chunk_n, n_total - pos + 1))
      ch <- readBin(con, what = "numeric", n = take, size = 8L,
                    endian = "little")
      if (length(ch) != take) {
        stop("Truncated .npy matrix payload: ", npy_path, call. = FALSE)
      }
      frob_sq <- frob_sq + sum(ch * ch)
      ch <- ch * scale
      values[seq.int(pos, length.out = take)] <- ch
      pos <- pos + take
    }
    rm(ch)
  } else if (n_total > max_single_read) {
    values <- numeric(n_total)
    chunk_n <- 67108864L
    pos <- 1
    while (pos <= n_total) {
      take <- as.integer(min(chunk_n, n_total - pos + 1))
      ch <- readBin(con, what = "numeric", n = take, size = 8L,
                    endian = "little")
      if (length(ch) != take) {
        stop("Truncated .npy matrix payload: ", npy_path, call. = FALSE)
      }
      values[seq.int(pos, length.out = take)] <- ch
      pos <- pos + take
    }
    rm(ch)
    frob_sq <- if (scale != 1 && has_c_fast) {
      .Call("gwasfm_frob_and_scale_inplace", values, scale,
            PACKAGE = "gwasfinemapping")
    } else {
      NULL
    }
  } else {
    values <- readBin(
      con,
      what = "numeric",
      n = as.integer(n_total),
      size = 8L,
      endian = "little"
    )
    if (length(values) != n_total) {
      stop("Truncated .npy matrix payload: ", npy_path, call. = FALSE)
    }
    frob_sq <- if (scale != 1 && has_c_fast) {
      .Call("gwasfm_frob_and_scale_inplace", values, scale,
            PACKAGE = "gwasfinemapping")
    } else {
      NULL
    }
  }

  if (fortran_order) {
    attr(values, "dim") <- c(dims[[1L]], dims[[2L]])
    if (!is.null(ids)) {
      attr(values, "dimnames") <- list(ids, ids)
    }
    if (isTRUE(prevalidated)) {
      attr(values, "gwasfm_prevalidated") <- TRUE
    }
    if (scale != 1) {
      attr(values, "ld_scale") <- scale
      attr(values, "frob_sq") <- frob_sq
      d_vec <- rep(scale, dims[[1L]])
      if (!is.null(ids)) names(d_vec) <- ids
      attr(values, "d") <- d_vec
      attr(values, "scaled:scale") <- rep(1, dims[[1L]])
    }
    values
  } else {
    mat <- matrix(values, nrow = dims[[1L]], ncol = dims[[2L]], byrow = TRUE)
    if (!is.null(ids)) {
      dimnames(mat) <- list(ids, ids)
    }
    if (isTRUE(prevalidated)) {
      attr(mat, "gwasfm_prevalidated") <- TRUE
    }
    if (scale != 1) {
      attr(mat, "ld_scale") <- scale
      attr(mat, "frob_sq") <- frob_sq
    }
    mat
  }
}

read_ld_matrix <- function(matrix_path, ids_path = NULL, scale = 1) {
  if (!is.character(matrix_path) || length(matrix_path) != 1L ||
      !file.exists(matrix_path)) {
    stop("LD matrix file does not exist: ", matrix_path, call. = FALSE)
  }

  if (grepl("\\.npy$", matrix_path, ignore.case = TRUE)) {
    resolved_ids <- if (!is.null(ids_path) && nzchar(ids_path)) {
      ids_path
    } else {
      paste0(sub("\\.npy$", "", matrix_path, ignore.case = TRUE), ".ids.txt")
    }
    if (!file.exists(resolved_ids)) {
      stop("Companion variant ID file missing for .npy LD matrix: ",
           resolved_ids, call. = FALSE)
    }
    ids <- readLines(resolved_ids, warn = FALSE, encoding = "UTF-8")
    ids <- ids[nzchar(ids)]
    prep_qc_file <- file.path(dirname(matrix_path), "locus_prep_qc.json")
    val_json_file <- file.path(dirname(matrix_path), "validation.json")
    is_preval <- file.exists(prep_qc_file) || file.exists(val_json_file)
    return(read_npy_float64_matrix(
      matrix_path,
      ids = ids,
      prevalidated = is_preval,
      scale = scale
    ))
  }

  matrix_frame <- utils::read.delim(
    matrix_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ids <- as.character(matrix_frame[[1L]])
  R <- as.matrix(matrix_frame[-1L])
  storage.mode(R) <- "double"
  if (scale != 1) {
    frob_sq <- sum(R * R)
    R <- R * as.numeric(scale)
    attr(R, "ld_scale") <- as.numeric(scale)
    attr(R, "frob_sq") <- frob_sq
  }
  rownames(R) <- ids
  colnames(R) <- names(matrix_frame)[-1L]
  R
}

scale_ld_matrix_inplace <- function(R, scale) {
  scale <- as.numeric(scale)
  if (length(scale) != 1L || is.na(scale) || !is.finite(scale) || scale <= 0) {
    stop("`scale` must be one positive finite number.", call. = FALSE)
  }
  prev_scale <- attr(R, "ld_scale")
  if (!is.null(prev_scale) && is.finite(prev_scale) && abs(prev_scale - scale) < 1e-12) {
    return(R)
  }
  eff_mult <- if (!is.null(prev_scale) && is.finite(prev_scale) && prev_scale > 0) {
    scale / prev_scale
  } else {
    scale
  }
  has_c_fast <- is.loaded("gwasfm_frob_and_scale_inplace", PACKAGE = "gwasfinemapping")
  if (has_c_fast) {
    raw_frob <- .Call("gwasfm_frob_and_scale_inplace", R, eff_mult,
                      PACKAGE = "gwasfinemapping")
    frob_unscaled <- if (!is.null(prev_scale) && is.finite(prev_scale) && prev_scale > 0) {
      raw_frob / (prev_scale * prev_scale)
    } else {
      raw_frob
    }
  } else {
    denom <- if (!is.null(prev_scale) && is.finite(prev_scale) && prev_scale > 0) {
      prev_scale * prev_scale
    } else {
      1
    }
    frob_unscaled <- sum(R * R) / denom
    R <- R * eff_mult
  }
  p <- nrow(R)
  ids <- colnames(R)
  d_vec <- rep(scale, p)
  if (!is.null(ids)) names(d_vec) <- ids
  attr(R, "ld_scale") <- scale
  attr(R, "frob_sq") <- frob_unscaled
  attr(R, "d") <- d_vec
  attr(R, "scaled:scale") <- rep(1, p)
  R
}

read_npy_float64_matrix <- function(npy_path) {
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

  n_total <- as.numeric(dims[[1L]]) * as.numeric(dims[[2L]])
  values <- readBin(
    con,
    what = "numeric",
    n = n_total,
    size = 8L,
    endian = "little"
  )
  if (length(values) != n_total) {
    stop("Truncated .npy matrix payload: ", npy_path, call. = FALSE)
  }

  if (fortran_order) {
    dim(values) <- c(dims[[1L]], dims[[2L]])
    values
  } else {
    matrix(values, nrow = dims[[1L]], ncol = dims[[2L]], byrow = TRUE)
  }
}

read_ld_matrix <- function(matrix_path, ids_path = NULL) {
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
    R <- read_npy_float64_matrix(matrix_path)
    if (nrow(R) != length(ids) || ncol(R) != length(ids)) {
      stop(".npy LD matrix dimensions (", nrow(R), "x", ncol(R),
           ") do not match companion ID count (", length(ids), ").",
           call. = FALSE)
    }
    rownames(R) <- ids
    colnames(R) <- ids
    return(R)
  }

  matrix_frame <- utils::read.delim(
    matrix_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ids <- as.character(matrix_frame[[1L]])
  R <- as.matrix(matrix_frame[-1L])
  storage.mode(R) <- "double"
  rownames(R) <- ids
  colnames(R) <- names(matrix_frame)[-1L]
  R
}

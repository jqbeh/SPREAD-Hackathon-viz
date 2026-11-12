# =====================================================================
# Plasmid Explorer
# R/Shiny version of the Python/Dash dashboard (app_1.py)
#
# Input: one or more recognised CSV/TSV/TXT outputs from supported plasmid
# analysis pipelines (plasmid metadata, AMRFinder, MOB-suite, ISFinder,
# MobileElementFinder, and a pairwise Mash-distance table). Put them in
# ./data (or the folder named in PLASMID_DATA_DIR), or upload them in the
# "Data import" tab. Files are recognised from their column headers and
# joined by plasmid accession. A separate "Consensus network" tab compares
# Pling typing/distances with mge-cluster results, and a "Genetic context" tab
# shows gene neighbourhoods from a Bakta table or a Clinker plot. Both tabs
# also pick up their files from the data folder (uploads take priority).
#
# Install packages once:
#   install.packages(c("shiny", "bslib", "plotly", "DT", "dplyr",
#                      "stringr", "readr", "igraph", "visNetwork", "ggplot2", "gggenes",
#                      "jsonlite"))
#   Optional (only for .gff3 / GenBank annotation input), from Bioconductor:
#   BiocManager::install(c("rtracklayer", "genbankr"))
#
# Run (save this file as PROJECT/app.R, with input files in PROJECT/data/):
#   shiny::runApp("path/to/PROJECT")
# or open app.R in RStudio and click "Run App".
# =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(dplyr)
  library(visNetwork)
  library(ggplot2)
  library(gggenes)
})

options(shiny.maxRequestSize = 200 * 1024^2)   # allow large uploads (200 MB)

`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------
# 1. INPUT FILE LOCATIONS AND COLUMN DEFINITIONS
# ---------------------------------------------------------------------

# Folder containing this script (the equivalent of Python's BASE_DIR). Works
# for runApp(), source() and Rscript; falls back to the working directory.
script_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f) && file.exists(f)) return(dirname(normalizePath(f)))
  }
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}

# Folder of this file when it is open in the RStudio editor (NULL otherwise)
rstudio_file_dir <- function() {
  if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) return(NULL)
  p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  if (length(p) == 1 && nzchar(p)) dirname(normalizePath(p)) else NULL
}

# Backup input folder, used only if no "data" folder with input files is found
# next to app.R. Edit this path if your files live somewhere else.
DEFAULT_DATA_DIR <- "/Users/jbbeh/Downloads/SPREAD Hackathon 2026"

# TRUE if a folder exists and directly contains at least one CSV/TSV/TXT file
has_input_files <- function(dir) {
  dir.exists(dir) &&
    length(list.files(dir, pattern = "\\.(csv|tsv|txt)$", ignore.case = TRUE)) > 0
}

# Data folder: the PLASMID_DATA_DIR environment variable if set; otherwise the
# FIRST of these folders that contains input files:
#   1. "data" next to this script (runApp(), source(), Rscript)
#   2. "data" next to the file open in the RStudio editor
#   3. "data" in the current working directory
#   4. DEFAULT_DATA_DIR, then DEFAULT_DATA_DIR/data
# The app does not change the working directory (getwd()); it only reads from DATA_DIR.
resolve_data_dir <- function() {
  env <- Sys.getenv("PLASMID_DATA_DIR", unset = "")       # name of a setting, not a path
  if (nzchar(env)) return(normalizePath(path.expand(env), mustWork = FALSE))
  
  roots <- unique(c(script_dir(), rstudio_file_dir(), getwd()))
  candidates <- c(file.path(roots, "data"), DEFAULT_DATA_DIR, file.path(DEFAULT_DATA_DIR, "data"))
  candidates <- unique(normalizePath(path.expand(candidates), mustWork = FALSE))
  
  with_files <- candidates[vapply(candidates, has_input_files, logical(1))]
  if (length(with_files)) return(with_files[1])
  
  warning("No folder with CSV/TSV input files found. Looked in:\n  ",
          paste(candidates, collapse = "\n  "),
          "\nPut your files in a 'data' folder next to app.R, or set ",
          "Sys.setenv(PLASMID_DATA_DIR = \"<path to data>\").", call. = FALSE)
  existing <- candidates[dir.exists(candidates)]
  if (length(existing)) existing[1] else candidates[1]
}
DATA_DIR      <- resolve_data_dir()
message("Reading input files from: ", DATA_DIR)
message("  files found: ", paste(list.files(DATA_DIR, pattern = "\\.(csv|tsv|txt)$",
                                            ignore.case = TRUE), collapse = ", "))
# All CSV/TSV/TXT files in DATA_DIR are read; file names do not matter.

STANDARD_METADATA_COLUMNS <- c(
  "plasmid_id", "isolate_id", "species", "mlst", "source_type", "country",
  "rep_type", "length", "Plasmid Size", "Plasmid GC content", "genes", "ARG genes",
  "AMR Class", "MOB type", "MPF present", "MPF type", "IS elements", "Transposons",
  "Virulence Factors", "collection_date"
)
DISTANCE_COLUMNS <- c("plasmid1", "plasmid2", "mash_distance")

# Fields whose values are delimited lists (";" "," or "|")
DELIMITED_FIELDS <- c("rep_type", "genes", "ARG genes", "AMR Class", "MOB type", "MPF type",
                      "IS elements", "Transposons", "Virulence Factors")

# Alternative header names recognised for each standard field
METADATA_ALIASES <- list(
  "plasmid_id"         = c("plasmid id", "plasmid", "plasmid accession", "accession", "contig id"),
  "isolate_id"         = c("isolate id", "isolate", "sample id", "strain id"),
  "species"            = c("organism", "organism name", "host species"),
  "mlst"               = c("st", "sequence type", "mlst st", "host st"),
  "source_type"        = c("source", "host source", "isolation source"),
  "country"            = c("country of origin", "location country"),
  "rep_type"           = c("rep type", "replicon", "replicon type", "replicon types", "rep_type(s)"),
  "length"             = c("length bp", "plasmid length", "plasmid size bp", "size bp", "size"),
  "Plasmid Size"       = c("plasmid size kb", "size kb", "plasmid length kb"),
  "Plasmid GC content" = c("gc", "gc content", "gc percent", "gc content percent"),
  "genes"              = c("amr genes", "resistance genes", "arg genes", "args", "arg"),
  "ARG genes"          = c("amr genes", "resistance genes", "arg genes", "args", "arg", "genes"),
  "AMR Class"          = c("amr class", "amr classes", "resistance class", "antibiotic class", "drug class"),
  "MOB type"           = c("mobility type", "relaxase type", "relaxase_type(s)", "mob"),
  "MPF present"        = c("mpf", "mpf presence", "conjugative machinery present"),
  "MPF type"           = c("mpf type", "mating pair formation type"),
  "IS elements"        = c("is element", "insertion sequences", "insertion sequence"),
  "Transposons"        = c("transposon", "tn elements", "tn"),
  "Virulence Factors"  = c("virulence factor", "virulence genes", "virulence"),
  "collection_date"    = c("date", "isolation date", "collection date")
)

# AMRFinder columns (compared after header_key(): lower case, no spaces/underscores)
AMR_GENE_KEYS  <- c("amrgenes", "amrgene", "arggenes", "arggene", "resistancegenes",
                    "elementsymbol", "genesymbol")
AMR_CLASS_KEYS <- c("drugclasses", "drugclass", "amrclasses", "amrclass", "resistanceclasses",
                    "resistanceclass", "antibioticclasses", "antibioticclass", "class")
# Columns that identify a one-row-per-plasmid AMRFinder table (kept narrow so a
# metadata table with an "ARG genes" or "AMR Class" column is not mistaken for it)
AMR_DETECT_KEYS <- c("amrgenes", "amrgene", "drugclasses", "drugclass")

DISTANCE_ALIASES <- list(
  "plasmid1"      = c("plasmid 1", "plasmid id 1", "source plasmid", "from", "source"),
  "plasmid2"      = c("plasmid 2", "plasmid id 2", "target plasmid", "to", "target"),
  "mash_distance" = c("mash distance", "distance", "mash", "similarity distance")
)

# What each plot needs in order to be drawn (shown in the import summary)
PLOT_REQUIREMENTS <- list(
  "Source type chart"                       = c("source_type"),
  "Country chart"                           = c("country"),
  "Size and GC distribution"                = c("Plasmid Size", "Plasmid GC content"),
  "Species and MLST distribution"           = c("species", "mlst"),
  "Replicon distribution"                   = c("rep_type"),
  "ARG gene frequency"                      = c("ARG genes"),
  "Species \u2192 replicon \u2192 AMR Class Sankey" = c("species", "rep_type", "AMR Class"),
  "MOB counts"                              = c("MOB type"),
  "MPF counts"                              = c("MPF type"),
  "IS element counts"                       = c("IS elements"),
  "Transposon counts"                       = c("Transposons"),
  "Virulence factor counts"                 = c("Virulence Factors")
)

PLOT_COLORS <- c("#167D8D", "#E28E48", "#5A6F9B", "#72A98F", "#C45B67",
                 "#9A7BB5", "#D4B64C", "#4D8AAB", "#D07C5A", "#7A9E45")
# Equivalent of plotly.express' qualitative "Safe" palette
SAFE_COLORS <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                 "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
EXTENDED_COLORS <- c(PLOT_COLORS, SAFE_COLORS)

ALL_VALUE <- "__ALL__"
NA_STRINGS <- c("", "#N/A", "#N/A N/A", "#NA", "-NaN", "-nan", "<NA>", "N/A", "NA",
                "NULL", "NaN", "None", "n/a", "nan", "null")


# ---------------------------------------------------------------------
# 2. SMALL HELPERS
# ---------------------------------------------------------------------

header_key <- function(x) gsub("[^a-z0-9]+", "", tolower(as.character(x)))

# C-locale (code point) ordering, like Python's sorted()
sorted_unique <- function(x) sort(unique(x[!is.na(x)]), method = "radix")

fmt_count <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)

has_usable_values <- function(frame, column) {
  if (is.null(frame) || !column %in% names(frame)) return(FALSE)
  v <- stringr::str_trim(as.character(frame[[column]]))
  any(!is.na(v) & v != "")
}

# Split delimited values on ; , or |  (list of character vectors)
split_multi_values <- function(x) {
  strsplit(coalesce(as.character(x), ""), "\\s*[;,|]\\s*")
}

# One row per value of a delimited column
explode_values <- function(df, col) {
  pieces <- split_multi_values(df[[col]])
  out <- df[rep(seq_len(nrow(df)), lengths(pieces)), , drop = FALSE]
  out[[col]] <- stringr::str_trim(unlist(pieces, use.names = FALSE) %||% character(0))
  out
}

count_delimited_values <- function(x) {
  v <- stringr::str_trim(unlist(split_multi_values(x), use.names = FALSE))
  length(unique(v[!is.na(v) & v != ""]))
}

distinct_count <- function(x) n_distinct(x, na.rm = TRUE)

# One row per plasmid / gene pair
expand_genes <- function(metadata) {
  genes <- explode_values(metadata, "genes")
  genes <- genes[!is.na(genes$genes) & genes$genes != "", ]
  distinct(genes, plasmid_id, genes, .keep_all = TRUE)
}

arg_gene_values <- function(metadata) {
  g <- stringr::str_trim(unlist(split_multi_values(metadata[["ARG genes"]]), use.names = FALSE))
  sorted_unique(g[!is.na(g) & g != ""])
}


# ---------------------------------------------------------------------
# 3. READING AND RECOGNISING PIPELINE FILES
# ---------------------------------------------------------------------

sniff_delimiter <- function(header_line) {
  candidates <- c(",", "\t", ";", "|")
  counts <- vapply(candidates, function(d) {
    lengths(regmatches(header_line, gregexpr(d, header_line, fixed = TRUE)))
  }, integer(1))
  if (max(counts) == 0) "," else candidates[which.max(counts)]
}

# Read CSV/TSV text. Leading "#" preamble lines (e.g. MobileElementFinder) are
# skipped, but a "#FILE," header line is kept. All columns are read as text.
read_table_file <- function(path, filename = path) {
  text <- readr::read_file(path)
  text <- sub("^\ufeff", "", text)
  lines <- strsplit(text, "\r\n|\n|\r")[[1]]
  is_preamble <- function(l) {
    l <- sub("^\\s+", "", l)
    startsWith(l, "#") && !startsWith(toupper(l), "#FILE,")
  }
  while (length(lines) && is_preamble(lines[1])) lines <- lines[-1]
  non_blank <- lines[nzchar(trimws(lines))]
  if (!length(non_blank)) stop("file is empty", call. = FALSE)
  delim <- if (tolower(tools::file_ext(filename)) == "tsv") "\t" else sniff_delimiter(non_blank[1])
  frame <- suppressWarnings(readr::read_delim(
    I(paste(lines, collapse = "\n")), delim = delim,
    col_types = readr::cols(.default = readr::col_character()),
    na = NA_STRINGS, name_repair = "minimal", progress = FALSE,
    show_col_types = FALSE, trim_ws = FALSE, comment = ""
  ))
  as.data.frame(frame, stringsAsFactors = FALSE, check.names = FALSE)
}

# Copy recognised columns into canonical names without dropping extras.
map_input_columns <- function(frame, aliases) {
  if (ncol(frame) == 0) return(list(frame = frame, mappings = character(0)))
  result <- frame
  keyed <- names(frame)
  names(keyed) <- header_key(names(frame))
  keyed <- keyed[!duplicated(names(keyed), fromLast = TRUE)]
  mappings <- character(0)
  for (canonical in names(aliases)) {
    keys <- header_key(c(canonical, aliases[[canonical]]))
    hit <- keys[keys %in% names(keyed)]
    if (length(hit)) {
      source <- keyed[[hit[1]]]
      result[[canonical]] <- frame[[source]]
      if (source != canonical) mappings <- c(mappings, paste0(source, " \u2192 ", canonical))
    }
  }
  list(frame = result, mappings = mappings)
}

# Identify a supported pipeline output from its headers, not its file name.
classify_pipeline_table <- function(frame) {
  keys <- unique(header_key(names(frame)))
  has <- function(...) all(c(...) %in% keys)
  any_of <- function(...) any(c(...) %in% keys)
  
  # Inputs of the "Consensus network" tab (read from the same folder there)
  if (has("plasmid1", "plasmid2", "containmentdistance")) return("consensus_input")
  if (has("samplename", "standardcluster")) return("consensus_input")
  if (has("plasmid", "type") && length(keys) <= 4 && !any_of("accession", "plasmidid")) return("consensus_input")
  
  if (has("plasmid1", "plasmid2") || has("source", "target", "distance")) return("distances")
  if (has("sequence", "start") && any_of("isfamily", "gene")) return("isfinder")
  if (has("contig", "name", "type") && any_of("mgeno", "prediction")) return("mobileelementfinder")
  # AMRFinder: one row per plasmid (amr_genes / drug_classes) ...
  if (has("accession") && any(AMR_DETECT_KEYS %in% keys)) return("amrfinder")
  # ... or native AMRFinderPlus output, one row per gene (Contig id / Element symbol / Class)
  if (any_of("elementsymbol", "genesymbol") && any_of("contigid", "name")) return("amrfinder")
  if (has("accession") && any_of("predictedmobility", "relaxasetypes", "mpftype")) return("mobsuite")
  if (has("accession") && any_of("hostspecies", "hostst", "country", "sourcetype")) return("metadata")
  
  # An already-aggregated metadata table with common headers
  mapped <- map_input_columns(frame, METADATA_ALIASES)$frame
  if ("plasmid_id" %in% names(mapped) &&
      any(vapply(setdiff(STANDARD_METADATA_COLUMNS, "plasmid_id"),
                 function(col) has_usable_values(mapped, col), logical(1)))) {
    return("metadata")
  }
  "unknown"
}

clean_pipeline_values <- function(x) {
  v <- stringr::str_trim(as.character(x))
  # placeholders for "no value", incl. "unresolved"/"unknown" from curated tables
  bad <- is.na(v) | v == "" |
    tolower(v) %in% c("-", "na", "nan", "none", "null", "unresolved", "unknown")
  v[bad] <- NA_character_
  v
}

# Unique values of a group joined with "; " (delimited fields are split first)
combine_group_values <- function(v, is_list) {
  v <- v[!is.na(v)]
  tokens <- if (is_list) {
    t <- stringr::str_trim(unlist(strsplit(v, "[;,|]"), use.names = FALSE))
    t[!is.na(t) & nzchar(t)]
  } else v
  u <- unique(tokens)
  if (length(u)) paste(u, collapse = "; ") else NA_character_
}

# Collapse hit-level output to one row per plasmid
collapse_pipeline_frame <- function(frame, list_columns = character(0)) {
  if (!"plasmid_id" %in% names(frame)) return(data.frame(plasmid_id = character(0)))
  frame$plasmid_id <- stringr::str_trim(as.character(frame$plasmid_id))
  frame <- frame[!is.na(frame$plasmid_id) & frame$plasmid_id != "", , drop = FALSE]
  if (nrow(frame) == 0) return(frame)
  
  ids <- sorted_unique(frame$plasmid_id)
  idx <- split(seq_len(nrow(frame)), factor(frame$plasmid_id, levels = ids))
  out <- data.frame(plasmid_id = ids, stringsAsFactors = FALSE)
  for (col in setdiff(names(frame), "plasmid_id")) {
    vals <- clean_pipeline_values(frame[[col]])
    is_list <- col %in% list_columns
    out[[col]] <- vapply(idx, function(i) combine_group_values(vals[i], is_list),
                         character(1), USE.NAMES = FALSE)
  }
  out
}

# Map one recognised pipeline table into plasmid-level fields
convert_pipeline_table <- function(frame, kind) {
  mapped <- map_input_columns(frame, METADATA_ALIASES)$frame
  keymap <- names(frame)
  names(keymap) <- header_key(names(frame))
  keymap <- keymap[!duplicated(names(keymap), fromLast = TRUE)]
  col_by_key <- function(key) if (key %in% names(keymap)) keymap[[key]] else NA_character_
  
  id_values <- if (kind == "isfinder") {
    c <- col_by_key("sequence"); if (is.na(c)) NULL else frame[[c]]
  } else if (kind == "mobileelementfinder") {
    c <- col_by_key("contig"); if (is.na(c)) NULL else frame[[c]]
  } else mapped[["plasmid_id"]]
  if (is.null(id_values)) {
    stop(kind, " output has no recognized plasmid accession column.", call. = FALSE)
  }
  out <- data.frame(plasmid_id = as.character(id_values), stringsAsFactors = FALSE)
  
  if (kind == "metadata") {
    for (column in STANDARD_METADATA_COLUMNS) {
      if (column %in% names(mapped)) out[[column]] <- as.character(mapped[[column]])
    }
  } else if (kind == "amrfinder") {
    first_col <- function(keys) {
      hit <- keys[keys %in% names(keymap)]
      if (length(hit)) keymap[[hit[1]]] else NA_character_
    }
    gene_col  <- first_col(AMR_GENE_KEYS)
    class_col <- first_col(AMR_CLASS_KEYS)
    type_col  <- first_col("type")
    # (Native AMRFinderPlus: the plasmid ID comes from "Contig id", a plasmid_id alias)
    genes   <- if (!is.na(gene_col))  as.character(frame[[gene_col]])  else rep(NA_character_, nrow(frame))
    classes <- if (!is.na(class_col)) as.character(frame[[class_col]]) else rep(NA_character_, nrow(frame))
    if (!is.na(type_col)) {
      # Per-gene output also lists STRESS and VIRULENCE hits: keep AMR genes as
      # ARG genes, and VIRULENCE genes as Virulence Factors.
      types <- toupper(trimws(as.character(frame[[type_col]])))
      is_amr <- is.na(types) | types == "AMR"
      if (any(types == "VIRULENCE", na.rm = TRUE)) {
        out[["Virulence Factors"]] <- ifelse(!is.na(types) & types == "VIRULENCE", genes, NA_character_)
      }
      genes[!is_amr]   <- NA_character_
      classes[!is_amr] <- NA_character_
    }
    if (!is.na(gene_col))  out[["ARG genes"]] <- genes
    if (!is.na(class_col)) out[["AMR Class"]] <- classes
    if ("ARG genes" %in% names(out)) out[["genes"]] <- out[["ARG genes"]]
  } else if (kind == "mobsuite") {
    for (column in c("rep_type", "length", "Plasmid GC content", "MOB type", "MPF type")) {
      if (column %in% names(mapped)) out[[column]] <- as.character(mapped[[column]])
    }
  } else if (kind == "isfinder") {
    # Use the IS element name (GENE, e.g. IS1216E) so it matches MobileElementFinder;
    # fall back to the IS family only if there is no GENE column.
    value_col <- if ("gene" %in% names(keymap)) col_by_key("gene") else col_by_key("isfamily")
    out[["IS elements"]] <- as.character(frame[[value_col]])
  } else if (kind == "mobileelementfinder") {
    names_v <- as.character(frame[[col_by_key("name")]])
    types   <- tolower(as.character(frame[[col_by_key("type")]]))
    is_mask <- (!is.na(types) & grepl("insertion sequence", types, fixed = TRUE)) |
      (!is.na(names_v) & grepl("^IS\\d", names_v, ignore.case = TRUE, perl = TRUE))
    tn_mask <- (!is.na(types) & grepl("transposon", types, fixed = TRUE)) |
      (!is.na(names_v) & grepl("^Tn", names_v, ignore.case = TRUE, perl = TRUE))
    out[["IS elements"]] <- ifelse(is_mask, names_v, NA_character_)
    # Composite transposons are named cn_<spacer length>_<IS>, which differs per
    # plasmid; group them by their flanking IS (e.g. "Composite (IS1216E)").
    tn_names <- sub("^cn_[0-9]+_(.+)$", "Composite (\\1)", names_v)
    out[["Transposons"]] <- ifelse(tn_mask, tn_names, NA_character_)
  }
  
  collapse_pipeline_frame(out, DELIMITED_FIELDS)
}

# Union of two delimited strings (de-duplicated, "-" ignored)
merge_delimited <- function(a, b) {
  tokens <- c(strsplit(coalesce(a, ""), "[;,|]")[[1]], strsplit(coalesce(b, ""), "[;,|]")[[1]])
  tokens <- stringr::str_trim(tokens)
  tokens <- unique(tokens[nzchar(tokens) & tokens != "-"])
  if (length(tokens)) paste(tokens, collapse = "; ") else NA_character_
}

normalize_metadata <- function(metadata) {
  if (nrow(metadata) == 0 || ncol(metadata) == 0) {
    stop("The metadata file has no rows or columns.", call. = FALSE)
  }
  m <- map_input_columns(metadata, METADATA_ALIASES)
  metadata <- m$frame
  mappings <- m$mappings
  
  for (column in STANDARD_METADATA_COLUMNS) {
    if (!column %in% names(metadata)) metadata[[column]] <- NA_character_
  }
  
  ids <- as.character(metadata$plasmid_id)
  missing_ids <- is.na(ids) | stringr::str_trim(ids) == ""
  if (all(missing_ids)) {
    metadata$plasmid_id <- sprintf("ROW_%05d", seq_len(nrow(metadata)))
    mappings <- c(mappings, "No plasmid identifier found; generated row IDs for this upload")
  } else if (any(missing_ids)) {
    ids[missing_ids] <- sprintf("ROW_%05d", which(missing_ids))
    metadata$plasmid_id <- ids
    mappings <- c(mappings, sprintf("Generated row IDs for %d rows without plasmid IDs",
                                    sum(missing_ids)))
  }
  
  text_columns <- c(
    "plasmid_id", "isolate_id", "species", "mlst", "source_type", "country",
    "rep_type", "genes", "ARG genes", "AMR Class", "MOB type", "MPF present",
    "MPF type", "IS elements", "Transposons", "Virulence Factors"
  )
  for (column in text_columns) metadata[[column]] <- stringr::str_trim(as.character(metadata[[column]]))
  
  for (column in c("length", "Plasmid Size", "Plasmid GC content")) {
    metadata[[column]] <- suppressWarnings(as.numeric(metadata[[column]]))
  }
  size <- metadata[["Plasmid Size"]]
  if (sum(!is.na(size)) == 0 && any(!is.na(metadata[["length"]]))) {
    metadata[["Plasmid Size"]] <- metadata[["length"]] / 1000
    mappings <- c(mappings, "length (bp) \u2192 Plasmid Size (kb)")
  } else if (any(size > 1000, na.rm = TRUE)) {
    metadata[["Plasmid Size"]] <- ifelse(!is.na(size) & size <= 1000, size, size / 1000)
    mappings <- c(mappings, "Large Plasmid Size values converted from bp to kb")
  }
  
  duplicates <- sum(duplicated(metadata$plasmid_id) | duplicated(metadata$plasmid_id, fromLast = TRUE))
  notes <- character(0)
  if (duplicates > 0) {
    notes <- c(notes, sprintf("%d rows share plasmid IDs; counts use distinct plasmids where applicable.",
                              duplicates))
  }
  list(metadata = metadata, mappings = mappings, notes = notes)
}

empty_distances <- function() {
  data.frame(plasmid1 = character(0), plasmid2 = character(0),
             mash_distance = numeric(0), stringsAsFactors = FALSE)
}

# An optional pairwise-distance table; incomplete tables disable the network.
normalize_distances <- function(distances) {
  if (is.null(distances) || nrow(distances) == 0) {
    return(list(distances = empty_distances(), mappings = character(0),
                notes = "Pairwise distance data not provided; similarity network is unavailable."))
  }
  m <- map_input_columns(distances, DISTANCE_ALIASES)
  distances <- m$frame
  if (!all(DISTANCE_COLUMNS %in% names(distances))) {
    return(list(distances = empty_distances(), mappings = m$mappings,
                notes = "Distance file is missing plasmid pair or Mash distance fields; similarity network is unavailable."))
  }
  distances <- distances[, DISTANCE_COLUMNS, drop = FALSE]
  distances$plasmid1 <- stringr::str_trim(as.character(distances$plasmid1))
  distances$plasmid2 <- stringr::str_trim(as.character(distances$plasmid2))
  distances$mash_distance <- suppressWarnings(as.numeric(distances$mash_distance))
  valid <- !is.na(distances$plasmid1) & !is.na(distances$plasmid2) &
    distances$plasmid1 != "" & distances$plasmid2 != "" &
    distances$plasmid1 != distances$plasmid2 &
    !is.na(distances$mash_distance) &
    distances$mash_distance >= 0 & distances$mash_distance <= 1
  removed <- sum(!valid)
  list(distances = distances[valid, , drop = FALSE], mappings = m$mappings,
       notes = if (removed > 0) sprintf("Removed %d invalid distance rows.", removed) else character(0))
}

# Recognise, normalise and outer-join any supplied set of pipeline outputs.
# `named_tables` is a list of list(name = <file name>, frame = <data.frame>).
combine_pipeline_tables <- function(named_tables, notes = character(0)) {
  pieces <- list()
  distance_frames <- list()
  recognized <- character(0)
  for (item in named_tables) {
    kind <- classify_pipeline_table(item$frame)
    if (kind == "unknown") {
      notes <- c(notes, sprintf("Skipped %s: its columns did not match a supported input format.",
                                item$name))
      next
    }
    if (kind == "consensus_input") {
      notes <- c(notes, sprintf("%s: Pling / mge-cluster file, used by the Consensus network tab.",
                                item$name))
      next
    }
    recognized <- c(recognized, sprintf("%s (%s)", item$name, kind))
    if (kind == "distances") {
      distance_frames[[length(distance_frames) + 1]] <- item$frame
      next
    }
    pieces[[length(pieces) + 1]] <- convert_pipeline_table(item$frame, kind)
  }
  if (!length(pieces)) {
    stop("No recognized plasmid metadata or pipeline output was found in the supplied files.",
         call. = FALSE)
  }
  
  combined <- pieces[[1]]
  for (piece in pieces[-1]) {
    combined <- merge(combined, piece, by = "plasmid_id", all = TRUE,
                      suffixes = c("", "__new"), sort = FALSE)
    for (column in grep("__new$", names(combined), value = TRUE)) {
      original <- sub("__new$", "", column)
      if (original %in% names(combined)) {
        if (original %in% DELIMITED_FIELDS) {
          a <- combined[[original]]; b <- combined[[column]]
          combined[[original]] <- vapply(seq_along(a), function(i) merge_delimited(a[i], b[i]),
                                         character(1))
        } else {
          combined[[original]] <- ifelse(is.na(combined[[original]]), combined[[column]],
                                         combined[[original]])
        }
        combined[[column]] <- NULL
      } else {
        names(combined)[names(combined) == column] <- original
      }
    }
  }
  combined <- combined[order(combined$plasmid_id, method = "radix"), , drop = FALSE]
  rownames(combined) <- NULL
  
  norm <- normalize_metadata(combined)
  
  if (length(distance_frames)) {
    dist_result <- normalize_distances(dplyr::bind_rows(distance_frames))
  } else {
    dist_result <- normalize_distances(NULL)
  }
  notes <- c(notes, norm$notes, dist_result$notes)
  
  list(
    metadata  = norm$metadata,
    distances = dist_result$distances,
    info = list(
      metadata_mappings = norm$mappings,
      distance_mappings = dist_result$mappings,
      notes = notes,
      recognized_files = recognized,
      metadata_source = if (length(recognized)) paste(recognized, collapse = ", ") else "Uploaded files",
      distance_source = if (length(distance_frames)) "Pairwise distance input" else "Not provided"
    )
  )
}

# Read a set of files (skipping unreadable ones with a note) and combine them.
combine_files <- function(paths, names, extra_tables = list()) {
  tables <- extra_tables
  notes <- character(0)
  for (i in seq_along(paths)) {
    if (is_bakta_table(paths[i])) {
      notes <- c(notes, sprintf("%s: Bakta annotation table, used by the Genetic context tab.", names[i]))
      next
    }
    frame <- tryCatch(read_table_file(paths[i], names[i]), error = function(e) e)
    if (inherits(frame, "error")) {
      notes <- c(notes, sprintf("Skipped %s: could not be read (%s).", names[i], conditionMessage(frame)))
    } else {
      tables[[length(tables) + 1]] <- list(name = names[i], frame = frame)
    }
  }
  combine_pipeline_tables(tables, notes)
}

# Auto-detect and combine all CSV/TSV/TXT files in the data folder.
load_data <- function() {
  if (!dir.exists(DATA_DIR)) {
    stop("Input data folder not found: ", DATA_DIR, call. = FALSE)
  }
  files <- list.files(DATA_DIR, full.names = FALSE)
  files <- sort(files[tolower(tools::file_ext(files)) %in% c("csv", "tsv", "txt") &
                        !dir.exists(file.path(DATA_DIR, files))], method = "radix")
  if (!length(files)) stop("No CSV/TSV input files found in: ", DATA_DIR, call. = FALSE)
  combine_files(file.path(DATA_DIR, files), files)
}

make_store <- function(result) {
  list(metadata = result$metadata, distances = result$distances,
       error = NULL, info = result$info)
}

FALLBACK_INFO <- list(metadata_mappings = character(0), distance_mappings = character(0),
                      notes = character(0), metadata_source = DATA_DIR,
                      distance_source = "Not provided")

load_store <- function() {
  tryCatch(make_store(load_data()), error = function(e) {
    list(metadata = NULL, distances = NULL, error = conditionMessage(e), info = FALLBACK_INFO)
  })
}



# ---------------------------------------------------------------------
# 4. IMPORT SUMMARY TEXT
# ---------------------------------------------------------------------

build_data_status <- function(metadata, distances, store) {
  lines <- c(
    "DATA IMPORT SUMMARY",
    "===================",
    paste0("Metadata rows: ", fmt_count(nrow(metadata))),
    paste0("Distinct plasmids: ", fmt_count(n_distinct(metadata$plasmid_id))),
    if (has_usable_values(metadata, "isolate_id"))
      paste0("Distinct isolates: ", fmt_count(n_distinct(metadata$isolate_id)))
    else "Distinct isolates: Not available",
    paste0("Valid pairwise distances: ", fmt_count(nrow(distances))),
    "",
    "PLOT AVAILABILITY",
    "================="
  )
  for (plot_name in names(PLOT_REQUIREMENTS)) {
    required <- PLOT_REQUIREMENTS[[plot_name]]
    missing <- unique(required[!vapply(required, function(c) has_usable_values(metadata, c), logical(1))])
    lines <- c(lines, paste0(plot_name, ": ",
                             if (!length(missing)) "Available" else paste0("Needs ", paste(missing, collapse = ", "))))
  }
  
  info <- (store %||% list())$info %||% list()
  detected <- STANDARD_METADATA_COLUMNS[vapply(STANDARD_METADATA_COLUMNS,
                                               function(c) has_usable_values(metadata, c), logical(1))]
  lines <- c(lines, "", "DETECTED METADATA FIELDS", "========================",
             if (length(detected)) detected else "No recognized fields found")
  if (length(info$recognized_files)) {
    lines <- c(lines, "", "INPUT FILES RECOGNIZED", "======================", info$recognized_files)
  }
  mappings <- c(info$metadata_mappings, info$distance_mappings)
  if (length(mappings)) lines <- c(lines, "", "COLUMN MAPPINGS", "===============", mappings)
  if (length(info$notes)) lines <- c(lines, "", "DATA NOTES", "==========", info$notes)
  paste(lines, collapse = "\n")
}


# ---------------------------------------------------------------------
# 5. FIGURE HELPERS
# ---------------------------------------------------------------------

# Shared, clean scientific figure styling.
polish_figure <- function(p, title = NULL, legend_title = NULL, show_legend = NULL) {
  p <- layout(
    p,
    paper_bgcolor = "white",
    plot_bgcolor  = "white",
    font = list(family = "Inter, Segoe UI, Arial, sans-serif", size = 12, color = "#263746"),
    margin = list(l = 24, r = 24, t = 48, b = 32),
    hoverlabel = list(bgcolor = "#17324D", bordercolor = "#17324D",
                      font = list(color = "white", size = 12)),
    legend = list(bgcolor = "rgba(255,255,255,0.85)", borderwidth = 0,
                  font = list(size = 11)),
    xaxis = list(showline = FALSE, zeroline = FALSE, gridcolor = "#E9EEF2",
                 ticks = "outside", tickcolor = "#C9D3DC"),
    yaxis = list(showline = FALSE, zeroline = FALSE, gridcolor = "#E9EEF2",
                 ticks = "outside", tickcolor = "#C9D3DC")
  )
  if (!is.null(title)) {
    p <- layout(p, title = list(text = title, font = list(size = 16, color = "#17324D"),
                                x = 0.02, xanchor = "left"))
  }
  if (!is.null(legend_title)) p <- layout(p, legend = list(title = list(text = legend_title)))
  if (!is.null(show_legend))  p <- layout(p, showlegend = show_legend)
  config(p, displaylogo = FALSE)
}

# Placeholder figure with a centred message and no axes.
unavailable_figure <- function(message) {
  plot_ly(type = "scatter", mode = "markers") %>%
    layout(
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(list(
        text = message, showarrow = FALSE, x = 0.5, y = 0.5,
        xref = "paper", yref = "paper", font = list(size = 13, color = "#526875")
      ))
    ) %>%
    polish_figure(show_legend = FALSE)
}

# Horizontal bar chart: one bar (and colour) per category, largest on top.
hbar_figure <- function(categories, values, colors, hover_label, y_title,
                        title = NULL, margin = list(l = 26, r = 24, t = 48, b = 30)) {
  plot_ly(
    x = values, y = categories, type = "bar", orientation = "h",
    marker = list(color = colors, line = list(color = "white", width = 0.8)),
    hovertemplate = paste0(hover_label, ": %{y}<br>Plasmids: %{x}<extra></extra>"),
    showlegend = FALSE
  ) %>%
    polish_figure(title = title, show_legend = FALSE) %>%
    layout(
      margin = margin,
      xaxis = list(title = "Number of plasmids"),
      yaxis = list(title = y_title, type = "category", automargin = TRUE,
                   categoryorder = "array", categoryarray = categories)
    )
}

# Count plasmids carrying each ARG gene.
make_arg_gene_bar_figure <- function(metadata, selected_gene = NULL) {
  if (!has_usable_values(metadata, "ARG genes")) {
    return(unavailable_figure("ARG gene counts need an ARG genes, ARGs, genes, or resistance genes column."))
  }
  genes <- explode_values(
    data.frame(plasmid_id = metadata$plasmid_id, arg = metadata[["ARG genes"]],
               stringsAsFactors = FALSE),
    "arg"
  )
  genes <- distinct(genes[!is.na(genes$arg) & genes$arg != "", ], plasmid_id, arg)
  filtered <- !is.null(selected_gene) && selected_gene != ALL_VALUE
  if (filtered) genes <- genes[genes$arg == selected_gene, ]
  
  counts <- genes %>%
    group_by(arg) %>%
    summarise(plasmids = n_distinct(plasmid_id), .groups = "drop") %>%
    arrange(plasmids, arg)
  
  if (nrow(counts) == 0) return(unavailable_figure("No ARG genes match these filters."))
  
  hbar_figure(
    counts$arg, counts$plasmids, "#167D8D",
    hover_label = "ARG gene", y_title = "ARG gene",
    title = if (filtered) paste0("ARG gene frequency: ", selected_gene) else "ARG gene frequency",
    margin = list(l = 30, r = 38, t = 56, b = 30)
  )
}

# Species -> replicon -> AMR class flows, counted by plasmid/class association.
make_amr_class_sankey <- function(metadata, selected_gene = NULL) {
  needed <- c("species", "rep_type", "AMR Class")
  missing <- needed[!vapply(needed, function(c) has_usable_values(metadata, c), logical(1))]
  if (length(missing)) {
    return(unavailable_figure(paste0("Sankey needs data for: ", paste(missing, collapse = ", "), ".")))
  }
  if (!is.null(selected_gene) && selected_gene != ALL_VALUE) {
    if (!has_usable_values(metadata, "ARG genes")) {
      return(unavailable_figure("Select All; ARG gene data is unavailable for gene filtering."))
    }
    has_gene <- vapply(split_multi_values(metadata[["ARG genes"]]),
                       function(v) selected_gene %in% v, logical(1))
    metadata <- metadata[has_gene, ]
  }
  
  frame <- data.frame(
    plasmid_id = metadata$plasmid_id,
    species    = metadata$species,
    rep_type   = metadata$rep_type,
    amr        = metadata[["AMR Class"]],
    stringsAsFactors = FALSE
  )
  for (col in c("species", "rep_type", "amr")) {
    x <- stringr::str_trim(coalesce(frame[[col]], "Unknown"))
    x[x == ""] <- "Unknown"
    frame[[col]] <- x
  }
  flows <- explode_values(explode_values(frame, "rep_type"), "amr")
  flows <- flows[!is.na(flows$rep_type) & flows$rep_type != "" &
                   !is.na(flows$amr) & flows$amr != "" & flows$amr != "-", ]
  flows <- distinct(flows, plasmid_id, species, rep_type, amr)
  
  if (nrow(flows) == 0) return(unavailable_figure("No AMR class data available."))
  
  species_values <- sorted_unique(flows$species)
  rep_values     <- sorted_unique(flows$rep_type)
  amr_values     <- sorted_unique(flows$amr)
  n_sp  <- length(species_values)
  n_rep <- length(rep_values)
  
  # Sankey node indices are 0-based in plotly.
  sp_idx  <- function(x) match(x, species_values) - 1L
  rep_idx <- function(x) n_sp + match(x, rep_values) - 1L
  amr_idx <- function(x) n_sp + n_rep + match(x, amr_values) - 1L
  
  links1 <- count(flows, species, rep_type, name = "value")
  links2 <- count(flows, rep_type, amr, name = "value")
  
  link_source <- c(sp_idx(links1$species), rep_idx(links2$rep_type))
  link_target <- c(rep_idx(links1$rep_type), amr_idx(links2$amr))
  link_value  <- c(links1$value, links2$value)
  link_color  <- c(rep("rgba(22,125,141,0.22)", nrow(links1)),
                   rep("rgba(226,142,72,0.24)", nrow(links2)))
  
  node_colors <- c(rep("#167D8D", n_sp), rep("#5A6F9B", n_rep),
                   rep("#E28E48", length(amr_values)))
  
  plot_ly(
    type = "sankey", orientation = "h", arrangement = "snap",
    node = list(
      label = c(species_values, rep_values, amr_values),
      color = node_colors, pad = 22, thickness = 18,
      line = list(color = "rgba(23,50,77,0.16)", width = 1),
      hovertemplate = "%{label}<extra></extra>"
    ),
    link = list(
      source = link_source, target = link_target, value = link_value,
      color = link_color,
      hovertemplate = paste0("%{source.label} \u2192 %{target.label}<br>",
                             "%{value} plasmid\u2013AMR class associations<extra></extra>")
    )
  ) %>%
    polish_figure(title = "Species \u2192 replicon type \u2192 AMR class", show_legend = FALSE) %>%
    layout(margin = list(l = 18, r = 18, t = 56, b = 18))
}

# Distinct plasmids per value of a metadata category.
make_category_count_figure <- function(metadata, column, title, category_label) {
  if (!has_usable_values(metadata, column)) {
    return(unavailable_figure(paste0(category_label, " data was not provided.")))
  }
  missing_label <- "Not present / blank"
  frame <- data.frame(plasmid_id = metadata$plasmid_id,
                      category = stringr::str_trim(as.character(metadata[[column]])),
                      stringsAsFactors = FALSE)
  if (column %in% c("rep_type", "MOB type", "MPF type", "IS elements", "Transposons",
                    "Virulence Factors")) {
    frame <- explode_values(frame, "category")
    frame <- frame[!is.na(frame$category) & frame$category != "" & frame$category != "-", ]
  }
  frame$category[is.na(frame$category) | frame$category == ""] <- missing_label
  
  counts <- frame %>%
    distinct(plasmid_id, category) %>%
    group_by(category) %>%
    summarise(plasmids = n_distinct(plasmid_id), .groups = "drop") %>%
    arrange(plasmids, category)
  
  if (nrow(counts) == 0) return(unavailable_figure("No data available."))
  
  hbar_figure(
    counts$category, counts$plasmids, rep_len(PLOT_COLORS, nrow(counts)),
    hover_label = category_label, y_title = category_label, title = title,
    margin = list(l = 26, r = 24, t = 48, b = 30)
  )
}

# Donut chart of distinct isolates (or plasmids when isolate IDs are absent).
# Values differing only by case/spacing are merged under the first spelling
# encountered; blanks become "Unknown".
make_donut_figure <- function(metadata, column, label) {
  if (!has_usable_values(metadata, column)) {
    return(unavailable_figure(paste0(label, " chart needs ", tolower(label), " data.")))
  }
  count_column <- if (has_usable_values(metadata, "isolate_id")) "isolate_id" else "plasmid_id"
  count_label  <- if (count_column == "isolate_id") "isolates" else "plasmids"
  
  data <- data.frame(unit = metadata[[count_column]], category = metadata[[column]],
                     stringsAsFactors = FALSE)
  cleaned <- stringr::str_squish(coalesce(as.character(data$category), ""))
  cleaned[cleaned == ""] <- "Unknown"
  first_spelling <- cleaned[!duplicated(tolower(cleaned))]
  names(first_spelling) <- tolower(first_spelling)
  data$category <- unname(first_spelling[tolower(cleaned)])
  
  counts <- data %>%
    distinct(unit, category) %>%
    group_by(category) %>%
    summarise(count = n_distinct(unit, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(count), category)
  
  total <- sum(counts$count)
  
  plot_ly(
    counts, labels = ~category, values = ~count, type = "pie",
    hole = 0.5, sort = FALSE, textinfo = "percent", textposition = "inside",
    insidetextorientation = "radial",
    marker = list(colors = rep_len(EXTENDED_COLORS, nrow(counts)),
                  line = list(color = "white", width = 2)),
    hovertemplate = paste0("%{label}<br>%{value} ", count_label, " (%{percent})<extra></extra>")
  ) %>%
    polish_figure(legend_title = label) %>%
    layout(
      margin = list(l = 12, r = 12, t = 28, b = 12),
      annotations = list(list(
        text = paste0(total, "<br>", count_label), x = 0.5, y = 0.5,
        xref = "paper", yref = "paper", showarrow = FALSE,
        font = list(size = 15, color = "#17324D")
      ))
    )
}

# Stacked bar: x = species, stack segment = MLST ST lineage.
make_species_st_figure <- function(metadata) {
  needed <- c("species", "mlst")
  missing <- needed[!vapply(needed, function(c) has_usable_values(metadata, c), logical(1))]
  if (length(missing)) {
    return(unavailable_figure(paste0("Species/ST chart needs data for: ", paste(missing, collapse = ", "), ".")))
  }
  count_column <- if (has_usable_values(metadata, "isolate_id")) "isolate_id" else "plasmid_id"
  count_label  <- if (count_column == "isolate_id") "Number of isolates" else "Number of plasmids"
  
  data <- data.frame(unit = metadata[[count_column]],
                     species = coalesce(as.character(metadata$species), "Unknown"),
                     mlst = coalesce(as.character(metadata$mlst), "Unknown"),
                     stringsAsFactors = FALSE)
  data$species[data$species == ""] <- "Unknown"
  data$mlst[data$mlst == ""] <- "Unknown"
  
  counts <- data %>%
    distinct(unit, species, mlst) %>%
    group_by(species, mlst) %>%
    summarise(count = n_distinct(unit, na.rm = TRUE), .groups = "drop") %>%
    arrange(species, mlst)
  
  species_order <- counts %>%
    group_by(species) %>%
    summarise(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    pull(species)
  
  st_levels <- unique(counts$mlst)
  counts$mlst <- factor(counts$mlst, levels = st_levels)
  color_map <- setNames(rep_len(EXTENDED_COLORS, length(st_levels)), st_levels)
  
  plot_ly(
    counts, x = ~species, y = ~count, color = ~mlst, colors = color_map,
    type = "bar",
    marker = list(line = list(color = "white", width = 0.7)),
    hovertemplate = paste0("Species: %{x}<br>ST: %{fullData.name}<br>",
                           count_label, ": %{y}<extra></extra>")
  ) %>%
    polish_figure(legend_title = "ST lineage") %>%
    layout(
      barmode = "stack",
      xaxis = list(title = "Species", categoryorder = "array", categoryarray = species_order),
      yaxis = list(title = count_label),
      legend = list(orientation = "h", yanchor = "bottom", y = 1.02,
                    xanchor = "left", x = 0)
    )
}

# Plasmid size and GC content distributions on separate axes.
make_size_gc_figure <- function(metadata) {
  needed <- c("Plasmid Size", "Plasmid GC content")
  missing <- needed[!vapply(needed, function(c) has_usable_values(metadata, c), logical(1))]
  if (length(missing)) {
    return(unavailable_figure(paste0("Size/GC plot needs data for: ", paste(missing, collapse = ", "), ".")))
  }
  size <- metadata[["Plasmid Size"]];        size <- size[!is.na(size)]
  gc   <- metadata[["Plasmid GC content"]];  gc   <- gc[!is.na(gc)]
  
  plot_ly() %>%
    add_trace(
      type = "box", x = rep("Plasmid size", length(size)), y = size,
      name = "Plasmid size", yaxis = "y",
      boxpoints = "all", jitter = 0.34, pointpos = 0,
      marker = list(color = "#167D8D", size = 7, opacity = 0.68,
                    line = list(color = "white", width = 0.7)),
      line = list(color = "#167D8D", width = 2),
      fillcolor = "rgba(22,125,141,0.16)",
      hovertemplate = "Plasmid size: %{y:.1f} kb<extra></extra>"
    ) %>%
    add_trace(
      type = "box", x = rep("GC content", length(gc)), y = gc,
      name = "GC content", yaxis = "y2",
      boxpoints = "all", jitter = 0.34, pointpos = 0,
      marker = list(color = "#E28E48", size = 7, opacity = 0.68,
                    line = list(color = "white", width = 0.7)),
      line = list(color = "#E28E48", width = 2),
      fillcolor = "rgba(226,142,72,0.18)",
      hovertemplate = "GC content: %{y:.2f}%<extra></extra>"
    ) %>%
    polish_figure(show_legend = FALSE) %>%
    layout(
      xaxis = list(title = "Metric", categoryorder = "array",
                   categoryarray = c("Plasmid size", "GC content"), showgrid = FALSE),
      yaxis = list(title = "Plasmid size (kb)", rangemode = "tozero", gridcolor = "#E9EEF2"),
      yaxis2 = list(title = "GC content (%)", overlaying = "y", side = "right",
                    showgrid = FALSE, ticksuffix = "%"),
      margin = list(l = 68, r = 72, t = 24, b = 50),
      hovermode = "closest"
    )
}


# ---------------------------------------------------------------------
# 6. CONSENSUS NETWORK MODULE (Pling + mge-cluster)
#    Self-contained Shiny module so its input/output ids cannot clash
#    with the rest of the dashboard.
# ---------------------------------------------------------------------

consensus_normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.(fna|fasta|fa)$", "", x, ignore.case = TRUE)
  # mge-cluster changes accession version dots to underscores
  x <- sub("^([A-Za-z]{1,5}[0-9]+)_([0-9]+)$", "\\1.\\2", x)
  x
}

consensus_read_input <- function(path, filename) {
  is_tab <- grepl("\\.(tabular|txt|tsv)$", filename, ignore.case = TRUE)
  read.csv(path, sep = if (is_tab) "\t" else ",", check.names = FALSE,
           stringsAsFactors = FALSE)
}

# First non-comment, non-blank line of a file, split into header names
consensus_header <- function(path) {
  con <- file(path, "r", encoding = "UTF-8-BOM"); on.exit(close(con))
  repeat {
    l <- readLines(con, n = 1, warn = FALSE)
    if (!length(l)) return(character(0))
    if (nzchar(trimws(l)) && !startsWith(trimws(l), "#")) break
  }
  gsub('^"|"$', "", trimws(strsplit(l, sniff_delimiter(l), fixed = TRUE)[[1]]))
}

# Find the Pling typing, Pling distance and mge-cluster files in a folder by header
consensus_scan_dir <- function(dir) {
  out <- list(typing = NULL, dist = NULL, mge = NULL)
  if (!dir.exists(dir)) return(out)
  files <- list.files(dir, full.names = TRUE)
  files <- files[!dir.exists(files) &
                   tolower(tools::file_ext(files)) %in% c("csv", "tsv", "txt", "tabular")]
  for (f in sort(files, method = "radix")) {
    h <- tryCatch(consensus_header(f), error = function(e) character(0))
    ref <- list(name = basename(f), datapath = f)
    if (all(c("plasmid_1", "plasmid_2", "containment_distance", "DCJI_distance") %in% h)) {
      if (is.null(out$dist)) out$dist <- ref
    } else if (all(c("Sample_Name", "Standard_Cluster") %in% h)) {
      if (is.null(out$mge)) out$mge <- ref
    } else if (all(c("plasmid", "type") %in% h) && !"accession" %in% h) {
      if (is.null(out$typing)) out$typing <- ref
    }
  }
  out
}

consensusUI <- function(id) {
  ns <- NS(id)
  layout_sidebar(
    sidebar = sidebar(
      width = 330,
      h5("Input files"),
      helpText("Files in the data folder are detected automatically; uploading a file replaces it."),
      verbatimTextOutput(ns("sources")),
      fileInput(ns("typing"), "Pling typing (.tabular / .tsv)", width = "100%"),
      fileInput(ns("dist"), "Pling pairwise distances (.txt / .tsv)", width = "100%"),
      fileInput(ns("mge"), "mge-cluster results (.csv)", width = "100%"),
      helpText("Needs all three files. No sequences are needed."),
      hr(),
      h5("Pling edge thresholds"),
      sliderInput(ns("contain"), "Maximum containment distance",
                  min = 0, max = 1, value = 0.15, step = 0.01, width = "100%"),
      sliderInput(ns("dcj"), "Maximum DCJ-Indel distance",
                  min = 0, max = 60, value = 5, step = 1, width = "100%"),
      checkboxInput(ns("samecommunity"), "Require same Pling community", TRUE),
      hr(),
      h5("Network display"),
      selectInput(ns("view"), "Show edges",
                  c("All supported relationships" = "all", "Consensus only" = "both",
                    "Pling only" = "pling", "mge-cluster only" = "mge"), width = "100%"),
      selectInput(ns("nodecolor"), "Colour nodes by",
                  c("mge-cluster" = "mge", "Pling community" = "community",
                    "Pling subcommunity" = "subcommunity"), width = "100%"),
      checkboxInput(ns("isolates"), "Include isolated nodes", TRUE),
      downloadButton(ns("download_edges"), "Download edge table"),
      downloadButton(ns("download_graph"), "Download GraphML"),
      downloadButton(ns("download_nodes"), "Download node table")
    ),
    navset_tab(
      nav_panel(
        "Network",
        div(class = "pt-3",
            verbatimTextOutput(ns("summary")),
            visNetworkOutput(ns("network"), height = "650px"),
            helpText("Green = both; orange = Pling only; purple = mge-cluster only. ",
                     "Hover for details. Edge definitions depend on selected thresholds."))
      ),
      nav_panel(
        "Cluster correspondence",
        div(class = "pt-3",
            plotOutput(ns("heatmap"), height = "520px"),
            h4("Cross-tabulation"),
            DT::DTOutput(ns("crosstab")))
      ),
      nav_panel("Pairwise relationships",
                div(class = "pt-3", DT::DTOutput(ns("edge_table")))),
      nav_panel("Input diagnostics",
                div(class = "pt-3",
                    verbatimTextOutput(ns("diagnostics")),
                    DT::DTOutput(ns("node_table"))))
    )
  )
}

consensusServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Files found in the data folder (re-checked every 5 s for changes)
    folder_files <- reactivePoll(
      5000, session,
      checkFunc = function() {
        f <- if (dir.exists(DATA_DIR)) list.files(DATA_DIR, full.names = TRUE) else character(0)
        paste(f, file.info(f)$mtime, file.info(f)$size, collapse = "|")
      },
      valueFunc = function() consensus_scan_dir(DATA_DIR)
    )
    # An uploaded file takes priority over the one found in the folder
    src <- reactive({
      ff <- folder_files()
      list(typing = input$typing %||% ff$typing,
           dist   = input$dist   %||% ff$dist,
           mge    = input$mge    %||% ff$mge)
    })
    output$sources <- renderText({
      s <- src(); ff <- folder_files()
      line <- function(label, x, up) paste0(label, ": ",
                                            if (is.null(x)) "not found" else paste0(x$name, if (!is.null(up)) " (uploaded)" else " (data folder)"))
      paste(c(line("Pling typing", s$typing, input$typing),
              line("Pling distances", s$dist, input$dist),
              line("mge-cluster", s$mge, input$mge)), collapse = "\n")
    })
    
    data <- reactive({
      s <- src()
      validate(need(!is.null(s$typing) && !is.null(s$dist) && !is.null(s$mge),
                    "Provide the Pling typing, Pling distance and mge-cluster files (upload them or put them in the data folder)."))
      ty <- consensus_read_input(s$typing$datapath, s$typing$name)
      pd <- consensus_read_input(s$dist$datapath,   s$dist$name)
      mg <- consensus_read_input(s$mge$datapath,    s$mge$name)
      validate(
        need(all(c("plasmid", "type") %in% names(ty)),
             "Pling typing needs plasmid and type columns."),
        need(all(c("plasmid_1", "plasmid_2", "containment_distance", "DCJI_distance") %in% names(pd)),
             "Pling distances need plasmid_1, plasmid_2, containment_distance, DCJI_distance."),
        need(all(c("Sample_Name", "Standard_Cluster") %in% names(mg)),
             "mge-cluster needs Sample_Name and Standard_Cluster.")
      )
      ty$id <- consensus_normalize_id(ty$plasmid)
      mg$id <- consensus_normalize_id(mg$Sample_Name)
      pd$a  <- consensus_normalize_id(pd$plasmid_1)
      pd$b  <- consensus_normalize_id(pd$plasmid_2)
      validate(
        need(!anyDuplicated(ty$id), "Duplicate Pling IDs after normalisation."),
        need(!anyDuplicated(mg$id), "Duplicate mge-cluster IDs after normalisation.")
      )
      ty$community    <- sub("_subcommunity_.*$", "", ty$type)
      ty$subcommunity <- ty$type
      mg$mge_cluster  <- as.character(mg$Standard_Cluster)
      mg$mge_cluster[is.na(mg$Standard_Cluster) | mg$mge_cluster %in% c("-1", "-1.0", "nan")] <- NA_character_
      if (!"Membership_Probability" %in% names(mg)) mg$Membership_Probability <- NA_real_
      
      nodes <- merge(ty[, c("id", "plasmid", "community", "subcommunity")],
                     mg[, c("id", "Sample_Name", "mge_cluster", "Membership_Probability")],
                     by = "id", all = TRUE)
      nodes$matched <- !is.na(nodes$plasmid) & !is.na(nodes$Sample_Name)
      ids <- nodes$id[nodes$matched]
      if (length(ids) < 2) {
        return(list(nodes = nodes, pairs = data.frame(), p = pd, typing = ty, mge = mg))
      }
      
      comb  <- t(combn(sort(ids), 2))
      pairs <- data.frame(a = comb[, 1], b = comb[, 2], stringsAsFactors = FALSE)
      pd$k1 <- pmin(pd$a, pd$b)
      pd$k2 <- pmax(pd$a, pd$b)
      pd <- pd[!duplicated(pd[, c("k1", "k2")]), ]
      pairs <- merge(pairs, pd[, c("k1", "k2", "containment_distance", "DCJI_distance")],
                     by.x = c("a", "b"), by.y = c("k1", "k2"), all.x = TRUE, sort = FALSE)
      ia <- match(pairs$a, nodes$id)
      ib <- match(pairs$b, nodes$id)
      pairs$pling_community_a    <- nodes$community[ia];    pairs$pling_community_b    <- nodes$community[ib]
      pairs$pling_subcommunity_a <- nodes$subcommunity[ia]; pairs$pling_subcommunity_b <- nodes$subcommunity[ib]
      pairs$mge_cluster_a        <- nodes$mge_cluster[ia];  pairs$mge_cluster_b        <- nodes$mge_cluster[ib]
      pairs$pairwise_available <- !is.na(pairs$containment_distance) & !is.na(pairs$DCJI_distance)
      list(nodes = nodes, pairs = pairs, p = pd, typing = ty, mge = mg)
    })
    
    edges <- reactive({
      d <- data()
      req(nrow(d$pairs) > 0)
      x <- d$pairs
      x$pling_close <- x$pairwise_available &
        x$containment_distance <= input$contain & x$DCJI_distance <= input$dcj
      if (input$samecommunity) {
        x$pling_close <- x$pling_close & x$pling_community_a == x$pling_community_b
      }
      x$mge_close <- !is.na(x$mge_cluster_a) & !is.na(x$mge_cluster_b) &
        x$mge_cluster_a == x$mge_cluster_b
      x$relationship <- ifelse(x$pling_close & x$mge_close, "Both",
                               ifelse(x$pling_close, "Pling only",
                                      ifelse(x$mge_close, "mge-cluster only", "Neither")))
      x
    })
    
    selected <- reactive({
      x <- edges()
      switch(input$view,
             all  = x[x$relationship != "Neither", ],
             both = x[x$relationship == "Both", ],
             pling = x[x$relationship == "Pling only", ],
             mge  = x[x$relationship == "mge-cluster only", ])
    })
    
    output$summary <- renderText({
      x <- edges(); d <- data()
      paste0("Matched plasmids: ", sum(d$nodes$matched), " / Pling ", nrow(d$typing),
             " / mge-cluster ", nrow(d$mge),
             "\nPairs: ", nrow(x), " | Pling distances available: ", sum(x$pairwise_available),
             "\nBoth: ", sum(x$relationship == "Both"),
             " | Pling only: ", sum(x$relationship == "Pling only"),
             " | mge-cluster only: ", sum(x$relationship == "mge-cluster only"),
             "\nNote: mge-cluster links mean shared cluster membership, not a measured pairwise similarity.")
    })
    
    output$network <- renderVisNetwork({
      d <- data(); e <- selected()
      n <- d$nodes[d$nodes$matched, ]
      if (!input$isolates) n <- n[n$id %in% c(e$a, e$b), ]
      validate(need(nrow(n) > 0, "No nodes to display at these settings."))
      group <- switch(input$nodecolor, mge = n$mge_cluster,
                      community = n$community, subcommunity = n$subcommunity)
      group[is.na(group)] <- "Unassigned"
      vn <- data.frame(
        id = n$id, label = n$id, group = group,
        title = paste0("<b>", n$id, "</b><br>Pling: ", n$subcommunity,
                       "<br>mge-cluster: ", n$mge_cluster,
                       "<br>Membership probability: ", round(n$Membership_Probability, 3)),
        stringsAsFactors = FALSE
      )
      ve <- if (nrow(e) > 0) {
        data.frame(
          from = e$a, to = e$b,
          color = ifelse(e$relationship == "Both", "#168A76",
                         ifelse(e$relationship == "Pling only", "#D97706", "#7C6EE6")),
          title = paste0(e$relationship, "<br>Containment: ", e$containment_distance,
                         "<br>DCJ-Indel: ", e$DCJI_distance,
                         "<br>Pairwise distance available: ", e$pairwise_available),
          stringsAsFactors = FALSE
        )
      } else NULL
      visNetwork(vn, ve, width = "100%", height = "650px") %>%
        visEdges(smooth = FALSE) %>% visNodes(size = 18) %>%
        visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visPhysics(stabilization = TRUE) %>% visInteraction(navigationButtons = TRUE)
    })
    
    cross <- reactive({
      d <- data()
      n <- d$nodes[d$nodes$matched, ]
      as.data.frame(table(Pling = n$subcommunity, MGE = n$mge_cluster, useNA = "ifany"))
    })
    
    output$heatmap <- renderPlot({
      ggplot2::ggplot(cross(), ggplot2::aes(x = MGE, y = Pling, fill = Freq)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(Freq > 0, Freq, "")), size = 3) +
        ggplot2::scale_fill_gradient(low = "white", high = "#246B9A") +
        ggplot2::labs(title = "Pling subcommunity vs mge-cluster",
                      x = "mge-cluster", y = "Pling subcommunity") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    })
    
    output$crosstab   <- DT::renderDT(DT::datatable(cross(), options = list(pageLength = 20)))
    output$edge_table <- DT::renderDT(DT::datatable(edges(), filter = "top",
                                                    options = list(pageLength = 15, scrollX = TRUE)))
    output$node_table <- DT::renderDT(DT::datatable(data()$nodes, filter = "top",
                                                    options = list(pageLength = 15, scrollX = TRUE)))
    
    output$diagnostics <- renderText({
      d <- data(); p <- d$pairs
      paste0("Unmatched Pling IDs: ", paste(setdiff(d$typing$id, d$mge$id), collapse = ", "),
             "\nUnmatched mge-cluster IDs: ", paste(setdiff(d$mge$id, d$typing$id), collapse = ", "),
             "\nMatched pairs lacking Pling distance: ",
             if (nrow(p)) sum(!p$pairwise_available) else "N/A",
             "\nImportant: missing Pling pairs are NOT interpreted as distant.")
    })
    
    output$download_edges <- downloadHandler(
      filename = function() "consensus_edges.csv",
      content = function(file) write.csv(edges(), file, row.names = FALSE))
    output$download_nodes <- downloadHandler(
      filename = function() "consensus_nodes.csv",
      content = function(file) write.csv(data()$nodes, file, row.names = FALSE))
    output$download_graph <- downloadHandler(
      filename = function() "consensus_network.graphml",
      content = function(file) {
        e <- selected(); n <- data()$nodes; n <- n[n$matched, ]
        if (!input$isolates) n <- n[n$id %in% c(e$a, e$b), ]
        graph <- igraph::graph_from_data_frame(
          data.frame(from = e$a, to = e$b, relationship = e$relationship, stringsAsFactors = FALSE),
          directed = FALSE,
          vertices = data.frame(name = n$id, pling_subcommunity = n$subcommunity,
                                mge_cluster = n$mge_cluster, stringsAsFactors = FALSE))
        igraph::write_graph(graph, file, format = "graphml")
      })
  })
}


# ---------------------------------------------------------------------
# 7. GENETIC CONTEXT MODULE (Bakta annotation + Clinker comparison)
#    Self-contained Shiny module. Optional extras: Bioconductor packages
#    rtracklayer / genbankr are only needed for .gff3 / GenBank input.
# ---------------------------------------------------------------------

# Reads a Bakta results table (.csv / .tsv / .txt). Handles:
#  - the "# Annotated with Bakta ..." comment block above the header
#  - the "#Sequence Id" header line (leading # is stripped)
#  - tab-, comma- or semicolon-delimited files (e.g. TSV re-saved as CSV
#    from Excel, including Danish-locale semicolon CSVs)
read_bakta_table <- function(path) {
  
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^\ufeff", "", lines)
  lines <- sub("\r$", "", lines)
  
  hdr <- which(
    grepl("^#?\\s*\"?Sequence[ _]?Id", lines, ignore.case = TRUE)
  )[1]
  
  if (is.na(hdr)) {
    stop("Could not find the 'Sequence Id' header row in the Bakta table.")
  }
  
  header_line <- lines[hdr]
  
  delim_counts <- c(
    "\t" = stringr::str_count(header_line, "\t"),
    "," = stringr::str_count(header_line, ","),
    ";" = stringr::str_count(header_line, ";")
  )
  delim <- names(which.max(delim_counts))
  
  body <- lines[hdr:length(lines)]
  body[1] <- sub("^#\\s*", "", body[1])
  
  # drop remaining comment lines and empty rows (e.g. ",,,,,,,,")
  body <- body[
    !grepl("^#", body) &
      nzchar(trimws(gsub("[,;\t\"]", "", body)))
  ]
  
  df <- utils::read.table(
    text = paste(body, collapse = "\n"),
    sep = delim,
    header = TRUE,
    quote = if (delim == "\t") "" else "\"",
    comment.char = "",
    colClasses = "character",
    na.strings = character(0),
    check.names = FALSE,
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  # "Sequence Id" -> sequence_id, "Locus Tag" -> locus_tag, etc.
  names(df) <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(names(df))))
  
  required <- c("sequence_id", "type", "start", "stop", "strand")
  missing_cols <- setdiff(required, names(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Bakta table is missing column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  col_or_blank <- function(name) {
    if (name %in% names(df)) df[[name]] else rep("", nrow(df))
  }
  
  data.frame(
    seqnames  = df$sequence_id,
    type      = df$type,
    start     = suppressWarnings(as.integer(df$start)),
    end       = suppressWarnings(as.integer(df$stop)),
    strand    = df$strand,
    locus_tag = col_or_blank("locus_tag"),
    gene      = col_or_blank("gene"),
    product   = col_or_blank("product"),
    dbxrefs   = col_or_blank("dbxrefs"),
    stringsAsFactors = FALSE
  )
}


read_bakta_cds <- function(path, filename) {
  
  ext <- tolower(tools::file_ext(filename))
  
  if (ext %in% c("csv", "tsv", "txt")) {
    
    ann <- read_bakta_table(path)
    
    # Keep coding sequences only (Bakta types: cds, tRNA, rRNA, sorf, ...).
    # Add e.g. "sorf" here if you also want small ORFs in the plot.
    ann <- ann %>%
      filter(tolower(type) %in% c("cds"))
    
  } else if (ext == "gff3") {
    
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop("Reading .gff3 files needs the Bioconductor package 'rtracklayer'.")
    }
    gr <- rtracklayer::import(path)
    ann <- as.data.frame(gr)
    
    if (!"type" %in% names(ann)) {
      stop("The GFF3 file does not contain a feature type column.")
    }
    
    ann <- ann %>%
      filter(toupper(as.character(type)) == "CDS")
    
  } else if (ext %in% c("gbff", "gb", "gbk")) {
    
    if (!requireNamespace("genbankr", quietly = TRUE)) {
      stop("Reading GenBank files needs the Bioconductor package 'genbankr'.")
    }
    gb <- genbankr::readGenBank(path)
    gr <- genbankr::cds(gb)
    ann <- as.data.frame(gr)
    ann$type <- "CDS"
    
  } else {
    
    stop("Please upload a Bakta .csv/.tsv table (or a .gff3, .gbff, .gb or .gbk file).")
    
  }
  
  if (nrow(ann) == 0) {
    stop("No CDS features were found in this annotation file.")
  }
  
  get_field <- function(data, candidates) {
    
    hit <- candidates[candidates %in% names(data)]
    
    if (length(hit) == 0) {
      return(rep("", nrow(data)))
    }
    
    x <- data[[hit[1]]]
    
    if (is.list(x)) {
      
      x <- vapply(x, function(z) {
        
        if (length(z) == 0 || all(is.na(z))) {
          ""
        } else {
          paste(as.character(z), collapse = "; ")
        }
        
      }, character(1))
      
    }
    
    x <- as.character(x)
    x[is.na(x)] <- ""
    
    x
  }
  
  seq_col <- if ("seqnames" %in% names(ann)) {
    "seqnames"
  } else {
    NA_character_
  }
  
  if (is.na(seq_col)) {
    stop("Could not identify sequence/contig names in the annotation.")
  }
  
  out <- tibble(
    
    feature_id = seq_len(nrow(ann)),
    
    seqid = as.character(ann[[seq_col]]),
    
    start = as.integer(ann$start),
    
    end = as.integer(ann$end),
    
    strand = as.character(ann$strand),
    
    gene = get_field(
      ann,
      c("gene", "gene_name", "Name")
    ),
    
    product = get_field(
      ann,
      c("product", "product_name", "description")
    ),
    
    locus_tag = get_field(
      ann,
      c("locus_tag", "locusTag")
    ),
    
    feature_name = get_field(
      ann,
      c("Name", "ID", "protein_id")
    ),
    
    type = "CDS"
    
  ) %>%
    
    mutate(
      
      strand = if_else(
        strand %in% c("+", "-"),
        strand,
        "?"
      ),
      
      start = pmin(start, end),
      
      end = pmax(start, end),
      
      display_name = case_when(
        
        gene != "" ~ gene,
        
        locus_tag != "" ~ locus_tag,
        
        product != "" ~ product,
        
        TRUE ~ paste0("CDS_", feature_id)
        
      ),
      
      search_label = paste(
        
        display_name,
        
        if_else(
          product != "",
          paste0("— ", product),
          ""
        ),
        
        paste0(
          "[",
          seqid,
          ":",
          start,
          "-",
          end,
          "; ",
          strand,
          "]"
        )
        
      )
      
    ) %>%
    
    filter(
      !is.na(start),
      !is.na(end),
      !is.na(seqid)
    )
  
  if (nrow(out) == 0) {
    stop("No usable CDS coordinates were found.")
  }
  
  out
}


# ---------------------------------------------------------
# 2.2 CLINKER OUTPUT HELPERS
# ---------------------------------------------------------
# Accepts the interactive HTML written by `clinker ... -p plot.html`.
# The HTML embeds all clusters, genes, links and gene groups as one
# JSON object (`const data = {...}`). To show only N genes upstream /
# downstream of an anchor gene, a small piece of JavaScript is
# injected right after that object; it trims clusters, links and
# groups before clinker draws the figure. The original HTML is never
# modified on disk.

# JavaScript injected into the Clinker HTML. __CLINKER_PARAMS__ is
# replaced by "null" (parameters then come from the URL query string,
# so sliders can change the view without rewriting the file) or by a
# JSON object (used for the downloadable, trimmed HTML).
clinker_filter_js <- r"---(
(function () {
  var P = __CLINKER_PARAMS__;
  if (!P) {
    var q = new URLSearchParams(window.location.search);
    P = {
      group: q.get('group') || '',
      up: parseInt(q.get('up'), 10),
      down: parseInt(q.get('down'), 10),
      drop: q.get('drop') !== '0'
    };
  }
  if (!P.group) return;
  if (isNaN(P.up)) P.up = 5;
  if (isNaN(P.down)) P.down = 5;
  var grp = (data.groups || []).find(function (g) { return g.uid === P.group; });
  if (!grp) return;
  var members = new Set(grp.genes);
  var keep = new Set();
  var clusters = [];
  data.clusters.forEach(function (cl) {
    var loci = [];
    cl.loci.forEach(function (loc) {
      var genes = loc.genes.slice().sort(function (a, b) { return a.start - b.start; });
      var i = genes.findIndex(function (g) { return members.has(g.uid); });
      if (i < 0) {
        if (!P.drop) {
          loc.genes.forEach(function (g) { keep.add(g.uid); });
          loci.push(loc);
        }
        return;
      }
      var rev = genes[i].strand === -1;
      var lo = Math.max(0, i - (rev ? P.down : P.up));
      var hi = Math.min(genes.length - 1, i + (rev ? P.up : P.down));
      var kept = genes.slice(lo, hi + 1);
      kept.forEach(function (g) { keep.add(g.uid); });
      loc.genes = kept;
      loc.start = Math.min.apply(null, kept.map(function (g) { return g.start; }));
      loc.end = Math.max.apply(null, kept.map(function (g) { return g.end; }));
      loci.push(loc);
    });
    if (loci.length) { cl.loci = loci; clusters.push(cl); }
  });
  data.clusters = clusters;
  data.links = (data.links || []).filter(function (l) {
    return keep.has(l.query.uid) && keep.has(l.target.uid);
  });
  data.groups = (data.groups || []).map(function (g) {
    g.genes = g.genes.filter(function (u) { return keep.has(u); });
    return g;
  }).filter(function (g) { return g.genes.length > 0; });
})();
)---"

# Read a Clinker HTML plot: split it around the embedded data (so the
# filter can be injected) and summarise the gene groups for the UI.
read_clinker_html <- function(path) {
  
  html <- paste(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  
  data_tag <- "const data="
  start <- regexpr(data_tag, html, fixed = TRUE)[1]
  
  if (start < 0) {
    stop(
      "No embedded clinker data found. Please upload the HTML file ",
      "written by 'clinker ... -p plot.html'."
    )
  }
  
  json_start <- start + nchar(data_tag)
  tail_html <- substr(html, json_start, nchar(html))
  end_rel <- regexpr(";function serialise", tail_html, fixed = TRUE)[1]
  
  if (end_rel < 0) {
    stop("Unrecognised clinker HTML layout (end of embedded data not found).")
  }
  
  json_txt <- substr(tail_html, 1, end_rel - 1)
  split_at <- json_start + end_rel - 1      # position of the ';'
  
  dat <- jsonlite::fromJSON(json_txt, simplifyVector = FALSE)
  
  one_string <- function(x) {
    if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]])
  }
  
  gene_tbl <- dplyr::bind_rows(lapply(dat$clusters, function(cl) {
    dplyr::bind_rows(lapply(cl$loci, function(loc) {
      dplyr::bind_rows(lapply(loc$genes, function(g) {
        tibble::tibble(
          uid     = g[["uid"]],
          cluster = cl[["name"]],
          gene    = one_string(g[["names"]][["gene"]]),
          product = one_string(g[["names"]][["product"]])
        )
      }))
    }))
  }))
  
  most_common <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x) == 0) NA_character_ else names(sort(table(x), decreasing = TRUE))[1]
  }
  
  groups <- dplyr::bind_rows(lapply(dat$groups, function(gr) {
    
    members <- unlist(gr$genes)
    sub_tbl <- gene_tbl[gene_tbl$uid %in% members, ]
    
    nm <- most_common(sub_tbl$gene)
    pr <- most_common(sub_tbl$product)
    desc <- paste(stats::na.omit(c(nm, pr)), collapse = " - ")
    
    n_clusters <- length(unique(sub_tbl$cluster))
    
    tibble::tibble(
      uid = gr[["uid"]],
      has_name = !is.na(nm),
      n_clusters = n_clusters,
      choice_label = paste0(
        gr[["label"]],
        if (nzchar(desc)) paste0(": ", desc) else "",
        " (", nrow(sub_tbl), " genes, ", n_clusters, " clusters)"
      )
    )
    
  }))
  
  if (nrow(groups) > 0) {
    groups <- groups %>%
      dplyr::arrange(dplyr::desc(has_name), dplyr::desc(n_clusters), choice_label)
  }
  
  max_genes <- max(
    1L,
    unlist(lapply(dat$clusters, function(cl) {
      vapply(cl$loci, function(loc) length(loc$genes), integer(1))
    }))
  )
  
  list(
    pre = substr(html, 1, split_at),
    post = substr(html, split_at + 1, nchar(html)),
    groups = groups,
    n_clusters = length(dat$clusters),
    max_genes = max_genes
  )
}


# Build the HTML text for a given set of filter parameters
# (params = NULL means "read them from the URL query string").
build_clinker_html <- function(parsed, params = NULL) {
  
  params_js <- if (is.null(params)) {
    "null"
  } else {
    as.character(jsonlite::toJSON(params, auto_unbox = TRUE))
  }
  
  js <- sub("__CLINKER_PARAMS__", params_js, clinker_filter_js, fixed = TRUE)
  
  paste0(parsed$pre, js, parsed$post)
}


# Does this file look like a Bakta results table? ("Sequence Id" header row)
is_bakta_table <- function(path) {
  if (!tolower(tools::file_ext(path)) %in% c("csv", "tsv", "txt")) return(FALSE)
  l <- tryCatch(readLines(path, n = 300, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  any(grepl('^#?\\s*"?Sequence[ _]?Id', sub("^\ufeff", "", l), ignore.case = TRUE))
}

# Does this file look like a Clinker plot (HTML with embedded data)?
is_clinker_html <- function(path) {
  if (!tolower(tools::file_ext(path)) %in% c("html", "htm")) return(FALSE)
  txt <- tryCatch(readChar(path, file.size(path), useBytes = TRUE), error = function(e) "")
  isTRUE(grepl("const data=", txt, fixed = TRUE))
}

# First Bakta table and first Clinker plot found in a folder
context_scan_dir <- function(dir) {
  out <- list(bakta = NULL, clinker = NULL)
  if (!dir.exists(dir)) return(out)
  files <- list.files(dir, full.names = TRUE)
  files <- sort(files[!dir.exists(files)], method = "radix")
  for (f in files) {
    ref <- list(name = basename(f), datapath = f)
    if (is.null(out$bakta) && is_bakta_table(f)) out$bakta <- ref
    else if (is.null(out$clinker) && is_clinker_html(f)) out$clinker <- ref
  }
  out
}

contextUI <- function(id) {
  ns <- NS(id)
  navset_tab(
    nav_panel(
      "Single plasmid",
      layout_sidebar(
        sidebar = sidebar(
          width = 330,
          verbatimTextOutput(ns("bakta_source")),
          fileInput(ns("bakta_file"), "Upload Bakta annotation table (.csv / .tsv)",
                    accept = c(".csv", ".tsv", ".txt", ".gff3", ".gbff", ".gb", ".gbk"),
                    width = "100%"),
          helpText("Upload the Bakta results table (.csv or .tsv), or put it in the data folder. ",
                   "The visualisation uses CDS features on the same contig as the selected gene."),
          selectInput(ns("context_gene"), "Select annotated gene/CDS:", choices = NULL, width = "100%"),
          sliderInput(ns("upstream_cds"), "Coding sequences upstream:",
                      min = 0, max = 20, value = 5, step = 1, width = "100%"),
          sliderInput(ns("downstream_cds"), "Coding sequences downstream:",
                      min = 0, max = 20, value = 5, step = 1, width = "100%")
        ),
        card(card_header("Genetic neighbourhood"), plotOutput(ns("context_plot"), height = "430px")),
        card(card_header("Features in selected region"), DT::DTOutput(ns("context_table"))),
        card(card_header("Selected feature details"), verbatimTextOutput(ns("context_details")))
      )
    ),
    nav_panel(
      "Multiple plasmids",
      layout_sidebar(
        sidebar = sidebar(
          width = 330,
          verbatimTextOutput(ns("clinker_source")),
          fileInput(ns("clinker_html"), "Upload Clinker output (plot.html)",
                    accept = c(".html", ".htm"), width = "100%"),
          selectInput(ns("clinker_focus"), "Anchor gene (homology group):",
                      choices = c("Upload a file first" = ""), width = "100%"),
          sliderInput(ns("clinker_up"), "Genes upstream:", min = 0, max = 20, value = 5, step = 1, width = "100%"),
          sliderInput(ns("clinker_down"), "Genes downstream:", min = 0, max = 20, value = 5, step = 1, width = "100%"),
          checkboxInput(ns("clinker_drop"), "Hide clusters that lack the anchor gene", value = TRUE),
          helpText("Choose an anchor gene group (e.g. optrA) to show only that many genes upstream and ",
                   "downstream of it in every cluster. Upstream/downstream follow the anchor gene's strand. ",
                   "With no anchor, the full Clinker plot is shown. Links to removed genes are dropped. ",
                   "You can still drag, recolour and export in the plot itself."),
          sliderInput(ns("clinker_height"), "Plot height (px):", min = 400, max = 1600, value = 800, step = 100,
                      width = "100%"),
          downloadButton(ns("download_clinker_trimmed"), "Download trimmed HTML")
        ),
        card(card_header("Clinker gene cluster comparison"),
             textOutput(ns("clinker_info")),
             uiOutput(ns("clinker_plot_ui")))
      )
    )
  )
}

contextServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Files found in the data folder (re-checked every 5 s for changes)
    folder_files <- reactivePoll(
      5000, session,
      checkFunc = function() {
        f <- if (dir.exists(DATA_DIR)) list.files(DATA_DIR, full.names = TRUE) else character(0)
        paste(f, file.info(f)$mtime, file.info(f)$size, collapse = "|")
      },
      valueFunc = function() context_scan_dir(DATA_DIR)
    )
    # An uploaded file takes priority over the one found in the folder
    bakta_src   <- reactive(input$bakta_file   %||% folder_files()$bakta)
    clinker_src <- reactive(input$clinker_html %||% folder_files()$clinker)
    
    output$bakta_source <- renderText({
      s <- bakta_src()
      paste0("Bakta table: ", if (is.null(s)) "not found" else
        paste0(s$name, if (!is.null(input$bakta_file)) " (uploaded)" else " (data folder)"))
    })
    output$clinker_source <- renderText({
      s <- clinker_src()
      paste0("Clinker plot: ", if (is.null(s)) "not found" else
        paste0(s$name, if (!is.null(input$clinker_html)) " (uploaded)" else " (data folder)"))
    })
    
    # Read the uploaded Bakta annotation and standardise CDS fields.
    
    bakta_features <- reactive({
      
      req(bakta_src())
      
      tryCatch(
        
        read_bakta_cds(
          bakta_src()$datapath,
          bakta_src()$name
        ),
        
        error = function(e) {
          
          # Show the problem to the user (and in the R console)
          # instead of failing silently.
          message("Bakta import error: ", conditionMessage(e))
          
          showNotification(
            paste("Could not read annotation:", conditionMessage(e)),
            type = "error",
            duration = NULL
          )
          
          NULL
          
        }
        
      )
      
    })
    
    
    # Populate the feature selector after a file is uploaded.
    
    observeEvent(
      
      bakta_features(),
      
      {
        
        req(bakta_features())
        
        x <- bakta_features()
        
        choices <- stats::setNames(
          as.character(x$feature_id),
          x$search_label
        )
        
        updateSelectInput(
          
          session,
          
          "context_gene",
          
          choices = choices,
          
          selected = as.character(x$feature_id[1])
          
        )
        
      },
      
      ignoreInit = FALSE
      
    )
    
    
    # Select the requested CDS plus its nearest neighbouring CDSs.
    
    context_region <- reactive({
      
      req(bakta_features(), input$context_gene)
      
      x <- bakta_features()
      
      target <- x %>%
        filter(feature_id == as.integer(input$context_gene))
      
      validate(
        need(
          nrow(target) == 1,
          "Please select a valid CDS feature."
        )
      )
      
      same_contig <- x %>%
        
        filter(
          seqid == target$seqid,
          feature_id != target$feature_id
        )
      
      
      # Upstream/downstream are defined relative to
      # the strand of the target gene.
      
      if (target$strand == "-") {
        
        upstream <- same_contig %>%
          
          filter(
            start > target$end
          ) %>%
          
          arrange(start) %>%
          
          slice_head(
            n = input$upstream_cds
          )
        
        
        downstream <- same_contig %>%
          
          filter(
            end < target$start
          ) %>%
          
          arrange(desc(end)) %>%
          
          slice_head(
            n = input$downstream_cds
          )
        
      } else {
        
        upstream <- same_contig %>%
          
          filter(
            end < target$start
          ) %>%
          
          arrange(desc(end)) %>%
          
          slice_head(
            n = input$upstream_cds
          )
        
        
        downstream <- same_contig %>%
          
          filter(
            start > target$end
          ) %>%
          
          arrange(start) %>%
          
          slice_head(
            n = input$downstream_cds
          )
        
      }
      
      
      region <- bind_rows(
        upstream,
        target,
        downstream
      ) %>%
        
        distinct(
          feature_id,
          .keep_all = TRUE
        ) %>%
        
        mutate(
          
          is_target = feature_id == target$feature_id,
          
          arrow_start = if_else(
            strand == "-",
            as.numeric(end),
            as.numeric(start)
          ),
          
          arrow_end = if_else(
            strand == "-",
            as.numeric(start),
            as.numeric(end)
          ),
          
          label = case_when(
            
            is_target ~ paste0(
              display_name,
              " (selected)"
            ),
            
            gene != "" ~ gene,
            
            locus_tag != "" ~ locus_tag,
            
            TRUE ~ display_name
            
          ),
          
          hover_text = paste0(
            
            "<b>", display_name, "</b>",
            
            if_else(
              gene != "",
              paste0("<br>Gene: ", gene),
              ""
            ),
            
            if_else(
              product != "",
              paste0("<br>Product: ", product),
              ""
            ),
            
            if_else(
              locus_tag != "",
              paste0("<br>Locus tag: ", locus_tag),
              ""
            ),
            
            "<br>Contig: ", seqid,
            
            "<br>Coordinates: ", start, "-", end,
            
            "<br>Strand: ", strand
            
          )
          
        ) %>%
        
        arrange(start)
      
      list(
        features = region,
        target = target
      )
      
    })
    
    
    # -------------------------------------------------------
    # 8.1 GENETIC CONTEXT PLOT USING GGENES
    # -------------------------------------------------------
    
    output$context_plot <- renderPlot({
      
      cr <- context_region()
      
      # Prepare the selected CDS neighbourhood
      # in the format expected by gggenes.
      
      plot_data <- cr$features %>%
        
        mutate(
          
          molecule = seqid,
          
          orientation = case_when(
            
            strand == "+" ~ TRUE,
            
            strand == "-" ~ FALSE,
            
            TRUE ~ TRUE
            
          ),
          
          gene_label = case_when(
            
            gene != "" ~ gene,
            
            locus_tag != "" ~ locus_tag,
            
            product != "" ~ product,
            
            TRUE ~ display_name
            
          ),
          
          # Highlight the selected gene;
          # show neighbouring CDSs in grey.
          
          colour = if_else(
            is_target,
            "red",
            "others"
          )
          
        ) %>%
        
        arrange(start)
      
      
      validate(
        need(
          nrow(plot_data) > 0,
          "No CDS features are available in this region."
        )
      )
      
      
      # Use the actual coordinate range of the selected
      # neighbourhood rather than a fixed x-axis limit.
      
      x_range <- range(
        c(plot_data$start, plot_data$end),
        na.rm = TRUE
      )
      
      x_padding <- max(
        100,
        diff(x_range) * 0.04
      )
      
      
      ggplot(
        
        plot_data,
        
        aes(
          
          xmin = start,
          
          xmax = end,
          
          y = molecule,
          
          fill = colour,
          
          forward = orientation
          
        )
        
      ) +
        
        geom_gene_arrow(
          
          arrowhead_width = grid::unit(6, "mm"),
          
          arrowhead_height = grid::unit(8, "mm"),
          
          arrow_body_height = grid::unit(6, "mm")
          
        ) +
        
        geom_gene_label(
          
          aes(
            label = gene_label
          ),
          
          align = "centre",
          
          padding.x = grid::unit(2, "mm"),
          
          padding.y = grid::unit(0.1, "lines"),
          
          min.size = 3,
          
          height = grid::unit(3, "mm")
          
        ) +
        
        facet_wrap(
          ~ molecule,
          scales = "free",
          ncol = 1
        ) +
        
        scale_fill_manual(
          
          values = c(
            "red" = "#d62828",
            "others" = "#bdbdbd"
          ),
          
          breaks = c(
            "red",
            "others"
          ),
          
          labels = c(
            "Selected gene",
            "Neighbouring CDSs"
          ),
          
          name = NULL
          
        ) +
        
        scale_x_continuous(
          
          limits = c(
            max(0, x_range[1] - x_padding),
            x_range[2] + x_padding
          ),
          
          name = "Genomic coordinate (bp)",
          
          expand = expansion(mult = 0)
          
        ) +
        
        labs(
          
          title = paste0(
            "Genetic context of ",
            cr$target$display_name
          ),
          
          subtitle = paste0(
            
            nrow(plot_data) - 1,
            " neighbouring CDS features; ",
            
            "selected feature highlighted in red"
            
          ),
          
          y = NULL
          
        ) +
        
        theme_genes() +
        
        theme(
          
          legend.position = "top",
          
          strip.text = element_text(
            face = "bold"
          ),
          
          axis.title.y = element_blank(),
          
          axis.text.y = element_blank(),
          
          axis.ticks.y = element_blank(),
          
          plot.title = element_text(
            face = "bold"
          )
          
        )
      
    })
    
    
    # -------------------------------------------------------
    # 8.2 GENETIC CONTEXT FEATURE TABLE
    # -------------------------------------------------------
    
    output$context_table <- DT::renderDT({
      
      cr <- context_region()
      
      cr$features %>%
        
        arrange(start) %>%
        
        select(
          
          seqid,
          start,
          end,
          strand,
          gene,
          product,
          locus_tag,
          is_target
          
        ) %>%
        
        DT::datatable(
          
          options = list(
            
            pageLength = 15,
            
            scrollX = TRUE,
            
            dom = "Bfrtip",
            
            buttons = c(
              "copy",
              "csv",
              "excel"
            )
            
          ),
          
          extensions = "Buttons",
          
          rownames = FALSE
          
        ) %>%
        
        DT::formatStyle(
          
          "is_target",
          
          target = "row",
          
          backgroundColor = DT::styleEqual(
            c(TRUE, FALSE),
            c("#FCE4E4", "white")
          )
          
        )
      
    })
    
    
    # -------------------------------------------------------
    # 8.3 SELECTED FEATURE DETAILS
    # -------------------------------------------------------
    
    output$context_details <- renderPrint({
      
      cr <- context_region()
      
      target <- cr$target
      
      cat("SELECTED CDS FEATURE\n")
      cat("====================\n")
      
      cat(
        "Gene/name: ",
        target$display_name,
        "\n",
        sep = ""
      )
      
      cat(
        "Product: ",
        ifelse(
          target$product == "",
          "Not annotated",
          target$product
        ),
        "\n",
        sep = ""
      )
      
      cat(
        "Locus tag: ",
        ifelse(
          target$locus_tag == "",
          "Not annotated",
          target$locus_tag
        ),
        "\n",
        sep = ""
      )
      
      cat(
        "Contig: ",
        target$seqid,
        "\n",
        sep = ""
      )
      
      cat(
        "Coordinates: ",
        target$start,
        "-",
        target$end,
        "\n",
        sep = ""
      )
      
      cat(
        "Strand: ",
        target$strand,
        "\n",
        sep = ""
      )
      
      cat(
        "CDS features displayed: ",
        nrow(cr$features),
        "\n",
        sep = ""
      )
      
      cat(
        "\nNote: upstream/downstream are defined relative to the selected gene's strand."
      )
      
    })
    
    
    # -------------------------------------------------------
    # 8.4 CLINKER OUTPUT VIEWER (upstream / downstream trimming)
    # -------------------------------------------------------
    
    # Parse the uploaded Clinker HTML and write a copy with the
    # trimming script injected (served to an iframe below).
    
    clinker_loaded <- reactive({
      
      req(clinker_src())
      
      tryCatch({
        
        parsed <- read_clinker_html(clinker_src()$datapath)
        
        out_dir <- tempfile("clinkerplot_")
        dir.create(out_dir)
        
        writeLines(
          build_clinker_html(parsed),
          file.path(out_dir, "clinker_plot.html"),
          useBytes = TRUE
        )
        
        parsed$prefix <- basename(out_dir)
        shiny::addResourcePath(parsed$prefix, out_dir)
        
        parsed
        
      },
      
      error = function(e) {
        
        message("Clinker import error: ", conditionMessage(e))
        
        showNotification(
          paste("Could not read Clinker file:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
        
        NULL
        
      })
      
    })
    
    
    # Fill the anchor selector and slider limits after upload.
    
    observeEvent(
      
      clinker_loaded(),
      
      {
        
        res <- clinker_loaded()
        req(res)
        
        choices <- c(
          "No anchor - show full clusters" = "",
          stats::setNames(res$groups$uid, res$groups$choice_label)
        )
        
        updateSelectInput(
          session,
          "clinker_focus",
          choices = choices,
          selected = ""
        )
        
        updateSliderInput(session, "clinker_up", max = res$max_genes)
        updateSliderInput(session, "clinker_down", max = res$max_genes)
        
      }
      
    )
    
    
    # Debounced so dragging a slider does not reload the plot each step.
    
    clinker_view <- shiny::debounce(
      
      reactive({
        
        focus <- input$clinker_focus
        if (is.null(focus)) focus <- ""
        
        list(
          group = focus,
          up = input$clinker_up,
          down = input$clinker_down,
          drop = isTRUE(input$clinker_drop)
        )
        
      }),
      
      600
      
    )
    
    
    output$clinker_info <- renderText({
      
      res <- clinker_loaded()
      req(res)
      
      v <- clinker_view()
      
      msg <- paste0(
        res$n_clusters, " clusters and ",
        nrow(res$groups), " gene groups loaded."
      )
      
      if (nzchar(v$group)) {
        
        row <- res$groups[res$groups$uid == v$group, ]
        
        if (nrow(row) == 1) {
          msg <- paste0(
            msg, " The selected anchor gene is present in ",
            row$n_clusters, " of ", res$n_clusters, " clusters."
          )
        }
        
      }
      
      msg
      
    })
    
    
    output$clinker_plot_ui <- renderUI({
      
      validate(
        need(
          !is.null(clinker_src()),
          "Upload the plot.html written by 'clinker ... -p plot.html' (or put it in the data folder)."
        )
      )
      
      res <- clinker_loaded()
      req(res)
      
      v <- clinker_view()
      
      url <- paste0(
        res$prefix, "/clinker_plot.html",
        "?group=", utils::URLencode(v$group, reserved = TRUE),
        "&up=", v$up,
        "&down=", v$down,
        "&drop=", as.integer(v$drop)
      )
      
      tagList(
        
        tags$p(
          tags$a(
            href = url,
            target = "_blank",
            "Open plot in a new tab"
          )
        ),
        
        tags$iframe(
          src = url,
          style = paste0(
            "width:100%; border:0; height:",
            input$clinker_height,
            "px;"
          )
        )
        
      )
      
    })
    
    
    output$download_clinker_trimmed <- downloadHandler(
      
      filename = function() "clinker_plot_trimmed.html",
      
      content = function(file) {
        
        res <- clinker_loaded()
        req(res)
        
        focus <- input$clinker_focus
        if (is.null(focus)) focus <- ""
        
        html <- build_clinker_html(
          res,
          params = list(
            group = focus,
            up = input$clinker_up,
            down = input$clinker_down,
            drop = isTRUE(input$clinker_drop)
          )
        )
        
        writeLines(html, file, useBytes = TRUE)
        
      }
      
    )
    
  })
}


# Initial load (once at startup; the Reload button re-reads the files)
INITIAL_STORE <- load_store()

# Report the result in the R console, so loading problems are easy to spot
if (!is.null(INITIAL_STORE$error)) {
  message("DATA NOT LOADED: ", INITIAL_STORE$error)
} else {
  message("Loaded ", nrow(INITIAL_STORE$metadata), " plasmids from:\n  ",
          paste(INITIAL_STORE$info$recognized_files, collapse = "\n  "))
  if (length(INITIAL_STORE$info$notes)) {
    message("Notes:\n  ", paste(INITIAL_STORE$info$notes, collapse = "\n  "))
  }
}


# ---------------------------------------------------------------------
# 8. UI
# ---------------------------------------------------------------------

dash_card <- function(title, ...) {
  card(class = "mb-4 shadow-sm", card_header(title), card_body(...))
}

metric_card <- function(title, output_id, icon_text) {
  card(
    class = "shadow-sm h-100",
    card_body(
      div(class = "text-muted fs-4", icon_text),
      div(class = "text-muted", title),
      textOutput(output_id, container = function(...) h3(class = "mb-0", ...))
    )
  )
}

# Compact metric used in the strip at the top of the AMR tab
overview_metric <- function(label, output_id, icon_text, is_last = FALSE) {
  div(
    class = "d-flex align-items-center gap-2 px-3 py-2 h-100",
    style = if (is_last) NULL else "border-right: 1px solid #D6E5EF;",
    span(icon_text, class = "fs-3",
         style = "color: #168BC4; width: 34px; text-align: center;"),
    div(
      div(toupper(label), class = "fw-bold",
          style = "font-size: 0.72rem; color: #42647D; letter-spacing: 0.03em;"),
      textOutput(output_id, container = function(...)
        div(class = "fw-bold fs-4 lh-1", style = "color: #17324D;", ...))
    )
  )
}

# Gene dropdown + reset button shown above a plot in the AMR tab
gene_filter <- function(select_id, reset_id) {
  tagList(
    tags$label("Resistance Genes", class = "fw-bold"),
    div(
      class = "d-flex align-items-center",
      style = "max-width: 340px; margin-bottom: 0.4rem;",
      div(style = "flex: 1; min-width: 0;",
          selectInput(select_id, label = NULL, choices = c("All" = ALL_VALUE),
                      selected = ALL_VALUE, width = "100%")),
      actionButton(reset_id, "\u21BB", title = "Reset filter",
                   class = "btn-link p-0 ms-2 mb-3",
                   style = "font-size: 1.5rem; line-height: 1;")
    )
  )
}

ui <- page_navbar(
  title = "Plasmid Explorer",
  id = "tabs",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  window_title = "Plasmid Explorer",
  
  # -------------------------------------------------------------------
  # TAB 1: OVERVIEW
  # -------------------------------------------------------------------
  nav_panel(
    "Overview",
    layout_columns(
      col_widths = c(4, 4, 4),
      metric_card("Genome isolates", "n_isolates", "\u25A6"),
      metric_card("Total plasmids",  "n_plasmids", "\u25CE"),
      metric_card("Rep types",       "n_reps",     "\u25A4")
    ),
    dash_card(
      "Isolate source type and country of origin",
      layout_columns(
        col_widths = c(6, 6),
        plotlyOutput("source_plot",  height = "430px"),
        plotlyOutput("country_plot", height = "430px")
      )
    ),
    dash_card("Plasmid size and GC content overview",
              plotlyOutput("size_gc_plot", height = "420px")),
    layout_columns(
      col_widths = c(6, 6),
      dash_card("Species and MLST (ST) lineage distribution",
                plotlyOutput("species_st_plot", height = "450px")),
      dash_card("Plasmid replicon distribution",
                plotlyOutput("overview_plot", height = "450px"))
    ),
    dash_card("Dataset information", verbatimTextOutput("data_status"))
  ),
  
  # -------------------------------------------------------------------
  # TAB 2: AMR GENE ASSOCIATIONS
  # -------------------------------------------------------------------
  nav_panel(
    "AMR gene associations",
    card(
      class = "mb-3 shadow-sm",
      style = "border: 1px solid #D6E5EF; border-radius: 12px;",
      card_body(
        layout_columns(
          col_widths = breakpoints(sm = 6, md = 4, lg = 2),
          gap = 0,
          overview_metric("Total plasmids",          "amr_total_plasmids",         "\u25C9"),
          overview_metric("Total ARG types",         "amr_total_arg_types",        "\u2301"),
          overview_metric("Total Rep types",         "amr_total_rep_types",        "\u232C"),
          overview_metric("Total IS Types",          "amr_total_is_elements",      "\u27F3"),
          overview_metric("Total Transposons Types", "amr_total_transposons",      "\u21C4"),
          overview_metric("Total Virulence Types",   "amr_total_virulence_factors", "\u2726",
                          is_last = TRUE)
        )
      )
    ),
    dash_card(
      "ARG gene frequency and AMR class flows",
      layout_columns(
        col_widths = c(6, 6),
        div(gene_filter("arg_bar_gene_filter", "arg_bar_reset"),
            plotlyOutput("arg_gene_bar", height = "520px")),
        div(gene_filter("amr_sankey_gene_filter", "amr_sankey_reset"),
            plotlyOutput("amr_class_sankey", height = "520px"),
            helpText("Flows count plasmid\u2013AMR class associations; plasmids with ",
                     "multiple classes contribute to each class."))
      )
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      dash_card("MOB type counts",      plotlyOutput("mob_count_plot",       height = "360px")),
      dash_card("MPF type counts",      plotlyOutput("mpf_count_plot",       height = "360px")),
      dash_card("Replicon type counts", plotlyOutput("rep_type_count_plot",  height = "360px"))
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      dash_card("IS element counts",       plotlyOutput("is_element_count_plot",       height = "360px")),
      dash_card("Transposon counts",       plotlyOutput("transposon_count_plot",       height = "360px")),
      dash_card("Virulence factor counts", plotlyOutput("virulence_factor_count_plot", height = "360px"))
    )
  ),
  
  # -------------------------------------------------------------------
  # TAB: GENETIC CONTEXT (Bakta neighbourhood + Clinker comparison)
  # -------------------------------------------------------------------
  nav_panel("Genetic context", contextUI("context")),
  
  # -------------------------------------------------------------------
  # TAB 3: CONSENSUS NETWORK (Pling + mge-cluster)
  # -------------------------------------------------------------------
  nav_panel("Consensus network", consensusUI("consensus")),
  
  # -------------------------------------------------------------------
  # TAB 4: DATA IMPORT / PREVIEW
  # -------------------------------------------------------------------
  nav_panel(
    "Data import",
    dash_card(
      "Input data location",
      p(class = "mb-2",
        "Upload any combination of supported pipeline CSV/TSV outputs. Files are recognized ",
        "from their columns and joined by plasmid accession."),
      fileInput("data_upload", label = NULL, multiple = TRUE,
                accept = c(".csv", ".tsv", ".txt"),
                buttonLabel = "Choose files", placeholder = "or drop files here",
                width = "100%"),
      helpText(class = "d-block mb-3",
               "Supports plasmid metadata, AMRFinder, MOB-suite, IS Finder and ",
               "MobileElementFinder outputs. Missing file types only limit the plots ",
               "that need their fields."),
      textOutput("data_directory", container = tags$p),
      textOutput("metadata_file",  container = tags$p),
      textOutput("distance_file",  container = tags$p),
      actionButton("refresh_data", "Reload input files", class = "btn-outline-primary"),
      span(class = "ms-3", textOutput("reload_message", inline = TRUE))
    ),
    dash_card("Plasmid metadata", DT::DTOutput("metadata_preview")),
    dash_card("Import summary and plot readiness",
              verbatimTextOutput("import_summary")),
    dash_card("Pairwise distance data", DT::DTOutput("distance_preview"))
  ),
  
  footer = div(class = "text-muted text-center py-4", "Plasmid Explorer")
)


# ---------------------------------------------------------------------
# 9. SERVER
# ---------------------------------------------------------------------

server <- function(input, output, session) {
  
  consensusServer("consensus")
  contextServer("context")
  
  # ---- data store ----------------------------------------------------
  store      <- reactiveVal(INITIAL_STORE)
  reload_msg <- reactiveVal("")
  
  # "Reload input files": re-read the data folder. On failure the current
  # data is kept and the error is shown next to the button.
  observeEvent(input$refresh_data, {
    result <- tryCatch(make_store(load_data()), error = function(e) e)
    if (inherits(result, "error")) {
      reload_msg(paste0("Load error: ", conditionMessage(result)))
    } else {
      store(result)
      reload_msg("Default input files loaded.")
    }
  })
  
  # Upload one or more pipeline files (replaces the current data). A set of
  # distance-only files is added to the data that is already loaded.
  observeEvent(input$data_upload, {
    up <- input$data_upload
    result <- tryCatch({
      if (is.null(up) || nrow(up) == 0) {
        stop("Choose one or more CSV/TSV pipeline output files to upload.", call. = FALSE)
      }
      tables <- list()
      notes <- character(0)
      for (i in seq_len(nrow(up))) {
        frame <- tryCatch(read_table_file(up$datapath[i], up$name[i]), error = function(e) e)
        if (inherits(frame, "error")) {
          notes <- c(notes, sprintf("Skipped %s: could not be read (%s).", up$name[i],
                                    conditionMessage(frame)))
        } else {
          tables[[length(tables) + 1]] <- list(name = up$name[i], frame = frame)
        }
      }
      kinds <- vapply(tables, function(t) classify_pipeline_table(t$frame), character(1))
      if (length(tables) && all(kinds == "distances")) {
        current <- store()
        if (is.null(current) || !is.null(current$error)) {
          stop("Upload metadata or pipeline output files before a distance file.", call. = FALSE)
        }
        tables <- c(list(list(name = "Previously loaded plasmid data", frame = current$metadata)),
                    tables)
      }
      make_store(combine_pipeline_tables(tables, notes))
    }, error = function(e) e)
    
    if (inherits(result, "error")) {
      reload_msg(paste0("Load error: ", conditionMessage(result)))
    } else {
      store(result)
      reload_msg("Uploaded data loaded.")
    }
  })
  
  # NULL when no data could be loaded
  frames <- reactive({
    s <- store()
    if (!is.null(s$error)) return(NULL)
    list(metadata = s$metadata, distances = s$distances,
         genes = expand_genes(s$metadata))
  })
  
  # ---- dropdown choices ---------------------------------------------
  observe({
    f <- frames()
    opts <- c("All" = ALL_VALUE)
    if (!is.null(f)) {
      vals <- arg_gene_values(f$metadata)
      opts <- c(opts, setNames(vals, vals))
    }
    for (id in c("arg_bar_gene_filter", "amr_sankey_gene_filter")) {
      current <- isolate(input[[id]])
      selected <- if (!is.null(current) && current %in% opts) current else ALL_VALUE
      updateSelectInput(session, id, choices = opts, selected = selected)
    }
  })
  
  observeEvent(input$arg_bar_reset, {
    updateSelectInput(session, "arg_bar_gene_filter", selected = ALL_VALUE)
  })
  observeEvent(input$amr_sankey_reset, {
    updateSelectInput(session, "amr_sankey_gene_filter", selected = ALL_VALUE)
  })
  
  # ---- overview tab ---------------------------------------------------
  output$n_isolates <- renderText({
    f <- frames(); if (is.null(f)) return("\u2014")
    if (has_usable_values(f$metadata, "isolate_id")) fmt_count(distinct_count(f$metadata$isolate_id)) else "N/A"
  })
  output$n_plasmids <- renderText({
    f <- frames(); if (is.null(f)) return("\u2014")
    fmt_count(nrow(f$metadata))
  })
  output$n_reps <- renderText({
    f <- frames(); if (is.null(f)) return("\u2014")
    if (has_usable_values(f$metadata, "rep_type")) fmt_count(distinct_count(f$metadata$rep_type)) else "N/A"
  })
  
  no_data_msg <- "Load valid input data to display the dashboard."
  
  output$source_plot <- renderPlotly({
    f <- frames(); if (is.null(f)) return(unavailable_figure(no_data_msg))
    make_donut_figure(f$metadata, "source_type", "Source type")
  })
  output$country_plot <- renderPlotly({
    f <- frames(); if (is.null(f)) return(unavailable_figure(no_data_msg))
    make_donut_figure(f$metadata, "country", "Country")
  })
  output$size_gc_plot <- renderPlotly({
    f <- frames(); if (is.null(f)) return(unavailable_figure(no_data_msg))
    make_size_gc_figure(f$metadata)
  })
  output$species_st_plot <- renderPlotly({
    f <- frames(); if (is.null(f)) return(unavailable_figure(no_data_msg))
    make_species_st_figure(f$metadata)
  })
  output$overview_plot <- renderPlotly({
    f <- frames(); if (is.null(f)) return(unavailable_figure(no_data_msg))
    md <- f$metadata
    if (!has_usable_values(md, "rep_type")) {
      return(unavailable_figure("Replicon distribution needs a Rep type column."))
    }
    rep_data <- explode_values(data.frame(plasmid_id = md$plasmid_id, rep_type = md$rep_type,
                                          stringsAsFactors = FALSE), "rep_type")
    rep_data <- rep_data[!is.na(rep_data$rep_type) & rep_data$rep_type != "" &
                           rep_data$rep_type != "-", ]
    counts <- rep_data %>%
      distinct(plasmid_id, rep_type) %>%
      group_by(rep_type) %>%
      summarise(plasmids = n_distinct(plasmid_id), .groups = "drop") %>%
      arrange(plasmids, rep_type)
    if (nrow(counts) == 0) return(unavailable_figure("No replicon types available."))
    hbar_figure(counts$rep_type, counts$plasmids, rep_len(PLOT_COLORS, nrow(counts)),
                hover_label = "Replicon", y_title = "Replicon type",
                margin = list(l = 24, r = 38, t = 24, b = 24))
  })
  
  # Shared by the overview "Dataset information" box and the import tab
  status_text <- reactive({
    f <- frames()
    if (is.null(f)) {
      err <- store()$error %||% "No data loaded. Place the two CSV files in the data directory and reload."
      return(paste0("DATA NOT LOADED\n================\n", err))
    }
    build_data_status(f$metadata, f$distances, store())
  })
  output$data_status    <- renderText(status_text())
  output$import_summary <- renderText(status_text())
  
  # ---- AMR gene association tab ---------------------------------------
  output$amr_total_plasmids <- renderText({
    f <- frames(); if (is.null(f)) return("\u2014")
    fmt_count(n_distinct(f$metadata$plasmid_id))
  })
  delimited_total <- function(column) {
    renderText({
      f <- frames(); if (is.null(f)) return("\u2014")
      if (has_usable_values(f$metadata, column)) fmt_count(count_delimited_values(f$metadata[[column]])) else "N/A"
    })
  }
  output$amr_total_arg_types  <- delimited_total("ARG genes")
  output$amr_total_is_elements <- delimited_total("IS elements")
  output$amr_total_transposons <- delimited_total("Transposons")
  output$amr_total_rep_types <- renderText({
    f <- frames(); if (is.null(f)) return("\u2014")
    if (!has_usable_values(f$metadata, "rep_type")) return("N/A")
    r <- stringr::str_trim(as.character(f$metadata$rep_type))
    fmt_count(length(unique(r[!is.na(r) & r != ""])))
  })
  output$amr_total_virulence_factors <- renderText({
    f <- frames(); if (is.null(f)) return("N/A")
    if (has_usable_values(f$metadata, "Virulence Factors"))
      fmt_count(count_delimited_values(f$metadata[["Virulence Factors"]])) else "N/A"
  })
  
  output$arg_gene_bar <- renderPlotly({
    f <- frames()
    if (is.null(f)) return(unavailable_figure("Load valid input data to display the ARG gene chart."))
    make_arg_gene_bar_figure(f$metadata, input$arg_bar_gene_filter)
  })
  
  output$amr_class_sankey <- renderPlotly({
    f <- frames()
    if (is.null(f)) return(unavailable_figure("Load valid input data to display the AMR class Sankey."))
    make_amr_class_sankey(f$metadata, input$amr_sankey_gene_filter)
  })
  
  count_plot <- function(column, title, label) {
    renderPlotly({
      f <- frames()
      if (is.null(f)) return(unavailable_figure("Load valid input data to display this chart."))
      make_category_count_figure(f$metadata, column, title, label)
    })
  }
  output$mob_count_plot             <- count_plot("MOB type",          "MOB type counts",         "MOB type")
  output$mpf_count_plot             <- count_plot("MPF type",          "MPF type counts",         "MPF type")
  output$rep_type_count_plot        <- count_plot("rep_type",          "Replicon type counts",    "Replicon type")
  output$is_element_count_plot      <- count_plot("IS elements",       "IS element counts",       "IS element")
  output$transposon_count_plot      <- count_plot("Transposons",       "Transposon counts",       "Transposon")
  output$virulence_factor_count_plot <- count_plot("Virulence Factors", "Virulence factor counts", "Virulence factor")
  
  # ---- data import tab --------------------------------------------------
  output$data_directory <- renderText(paste0("Data directory: ", DATA_DIR))
  output$metadata_file  <- renderText({
    s <- store()
    paste0("Input files: ", s$info$metadata_source %||% DATA_DIR)
  })
  output$distance_file  <- renderText({
    s <- store()
    paste0("Distance file: ", s$info$distance_source %||% "Not provided")
  })
  output$reload_message <- renderText(reload_msg())
  
  preview_table <- function(df) {
    DT::datatable(utils::head(df, 100), rownames = FALSE, filter = "top",
                  options = list(pageLength = 10, scrollX = TRUE))
  }
  output$metadata_preview <- DT::renderDT({
    f <- frames(); if (is.null(f)) return(NULL)
    preview_table(f$metadata)
  })
  output$distance_preview <- DT::renderDT({
    f <- frames(); if (is.null(f)) return(NULL)
    preview_table(f$distances)
  })
}

cat("Data directory:", DATA_DIR, "\n")

shinyApp(ui, server)
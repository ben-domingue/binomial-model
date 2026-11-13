# Build data_refs.bib: one BibTeX entry per source study for every IRW table
# used in the empirical evaluation (same table list as make_si_dataset_table.R).
# Cited in si.tex via multibib (\newcites{data}{Data References}).
#
# Matching of analysis table names to the IRW biblio table:
#   1. exact name; 2. case-insensitive; 3. sibling table from the same study
#   (same name minus the final _suffix, case-insensitive). Unresolved tables
#   are written to data_refs_unmatched.txt for manual follow-up.

suppressPackageStartupMessages({ library(irw); library(httr) })

OUT      <- "/home/ben/Dropbox/Apps/Overleaf/binomial_model"
poly_dir <- Sys.getenv("FIG_POLY", "/home/ben/Dropbox/projects/binomial_model/src/irw/poly_compare_data_imv040")

fs  <- list.files(poly_dir, pattern = "^multi_.+\\.rds$", full.names = TRUE)
fs  <- fs[!grepl("multi_(results|summary)\\.rds$", fs)]
raw <- lapply(fs, readRDS)
names(raw) <- sub("^multi_", "", sub("\\.rds$", "", basename(fs)))
raw <- Filter(Negate(is.null), raw)
raw <- lapply(raw, function(r) { if (is.null(r[["K"]])) r$K <- r$K_range + 1L; r })
tables <- sort(names(Filter(function(r) r$K >= 3, raw)))
cat("tables used:", length(tables), "\n")

biblio <- suppressWarnings(irw:::.fetch_biblio_table())
lo <- tolower(biblio$table)

# Tables renamed in IRW since the analysis snapshot (verified by matching N/items)
ALIAS <- c("2024_online_addiction_bsmas" = "ma2026_bsmas",
           "2024_online_addiction_sabas" = "ma2026_sabas",
           "MotAcademica_Ribeiro_2019"   = "ribeiro_2019_academic_motivation",
           "test_taking_much_2025_ef"    = "much_tte_2025_effort")
# Tables with no IRW biblio row: source DOI supplied by hand
MANUAL_DOI <- c("joreskog_moustaki_2001" = "10.1207/S15327906347-387")
# IRW biblio rows whose BibTeX/DOI points at the wrong work: override by hand
OVERRIDE <- c(
  # IRW's DOI is Wyse (2010), a Psychometrika *review* of this book
  "mpsycho_condom" = "@book{deAyala2009,
  author={de Ayala, R. J.},
  title={The Theory and Practice of Item Response Theory},
  publisher={Guilford Press},
  address={New York},
  year={2009}
}")

match_row <- function(t) {
  if (t %in% names(ALIAS)) t <- ALIAS[[t]]
  i <- which(biblio$table == t)
  if (!length(i)) i <- which(lo == tolower(t))
  if (!length(i)) {
    stem <- tolower(sub("_[^_]+$", "", t))
    i <- which(startsWith(lo, paste0(stem, "_")) | lo == stem)
  }
  if (length(i)) i[1] else NA_integer_
}
rows <- vapply(tables, match_row, integer(1))

doi_bib <- function(doi) {
  doi <- sub("^https?://(dx\\.)?doi\\.org/", "", trimws(doi))
  r <- tryCatch(GET(paste0("https://doi.org/", doi),
                    add_headers(Accept = "application/x-bibtex")),
                error = function(e) NULL)
  if (is.null(r) || status_code(r) != 200) return(NA_character_)
  txt <- trimws(content(r, as = "text", encoding = "UTF-8"))
  if (grepl("@", txt)) txt else NA_character_
}

nz <- function(x) !is.na(x) && nzchar(trimws(x)) && !(trimws(x) %in% c("NA", "N/A"))
get_bib <- function(i) {
  if (biblio$table[i] %in% names(OVERRIDE)) return(OVERRIDE[[biblio$table[i]]])
  complete <- function(b) grepl("\\btitle\\s*=", b) && grepl("\\bauthor\\s*=", b)
  bt <- trimws(biblio$BibTex[i])
  if (nz(bt) && grepl("@", bt) && complete(bt)) return(bt)
  doi <- biblio$DOI__for_paper_[i]
  if (!nz(doi)) doi <- biblio$DOI__for_data_[i]
  stem <- sub("_[^_]+$", "", biblio$table[i])
  if (nz(doi)) {
    b <- doi_bib(doi)
    # doi.org metadata is sometimes skeletal (no author/title, e.g. dataset or
    # figure DOIs); then prefer IRW's free-text reference below.
    if (!is.na(b) && complete(b)) return(b)
  }
  ref <- biblio$Reference_x[i]
  if (nz(ref) && !grepl("^https?://\\S+$", trimws(ref)))
    return(sprintf("@misc{%s,\n  note={%s}\n}", stem, trimws(ref)))
  # Last resort: a @misc entry pointing at the data URL
  url <- biblio$URL__for_data_[i]; if (!nz(url)) url <- ref
  if (!nz(url)) return(NA_character_)
  sprintf("@misc{%s,\n  title={%s},\n  howpublished={\\url{%s}},\n  note={Data obtained via the Item Response Warehouse}\n}",
          stem, trimws(biblio$Description[i]), trimws(url))
}

# Clean up common problems in IRW-supplied / doi.org BibTeX
STD_MONTHS <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
sanitize <- function(b) {
  # split off "@type{key," so the key is never touched by the text fixes
  hdr  <- regmatches(b, regexpr("^\\s*@\\w+\\s*\\{[^,]*,", b, perl = TRUE))
  b    <- substring(b, nchar(hdr) + 1)
  key  <- sub("^\\s*@\\w+\\s*\\{\\s*([^,]*),$", "\\1", hdr, perl = TRUE)
  type <- tolower(sub("^\\s*@(\\w+).*$", "\\1", hdr, perl = TRUE))
  key  <- gsub("[^A-Za-z0-9_:.-]", "_", sub("^https?://(dx\\.)?(doi\\.org/)?", "", key))
  if (type == "dataset") type <- "misc"
  # @article without a journal (preprints, OSF records) prints badly: use @misc
  if (type == "article" && !grepl("\\bjournal\\s*=", b)) {
    type <- "misc"
    pub <- regmatches(b, regexpr("\\bpublisher\\s*=\\s*\\{[^}]*\\}", b, perl = TRUE))
    if (length(pub)) b <- paste0(" howpublished=", sub("^[^=]*=\\s*", "", pub), ",", b)
  }
  b <- gsub("\\\\\\\\url", "\\\\url", b)                         # \\url -> \url
  b <- gsub("(?<!\\\\)&", "\\\\&", b, perl = TRUE)               # bare & -> \&
  b <- gsub("(?<!\\\\)\\$", "\\\\$", b, perl = TRUE)             # bare $ -> \$
  # d-stroke has no OT1 glyph (si.tex uses default font encoding): transliterate
  b <- gsub("\\{?\\\\dj\\}?", "dj", gsub("\\{?\\\\DJ\\}?", "Dj", b, perl = TRUE), perl = TRUE)
  b <- gsub("đ", "dj", gsub("Đ", "Dj", b))
  b <- gsub("(?<![{])(https?://[^\\s},]+)", "\\\\url{\\1}", b, perl = TRUE)  # bare URLs
  m <- regmatches(b, regexpr("month\\s*=\\s*[A-Za-z]+", b, perl = TRUE))
  if (length(m)) {
    v <- sub("month\\s*=\\s*", "", m)
    if (!(v %in% STD_MONTHS)) b <- sub(m, sprintf("month={%s}", v), b, fixed = TRUE)
  }
  sprintf("@%s{%s,%s", type, key, b)
}

ok    <- !is.na(rows)
bibs  <- vapply(unique(rows[ok]), get_bib, character(1))
names(bibs) <- unique(rows[ok])
nobib <- tables[ok][is.na(bibs[as.character(rows[ok])])]
man   <- intersect(tables[!ok], names(MANUAL_DOI))
bibs  <- c(bibs, vapply(MANUAL_DOI[man], doi_bib, character(1)))
ok[tables %in% man] <- TRUE

# Re-key every entry as data:<key> so labels cannot collide with refs.bib;
# dedupe on the normalized entry body (several tables share a paper).
bibs <- bibs[!is.na(bibs)]
rekey <- function(b) sub("^(\\s*@\\w+\\s*\\{)\\s*([^,]*),",
                         "\\1data:\\2,", b, perl = TRUE)
bibs <- vapply(bibs, function(b) rekey(sanitize(b)), character(1))
keys <- sub("(?s)^\\s*@\\w+\\s*\\{\\s*([^,]*),.*$", "\\1", bibs, perl = TRUE)
norm <- gsub("\\s+", "", sub("^[^,]*,", "", bibs, perl = TRUE))
# URL-only fallback entries: one per data URL (several tables share a survey)
fb   <- grepl("Data obtained via the Item Response Warehouse", bibs, fixed = TRUE)
norm[fb] <- sub("(?s).*howpublished=\\{(.*?)\\}\\},.*", "\\1", bibs[fb], perl = TRUE)
keep <- !duplicated(norm)
bibs <- bibs[keep]; keys <- keys[keep]
# Distinct entries that happen to share a key: suffix a, b, ...
dup <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
for (k in unique(keys[dup])) {
  j <- which(keys == k)
  new <- paste0(k, letters[seq_along(j)])
  for (m in seq_along(j))
    bibs[j[m]] <- sub(k, new[m], bibs[j[m]], fixed = TRUE)
  keys[j] <- new
}

writeLines(c("% Generated by src/make_si_data_bib.R -- do not edit by hand.",
             "", paste0(bibs, "\n")),
           file.path(OUT, "data_refs.bib"))
cat("wrote", length(bibs), "entries to data_refs.bib\n")

unmatched <- c(if (any(!ok)) paste("NO_BIBLIO_ROW", tables[!ok]),
               if (length(nobib)) paste("NO_BIBTEX_OR_DOI", nobib))
writeLines(as.character(unmatched), file.path(dirname(poly_dir), "..", "data_refs_unmatched.txt"))
cat("unmatched:", length(unmatched), "\n"); if (length(unmatched)) print(unmatched)

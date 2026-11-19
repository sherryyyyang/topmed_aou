#!/usr/bin/env Rscript
# Post-run sanity check on the PC_SETTING output folder.
#
# Two checks, both recorded row by row in
# final_paper_draft_pc<PC>/pipeline_verification.tsv:
#
#   inventory   every published deliverable in final_paper_draft_pc55/ has a
#               counterpart here, or a recorded reason for not having one;
#               every file present only here is classified rather than waved
#               through.
#   png_render  every PNG this run produced decodes, is not truncated and is
#               not blank.
#
# The script stops with a non-zero status if either check finds a problem, so
# the pipeline cannot report success over a missing or broken deliverable.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble); library(png)
})

ROOT <- "/net/topmed11/working/porchard/for-others/julie-and-esther/2025-11-17/sheyang/pre_Watershed_scripts"
source(file.path(ROOT, "scripts/helpers/_paths.R"))

REF <- file.path(ROOT, "final_paper_draft_pc55")
OUT <- FINAL_DRAFT
MIN_BYTES <- 20000L

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
msg("verifying ", OUT, " against ", REF)

rel <- function(dir) {
  if (!dir.exists(dir)) return(character(0))
  sub(paste0("^", dir, "/"), "", list.files(dir, recursive = TRUE, full.names = TRUE))
}
ref_files <- sort(rel(REF))
out_files <- sort(rel(OUT))

# --------------------------------------------------------------------------
# What counts as a deliverable
# --------------------------------------------------------------------------
# Every published deliverable counts, not just figures and tables: the docx,
# the pptx, the hand-authored prose documents and source_panels/as_in_docx/ are
# all files a reader of final_paper_draft_pc55/ would expect to find here too.
# This check used to look at png/pdf/tsv/xlsx only, which is why it once
# reported a clean run while the manuscript, the panel deck and the four prose
# documents had no 2005-PC counterpart at all.
#
# The opposite mistake is just as bad. Editor backups, crash droppings and
# per-run working state are *not* published files: they appear in whichever
# tree someone happened to be working in, and demanding a counterpart for them
# fails the run for reasons that have nothing to do with the analysis. The
# 2026-09-23 01:38 run failed exactly this way, on a
# Supplementary_Figures_methods_BACKUP_<stamp>.docx that an unrelated editing
# session had written into the 55-PC tree one minute after the run started.
#
# The rules below are applied to BOTH trees so the diff stays symmetric, and
# nothing is dropped silently: every skipped file is written to
# pipeline_verification.tsv as check = "ignored" with the reason that skipped
# it. The list is deliberately explicit -- anything a future run leaves behind
# that is not matched here is treated as a deliverable and must be accounted
# for, which is the failure mode we want.
IGNORE_RULES <- tribble(
  ~pattern,                                  ~reason,
  # -- backup copies left by an editor, a human or an agent ------------------
  "(^|/)[^/]*_BACKUP_[^/]*$",                "timestamped backup copy (rollback for a hand edit)",
  "\\.bak([._-][^/]*)?$",                    "'.bak' backup copy",
  "\\.orig$",                                "merge/patch backup copy",
  "~$",                                      "editor backup file",
  "(^|/)#[^/]*#$",                           "Emacs autosave file",
  "\\.sw[op]$",                              "vim swap file",
  # -- per-run working state, not something a reader of the folder opens -----
  "^logs/",                                  "per-run pipeline log",
  "\\.log$",                                 "stray run log outside logs/",
  "\\.rds$",                                 "cached intermediate (no .rds is a published deliverable)",
  # -- partial writes left by an interrupted run -----------------------------
  "\\.(tmp|temp|part|partial|crdownload)$",  "partial or interrupted write",
  "(^|/)\\.nfs[0-9a-f]+$",                   "NFS silly-rename stub for a deleted-but-open file",
  # -- toolchain and operating-system droppings ------------------------------
  "(^|/)__pycache__/",                       "Python bytecode cache",
  "\\.pyc$",                                 "Python bytecode",
  "(^|/)\\.ipynb_checkpoints/",              "Jupyter checkpoint",
  "(^|/)Rplots\\.pdf$",                      "stray default graphics-device output",
  "(^|/)\\.DS_Store$",                       "macOS Finder metadata",
  "(^|/)Thumbs\\.db$",                       "Windows Explorer thumbnail cache",
  # -- package reorganization / non-delivery trees -----------------------------
  "^_package_archive/",                      "archived pre-reorganization copies",
  "^scripts/",                               "figure-generation scripts shipped in package",
  "^docs/",                                  "package documentation",
  "^manuscript/",                            "manuscript subtree (docx also at package root)"
)

# First matching rule wins, so the reason recorded is the most specific one.
ignore_reason <- function(x) {
  out <- rep(NA_character_, length(x))
  for (i in seq_len(nrow(IGNORE_RULES))) {
    hit <- is.na(out) & grepl(IGNORE_RULES$pattern[i], x)
    out[hit] <- IGNORE_RULES$reason[i]
  }
  out
}

ref_skip <- ignore_reason(ref_files)
out_skip <- ignore_reason(out_files)
ref_art  <- ref_files[is.na(ref_skip)]
out_art  <- out_files[is.na(out_skip)]

# The S6 enrichment slice is named after the PC setting it was taken at, so the
# 55-PC name is expected to be absent here and its counterpart expected to be
# extra. Normalize both sides so the diff shows only real gaps.
unpc <- function(x) sub("S6_enrichment_or_[0-9]+pcs\\.tsv$", "S6_enrichment_or_PCS.tsv", x)

# Published files with no counterpart here, each with the reason recorded.
# Anything missing and *not* on this list is a real gap and stops the run.
INTENTIONALLY_ABSENT <- tribble(
  ~file, ~reason
  # (none: every published deliverable now has a PC_SETTING counterpart)
)

# Files that exist only here, by design rather than by accident. These do not
# fail the run, but an extra that is not on this list is called out so a stray
# artifact cannot accumulate unnoticed.
EXPECTED_EXTRA <- tribble(
  ~pattern,                        ~reason,
  "^README\\.md$",                 "this folder's own README (scripts/pc_swap/write_readme.R)",
  "^pipeline_verification\\.tsv$", "this verification report",
  "^figures/",                     "reorganized figures/{main,supplementary}/ delivery tree",
  "^supporting_data/",             "supporting tables / gene lists",
  "^source_panels/",               "individual panel exports",
  "^supp_figs/",                   "supplementary composites and tables",
  "^supplemental_data/",           "source-data tables for figures",
  "^docs/",                        "package documentation",
  "^CHANGELOG",                    "package changelog docs",
  "^PATH_TRACEABILITY",            "path documentation",
  "^CANONICAL",                    "figure index docs",
  "^LAYOUT",                       "layout notes",
  "^WHERE_ARE",                    "figure location notes",
  "pc55_vs_pc2005|PC55_vs_PC2005", "direct 55-vs-2005 comparison; exists only because there are two PC settings to compare"
)
expected_extra_reason <- function(x) {
  out <- rep(NA_character_, length(x))
  for (i in seq_len(nrow(EXPECTED_EXTRA))) {
    hit <- is.na(out) & grepl(EXPECTED_EXTRA$pattern[i], x)
    out[hit] <- EXPECTED_EXTRA$reason[i]
  }
  out
}

msg("published files: ", length(ref_files), " (", sum(!is.na(ref_skip)), " ignored, ",
    length(ref_art), " deliverables)")
msg("produced here:   ", length(out_files), " (", sum(!is.na(out_skip)), " ignored, ",
    length(out_art), " deliverables)")

# --- inventory ------------------------------------------------------------
missing <- ref_art[!unpc(ref_art) %in% unpc(out_art)]
extra   <- out_art[!unpc(out_art) %in% unpc(ref_art)]

missing_reason <- INTENTIONALLY_ABSENT$reason[match(missing, INTENTIONALLY_ABSENT$file)]
unexplained    <- missing[is.na(missing_reason)]
extra_reason   <- expected_extra_reason(extra)
unexpected_extra <- extra[is.na(extra_reason)]

# ifelse() on a zero-length test returns logical(0), which would make `status`
# a logical column on a clean run and break the bind_rows() below.
chr_ifelse <- function(test, yes, no) as.character(ifelse(test, yes, no))

inventory <- bind_rows(
  tibble(check = "inventory", tree = "both",
         file = setdiff(ref_art, missing), status = "MATCHED",
         detail = "published deliverable has a counterpart here"),
  tibble(check = "inventory", tree = "pc55_only",
         file = missing,
         status = chr_ifelse(is.na(missing_reason), "MISSING_UNEXPLAINED", "MISSING_EXPLAINED"),
         detail = chr_ifelse(is.na(missing_reason),
                             "published deliverable with no counterpart here and no recorded reason",
                             missing_reason)),
  tibble(check = "inventory", tree = paste0("pc", PC_SETTING, "_only"),
         file = extra,
         status = chr_ifelse(is.na(extra_reason), "EXTRA_UNEXPECTED", "EXTRA_EXPECTED"),
         detail = chr_ifelse(is.na(extra_reason),
                             "present only here and not on the expected-extra list",
                             extra_reason)),
  tibble(check = "ignored", tree = "pc55",
         file = ref_files[!is.na(ref_skip)], status = "SKIPPED",
         detail = ref_skip[!is.na(ref_skip)]),
  tibble(check = "ignored", tree = paste0("pc", PC_SETTING),
         file = out_files[!is.na(out_skip)], status = "SKIPPED",
         detail = out_skip[!is.na(out_skip)])
)

msg("matched: ", sum(inventory$status == "MATCHED"))
msg("missing vs published: ", length(missing),
    " (", length(unexplained), " unexplained)")
if (length(missing)) cat(paste0("    - ", missing, collapse = "\n"), "\n")
msg("present only here: ", length(extra),
    " (", length(unexpected_extra), " not on the expected-extra list)")
if (length(extra)) cat(paste0("    + ", head(extra, 40), collapse = "\n"), "\n")
msg("ignored as non-deliverable: ", sum(!is.na(ref_skip)), " published / ",
    sum(!is.na(out_skip)), " here")
ign <- inventory |> filter(check == "ignored") |> count(detail, name = "n") |> arrange(desc(n))
if (nrow(ign)) cat(paste0("    ~ ", ign$n, "  ", ign$detail, collapse = "\n"), "\n")

# --- size + decode check on every PNG we produced -------------------------
pngs <- out_files[grepl("\\.png$", out_files) & is.na(out_skip)]
rows <- lapply(pngs, function(p) {
  f <- file.path(OUT, p)
  sz <- file.size(f)
  d <- tryCatch({
    img <- png::readPNG(f)
    # A blank panel has (almost) no variation across pixels.
    list(ok = TRUE, w = dim(img)[2], h = dim(img)[1], sd = stats::sd(as.numeric(img)))
  }, error = function(e) list(ok = FALSE, w = NA, h = NA, sd = NA_real_))
  tibble(check = "png_render", tree = paste0("pc", PC_SETTING),
         file = p, bytes = sz, decodes = d$ok, width = d$w, height = d$h,
         pixel_sd = d$sd,
         status = dplyr::case_when(
           !d$ok ~ "UNREADABLE",
           is.na(sz) | sz < MIN_BYTES ~ "SUSPICIOUSLY_SMALL",
           !is.na(d$sd) & d$sd < 0.01 ~ "LOOKS_BLANK",
           TRUE ~ "OK"),
         detail = "")
})
chk <- bind_rows(rows)

report <- bind_rows(inventory, chk) |>
  select(check, tree, file, status, detail, bytes, decodes, width, height, pixel_sd)
write_tsv(report, file.path(OUT, "pipeline_verification.tsv"))

bad <- filter(chk, status != "OK")
msg("PNGs checked: ", nrow(chk), "; problems: ", nrow(bad))
if (nrow(bad)) print(as.data.frame(bad))
msg("wrote ", file.path(OUT, "pipeline_verification.tsv"))

fail <- character(0)
if (length(unexplained) > 0) {
  fail <- c(fail, paste0(length(unexplained),
    " published deliverable(s) have no counterpart here and no recorded reason:\n  ",
    paste(unexplained, collapse = "\n  ")))
}
if (nrow(bad) > 0) {
  fail <- c(fail, paste0(nrow(bad), " PNG(s) failed the render check:\n  ",
    paste(bad$file, bad$status, collapse = "\n  ")))
}
# An unexpected extra is reported but not fatal: a file present only here
# cannot invalidate a published deliverable, and forcing a listing for every
# scratch file would just recreate the failure this exclusion list fixed.
if (length(unexpected_extra) > 0) {
  warning("present only here and not on the expected-extra list:\n  ",
          paste(unexpected_extra, collapse = "\n  "))
}
if (length(fail) > 0) stop("verification failed: ", paste(fail, collapse = "\n"))

msg("DONE")

# 04_haplogroup_concordance.R
# Agreement between Haplogrep 3 and MToolBox haplogroup assignments.
# Run from project root:  Rscript R/04_haplogroup_concordance.R

source(here::here("R", "00_setup.R"))

# Sheet "Concordance chec_": one row per sample with a call from either tool
conc_raw <- read_excel(data_file, sheet = "Concordance chec_", na = "")

# ---- Major haplogroup from a full haplogroup string ----
# Takes the leading letter plus digits, e.g. "L2a1b1a" -> "L2", "H2a2a1" -> "H2".
# MToolBox may return several candidates separated by ";" or "'"; the first
# candidate is used, and a call is counted only if all candidates share a
# major haplogroup.
major <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "nan")] <- NA
  sub("^([A-Z][0-9]*).*$", "\\1", x)
}
major_set <- function(x) {
  sapply(strsplit(as.character(x), "[;]"), function(v) {
    m <- unique(major(trimws(v)))
    m <- m[!is.na(m)]
    if (length(m) == 1) m else NA_character_
  })
}

conc <- conc_raw %>%
  rename(sample_name = ID,
         haplogrep = Haplogroup_HaploGrep,
         mtoolbox  = `Best predicted haplogroup(s)_MToolBox`) %>%
  mutate(hg_major = major(haplogrep),
         mt_major = major_set(mtoolbox),
         both     = !is.na(hg_major) & !is.na(mt_major),
         agree    = both & hg_major == mt_major) %>%
  left_join(select(dat, sample_name, Sequencing, Country, Haplogroup),
            by = "sample_name")

header("1. Calls returned by each tool")
cat("Samples in concordance set:", nrow(conc), "\n")
cat("Haplogrep 3 call:", sum(!is.na(conc$hg_major)), "\n")
cat("MToolBox call:   ", sum(!is.na(conc$mt_major)), "\n")
cat("Both tools:      ", sum(conc$both), "\n")
cat("Haplogrep only:  ", sum(!is.na(conc$hg_major) & is.na(conc$mt_major)), "\n")
cat("MToolBox only:   ", sum(is.na(conc$hg_major) & !is.na(conc$mt_major)), "\n")

header("2. Agreement at major haplogroup level")
both <- filter(conc, both)
cat("Samples compared:", nrow(both), "\n")
cat("Agree:", sum(both$agree),
    sprintf("(%.1f%%)", 100 * mean(both$agree)), "\n")
cat("Disagree:", sum(!both$agree),
    sprintf("(%.1f%%)", 100 * mean(!both$agree)), "\n")

cat("\nAgreement by sequencing type:\n")
print(both %>% group_by(Sequencing) %>%
        summarise(n = n(), agree = sum(agree),
                  percent = round(100 * sum(agree) / n(), 1), .groups = "drop") %>%
        as.data.frame())

# Cohen's kappa: agreement beyond chance
tab <- table(both$hg_major, both$mt_major)
k <- length(union(rownames(tab), colnames(tab)))
sq <- matrix(0, k, k, dimnames = list(sort(union(rownames(tab), colnames(tab))),
                                      sort(union(rownames(tab), colnames(tab)))))
sq[rownames(tab), colnames(tab)] <- tab
po <- sum(diag(sq)) / sum(sq)
pe <- sum(rowSums(sq) * colSums(sq)) / sum(sq)^2
cat("\nObserved agreement:", round(po, 4),
    "| Expected by chance:", round(pe, 4),
    "| Cohen's kappa:", round((po - pe) / (1 - pe), 4), "\n")

header("3. Agreement for L2 specifically")
l2 <- both %>%
  mutate(hg_L2 = hg_major == "L2", mt_L2 = mt_major == "L2")
print(table(Haplogrep = l2$hg_L2, MToolBox = l2$mt_L2))
cat("\nL2 calls discordant between tools:",
    sum(l2$hg_L2 != l2$mt_L2), "\n")

header("4. Discordant samples")
disc <- filter(both, !agree)
if (nrow(disc) > 0) {
  print(as.data.frame(select(disc, sample_name, Sequencing, Country,
                             haplogrep, mtoolbox, hg_major, mt_major)))
  cat("\nDirection of discordance (Haplogrep -> MToolBox):\n")
  print(table(paste(disc$hg_major, "->", disc$mt_major)))
} else {
  cat("None.\n")
}

header("5. Samples without an MToolBox call")
miss <- filter(conc, !is.na(hg_major), is.na(mt_major))
cat("n =", nrow(miss), "\n")
if (nrow(miss) > 0) {
  cat("\nBy Haplogrep major haplogroup:\n")
  print(table(miss$hg_major))
  cat("\nBy sequencing type:\n")
  print(table(miss$Sequencing, useNA = "ifany"))
}

header("Session info")
print(sessionInfo())

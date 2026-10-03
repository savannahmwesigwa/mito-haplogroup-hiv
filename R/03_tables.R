# 03_tables.R
# Builds Tables 1-2 and Supplementary Tables S1-S7 as Word documents.
# Layout and tests follow the submitted manuscript:
#   Table 1: n (%), Pearson's Chi-squared test (no continuity correction)
#   Table 2: n (%) and median (Q1, Q3); Fisher's exact and Wilcoxon rank sum
#   Supplementary S1-S3: n (%), Fisher's exact test
# Run from project root:  Rscript R/03_tables.R
# Output: output/Table_1.docx, output/Table_2.docx, output/Supplementary_Tables.docx

source(here::here("R", "00_setup.R"))
library(flextable)   # Word tables
library(officer)     # Word document assembly

set_flextable_defaults(font.family = "Calibri", font.size = 10, padding = 3)


# ---- Formatting helpers ----

# p-value: 2 significant digits, "<0.001" below that
fmt_p <- function(p) if (p < 0.001) "<0.001" else format(signif(p, 2), scientific = FALSE)

# Quantile definition: type 6 matches the quartiles in the submitted tables
# (R's default is type 7, which can differ by one unit)
quantile_type <- 6
q1 <- function(v) quantile(v, 0.25, type = quantile_type)
q3 <- function(v) quantile(v, 0.75, type = quantile_type)

p_chisq  <- function(tab) suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
p_fisher <- function(tab) fisher.test(tab)$p.value

# Categorical variable: label row carrying the p-value, then one row per level
cat_rows <- function(df, var, by, var_label, test = p_chisq) {
  tab <- table(df[[var]], df[[by]])
  pct <- round(100 * prop.table(tab, margin = 2))
  cells <- matrix(paste0(tab, " (", pct, "%)"), nrow = nrow(tab),
                  dimnames = dimnames(tab))
  out <- data.frame(Characteristic = c(var_label, paste0("   ", rownames(tab))),
                    stringsAsFactors = FALSE)
  for (g in colnames(tab)) out[[g]] <- c("", cells[, g])
  out$p <- c(fmt_p(test(tab)), rep("", nrow(tab)))
  out
}

# Continuous variable: one row, median (Q1, Q3) per group, Wilcoxon p
num_rows <- function(df, var, by, var_label) {
  d <- df[!is.na(df[[var]]), ]
  out <- data.frame(Characteristic = var_label, stringsAsFactors = FALSE)
  for (g in levels(droplevels(d[[by]]))) {
    v <- d[[var]][d[[by]] == g]
    out[[g]] <- sprintf("%.0f (%.0f, %.0f)", median(v), q1(v), q3(v))
  }
  out$p <- fmt_p(wilcox.test(as.formula(paste(var, "~", by)), data = d,
                             exact = FALSE)$p.value)
  out
}

# Overall column, matching add_overall()
overall_cat <- function(df, var) {
  tab <- table(df[[var]])
  c("", paste0(tab, " (", round(100 * prop.table(tab)), "%)"))
}
overall_num <- function(df, var) {
  v <- df[[var]][!is.na(df[[var]])]
  sprintf("%.0f (%.0f, %.0f)", median(v), q1(v), q3(v))
}

# Group label with sample size, e.g. "LTNP (N = 531)"
n_lab <- function(df, ph) paste0(ph, " (N = ", sum(df$Phenotype == ph), ")")

# flextable: caption above the header, footnotes below
make_ft <- function(df, caption, footnotes, header_labels = NULL) {
  ft <- flextable(df)
  if (!is.null(header_labels)) ft <- set_header_labels(ft, values = header_labels)
  ft <- ft |>
    add_header_lines(values = caption) |>
    bold(part = "header") |>
    bold(i = 1, bold = FALSE, part = "header") |>
    italic(i = 1, part = "header") |>
    fontsize(i = 1, size = 9, part = "header") |>
    align(align = "center", part = "all") |>
    align(j = 1, align = "left", part = "all") |>
    align(i = 1, align = "left", part = "header")
  if (!identical(footnotes, "")) {
    ft <- ft |> add_footer_lines(values = footnotes) |> fontsize(size = 8, part = "footer")
  }
  autofit(ft)
}


# =====================================================================
# Table 1: sex and sequencing by phenotype, stratified by country
# =====================================================================
header("Table 1")

t1_bwr <- bind_rows(
  cat_rows(bwr, "Sex", "Phenotype", "Sex"),
  cat_rows(bwr, "Sequencing", "Phenotype", "Sequencing")
) |> rename(BWR_LTNP = LTNP, BWR_RP = RP, BWR_p = p)

t1_ugr <- bind_rows(
  cat_rows(ugr, "Sex", "Phenotype", "Sex"),
  cat_rows(ugr, "Sequencing", "Phenotype", "Sequencing")
) |> rename(UGR_LTNP = LTNP, UGR_RP = RP, UGR_p = p)

table1 <- bind_cols(t1_bwr, select(t1_ugr, -Characteristic))
print(table1)

# Country spanner sits above the column labels, below the caption
ft1 <- flextable(table1) |>
  set_header_labels(values = list(
    Characteristic = "Characteristic",
    BWR_LTNP = n_lab(bwr, "LTNP"), BWR_RP = n_lab(bwr, "RP"), BWR_p = "p-value",
    UGR_LTNP = n_lab(ugr, "LTNP"), UGR_RP = n_lab(ugr, "RP"), UGR_p = "p-value"
  )) |>
  add_header_row(values = c("", "Botswana (BWR)", "Uganda (UGR)"),
                 colwidths = c(1, 3, 3), top = TRUE) |>
  add_header_lines(values = paste0(
    "Table 1. Demographic and sequencing characteristics of LTNPs and RPs ",
    "stratified by country: Botswana (BWR) and Uganda (UGR).")) |>
  add_footer_lines(values = c("n (%)", "Pearson's Chi-squared test")) |>
  bold(part = "header") |>
  bold(i = 1, bold = FALSE, part = "header") |>
  italic(i = 1, part = "header") |>
  fontsize(i = 1, size = 9, part = "header") |>
  fontsize(size = 8, part = "footer") |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "all") |>
  align(i = 1, align = "left", part = "header") |>
  autofit()

read_docx() |>
  body_add_flextable(ft1) |>
  print(target = file.path(output_dir, "Table_1.docx"))


# =====================================================================
# Table 2: WES subset characteristics
# =====================================================================
header("Table 2")

# Country rows listed BWR first, matching the submitted table
wes <- dat |>
  filter(Sequencing == "WES") |>
  mutate(Country = factor(Country, levels = c("BWR", "UGR")))

table2 <- bind_rows(
  mutate(cat_rows(wes, "Country", "Phenotype", "Country", test = p_fisher),
         Overall = overall_cat(wes, "Country")),
  mutate(cat_rows(wes, "Sex", "Phenotype", "Gender", test = p_fisher),
         Overall = overall_cat(wes, "Sex")),
  mutate(num_rows(wes, "time", "Phenotype", "Time to progression (Months)"),
         Overall = overall_num(wes, "time")),
  mutate(num_rows(wes, "Age_at_sampleCollection_Months", "Phenotype",
                  "Age at enrollment (Months)"),
         Overall = overall_num(wes, "Age_at_sampleCollection_Months")),
  mutate(num_rows(wes, "ART_duration_at_Collection_recalculated", "Phenotype",
                  "Duration of HAART at enrollment (Months)"),
         Overall = overall_num(wes, "ART_duration_at_Collection_recalculated"))
) |>
  select(Characteristic, Overall, LTNP, RP, p)
print(table2)

ft2 <- make_ft(
  table2,
  paste0("Table 2. Demographic and clinical characteristics of the LTNP and RP subset ",
         "from the whole-exome sequencing (WES) batch (N = ", nrow(wes), ")."),
  c("n (%); Median (Q1, Q3)",
    "Fisher's exact test; Wilcoxon rank sum test",
    paste0("Time to progression is time to ART initiation, or time to last ",
           "follow-up for ", sum(wes$event == 0), " censored participants."),
    paste0("Age at enrollment and duration of HAART were available for ",
           sum(!is.na(wes$Age_at_sampleCollection_Months)), " participants.")),
  header_labels = list(
    Characteristic = "Characteristic",
    Overall = paste0("Overall (N = ", nrow(wes), ")"),
    LTNP = n_lab(wes, "LTNP"), RP = n_lab(wes, "RP"), p = "p-value"
  )
)

read_docx() |>
  body_add_flextable(ft2) |>
  print(target = file.path(output_dir, "Table_2.docx"))


# =====================================================================
# Supplementary tables
# =====================================================================
header("Supplementary tables")

supp <- list()

supp$S1 <- list(
  data = cat_rows(dat, "Country", "Phenotype", "Country", test = p_fisher),
  caption = "Supplementary Table 1. Country distribution of LTNPs and RPs.",
  foot = c("n (%)", "Fisher's exact test"),
  labels = list(Characteristic = "Characteristic",
                LTNP = n_lab(dat, "LTNP"), RP = n_lab(dat, "RP"), p = "p-value")
)
supp$S2 <- list(
  data = cat_rows(dat, "Sex", "Phenotype", "Gender", test = p_fisher),
  caption = "Supplementary Table 2. Gender distribution of LTNPs and RPs.",
  foot = c("n (%)", "Fisher's exact test"),
  labels = list(Characteristic = "Characteristic",
                LTNP = n_lab(dat, "LTNP"), RP = n_lab(dat, "RP"), p = "p-value")
)

supp$S3 <- list(
  data = bind_rows(
    cat_rows(ugr, "Phenotype", "L2", "Phenotype", test = p_fisher),
    cat_rows(ugr, "Sex", "L2", "Gender", test = p_fisher),
    cat_rows(ugr, "Sequencing", "L2", "Sequencing", test = p_fisher)
  ) |> select(Characteristic, L2, `Non-L2`, p),
  caption = "Supplementary Table 3. Characteristics of L2 and non-L2 groups, Uganda cohort.",
  foot = c("n (%)", "Fisher's exact test"),
  labels = list(Characteristic = "Characteristic",
                L2 = paste0("L2 (N = ", sum(ugr$L2 == "L2"), ")"),
                `Non-L2` = paste0("non-L2 (N = ", sum(ugr$L2 == "Non-L2"), ")"),
                p = "p-value")
)

hap_table <- function(country_code) {
  filter(dat, Country == country_code) |>
    count(Haplogroup, name = "Frequency") |>
    mutate(`Percentage (%)` = sprintf("%.2f", 100 * Frequency / sum(Frequency))) |>
    arrange(desc(Frequency))
}
supp$S4 <- list(
  data = hap_table("BWR"),
  caption = "Supplementary Table 4. Distribution of mitochondrial haplogroups, Botswana cohort.",
  foot = "", labels = NULL
)
supp$S5 <- list(
  data = hap_table("UGR"),
  caption = "Supplementary Table 5. Distribution of mitochondrial haplogroups, Uganda cohort.",
  foot = "", labels = NULL
)

# S6: Ugandan RPs by birth period, following ART eligibility guideline changes
rp_era <- ugr_rp |>
  mutate(Era = cut(Birth_year, breaks = c(-Inf, 2005, 2010, Inf),
                   labels = c("pre-2006", "2006-2010", "post-2010")))
s6 <- rp_era |>
  count(Era, L2) |>
  tidyr::pivot_wider(names_from = L2, values_from = n, values_fill = 0) |>
  mutate(Total = L2 + `Non-L2`) |>
  select(Era, L2, `Non-L2`, Total) |>
  rename(`Birth period` = Era, `L2 (N)` = L2, `Non-L2 (N)` = `Non-L2`,
         `Total RPs (N)` = Total)
s6 <- bind_rows(s6, summarise(s6, `Birth period` = "Total", across(where(is.numeric), sum)))
supp$S6 <- list(
  data = s6,
  caption = paste0("Supplementary Table 6. Ugandan rapid progressors by birth period and L2 ",
                   "haplogroup status (N = ", nrow(ugr_rp), ")."),
  foot = "Birth periods follow national and WHO ART eligibility guideline changes.",
  labels = NULL
)

# S7: sensitivity analyses of the L2 association among Ugandan RPs
cox_row <- function(fit, label) {
  s <- summary(fit)
  data.frame(
    Analysis = label,
    N = fit$n,
    `L2 HR (95% CI)` = sprintf("%.3f (%.3f-%.3f)",
                               s$conf.int["L2L2", "exp(coef)"],
                               s$conf.int["L2L2", "lower .95"],
                               s$conf.int["L2L2", "upper .95"]),
    `p-value` = sprintf("%.3f", s$coefficients["L2L2", "Pr(>|z|)"]),
    check.names = FALSE
  )
}
s7 <- bind_rows(
  cox_row(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_rp), "All Ugandan RPs"),
  cox_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(ugr_rp, Birth_year < 2010)),
          "Born before 2010"),
  cox_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(ugr_rp, Birth_year <= 2010)),
          "Born 2010 or earlier"),
  cox_row(coxph(Surv(time, event) ~ L2 + Era + Sex, data = rp_era),
          "Adjusted for birth period"),
  cox_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == "pre-2006")),
          "Within pre-2006"),
  cox_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == "2006-2010")),
          "Within 2006-2010"),
  cox_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == "post-2010")),
          "Within post-2010")
)
supp$S7 <- list(
  data = s7,
  caption = paste0("Supplementary Table 7. Sensitivity analyses of the L2 association with ",
                   "time to progression among Ugandan rapid progressors."),
  foot = "Cox proportional hazards models adjusted for sex. Reference: non-L2.",
  labels = NULL
)

doc <- read_docx()
for (nm in names(supp)) {
  print(supp[[nm]]$data)
  cat("\n")
  doc <- doc |>
    body_add_flextable(make_ft(supp[[nm]]$data, supp[[nm]]$caption,
                               supp[[nm]]$foot, supp[[nm]]$labels)) |>
    body_add_par("")
}
print(doc, target = file.path(output_dir, "Supplementary_Tables.docx"))

cat("\nWritten to output/: Table_1.docx, Table_2.docx, Supplementary_Tables.docx\n")

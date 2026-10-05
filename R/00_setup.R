# 00_setup.R
# Packages, paths, data import, derived variables, helper functions.
# Sourced by 01_manuscript_analysis.R and 02_revision_analyses.R.

# ---- Packages ----
library(here)        # project-relative paths
library(readxl)      # read .xlsx
library(dplyr)       # data manipulation
library(tidyr)       # reshaping
library(ggplot2)     # plots
library(survival)    # KM, Cox
library(survminer)   # KM plots
library(waffle)      # Figure 1

# ---- Paths ----
data_file  <- here("data", "Mito_haplo_data.xlsx")
output_dir <- here("output")
dir.create(output_dir, showWarnings = FALSE)

# ---- Data import ----
# Sheet "Combined": one row per participant
raw <- read_excel(data_file, sheet = "Combined", na = "")

# Remove duplicate rows (one participant entered twice)
cat("Rows:", nrow(raw), "| unique sample_name:", n_distinct(raw$sample_name), "\n")
raw <- distinct(raw, sample_name, .keep_all = TRUE)

# ---- Clinical record corrections ----
# Six participants have Time_to_progression missing in the source file. In the
# clinical records their date of progression is a text code, not a date, and
# their time is recorded as 0. All six are LTNPs. Handling:
#   4 with no ART record  -> right-censored at age at enrolment (event = 0)
#   BWR0375, progression date flagged "deleted" but an ART start date present
#     -> ambiguous; held out of survival analyses, see sensitivity below
#   BWR0280, treatment field records HIV-negative status
#     -> ambiguous; held out of survival analyses
# Times are months from birth to enrolment, from the clinical records.
censored_ltnp <- data.frame(
  sample_name = c("BWR0281", "UGR0103", "UGR0285", "UGR0335"),
  cens_time   = c(190.0, 200.8, 242.2, 194.3)
)
# Alternative handling for BWR0375, used in the sensitivity analysis only
bwr0375_event_time <- 141.1

# ---- Derived variables ----
# time:  months from birth to date of progression (Time_to_progression)
# event: all participants reached the endpoint (no censoring)
dat <- raw %>%
  mutate(
    Country    = factor(Country, levels = c("UGR", "BWR")),
    Phenotype  = factor(RP_or_LTNP, levels = c("LTNP", "RP")),
    Sex        = factor(Gender, levels = c("F", "M")),
    Sequencing = factor(Sequencing, levels = c("WES", "WGS")),
    Haplogroup = Haplogroup_2_subclades,
    Subclade   = Haplogroup_3_subclades,
    L2         = factor(if_else(Haplogroup == "L2", "L2", "Non-L2"),
                        levels = c("Non-L2", "L2")),
    Birth_year = as.numeric(format(Birth_date, "%Y")),
    time       = Time_to_progression,
    event      = 1
  ) %>%
  left_join(censored_ltnp, by = "sample_name") %>%
  mutate(
    event = if_else(!is.na(cens_time), 0, event),
    time  = if_else(!is.na(cens_time), cens_time, time)
  ) %>%
  select(-cens_time)

# ---- Analysis subsets ----
ugr <- filter(dat, Country == "UGR")
bwr <- filter(dat, Country == "BWR")

ugr_rp   <- filter(ugr, Phenotype == "RP",   !is.na(time))
ugr_ltnp <- filter(ugr, Phenotype == "LTNP", !is.na(time))
bwr_rp   <- filter(bwr, Phenotype == "RP",   !is.na(time))
bwr_ltnp <- filter(bwr, Phenotype == "LTNP", !is.na(time))

# ---- Helper functions ----

# Section header in console output
header <- function(x) cat("\n\n==========", x, "==========\n")

# Fisher test: counts, p-value, OR with 95% CI (2x2 only)
fisher_report <- function(df, row_var, col_var, label) {
  tab <- table(df[[row_var]], df[[col_var]])
  ft  <- fisher.test(tab)
  cat("\n--", label, "--\n")
  print(tab)
  cat("Fisher p =", signif(ft$p.value, 3))
  if (all(dim(tab) == 2)) {
    cat("; OR =", round(ft$estimate, 3),
        "(95% CI", round(ft$conf.int[1], 3), "-", round(ft$conf.int[2], 3), ")")
  }
  cat("\n")
  invisible(ft$p.value)
}

# Log-rank test: group sizes, p-value
logrank_report <- function(df, group_var, label) {
  f  <- as.formula(paste("Surv(time, event) ~", group_var))
  sd <- survdiff(f, data = df)
  p  <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  cat("\n--", label, "--\n")
  print(sd$n)
  cat("Log-rank p =", signif(p, 3), "\n")
  invisible(p)
}

# Figure output size (inches) and resolution, matching submitted figures
fig2_dpi <- 330
fig2_w   <- 1950 / fig2_dpi          # Figure 2: 1950 x 1696 px
fig2_h   <- 1696 / fig2_dpi

fig_dpi <- 300
fig1_w  <- 1950 / fig_dpi            # Figure 1: 1950 x 1577 px
fig1_h  <- 1577 / fig_dpi
fig3_w  <- 1950 / fig_dpi            # Figure 3: 1950 x 695 px
fig3_h  <-  695 / fig_dpi

# Stack one KM plot above its risk table, for multi-panel figures
# base_size: point size at the final figure dimensions
km_panel <- function(p, heights = c(1, 0.7), base_size = 7) {
  plt <- p$plot +
    theme_bw(base_size = base_size) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          panel.grid = element_blank(),
          plot.margin = margin(2, 4, 2, 2))
  tab <- p$table +
    theme_bw(base_size = base_size) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(size = base_size),
          plot.margin = margin(2, 4, 2, 2))
  ggarrange(plt, tab, ncol = 1, nrow = 2, heights = heights)
}

# Cox model: n, events, concordance, HR, 95% CI, p-value
cox_report <- function(fit, label) {
  s <- summary(fit)
  out <- data.frame(
    term    = rownames(s$coefficients),
    HR      = round(s$conf.int[, "exp(coef)"], 3),
    lower95 = round(s$conf.int[, "lower .95"], 3),
    upper95 = round(s$conf.int[, "upper .95"], 3),
    p       = signif(s$coefficients[, "Pr(>|z|)"], 3),
    row.names = NULL
  )
  cat("\n--", label, "--\n")
  cat("n =", fit$n, "; events =", fit$nevent,
      "; concordance =", round(s$concordance[1], 3),
      "(SE", round(s$concordance[2], 3), ")\n")
  print(out)
  invisible(out)
}

# Pull one term from a cox_report() table as "HR (lower-upper), p"
hr_text <- function(tab, term) {
  r <- tab[tab$term == term, ]
  sprintf("%.3f (%.3f-%.3f), p = %.3f", r$HR, r$lower95, r$upper95, r$p)
}

# KM plot: saved as single PDF, returned for multi-panel assembly
make_km <- function(df, group_var, labels, xlab, file = NULL,
                    palette = c("#E26513", "#406CBF"), break_by = NULL,
                    pval_size = 3, table_fontsize = 2.6) {
  f   <- as.formula(paste("Surv(time, event) ~", group_var))
  fit <- surv_fit(f, data = df)
  p <- ggsurvplot(
    fit, data = df,
    size = 0.3, conf.int = TRUE, conf.int.alpha = 0.3,
    censor.shape = "|", censor.size = 2,
    palette = palette, break.time.by = break_by,
    pval = TRUE, pval.size = pval_size, surv.median.line = "hv",
    risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.4,
    fontsize = table_fontsize,
    legend.labs = labels, xlab = xlab
  )
  if (!is.null(file)) {
    pdf(file.path(output_dir, file), width = 5, height = 4.5, onefile = FALSE)
    print(p)
    invisible(dev.off())
  }
  invisible(p)
}

# Figure output size (inches) and resolution, matching submitted figures
fig2_dpi <- 330
fig2_w   <- 1950 / fig2_dpi          # Figure 2: 1950 x 1696 px
fig2_h   <- 1696 / fig2_dpi

fig_dpi <- 300
fig1_w  <- 1950 / fig_dpi            # Figure 1: 1950 x 1577 px
fig1_h  <- 1577 / fig_dpi
fig3_w  <- 1950 / fig_dpi            # Figure 3: 1950 x 695 px
fig3_h  <-  695 / fig_dpi

# Stack one KM plot above its risk table, for multi-panel figures
# base_size: point size at the final figure dimensions
km_panel <- function(p, heights = c(1, 0.7), base_size = 7) {
  plt <- p$plot +
    theme_bw(base_size = base_size) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          panel.grid = element_blank(),
          plot.margin = margin(2, 4, 2, 2))
  tab <- p$table +
    theme_bw(base_size = base_size) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(size = base_size),
          plot.margin = margin(2, 4, 2, 2))
  ggarrange(plt, tab, ncol = 1, nrow = 2, heights = heights)
}

# Manuscript vs re-run comparison table, printed at end of script
check <- data.frame(item = character(), manuscript = character(),
                    rerun = character(), note = character())
add_check <- function(item, manuscript, rerun, note = "") {
  check <<- rbind(check, data.frame(item = item, manuscript = manuscript,
                                    rerun = as.character(rerun), note = note))
}

# 02_revision_analyses.R
# Additional analyses for revision.
# Run from project root:  Rscript R/02_revision_analyses.R

source(here::here("R", "00_setup.R"))
library(coxphf)      # Firth-penalized Cox regression


# =====================================================================
# 1. Endpoint definition check
# =====================================================================
# Time_to_progression = months from birth to date of progression
# Check: how often progression date differs from ART initiation date
header("1. Endpoint check: time vs birth-to-ART interval")

endpoint <- dat %>%
  filter(!is.na(time), !is.na(Birth_date), !is.na(ART_inititation_date)) %>%
  mutate(birth_to_art = as.numeric(difftime(ART_inititation_date, Birth_date,
                                            units = "days")) / 30.44,
         diff = time - birth_to_art)

cat("Participants with time, birth date and ART date:", nrow(endpoint), "\n")
cat("|time - birth-to-ART| <= 1 month:", sum(abs(endpoint$diff) <= 1), "\n")
cat("|time - birth-to-ART| >  1 month:", sum(abs(endpoint$diff) > 1), "\n")
cat("\nMismatches > 1 month, by country and phenotype:\n")
print(with(filter(endpoint, abs(diff) > 1), table(Country, Phenotype)))
cat("\nSize of mismatch (months), summary:\n")
print(summary(endpoint$diff[abs(endpoint$diff) > 1]))


# =====================================================================
# 1b. Censoring and handling of incomplete clinical records
# =====================================================================
header("1b. Censoring")

cat("Participants with time data:", sum(!is.na(dat$time)), "\n")
cat("Events (ART initiated):", sum(dat$event[!is.na(dat$time)]), "\n")
cat("Right-censored (LTNPs not yet on ART):",
    sum(dat$event[!is.na(dat$time)] == 0), "\n")
cat("Held out of survival analyses (uninterpretable records): 2",
    "(BWR0375, BWR0280)\n\n")
print(filter(dat, event == 0) |>
        select(sample_name, Country, Phenotype, L2, time))

# Sensitivity: three handlings of the incomplete records
sens <- function(df, lab) {
  df <- filter(df, !is.na(time))
  m <- coxph(Surv(time, event) ~ L2 + Sex, data = df)
  s <- summary(m)
  data.frame(Handling = lab, n = m$n, events = m$nevent,
             HR = round(s$conf.int["L2L2", "exp(coef)"], 3),
             lower95 = round(s$conf.int["L2L2", "lower .95"], 3),
             upper95 = round(s$conf.int["L2L2", "upper .95"], 3),
             p = round(s$coefficients["L2L2", "Pr(>|z|)"], 3))
}
for (grp in list(list(d = ugr_ltnp, lab = "Uganda LTNP"),
                 list(d = bwr_ltnp, lab = "Botswana LTNP"),
                 list(d = ugr_rp,   lab = "Uganda RP"))) {
  base <- grp$d
  # (a) as submitted: incomplete records dropped
  a <- sens(filter(base, event == 1), "Records dropped (as submitted)")
  # (b) primary: four LTNPs censored
  b <- sens(base, "Four LTNPs censored (primary)")
  # (c) also treat BWR0375 as an event at its recorded ART start
  alt <- base
  if ("BWR0375" %in% dat$sample_name && grp$lab == "Botswana LTNP") {
    extra <- filter(dat, sample_name == "BWR0375") |>
      mutate(time = bwr0375_event_time, event = 1)
    alt <- bind_rows(base, extra)
  }
  c_ <- sens(alt, "Plus BWR0375 as an event")
  cat("\n--", grp$lab, "--\n")
  print(bind_rows(a, b, c_), row.names = FALSE)
}


# =====================================================================
# 2. Participant flow
# =====================================================================
header("2. Participant flow")

cat("Haplogroup assigned:", nrow(dat), "\n")
print(table(Country = dat$Country, Sequencing = dat$Sequencing))

cat("\nWith time-to-event data:", sum(!is.na(dat$time)), "\n")
print(table(Country = dat$Country, Phenotype = dat$Phenotype,
            has_time = !is.na(dat$time)))

cat("\nWith age at enrollment and HAART duration (Table 1 variables):\n")
print(table(Country = dat$Country, Phenotype = dat$Phenotype,
            has_age = !is.na(dat$Age_at_sampleCollection_Months)))

cat("\nAnalytic samples:\n")
cat("Uganda RP:", nrow(ugr_rp), "| L2:", sum(ugr_rp$L2 == "L2"), "\n")
cat("Uganda LTNP:", nrow(ugr_ltnp), "| L2:", sum(ugr_ltnp$L2 == "L2"), "\n")
cat("Botswana RP:", nrow(bwr_rp), "| L2:", sum(bwr_rp$L2 == "L2"), "\n")
cat("Botswana LTNP:", nrow(bwr_ltnp), "| L2:", sum(bwr_ltnp$L2 == "L2"), "\n")


# =====================================================================
# 3. Included vs excluded from survival analysis
# =====================================================================
header("3. With vs without time-to-event data")

dat_flag <- mutate(dat, has_time = factor(!is.na(time), labels = c("No", "Yes")))
for (v in c("Country", "Phenotype", "Sex", "Sequencing", "L2")) {
  fisher_report(dat_flag, v, "has_time", paste(v, "vs has time data"))
}


# =====================================================================
# 4. Effect sizes: L2 vs RP phenotype (OR, prevalence ratio)
# =====================================================================
# PR: proportion RP in L2 / proportion RP in non-L2, log-scale Wald CI
header("4. Effect sizes: L2 vs RP phenotype")

pr_report <- function(df, label) {
  a <- sum(df$L2 == "L2" & df$Phenotype == "RP");     n1 <- sum(df$L2 == "L2")
  b <- sum(df$L2 == "Non-L2" & df$Phenotype == "RP"); n0 <- sum(df$L2 == "Non-L2")
  pr <- (a / n1) / (b / n0)
  se <- sqrt(1 / a - 1 / n1 + 1 / b - 1 / n0)
  cat("\n--", label, "--\n")
  cat("RP in L2:", a, "/", n1, "; RP in non-L2:", b, "/", n0, "\n")
  cat("PR =", round(pr, 3), "(95% CI", round(exp(log(pr) - 1.96 * se), 3), "-",
      round(exp(log(pr) + 1.96 * se), 3), ")\n")
}

fisher_report(ugr, "L2", "Phenotype", "Uganda: OR")
pr_report(ugr, "Uganda: PR")
fisher_report(bwr, "L2", "Phenotype", "Botswana: OR")
pr_report(bwr, "Botswana: PR")


# =====================================================================
# 5. Country x L2 interaction
# =====================================================================
header("5. Country x L2 interaction")

glm_country <- glm(Phenotype ~ L2 * Country + Sex, family = binomial, data = dat)
cat("\n-- Logistic: RP ~ L2 * Country + Sex --\n")
print(round(cbind(OR = exp(coef(glm_country)),
                  exp(confint.default(glm_country)),
                  p = summary(glm_country)$coefficients[, 4]), 3))

rp_all <- filter(dat, Phenotype == "RP", !is.na(time))
cox_report(coxph(Surv(time, event) ~ L2 * Country + Sex, data = rp_all),
           "Cox: RP, L2 * Country + Sex")


# =====================================================================
# 6. Sequencing sensitivity
# =====================================================================
header("6. Sequencing sensitivity")

cat("Sequencing in Uganda RP survival sample:\n")
print(table(ugr_rp$Sequencing))

fisher_report(filter(ugr, Sequencing == "WES"), "L2", "Phenotype", "Uganda WES: L2 vs phenotype")
fisher_report(filter(ugr, Sequencing == "WGS"), "L2", "Phenotype", "Uganda WGS: L2 vs phenotype")

glm_seq <- glm(Phenotype ~ L2 + Sex + Sequencing, family = binomial, data = ugr)
cat("\n-- Logistic, Uganda: RP ~ L2 + Sex + Sequencing --\n")
print(round(cbind(OR = exp(coef(glm_seq)),
                  exp(confint.default(glm_seq)),
                  p = summary(glm_seq)$coefficients[, 4]), 3))


# =====================================================================
# 7. Proportional hazards assumption
# =====================================================================
header("7. Proportional hazards (cox.zph)")

cat("\n-- Uganda RP: L2 + Sex --\n")
print(cox.zph(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_rp)))
cat("\n-- Uganda LTNP: L2 + Sex --\n")
print(cox.zph(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_ltnp)))


# =====================================================================
# 8. Birth year sensitivity (universal ART for children <3 from 2010)
# =====================================================================
header("8. Birth year sensitivity, Uganda RP")

cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_rp),
           "All Uganda RP")
cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = filter(ugr_rp, Birth_year < 2010)),
           "Born before 2010")
cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = filter(ugr_rp, Birth_year <= 2010)),
           "Born 2010 or earlier")


# =====================================================================
# 9. Birth cohort model, RP only, ART guideline periods
# =====================================================================
# Periods follow ART eligibility changes: CD4 <200 pre-2006,
# CD4 <=350 2006-2010, universal ART for children <3 from 2010
header("9. Birth cohort model, Uganda RP only")

rp_era <- ugr_rp %>%
  mutate(Era = cut(Birth_year,
                   breaks = c(-Inf, 2005, 2010, Inf),
                   labels = c("pre-2006", "2006-2010", "post-2010")))

cat("Uganda RP by birth period and L2:\n")
print(table(Era = rp_era$Era, L2 = rp_era$L2))

cat("\nMedian time (months) by birth period and L2:\n")
print(as.data.frame(rp_era %>% group_by(Era, L2) %>%
                      summarise(n = n(), median_time = median(time), .groups = "drop")))

cox_report(coxph(Surv(time, event) ~ L2 + Era + Sex, data = rp_era),
           "Cox: L2 + Era + Sex")
cox_report(coxph(Surv(time, event) ~ L2 * Era + Sex, data = rp_era),
           "Cox: L2 * Era + Sex")

# L2 HR within each period
for (e in levels(rp_era$Era)) {
  cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == e)),
             paste("Cox within period:", e))
}

# Firth-penalized version of interaction model
cat("\n-- Firth-penalized Cox: L2 * Era + Sex --\n")
firth_fit <- coxphf(Surv(time, event) ~ L2 * Era + Sex, data = rp_era)
print(round(cbind(HR = exp(coef(firth_fit)),
                  lower95 = firth_fit$ci.lower,
                  upper95 = firth_fit$ci.upper,
                  p = firth_fit$prob), 3))


# =====================================================================
# 10. Covariate availability, Uganda RP
# =====================================================================
header("10. Covariate availability, Uganda RP")

ugr_rp %>%
  summarise(n = n(),
            age_at_enrollment = sum(!is.na(Age_at_sampleCollection_Months)),
            HAART_duration = sum(!is.na(ART_duration_at_Collection_recalculated)),
            birth_year = sum(!is.na(Birth_year)),
            sex = sum(!is.na(Sex))) %>%
  as.data.frame() %>%
  print()


# =====================================================================
# 11. Multiple testing: L2 survival tests (Benjamini-Hochberg)
# =====================================================================
header("11. FDR adjustment, L2 survival tests")

lr_p <- function(df, g) {
  sd <- survdiff(as.formula(paste("Surv(time, event) ~", g)), data = df)
  1 - pchisq(sd$chisq, length(sd$n) - 1)
}
l2_rp <- filter(ugr_rp, L2 == "L2") %>%
  mutate(L2b = if_else(Subclade == "L2b", "L2b", "Non L2b"))
l2_ltnp <- filter(ugr_ltnp, L2 == "L2") %>%
  mutate(L2b = if_else(Subclade == "L2b", "L2b", "Non L2b"))

fdr <- data.frame(
  test = c("Uganda RP: L2", "Uganda LTNP: L2", "Botswana RP: L2",
           "Botswana LTNP: L2", "Uganda RP: L2b", "Uganda LTNP: L2b"),
  p_raw = c(lr_p(ugr_rp, "L2"), lr_p(ugr_ltnp, "L2"), lr_p(bwr_rp, "L2"),
            lr_p(bwr_ltnp, "L2"), lr_p(l2_rp, "L2b"), lr_p(l2_ltnp, "L2b"))
)
fdr$p_BH <- p.adjust(fdr$p_raw, method = "BH")
print(mutate(fdr, across(where(is.numeric), ~ signif(.x, 3))))

header("Session info")
print(sessionInfo())


# =====================================================================
# 12. L2b subclade: composition and the discordant haplogroup call
# =====================================================================
# One sample, UGR0185, was discordant between the two haplogroup tools
# (Haplogrep 3: L2b2a; MToolBox: L3i1). Its position in the cohort is
# checked here, together with the composition of the L2b group.
header("12. L2b subclade composition")

cat("UGR0185 in the final dataset:\n")
print(as.data.frame(filter(dat, sample_name == "UGR0185") %>%
                      select(sample_name, Country, Phenotype, Sequencing,
                             Haplogroup, Subclade)))
cat("\nIn the Ugandan RP survival set:", "UGR0185" %in% ugr_rp$sample_name, "\n")

l2_rp <- ugr_rp %>%
  filter(L2 == "L2") %>%
  mutate(L2b = factor(if_else(Subclade == "L2b", "L2b", "Non L2b"),
                      levels = c("L2b", "Non L2b")),
         Subclade = factor(Subclade))

cat("\nL2 rapid progressors by subclade:\n")
print(as.data.frame(l2_rp %>% group_by(Subclade) %>%
                      summarise(n = n(), median_time = median(time),
                                .groups = "drop")))

# Worst case: drop one L2b participant and refit
cat("\nL2b group with each participant removed in turn:\n")
l2b_ids <- filter(l2_rp, L2b == "L2b")$sample_name
loo <- lapply(l2b_ids, function(id) {
  d <- filter(l2_rp, sample_name != id)
  sd <- survdiff(Surv(time, event) ~ L2b, data = d)
  m <- coxph(Surv(time, event) ~ Subclade + Sex, data = d)
  s <- summary(m)
  data.frame(removed = id, L2b_n = sum(d$L2b == "L2b"),
             logrank_p = round(1 - pchisq(sd$chisq, 1), 3),
             HR = round(s$conf.int["SubcladeL2b", "exp(coef)"], 3),
             lower95 = round(s$conf.int["SubcladeL2b", "lower .95"], 3),
             upper95 = round(s$conf.int["SubcladeL2b", "upper .95"], 3),
             p = round(s$coefficients["SubcladeL2b", "Pr(>|z|)"], 3))
})
print(bind_rows(loo), row.names = FALSE)


# =====================================================================
# 13. Haplogroup screen and correction for selection
# =====================================================================
# Each haplogroup with at least 10 carriers is tested against all others
# for association with progression phenotype. The tests are not
# independent, so significance across the screen is assessed by
# permutation of the phenotype labels.
header("13. Haplogroup screen")

min_carriers <- 10
n_perm       <- 10000

screen_one <- function(df, hap_col = "Haplogroup", min_n = min_carriers) {
  haps <- names(which(table(df[[hap_col]]) >= min_n))
  out <- lapply(haps, function(h) {
    ft <- fisher.test(table(df[[hap_col]] == h, df$Phenotype))
    data.frame(haplogroup = h,
               n       = sum(df[[hap_col]] == h, na.rm = TRUE),
               n_RP    = sum(df[[hap_col]] == h & df$Phenotype == "RP", na.rm = TRUE),
               OR      = round(unname(ft$estimate), 3),
               lower95 = round(ft$conf.int[1], 3),
               upper95 = round(ft$conf.int[2], 3),
               p_raw   = signif(ft$p.value, 3))
  })
  out <- do.call(rbind, out)
  out <- out[order(out$p_raw), ]
  out$p_BH <- signif(p.adjust(out$p_raw, "BH"), 3)
  rownames(out) <- NULL
  out
}

# Smallest p-value across the screen, for a given phenotype vector
min_p <- function(df, phenotype, hap_col = "Haplogroup", min_n = min_carriers) {
  haps <- names(which(table(df[[hap_col]]) >= min_n))
  min(sapply(haps, function(h)
    fisher.test(table(df[[hap_col]] == h, phenotype))$p.value))
}

for (cc in c("UGR", "BWR")) {
  d <- filter(dat, Country == cc)
  cat("\n--", cc, "--\n")
  print(screen_one(d))
}

header("14. Permutation test for the Ugandan screen")

set.seed(1)
obs      <- min_p(ugr, ugr$Phenotype)
null_min <- replicate(n_perm, min_p(ugr, sample(ugr$Phenotype)))

cat("Smallest observed p across the screen:", signif(obs, 3), "\n")
cat("Permutation-adjusted p:",
    signif((1 + sum(null_min <= obs)) / (1 + n_perm), 3), "\n")
cat("Permutations with at least one haplogroup at p < 0.05:",
    sprintf("%.0f%%", 100 * mean(null_min < 0.05)), "\n")

# Global test: does phenotype vary across haplogroups at all?
keep <- names(which(table(ugr$Haplogroup) >= min_carriers))
gtab <- table(filter(ugr, Haplogroup %in% keep)$Haplogroup,
              filter(ugr, Haplogroup %in% keep)$Phenotype)
cat("\nGlobal test across haplogroups (Uganda):\n")
print(gtab)
cat("Fisher, simulated p:",
    signif(fisher.test(gtab, simulate.p.value = TRUE, B = 1e5)$p.value, 3), "\n")

header("15. Sensitivity to the carrier threshold")
# The threshold changes the number of tests and therefore the adjustment.
# Reported so the choice is visible rather than implicit.
set.seed(1)
thr <- lapply(c(1, 5, 10, 20, 30), function(m) {
  haps <- names(which(table(ugr$Haplogroup) >= m))
  o    <- min_p(ugr, ugr$Phenotype, min_n = m)
  nm   <- replicate(2000, min_p(ugr, sample(ugr$Phenotype), min_n = m))
  data.frame(min_carriers = m, n_tests = length(haps),
             p_raw = signif(o, 3),
             p_BH  = signif(min(1, o * length(haps)), 3),
             p_perm = signif((1 + sum(nm <= o)) / 2001, 3))
})
print(bind_rows(thr), row.names = FALSE)


# =====================================================================
# 16. Population structure and influential observations
# =====================================================================
# Tribal affiliation is the only marker of population structure in the
# dataset. One Ugandan rapid progressor has a markedly longer time than
# the rest and is checked for influence on the primary estimate.
header("16. Tribe and influential observations")

l2_hr <- function(fit, label) {
  s <- summary(fit)
  cat(sprintf("%-34s n=%3d  HR=%.3f (%.3f-%.3f)  p=%.3f\n", label, fit$n,
              s$conf.int["L2L2", "exp(coef)"], s$conf.int["L2L2", "lower .95"],
              s$conf.int["L2L2", "upper .95"], s$coefficients["L2L2", "Pr(>|z|)"]))
}

rp_tribe <- ugr_rp %>%
  mutate(Tribe_grp = factor(if_else(Tribe == "BAGANDA", "Baganda", "Other"),
                            levels = c("Other", "Baganda")))

cat("Tribal affiliation among Ugandan RPs:\n")
print(table(Tribe = rp_tribe$Tribe_grp, L2 = rp_tribe$L2))
cat("Fisher, L2 vs tribe: p =",
    signif(fisher.test(table(rp_tribe$Tribe_grp, rp_tribe$L2))$p.value, 3), "\n\n")

l2_hr(coxph(Surv(time, event) ~ L2 + Sex, data = rp_tribe), "L2 + sex")
l2_hr(coxph(Surv(time, event) ~ L2 + Sex + Tribe_grp, data = rp_tribe),
      "L2 + sex + tribe")
l2_hr(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_tribe, Tribe_grp == "Baganda")),
      "Baganda only")

cat("\nUgandan RPs with time > 36 months:\n")
print(as.data.frame(filter(ugr_rp, time > 36) %>%
                      select(sample_name, Subclade, Sex, time)))
cat("\n")
l2_hr(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_rp), "All Ugandan RPs")
for (id in filter(ugr_rp, time > 36)$sample_name) {
  l2_hr(coxph(Surv(time, event) ~ L2 + Sex, data = filter(ugr_rp, sample_name != id)),
        paste("Excluding", id))
}

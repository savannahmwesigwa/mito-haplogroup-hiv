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

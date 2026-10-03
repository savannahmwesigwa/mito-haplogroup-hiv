# 01_manuscript_analysis.R
# Reproduces tables, figures and statistics reported in the manuscript.
# Run from project root:  Rscript R/01_manuscript_analysis.R
# Console output: results + manuscript vs re-run check table
# Figures: output/

source(here::here("R", "00_setup.R"))


# =====================================================================
# 1. Cohort description (Supp Tables S1-S2)
# =====================================================================
header("1. Cohort description")

cat("Total participants with haplogroup:", nrow(dat), "\n")
print(table(Country = dat$Country, Phenotype = dat$Phenotype))

p <- fisher_report(dat, "Country", "Phenotype", "Supp S1: country vs phenotype")
add_check("Supp S1 country vs phenotype p", "0.002", signif(p, 2))

p <- fisher_report(dat, "Sex", "Phenotype", "Supp S2: sex vs phenotype")
add_check("Supp S2 sex vs phenotype p", "<0.001", signif(p, 2))

# Minor descriptive p-values below: ms values not reproduced exactly
# by Fisher or chi-square; update ms with re-run values
p <- fisher_report(dat, "Sex", "Country", "Sex vs country")
add_check("Sex vs country p", "0.20", signif(p, 2), "minor; update ms")

p <- fisher_report(dat, "Sequencing", "Country", "Sequencing vs country")
add_check("Sequencing vs country p", "0.30", signif(p, 2), "minor; update ms")

p <- fisher_report(bwr, "Sex", "Phenotype", "Botswana: sex vs phenotype")
add_check("Botswana sex vs phenotype p", "0.30", signif(p, 2), "minor; update ms")

p <- fisher_report(bwr, "Sequencing", "Phenotype", "Botswana: sequencing vs phenotype")
add_check("Botswana sequencing vs phenotype p", "0.20", signif(p, 2), "minor; update ms")

p <- fisher_report(ugr, "Sex", "Phenotype", "Uganda: sex vs phenotype")
add_check("Uganda sex vs phenotype p", "<0.001", signif(p, 2))

p <- fisher_report(ugr, "Sequencing", "Phenotype", "Uganda: sequencing vs phenotype")
add_check("Uganda sequencing vs phenotype p", "0.13", signif(p, 2), "minor; update ms")


# =====================================================================
# 2. Haplogroup distribution (Supp Tables S4-S5, Figure 1)
# =====================================================================
header("2. Haplogroup distribution")

hap_by_country <- dat %>%
  count(Country, Haplogroup) %>%
  group_by(Country) %>%
  mutate(Percent = round(100 * n / sum(n), 2)) %>%
  arrange(Country, desc(n)) %>%
  ungroup()
print(as.data.frame(hap_by_country))

add_check("Supp S4 Botswana L0 n",
          "363", hap_by_country$n[hap_by_country$Country == "BWR" &
                                    hap_by_country$Haplogroup == "L0"])
add_check("Supp S5 Uganda L3 n",
          "212", hap_by_country$n[hap_by_country$Country == "UGR" &
                                    hap_by_country$Haplogroup == "L3"])

# Figure 1: waffle chart by country and phenotype
c25 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black",
         "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon",
         "orchid1", "darkturquoise", "darkorange4", "brown")

hap_counts <- dat %>%
  count(Country, Phenotype, Haplogroup, name = "Frequency") %>%
  mutate(Country = factor(recode(Country, UGR = "Uganda", BWR = "Botswana"),
                          levels = c("Botswana", "Uganda")))

fig1 <- ggplot(hap_counts, aes(fill = Haplogroup, values = Frequency)) +
  geom_waffle(color = "white", size = 0.25, n_rows = 15, flip = TRUE) +
  facet_grid(Phenotype ~ Country) +
  scale_x_discrete() +
  scale_y_continuous(labels = function(x) x * 10, expand = c(0, 0)) +
  scale_fill_manual(values = c25) +
  coord_equal() +
  labs(x = "Country", y = "Frequency") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.ticks.y = element_line(),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.45, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  guides(fill = guide_legend(ncol = 2, title = "Haplogroup"))

for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("Figure1_haplogroup_distribution.", ext)),
         fig1, width = fig1_w, height = fig1_h, units = "in", dpi = fig_dpi)
}


# =====================================================================
# 3. Table 1: WES subset characteristics
# =====================================================================
# Continuous variables: median (IQR), Wilcoxon rank-sum test
# Note: manuscript labels IQR as "95% CI"
header("3. Table 1 (WES subset)")

wes <- filter(dat, Sequencing == "WES")
cat("WES participants:", nrow(wes), "\n")

fisher_report(wes, "Country", "Phenotype", "Table 1: country")
fisher_report(wes, "Sex", "Phenotype", "Table 1: sex")

table1_vars <- c(time       = "Time to progression (months)",
                 Age_at_sampleCollection_Months          = "Age at enrollment (months)",
                 ART_duration_at_Collection_recalculated = "HAART duration at enrollment (months)")

for (v in names(table1_vars)) {
  summ <- wes %>%
    filter(!is.na(.data[[v]])) %>%
    group_by(Phenotype) %>%
    summarise(n = n(),
              median = median(.data[[v]]),
              q1 = quantile(.data[[v]], 0.25),
              q3 = quantile(.data[[v]], 0.75))
  p <- wilcox.test(as.formula(paste(v, "~ Phenotype")), data = wes, exact = FALSE)$p.value
  cat("\n--", table1_vars[v], "--\n")
  print(as.data.frame(summ))
  cat("Wilcoxon p =", signif(p, 3), "\n")
}

med <- function(v, ph) median(wes[[v]][wes$Phenotype == ph], na.rm = TRUE)
add_check("Table 1 time, LTNP median", "151", med("time", "LTNP"),
          "ms used Time_to_Progression_2 (n=520); re-run uses Time_to_progression (n=805)")
add_check("Table 1 time, RP median", "17", med("time", "RP"))
add_check("Table 1 age, LTNP / RP median", "220 / 130",
          paste(med("Age_at_sampleCollection_Months", "LTNP"),
                med("Age_at_sampleCollection_Months", "RP"), sep = " / "))
add_check("Table 1 HAART, LTNP / RP median", "51 / 110",
          paste(med("ART_duration_at_Collection_recalculated", "LTNP"),
                med("ART_duration_at_Collection_recalculated", "RP"), sep = " / "))


# =====================================================================
# 4. L2 vs phenotype (Supp Table S3)
# =====================================================================
header("4. L2 vs phenotype")

p <- fisher_report(ugr, "L2", "Phenotype", "Uganda: L2 vs phenotype")
add_check("Supp S3 Uganda L2 vs phenotype p", "0.018", round(p, 3))

pct_rp <- ugr %>% group_by(L2) %>% summarise(pct_RP = round(100 * mean(Phenotype == "RP")))
print(as.data.frame(pct_rp))
add_check("Supp S3 % RP, L2 / non-L2", "59 / 44",
          paste(pct_rp$pct_RP[pct_rp$L2 == "L2"], pct_rp$pct_RP[pct_rp$L2 == "Non-L2"], sep = " / "))

p <- fisher_report(ugr, "L2", "Sex", "Uganda: L2 vs sex")
add_check("Supp S3 L2 vs sex p", "0.813", round(p, 3))
p <- fisher_report(ugr, "L2", "Sequencing", "Uganda: L2 vs sequencing")
add_check("Supp S3 L2 vs sequencing p", "0.896", round(p, 3))

p <- fisher_report(bwr, "L2", "Phenotype", "Botswana: L2 vs phenotype")
add_check("Botswana L2 vs phenotype p", "0.183", round(p, 3))


# =====================================================================
# 5. Survival: L2 vs non-L2 (Figure 2A-B, Supp Figure S1)
# =====================================================================
# Time origin: birth; endpoint: ART initiation; no censoring
header("5. Survival: L2 vs non-L2")

# Uganda RP (Figure 2A)
p <- logrank_report(ugr_rp, "L2", "Figure 2A: Uganda RP, log-rank")
add_check("Fig 2A Uganda RP log-rank p", "0.009", round(p, 4))

cox_ugr_rp <- cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_rp),
                         "Cox: Uganda RP")
add_check("Cox Uganda RP, L2", "0.594 (0.398-0.886), p = 0.011", hr_text(cox_ugr_rp, "L2L2"))
add_check("Cox Uganda RP, male", "0.961, p = 0.813", hr_text(cox_ugr_rp, "SexM"))

# Uganda LTNP (Figure 2B)
p <- logrank_report(ugr_ltnp, "L2", "Figure 2B: Uganda LTNP, log-rank")
add_check("Fig 2B Uganda LTNP log-rank p", "0.650", round(p, 3),
          "ms KM used Time_to_Progression_2 (44 of 198 LTNPs)")

cox_ugr_ltnp <- cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = ugr_ltnp),
                           "Cox: Uganda LTNP")
add_check("Cox Uganda LTNP, L2", "1.138, p = 0.535", hr_text(cox_ugr_ltnp, "L2L2"))
add_check("Cox Uganda LTNP, male", "1.233 (0.916-1.658), p = 0.167", hr_text(cox_ugr_ltnp, "SexM"))

# Botswana (Supp Figure S1)
p <- logrank_report(bwr_rp, "L2", "Supp S1: Botswana RP, log-rank")
add_check("Supp Fig S1 Botswana RP log-rank p", "0.25", round(p, 3), "not reproduced; regenerate Supp Fig S1")
cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = bwr_rp), "Cox: Botswana RP")

p <- logrank_report(bwr_ltnp, "L2", "Supp S1: Botswana LTNP, log-rank")
add_check("Supp Fig S1 Botswana LTNP log-rank p", "0.29", round(p, 3))
cox_report(coxph(Surv(time, event) ~ L2 + Sex, data = bwr_ltnp), "Cox: Botswana LTNP")

# KM plots: panels A and B of Figure 2, and Supplementary Figure S1
km_2A <- make_km(ugr_rp,   "L2", c("Non L2", "L2"), "RP time to progression (months)",
                 "Figure2A_UGR_RP_KM.pdf", break_by = 10)
km_2B <- make_km(ugr_ltnp, "L2", c("Non L2", "L2"), "LTNP time to progression (months)",
                 "Figure2B_UGR_LTNP_KM.pdf", break_by = 30)

km_S1A <- make_km(bwr_rp,   "L2", c("Non L2", "L2"), "RP time to progression (months)",
                  "SuppS1_BWR_RP_KM.pdf", break_by = 10)
km_S1B <- make_km(bwr_ltnp, "L2", c("Non L2", "L2"), "LTNP time to progression (months)",
                  "SuppS1_BWR_LTNP_KM.pdf", break_by = 30)

# Supplementary Figure S1: two panels
figS1 <- ggarrange(km_panel(km_S1A), km_panel(km_S1B),
                   ncol = 2, nrow = 1, labels = c("A", "B"),
                   font.label = list(size = 11))
for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("SuppFigureS1_BWR_KM.", ext)), figS1,
         width = fig2_w, height = fig2_h / 2, units = "in", dpi = fig2_dpi)
}


# =====================================================================
# 6. L2 subclades (Figure 2C-D)
# =====================================================================
header("6. L2 subclades")

l2_rp <- ugr_rp %>%
  filter(L2 == "L2") %>%
  mutate(L2b = factor(if_else(Subclade == "L2b", "L2b", "Non L2b"),
                      levels = c("L2b", "Non L2b")),
         Subclade = factor(Subclade))
l2_ltnp <- ugr_ltnp %>%
  filter(L2 == "L2") %>%
  mutate(L2b = factor(if_else(Subclade == "L2b", "L2b", "Non L2b"),
                      levels = c("L2b", "Non L2b")))

cat("\nL2 subclades, Uganda RP with time data:\n")
print(table(l2_rp$Subclade))

p <- logrank_report(l2_rp, "L2b", "Figure 2C: L2b vs other L2, Uganda RP")
add_check("Fig 2C L2b log-rank p", "0.047", round(p, 3))

p <- logrank_report(l2_ltnp, "L2b", "Figure 2D: L2b vs other L2, Uganda LTNP")
add_check("Fig 2D L2b LTNP log-rank p", "0.620", round(p, 3))

# Cox: subclade + sex, reference L2a
cox_sub <- cox_report(coxph(Surv(time, event) ~ Subclade + Sex, data = l2_rp),
                      "Cox: L2 subclades, Uganda RP (ref = L2a)")
add_check("Cox subclade L2b", "0.353, p = 0.040", hr_text(cox_sub, "SubcladeL2b"))
add_check("Cox subclade L2c", "0.646, p = 0.559", hr_text(cox_sub, "SubcladeL2c"))
add_check("Cox subclade L2d", "1.854 (0.231-14.890), p = 0.562", hr_text(cox_sub, "SubcladeL2d"))
add_check("Cox subclade male", "2.163, p = 0.045", hr_text(cox_sub, "SexM"))

# Panels C and D of Figure 2; different palette from A-B to mark the subclade analysis
km_2C <- make_km(l2_rp,   "L2b", c("L2b", "Non L2b"), "Time to progression (months)",
                 "Figure2C_L2b_RP_KM.pdf", palette = c("#F8766D", "#00BFC4"), break_by = 10)
km_2D <- make_km(l2_ltnp, "L2b", c("L2b", "Non L2b"), "Time to progression (months)",
                 "Figure2D_L2b_LTNP_KM.pdf", palette = c("#F8766D", "#00BFC4"), break_by = 30)

# Figure 2: four panels, labelled A-D
fig2 <- ggarrange(km_panel(km_2A), km_panel(km_2B),
                  km_panel(km_2C), km_panel(km_2D),
                  ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
                  font.label = list(size = 11))
# 1950 x 1696 px at 330 dpi
for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("Figure2_KM_panels.", ext)), fig2,
         width = fig2_w, height = fig2_h, units = "in", dpi = fig2_dpi)
}


# =====================================================================
# 7. Birth cohort x L2 interaction (Figure 3, Supp Table S6), as submitted
# =====================================================================
# Model data: all Ugandan participants with time data (RP and LTNP)
# Note: manuscript describes this model as RP only
header("7. Birth cohort x L2 interaction (as submitted)")

# Restricted to participants with an observed event, so the figures match
# the submitted model exactly (censored records were absent from it)
cohort_df <- ugr %>%
  filter(!is.na(time), event == 1) %>%
  mutate(Year_group = cut(Birth_year,
                          breaks = seq(1993, 2013, by = 4),
                          right = FALSE,
                          labels = c("1993-1996", "1997-2000", "2001-2004",
                                     "2005-2008", "2009-2012"))) %>%
  filter(!is.na(Year_group))

cat("Supp Table S6 counts (all phenotypes):\n")
print(table(Year_group = cohort_df$Year_group, L2 = cohort_df$L2))
cat("\nSame, split by phenotype:\n")
print(table(Year_group = cohort_df$Year_group, L2 = cohort_df$L2,
            Phenotype = cohort_df$Phenotype))

cox_int_fit <- coxph(Surv(time, event) ~ L2 * Year_group, data = cohort_df)
cox_int <- cox_report(cox_int_fit, "Cox: L2 x birth cohort")

add_check("Fig 3 non-L2 1997-2000 HR", "5.35", round(cox_int$HR[cox_int$term == "Year_group1997-2000"], 2))
add_check("Fig 3 non-L2 2001-2004 HR", "103.27", round(cox_int$HR[cox_int$term == "Year_group2001-2004"], 2))
add_check("Fig 3 non-L2 2005-2008 HR", "884.31", round(cox_int$HR[cox_int$term == "Year_group2005-2008"], 2))
add_check("Fig 3 non-L2 2009-2012 HR", "1161.66", round(cox_int$HR[cox_int$term == "Year_group2009-2012"], 2))
add_check("Fig 3 L2 main effect HR", "2.02, p = 0.155", hr_text(cox_int, "L2L2"))
add_check("Fig 3 L2 x 2001-2004 HR", "0.334, p = 0.047", hr_text(cox_int, "L2L2:Year_group2001-2004"))
add_check("Fig 3 concordance", "0.838 (SE 0.008)",
          sprintf("%.3f (SE %.3f)", summary(cox_int_fit)$concordance[1],
                  summary(cox_int_fit)$concordance[2]))

# Figure 3 as submitted: forest plot of the interaction model above
forest_df <- cox_int |>
  mutate(term = gsub("L2L2", "L2", term),
         term = gsub("Year_group", "", term),
         term = gsub("L2:", "L2 x ", term),
         sig  = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                          p < 0.05 ~ "*", p < 0.1 ~ "\u2022", TRUE ~ "ns"),
         label = sprintf("HR: %.2f (%.2f-%.2f) %s", HR, lower95, upper95, sig))

fig3_old <- ggplot(forest_df, aes(y = term, x = HR, xmin = lower95, xmax = upper95)) +
  geom_point(size = 2) +
  geom_errorbar(orientation = "y", width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "firebrick") +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000, 5000),
                labels = c("0.1", "0.5", "1", "5", "10", "50", "100", "500", "1000", "5000")) +
  geom_text(aes(label = label, x = max(upper95) * 1e4), hjust = 1, size = 3) +
  coord_cartesian(xlim = c(0.1, max(forest_df$upper95) * 1.1e4)) +
  labs(x = "Hazard ratio (log scale)", y = NULL,
       caption = "*** p<0.001, ** p<0.01, * p<0.05, \u2022 p<0.1, ns p>=0.1") +
  theme_minimal()

ggsave(file.path(output_dir, "Figure3_forest_as_submitted.png"),
       fig3_old, width = 8.5, height = 3.54, dpi = 350)


# =====================================================================
# 7b. Birth period model, Ugandan RPs only (revised Figure 3)
# =====================================================================
# Birth periods follow ART eligibility guideline changes:
#   CD4 <200 before 2006; CD4 <=350 from 2006; universal ART for children
#   under 3 from 2010
header("7b. Birth period model, Ugandan RPs only")

rp_era <- ugr_rp |>
  mutate(Era = cut(Birth_year, breaks = c(-Inf, 2005, 2010, Inf),
                   labels = c("pre-2006", "2006-2010", "post-2010")))
print(table(Era = rp_era$Era, L2 = rp_era$L2))

cox_era <- cox_report(coxph(Surv(time, event) ~ L2 + Era + Sex, data = rp_era),
                      "Cox: L2 + birth period + sex")
cox_report(coxph(Surv(time, event) ~ L2 * Era + Sex, data = rp_era),
           "Cox: L2 x birth period + sex")

# L2 hazard ratio within each birth period, plus the adjusted overall estimate
l2_row <- function(fit, label) {
  s <- summary(fit)
  data.frame(label = label, n = fit$n,
             HR = s$conf.int["L2L2", "exp(coef)"],
             lower95 = s$conf.int["L2L2", "lower .95"],
             upper95 = s$conf.int["L2L2", "upper .95"],
             p = s$coefficients["L2L2", "Pr(>|z|)"])
}
fig3_data <- bind_rows(
  l2_row(coxph(Surv(time, event) ~ L2 + Era + Sex, data = rp_era),
         "All RPs, adjusted\nfor birth period"),
  l2_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == "post-2010")),
         "Born 2011-2012"),
  l2_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == "2006-2010")),
         "Born 2006-2010"),
  l2_row(coxph(Surv(time, event) ~ L2 + Sex, data = filter(rp_era, Era == "pre-2006")),
         "Born before 2006")
) |>
  mutate(label  = paste0(label, " (n = ", n, ")"),
         label  = factor(label, levels = label),
         hr_txt = sprintf("%.2f (%.2f-%.2f)", HR, lower95, upper95),
         p_txt  = sprintf("%.3f", p))
print(fig3_data)

# x positions for the right-hand text columns, on the log scale
# Forest panel: estimates and confidence intervals
p_forest <- ggplot(fig3_data, aes(y = label, x = HR, xmin = lower95, xmax = upper95)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "firebrick",
             linewidth = 0.3) +
  geom_hline(yintercept = 1.5, linewidth = 0.3, colour = "grey70") +
  geom_errorbar(orientation = "y", width = 0.12, linewidth = 0.4) +
  geom_point(size = 1.8) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
                labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10")) +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.9))) +
  coord_cartesian(xlim = c(0.1, 10)) +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.title.x = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(linewidth = 0.3),
        axis.ticks.y = element_blank(),
        plot.margin = margin(4, 2, 4, 4))

# Text panel: hazard ratio and p-value columns, aligned to the same rows
p_text <- ggplot(fig3_data, aes(y = label)) +
  geom_hline(yintercept = 1.5, linewidth = 0.3, colour = "grey70") +
  geom_text(aes(x = 0, label = hr_txt), hjust = 0, size = 2.9) +
  geom_text(aes(x = 1, label = p_txt), hjust = 0, size = 2.9) +
  annotate("text", x = c(0, 1), y = 4.7, hjust = 0, size = 2.9,
           fontface = "bold", label = c("HR (95% CI)", "p-value")) +
  scale_x_continuous(limits = c(-0.05, 1.45)) +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.9))) +
  # invisible axis reserves the same vertical space as the forest panel,
  # so the rows in the two panels line up
  labs(x = "Hazard ratio", y = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "white"),
        axis.title.x = element_text(colour = "white", size = 8),
        plot.margin = margin(4, 4, 4, 0))

fig3 <- ggarrange(p_forest, p_text, ncol = 2, widths = c(1.9, 1))

# 1300 x 620 px at 300 dpi
fig3_w <- 1300 / fig_dpi
fig3_h <-  620 / fig_dpi
for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("Figure3_L2_by_birth_period.", ext)),
         fig3, width = fig3_w, height = fig3_h, units = "in", dpi = fig_dpi)
}


# Supplementary Figure S2: survival curves for the two informative periods
# post-2010 omitted (n = 13, one L2 carrier)
km_pre06 <- make_km(filter(rp_era, Era == "pre-2006"), "L2", c("Non L2", "L2"),
                    "Time to progression (months)", break_by = 10)
km_0610  <- make_km(filter(rp_era, Era == "2006-2010"), "L2", c("Non L2", "L2"),
                    "Time to progression (months)", break_by = 10)

figS2 <- ggarrange(km_panel(km_pre06), km_panel(km_0610),
                   ncol = 2, nrow = 1, labels = c("A", "B"),
                   font.label = list(size = 11))
for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("SuppFigureS2_RP_by_birth_period_KM.", ext)),
         figS2, width = fig2_w, height = fig2_h / 2, units = "in", dpi = fig2_dpi)
}


# =====================================================================
# 7c. Participant flow diagram (Figure 4)
# =====================================================================
# Two analysis streams:
#   all participants -> haplogroup vs phenotype association
#   participants with time-to-event data -> survival analyses
# Eligible counts at the two centres come from the parent CAfGEN study;
# all later counts are computed from the data.
header("7c. Participant flow")

n_hap  <- nrow(dat)
n_wes  <- sum(dat$Sequencing == "WES")
n_wgs  <- sum(dat$Sequencing == "WGS")
n_time <- sum(!is.na(dat$time))
n_wes_no_time <- sum(dat$Sequencing == "WES" & is.na(dat$time))
n_l2_rp <- sum(ugr_rp$L2 == "L2")
sub_n   <- table(filter(ugr_rp, L2 == "L2")$Subclade)

cat("Haplogroup assigned:", n_hap, " (WES", n_wes, ", WGS", n_wgs, ")\n")
cat("Time-to-event data:", n_time, " (", sum(dat$event[!is.na(dat$time)] == 0),
    "censored )\n")
cat("Uganda:", nrow(ugr_ltnp), "LTNP,", nrow(ugr_rp), "RP |  Botswana:",
    nrow(bwr_ltnp), "LTNP,", nrow(bwr_rp), "RP\n")
cat("L2 subclades among Ugandan RPs with L2:",
    paste(names(sub_n), sub_n, collapse = ", "), "\n")

box <- function(x, y, label) data.frame(x = x, y = y, label = label,
                                        stringsAsFactors = FALSE)

main <- bind_rows(
  box(1.5, 5, paste0(
    "Eligible at the two clinical centres (parent CAfGEN study)\n",
    "Uganda: 440 LTNPs, 275 RPs      Botswana: 126 LTNPs, 500 RPs")),
  box(1.5, 4, paste0(
    "Recruited, sequenced and assigned a mitochondrial haplogroup\n",
    "N = ", n_hap, "   (WES ", n_wes, ", WGS ", n_wgs, ")"))
)

left <- bind_rows(
  box(0.72, 3, paste0(
    "Haplogroup and phenotype analysis\n",
    "All participants, N = ", n_hap, "\n",
    "Uganda ", nrow(ugr), ", Botswana ", nrow(bwr))),
  box(0.72, 2, paste0(
    "L2 vs non-L2 by phenotype\n",
    "Uganda and Botswana, WES and WGS"))
)

right <- bind_rows(
  box(2.28, 3, paste0(
    "Survival analyses\n",
    "Time-to-event data available, N = ", n_time,
    " (", sum(dat$event[!is.na(dat$time)] == 0), " censored)\n",
    "Uganda ", nrow(ugr_ltnp) + nrow(ugr_rp), ", Botswana ",
    nrow(bwr_ltnp) + nrow(bwr_rp))),
  box(2.28, 2, paste0(
    "Uganda: ", nrow(ugr_rp), " RPs, ", nrow(ugr_ltnp), " LTNPs\n",
    "Botswana: ", nrow(bwr_rp), " RPs, ", nrow(bwr_ltnp), " LTNPs")),
  box(2.28, 1, paste0(
    "L2 subclade analysis, Ugandan RPs with L2\n",
    "N = ", n_l2_rp, "   (L2a ", sub_n[["L2a"]], ", L2b ", sub_n[["L2b"]],
    ", L2c ", sub_n[["L2c"]], ", L2d ", sub_n[["L2d"]], ")"))
)

note <- box(3.6, 3, paste0(
  "No linked time-to-event data\n",
  "n = ", n_wgs, " (WGS second batch)\n",
  "Records uninterpretable: n = 2"))

arrows <- bind_rows(
  data.frame(x = 1.5,  xend = 1.5,  y = 4.72, yend = 4.3),   # box 1 to box 2
  data.frame(x = 1.5,  xend = 0.72, y = 3.68, yend = 3.42),  # split, left
  data.frame(x = 1.5,  xend = 2.28, y = 3.68, yend = 3.42),  # split, right
  data.frame(x = 0.72, xend = 0.72, y = 2.6,  yend = 2.3),
  data.frame(x = 2.28, xend = 2.28, y = 2.6,  yend = 2.35),
  data.frame(x = 2.28, xend = 2.28, y = 1.65, yend = 1.35),
  data.frame(x = 2.82, xend = 3.0,  y = 3,    yend = 3)      # to side note
)

fig4 <- ggplot() +
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
               linewidth = 0.3) +
  geom_label(data = bind_rows(main, left, right),
             aes(x = x, y = y, label = label),
             size = 2.3, lineheight = 1.2,
             label.padding = unit(0.22, "lines"), fill = "white") +
  geom_label(data = note, aes(x = x, y = y, label = label),
             size = 2.1, lineheight = 1.2,
             label.padding = unit(0.22, "lines"), fill = "grey96") +
  scale_x_continuous(limits = c(0, 4.2)) +
  scale_y_continuous(limits = c(0.6, 5.4)) +
  theme_void()

for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("Figure4_participant_flow.", ext)),
         fig4, width = 1950 / fig_dpi, height = 1300 / fig_dpi,
         units = "in", dpi = fig_dpi)
}


# =====================================================================
# 7d. Causal diagram (Figure 5)
# =====================================================================
# Assumed relationships between the exposure (L2 haplogroup), the outcome
# (time to progression) and other variables. Solid nodes are measured;
# dashed nodes are unmeasured in this cohort.
header("7d. Causal diagram")

nodes <- data.frame(
  name = c("L2", "Outcome", "Sex", "Birth", "Country", "Select", "Unmeas"),
  x    = c(1.0,  5.0,       3.0,   4.0,     1.0,       2.3,      5.0),
  y    = c(2.1,  2.1,       3.4,   1.0,     3.4,       0.6,      3.4),
  label = c("Mitochondrial\nL2 haplogroup",
            "Time to\nprogression",
            "Sex",
            "Birth period\n(ART eligibility)",
            "Country /\npopulation\nancestry",
            "Selection into\nextreme-phenotype\ncohort [conditioned]",
            "Nutrition,\nco-infections,\nsocioeconomic\nconditions"),
  measured = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
pos <- function(n, what) nodes[[what]][nodes$name == n]

edges <- data.frame(
  from = c("L2", "Sex", "Birth", "Country", "Country", "Unmeas", "L2", "Birth"),
  to   = c("Outcome", "Outcome", "Outcome", "L2", "Outcome", "Outcome",
           "Select", "Select"),
  stringsAsFactors = FALSE
) |>
  mutate(x = sapply(from, pos, "x"), y = sapply(from, pos, "y"),
         xend = sapply(to, pos, "x"), yend = sapply(to, pos, "y"))

# shorten each arrow so it stops short of the node labels
shrink <- function(d, f = 0.72) {
  dx <- d$xend - d$x; dy <- d$yend - d$y
  d$x <- d$x + dx * (1 - f) / 2; d$y <- d$y + dy * (1 - f) / 2
  d$xend <- d$xend - dx * (1 - f) / 2; d$yend <- d$yend - dy * (1 - f) / 2
  d
}
edges <- shrink(edges)

fig5 <- ggplot() +
  geom_segment(data = filter(edges, from != "Unmeas"),
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               linewidth = 0.3) +
  geom_segment(data = filter(edges, from == "Unmeas"),
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               linewidth = 0.3, linetype = "dashed") +
  geom_label(data = filter(nodes, measured),
             aes(x = x, y = y, label = label),
             size = 2.4, lineheight = 1.15, fill = "white") +
  geom_label(data = filter(nodes, !measured),
             aes(x = x, y = y, label = label),
             size = 2.4, lineheight = 1.15, fill = "grey95",
             linetype = "dashed") +
  labs(caption = paste("Arrows show assumed causal relationships.",
                       "Grey dashed node: unmeasured in this cohort.")) +
  scale_x_continuous(limits = c(0.1, 5.9)) +
  scale_y_continuous(limits = c(0.1, 4.0)) +
  theme_void(base_size = 8) +
  theme(plot.caption = element_text(hjust = 0, size = 6, colour = "grey30"))

for (ext in c("pdf", "tiff")) {
  ggsave(file.path(output_dir, paste0("Figure5_causal_diagram.", ext)),
         fig5, width = 1600 / fig_dpi, height = 1100 / fig_dpi,
         units = "in", dpi = fig_dpi)
}


# =====================================================================
# 8. Manuscript vs re-run check table
# =====================================================================
header("8. Manuscript vs re-run")
options(width = 200)
print(check, right = FALSE, row.names = FALSE)

header("Session info")
print(sessionInfo())

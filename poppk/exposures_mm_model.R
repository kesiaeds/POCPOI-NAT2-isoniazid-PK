## Compute model-derived AUC0-24, Cmax, C24 per subject-occasion from the
## final MM popPK fit (fit_Km10.rds), using individual EBE parameters and
## dense rxode2 simulation. Then re-run paired Wilcoxon and McNemar.
##
## Inputs:
##   /Users/kesiaeds/Documents/POCPOI_dose_finding/popPK_model/fit_Km10.rds
##   /Users/kesiaeds/Documents/POCPOI_dose_finding/data/pk_nonmem_combined.csv
## Outputs (do not overwrite NCA versions):
##   /tmp/mm_exposures/exposures_mm_model.csv
##   /tmp/mm_exposures/paired_stats_mm_model.csv
##   /tmp/mm_exposures/summary_phenotype_dosing_mm_model.csv
##   /tmp/mm_exposures/mcnemar_supra_C24_mm_model.txt
##   /tmp/mm_exposures/Figure_exposures_MM.pdf  (3-panel boxplot)

suppressPackageStartupMessages({
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

OUT <- "/tmp/mm_exposures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SUPRA_C24 <- 0.15

# ----------------------------------------------------------------------------
# 1. Load fit + raw data, extract per-subject-occasion individual parameters
# ----------------------------------------------------------------------------
fit <- readRDS("/Users/kesiaeds/Documents/POCPOI_dose_finding/popPK_model/fit_Km10.rds")
fd  <- as.data.frame(fit)

# Map fit-internal ID (1..104) back to original subject IDs.
raw <- read_csv("/Users/kesiaeds/Documents/POCPOI_dose_finding/data/pk_nonmem_combined.csv",
                show_col_types = FALSE) |>
  filter(!is.na(ACE_PROFILE), WT > 0) |>
  mutate(SUBJID = ID, ID_int = as.integer(factor(ID)))

# One row per (ID, dosenum) holding the individual MM parameter vector.
ind_pars <- fd |>
  mutate(ID = as.integer(as.character(ID))) |>
  group_by(ID, dosenum) |>
  summarise(Ka = first(Ka), Vmax = first(Vmax), Km = first(Km),
            V2 = first(V2), Q = first(Q), V3 = first(V3),
            WT = first(WT), ACE_PROFILE = first(ACE_PROFILE),
            eta_vmax = first(eta.vmax), .groups = "drop")
cat("Individual parameter sets:", nrow(ind_pars), "subject-occasions\n")

# Get the dose amount in mg per (ID_int, OCC) from raw data.
dose_lookup <- raw |>
  filter(EVID == 1) |>
  distinct(ID_int, OCC, AMT, DAY, DOSE) |>
  rename(ID = ID_int, dosenum = OCC, AMT_MG = AMT)

ind_pars <- ind_pars |> left_join(dose_lookup, by = c("ID", "dosenum"))
stopifnot(all(!is.na(ind_pars$AMT_MG)))

# Map fit-internal ID back to original SUBJID for traceability
id_map <- raw |> distinct(ID_int, SUBJID) |> rename(ID = ID_int)
ind_pars <- ind_pars |> left_join(id_map, by = "ID")

# ----------------------------------------------------------------------------
# 2. rxode2 simulation model (parameters fixed per subject; no eta sampling)
# ----------------------------------------------------------------------------
sim_mod <- rxode2::rxode2({
  Cp = central / V2
  d/dt(depot)      = -Ka * depot
  d/dt(central)    =  Ka * depot - Vmax * Cp / (Km + Cp) -
                       Q * (Cp - peripheral / V3)
  d/dt(peripheral) =  Q * (Cp - peripheral / V3)
})

# Dense observation times: 0, 0.05, 0.10, ..., 24 h (1-min grid).
sample_times <- seq(0, 24, by = 1/60)

# Simulate each subject-occasion individually, then collect.
sim_one <- function(row) {
  pars <- c(Ka = row$Ka, Vmax = row$Vmax, Km = row$Km,
            V2 = row$V2, Q = row$Q, V3 = row$V3)
  # Build event data frame: one dose row + sampling rows
  ev <- data.frame(
    id   = 1L,
    time = c(0, sample_times),
    amt  = c(row$AMT_MG, rep(0, length(sample_times))),
    cmt  = c("depot", rep("(obs)", length(sample_times))),
    evid = c(1L, rep(0L, length(sample_times)))
  )
  s <- rxode2::rxSolve(sim_mod, params = pars, events = ev,
                       returnType = "data.frame")
  data.frame(ID = row$ID, dosenum = row$dosenum,
             SUBJID = row$SUBJID, WT = row$WT,
             ACE_PROFILE = row$ACE_PROFILE, AMT_MG = row$AMT_MG,
             time = s$time, Cp = s$Cp)
}

cat("Simulating", nrow(ind_pars), "dense profiles...\n")
sim_all <- do.call(rbind, lapply(seq_len(nrow(ind_pars)),
                                  function(i) sim_one(ind_pars[i, ])))
cat("Dense rows:", nrow(sim_all), "\n")

# ----------------------------------------------------------------------------
# 3. Per-subject-occasion AUC0-24, Cmax, C24 from dense IPRED
# ----------------------------------------------------------------------------
exposures <- sim_all |>
  group_by(ID, dosenum, SUBJID, WT, ACE_PROFILE, AMT_MG) |>
  summarise(
    AUC24 = {
      o <- order(time); tt <- time[o]; cc <- Cp[o]
      sum(diff(tt) * (head(cc, -1) + tail(cc, -1)) / 2)
    },
    Cmax  = max(Cp),
    C24   = approx(time, Cp, xout = 24)$y,
    .groups = "drop")

PHEN  <- c("1" = "Slow", "2" = "Intermediate", "3" = "Rapid")
OCCN  <- c("1" = "Genotype-guided", "2" = "Flat 900 mg")
exposures$Phenotype <- PHEN[as.character(exposures$ACE_PROFILE)]
exposures$Dosing    <- OCCN[as.character(exposures$dosenum)]
exposures$Dosing    <- factor(exposures$Dosing,
                              levels = c("Genotype-guided", "Flat 900 mg"))
exposures$Phenotype <- factor(exposures$Phenotype,
                              levels = c("Slow", "Intermediate", "Rapid"))

write.csv(exposures, file.path(OUT, "exposures_mm_model.csv"), row.names = FALSE)
cat("Wrote exposures_mm_model.csv (", nrow(exposures), "rows )\n")

# ----------------------------------------------------------------------------
# 4. Summary by phenotype x dosing strategy (medians + IQR)
# ----------------------------------------------------------------------------
med_iqr <- function(x) sprintf("%.2f (%.2f-%.2f)",
                               median(x, na.rm = TRUE),
                               quantile(x, .25, na.rm = TRUE),
                               quantile(x, .75, na.rm = TRUE))
summ <- exposures |>
  group_by(Phenotype, Dosing) |>
  summarise(N = n(),
            AUC24_med_IQR = med_iqr(AUC24),
            Cmax_med_IQR  = med_iqr(Cmax),
            C24_med_IQR   = med_iqr(C24),
            pct_C24_supra = sprintf("%.1f%%", 100 * mean(C24 > SUPRA_C24, na.rm = TRUE)),
            .groups = "drop")
write.csv(summ, file.path(OUT, "summary_phenotype_dosing_mm_model.csv"),
          row.names = FALSE)
cat("\n=== Summary by phenotype x dosing (MM model-derived) ===\n")
print(as.data.frame(summ))

# ----------------------------------------------------------------------------
# 5. Within-subject paired Wilcoxon (Day 7 guided vs Day 14 flat 900 mg)
# ----------------------------------------------------------------------------
wide <- exposures |>
  pivot_wider(id_cols = c(SUBJID, Phenotype),
              names_from = Dosing,
              values_from = c(AUC24, Cmax, C24)) |>
  drop_na()

ps <- function(label, a, b) {
  n <- length(a)
  if (n < 2) return(data.frame(Metric = label, N_pairs = n,
                               guided_med = NA, flat_med = NA, p = NA))
  pv <- tryCatch(wilcox.test(a, b, paired = TRUE, exact = FALSE)$p.value,
                 error = function(e) NA)
  data.frame(Metric = label, N_pairs = n,
             guided_med = sprintf("%.2f (%.2f-%.2f)", median(a),
                                  quantile(a, .25), quantile(a, .75)),
             flat_med   = sprintf("%.2f (%.2f-%.2f)", median(b),
                                  quantile(b, .25), quantile(b, .75)),
             p = signif(pv, 3))
}

paired <- bind_rows(
  ps("AUC24 (all)",
     wide$`AUC24_Genotype-guided`, wide$`AUC24_Flat 900 mg`),
  ps("AUC24 (Slow)",
     wide$`AUC24_Genotype-guided`[wide$Phenotype == "Slow"],
     wide$`AUC24_Flat 900 mg`[wide$Phenotype == "Slow"]),
  ps("AUC24 (Intermediate)",
     wide$`AUC24_Genotype-guided`[wide$Phenotype == "Intermediate"],
     wide$`AUC24_Flat 900 mg`[wide$Phenotype == "Intermediate"]),
  ps("AUC24 (Rapid)",
     wide$`AUC24_Genotype-guided`[wide$Phenotype == "Rapid"],
     wide$`AUC24_Flat 900 mg`[wide$Phenotype == "Rapid"]),
  ps("Cmax (Slow)",
     wide$`Cmax_Genotype-guided`[wide$Phenotype == "Slow"],
     wide$`Cmax_Flat 900 mg`[wide$Phenotype == "Slow"]),
  ps("Cmax (Intermediate)",
     wide$`Cmax_Genotype-guided`[wide$Phenotype == "Intermediate"],
     wide$`Cmax_Flat 900 mg`[wide$Phenotype == "Intermediate"]),
  ps("Cmax (Rapid)",
     wide$`Cmax_Genotype-guided`[wide$Phenotype == "Rapid"],
     wide$`Cmax_Flat 900 mg`[wide$Phenotype == "Rapid"]),
  ps("C24 (Slow)",
     wide$`C24_Genotype-guided`[wide$Phenotype == "Slow"],
     wide$`C24_Flat 900 mg`[wide$Phenotype == "Slow"]),
  ps("C24 (Intermediate)",
     wide$`C24_Genotype-guided`[wide$Phenotype == "Intermediate"],
     wide$`C24_Flat 900 mg`[wide$Phenotype == "Intermediate"]),
  ps("C24 (Rapid)",
     wide$`C24_Genotype-guided`[wide$Phenotype == "Rapid"],
     wide$`C24_Flat 900 mg`[wide$Phenotype == "Rapid"])
)
write.csv(paired, file.path(OUT, "paired_stats_mm_model.csv"), row.names = FALSE)
cat("\n=== Paired Wilcoxon (MM model-derived) ===\n")
print(paired)

# ----------------------------------------------------------------------------
# 6. McNemar on C24 > 0.15 µg/mL (paired binary supratherapeutic flag)
# ----------------------------------------------------------------------------
mc <- wide |>
  transmute(SUBJID, Phenotype,
            supra_guided = `C24_Genotype-guided` > SUPRA_C24,
            supra_flat   = `C24_Flat 900 mg`     > SUPRA_C24)
tab <- table(guided = mc$supra_guided, flat = mc$supra_flat)
mc_test <- mcnemar.test(tab)
sink(file.path(OUT, "mcnemar_supra_C24_mm_model.txt"))
cat("McNemar (paired) on C24 > ", SUPRA_C24, " µg/mL (MM model-derived):\n", sep="")
cat("n pairs =", nrow(mc), "\n\n")
cat("Contingency:\n")
print(tab)
cat("\n")
print(mc_test)
cat("\n% supratherapeutic by phenotype × dosing:\n")
by_phen <- mc |>
  group_by(Phenotype) |>
  summarise(n = n(),
            pct_supra_guided = sprintf("%.1f%%", 100 * mean(supra_guided)),
            pct_supra_flat   = sprintf("%.1f%%", 100 * mean(supra_flat)))
print(as.data.frame(by_phen))
sink()
cat("\n=== McNemar (MM model-derived) ===\n")
print(tab)
print(mc_test)

# ----------------------------------------------------------------------------
# 7. Figure: 3-panel boxplot of exposures by dosing strategy, fill by phenotype
# ----------------------------------------------------------------------------
long <- exposures |>
  select(SUBJID, Phenotype, Dosing, AUC24, Cmax, C24) |>
  pivot_longer(c(AUC24, Cmax, C24),
               names_to = "Metric", values_to = "Value") |>
  mutate(Metric = factor(Metric,
                         levels = c("AUC24", "Cmax", "C24"),
                         labels = c("AUC[0-24]~'(mg%.%h/L)'",
                                    "C[max]~'('*mu*g/mL*')'",
                                    "C[24]~'('*mu*g/mL*')'")))

gg <- ggplot(long, aes(Dosing, Value, fill = Phenotype)) +
  geom_boxplot(outlier.shape = 1, outlier.alpha = .5) +
  facet_wrap(~ Metric, scales = "free_y",
             labeller = label_parsed) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = NULL,
       title = "Model-derived isoniazid exposures (MM popPK, dense IPRED simulation)",
       subtitle = "Per-subject AUC, Cmax, C24 from individual EBE parameters; n = 104 subjects, 2 occasions each") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9))

ggsave(file.path(OUT, "Figure_exposures_MM.pdf"), gg, width = 11, height = 4)
ggsave(file.path(OUT, "Figure_exposures_MM.png"), gg, width = 11, height = 4, dpi = 200)

cat("\nAll outputs in", OUT, "\n")

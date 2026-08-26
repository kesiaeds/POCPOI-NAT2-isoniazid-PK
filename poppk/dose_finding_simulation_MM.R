# ============================================================================
#  POCPOI dose-finding Monte Carlo simulation -- PRIMARY (MM popPK) model.
#  Reconstructed from briefing/dose_finding_MM_briefing.html (2026-06-14).
#  Same logical framework as scripts/dose_finding_simulation.py (linear popPK),
#  but the structural model is the 2-cmt + Michaelis-Menten ODE, solved per
#  subject with rxode2, using the parameters from fit_Km10.rds.
#
#  Outputs (results/):
#    mm_dose_sim_results.csv      per-phenotype x dose summary stats
#    mm_dose_sim_recommended.csv  recommended dose per phenotype
#  Run:  Rscript scripts/dose_finding_simulation_MM.R
# ============================================================================
suppressPackageStartupMessages({
  library(rxode2)
  library(dplyr)
  library(readr)
})

here <- tryCatch(dirname(normalizePath(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)[1]))),
          error = function(e) getwd())
proj <- normalizePath(file.path(here, ".."))
RES  <- file.path(proj, "results");  dir.create(RES, showWarnings = FALSE)

# ---------- MM parameters from the regenerated fit --------------------------
fit  <- readRDS(file.path(proj, "popPK_model", "fit_Km10.rds"))
est  <- fit$parFixedDf[, "Estimate"]; names(est) <- rownames(fit$parFixedDf)
g <- function(nm) unname(est[nm])
tKA   <- g("tKA");   tVmax <- g("tVmax"); tKm <- g("tKm")
tV2   <- g("tV2");   tQ    <- g("tQ");    tV3 <- g("tV3")
dVm   <- c(Slow = g("dVmdACEPROFILE1"), Intermediate = 0,
           Rapid = g("dVmdACEPROFILE3"))
omega <- sqrt(fit$omega["eta.vmax", "eta.vmax"])   # SD of eta on Vmax (log scale)

EXP_VMAX_WT <- 0.75   # allometric exponent on Vmax / Q (fixed)
EXP_V_WT    <- 1.0    # allometric exponent on V2 / V3   (fixed)
REF_WT      <- 70.0
SUPRA       <- 0.15   # ug/mL supratherapeutic C24 threshold
SUPRA_CAP   <- 15.0   # % allowed above SUPRA (ceiling-dose logic)

PHEN_ORDER <- c("Slow", "Intermediate", "Rapid")
PHEN_ID    <- c(Slow = 1, Intermediate = 2, Rapid = 3)

# ---------- reference target: Intermediate x flat 900 mg (observed NCA) ------
exps <- read_csv(file.path(proj, "data", "exposures_subject_occasion.csv"),
                 show_col_types = FALSE)
ref  <- exps %>% filter(Phenotype == "Intermediate", Dosing == "Flat 900 mg")
AUC_LO  <- quantile(ref$AUC24, .25, na.rm = TRUE); AUC_MED <- median(ref$AUC24, na.rm = TRUE)
AUC_HI  <- quantile(ref$AUC24, .75, na.rm = TRUE)
C24_LO  <- quantile(ref$C24, .25, na.rm = TRUE);  C24_MED <- median(ref$C24, na.rm = TRUE)
C24_HI  <- quantile(ref$C24, .75, na.rm = TRUE)
cat(sprintf("Reference (Inter x flat 900): AUC %.2f (IQR %.2f-%.2f) | C24 %.4f\n",
            AUC_MED, AUC_LO, AUC_HI, C24_MED))

# ---------- weight pool (PK cohort, per phenotype) --------------------------
nm <- read_csv(file.path(proj, "data", "pk_nonmem_combined.csv"),
               show_col_types = FALSE)
wt_pool <- nm %>% distinct(ID, .keep_all = TRUE) %>% select(ACE_PROFILE, WT)

# ---------- MM structural model (rxode2) ------------------------------------
mm <- rxode2({
  Cp = central / V2
  d/dt(depot)      = -Ka * depot
  d/dt(central)    =  Ka * depot - Vmax * Cp/(Km + Cp) - Q*(Cp - peripheral/V3)
  d/dt(peripheral) =  Q*(Cp - peripheral/V3)
})

samp <- sort(unique(c(0, seq(0.1, 4, 0.1), seq(4.5, 24, 0.5))))   # dense grid
i24  <- which(samp == 24)
N    <- 2000
DOSES <- seq(200, 1800, 50)
set.seed(20260601)

records <- list()
for (phen in PHEN_ORDER) {
  pool <- wt_pool$WT[wt_pool$ACE_PROFILE == PHEN_ID[phen]]
  wts  <- sample(pool, N, replace = TRUE)              # resampled weights
  eta  <- rnorm(N, 0, omega)                           # BSV on Vmax
  vmax_i <- tVmax * (wts/REF_WT)^EXP_VMAX_WT * exp(dVm[[phen]] + eta)
  v2_i   <- tV2   * (wts/REF_WT)^EXP_V_WT
  v3_i   <- tV3   * (wts/REF_WT)^EXP_V_WT
  q_i    <- tQ    * (wts/REF_WT)^EXP_VMAX_WT
  pars   <- data.frame(Ka = tKA, Vmax = vmax_i, Km = tKm,
                       V2 = v2_i, Q = q_i, V3 = v3_i)
  for (dose in DOSES) {
    ev  <- et(amt = dose, cmt = "depot") %>% et(samp)
    sol <- rxSolve(mm, pars, ev, returnType = "data.frame",
                   cores = 4, addDosing = FALSE)
    # rows ordered by sim.id then time -> reshape Cp to [time x subject]
    M    <- matrix(sol$Cp, nrow = length(samp))
    auc  <- colSums((head(M, -1) + tail(M, -1)) / 2 * diff(samp))   # trapezoid 0-24
    c24  <- M[i24, ]
    records[[length(records) + 1]] <- data.frame(
      Phenotype = phen, Dose_mg = as.integer(dose),
      AUC_med = median(auc), AUC_q1 = quantile(auc, .25), AUC_q3 = quantile(auc, .75),
      C24_med = median(c24), C24_q1 = quantile(c24, .25), C24_q3 = quantile(c24, .75),
      pct_in_auc_target  = mean(auc >= AUC_LO & auc <= AUC_HI) * 100,
      pct_in_c24_target  = mean(c24 >= C24_LO & c24 <= C24_HI) * 100,
      pct_in_dual_target = mean(auc >= AUC_LO & auc <= AUC_HI &
                                c24 >= C24_LO & c24 <= C24_HI) * 100,
      pct_supra = mean(c24 > SUPRA) * 100, row.names = NULL)
  }
  cat("  done:", phen, "\n")
}
sim <- bind_rows(records)
write_csv(sim, file.path(RES, "mm_dose_sim_results.csv"))

# ---------- recommended dose per phenotype ----------------------------------
pick <- function(phen, col, target) {
  sub <- sim[sim$Phenotype == phen, ]
  sub <- sub[order(abs(sub[[col]] - target), sub$pct_supra), ]
  sub[1, ]
}
pick_ceiling <- function(phen) {
  sub <- sim[sim$Phenotype == phen, ]; sub <- sub[order(sub$Dose_mg), ]
  ok  <- sub[sub$pct_supra <= SUPRA_CAP, ]
  if (nrow(ok)) ok[nrow(ok), ] else sub[1, ]
}
rec <- do.call(rbind, lapply(PHEN_ORDER, function(p) {
  a <- pick(p, "AUC_med", AUC_MED); c <- pick(p, "C24_med", C24_MED); k <- pick_ceiling(p)
  data.frame(Phenotype = p,
             AUC_match_dose_mg = a$Dose_mg, AUC_med = a$AUC_med, C24_med = a$C24_med,
             pct_in_auc_target = a$pct_in_auc_target, pct_supra = a$pct_supra,
             C24_match_dose_mg = c$Dose_mg, pct_in_c24_target = c$pct_in_c24_target,
             Ceiling_dose_mg = k$Dose_mg, pct_supra_at_ceiling = k$pct_supra,
             row.names = NULL)
}))
write_csv(rec, file.path(RES, "mm_dose_sim_recommended.csv"))

cat("\nAUC-match dose per phenotype:\n")
print(rec[, c("Phenotype", "AUC_match_dose_mg", "AUC_med",
              "pct_in_auc_target", "pct_supra")], row.names = FALSE)

# ---------- check against final tablet-count recommendations ----------------
cat("\nAt final recommended doses (600 / 900 / 1200 mg):\n")
final <- c(Slow = 600, Intermediate = 900, Rapid = 1200)
for (p in PHEN_ORDER) {
  r <- sim[sim$Phenotype == p & sim$Dose_mg == final[[p]], ]
  cat(sprintf("  %-12s %4d mg : AUC %.1f (IQR %.1f-%.1f), %%inIQR %.1f, %%supra %.1f\n",
              p, final[[p]], r$AUC_med, r$AUC_q1, r$AUC_q3,
              r$pct_in_auc_target, r$pct_supra))
}
cat("\nWrote results/mm_dose_sim_results.csv and mm_dose_sim_recommended.csv\n")

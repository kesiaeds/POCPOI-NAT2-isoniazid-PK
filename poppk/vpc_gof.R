# ============================================================================
#  POCPOI popPK model diagnostics: GOF panels (always) + manual VPC (rxode2).
#  Run AFTER pk_analysis_combined.R has produced output/fit_combined_focei.rds.
# ============================================================================
LIB <- Sys.getenv("R_LIBS_USER",
                  file.path("~", "R", "x86_64-pc-linux-gnu-library", "4.3"))
.libPaths(c(LIB, .libPaths()))

suppressPackageStartupMessages({
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

setwd(Sys.getenv("POCPOI_WORKDIR",
                 file.path("/scratch/users", Sys.getenv("USER"), "POCPOI")))
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

fit <- readRDS("output/fit_combined_focei.rds")
cat("Loaded fit. OFV =", fit$objf, "\n")

fit_df <- as.data.frame(fit)
if (!"TAD" %in% names(fit_df)) fit_df$TAD <- fit_df$TIME %% 168
fit_df <- dplyr::mutate(fit_df,
  Phenotype = factor(ACE_PROFILE, 1:3, c("Slow", "Intermediate", "Rapid")))

# ============================================================================
#  Goodness-of-fit panels (no tidyvpc dependency)
# ============================================================================
gof <- fit_df |> filter(EVID == 0, !is.na(DV))

p1 <- ggplot(gof, aes(PRED, DV)) +
  geom_point(alpha = .4) + geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = .6) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Population prediction (µg/mL)",
       y = "Observed conc. (µg/mL)", title = "DV vs PRED") + theme_bw()

p2 <- ggplot(gof, aes(IPRED, DV)) +
  geom_point(alpha = .4) + geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = .6) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Individual prediction (µg/mL)",
       y = "Observed conc. (µg/mL)", title = "DV vs IPRED") + theme_bw()

p3 <- ggplot(gof, aes(TAD, CWRES)) +
  geom_point(alpha = .4) + geom_hline(yintercept = 0) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = .6) +
  labs(x = "Time after dose (h)", y = "CWRES",
       title = "CWRES vs time after dose") + theme_bw()

p4 <- ggplot(gof, aes(PRED, CWRES)) +
  geom_point(alpha = .4) + geom_hline(yintercept = 0) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = .6) +
  scale_x_log10() +
  labs(x = "Population prediction (µg/mL)", y = "CWRES",
       title = "CWRES vs PRED") + theme_bw()

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  gg_gof <- (p1 | p2) / (p3 | p4)
  ggsave("output/figures/gof_panels.pdf", gg_gof, width = 10, height = 8)
} else {
  ggsave("output/figures/gof_dv_pred.pdf",   p1, width = 5, height = 4)
  ggsave("output/figures/gof_dv_ipred.pdf",  p2, width = 5, height = 4)
  ggsave("output/figures/gof_cwres_tad.pdf", p3, width = 5, height = 4)
  ggsave("output/figures/gof_cwres_pred.pdf",p4, width = 5, height = 4)
}
cat("GOF panels written.\n")

# ============================================================================
#  Visual predictive check (manual: rxode2 simulation + ggplot)
# ============================================================================
vpc_ok <- tryCatch({

  pf <- setNames(fit$parFixedDf$Estimate, rownames(fit$parFixedDf))
  om_cl_sd <- sqrt(as.matrix(fit$omega)["eta.cl", "eta.cl"])
  prop_sd  <- unname(pf["prop.err.p"])

  sim_mod <- rxode2::rxode2({
    Ka <- tka
    CL <- tcl * (WT/70)^0.75 *
            exp(d1 * (ACE_PROFILE == 1)) *
            exp(d3 * (ACE_PROFILE == 3)) *
            exp(eta_cl)
    V  <- tv2 * (WT/70)
    Q  <- tq  * (WT/70)^0.75
    V3 <- tv3 * (WT/70)
    cp <- linCmt()
  })

  # Build events from the source CSV (fit's as.data.frame drops AMT).
  # Reproduce the ID remapping pk_analysis_combined.R applied.
  dat_raw <- readr::read_csv("data/pk_nonmem_combined.csv",
                             show_col_types = FALSE) |>
    mutate(SUBJID = ID, ID = as.integer(factor(ID))) |>
    filter(WT > 0, !is.na(ACE_PROFILE))
  if (!"TAD" %in% names(dat_raw)) dat_raw$TAD <- dat_raw$TIME %% 168

  ev_obs <- dat_raw |> filter(EVID == 0) |>
    transmute(id = ID, time = TIME, evid = 0L, amt = 0,
              WT = WT, ACE_PROFILE = as.integer(ACE_PROFILE), TAD = TAD)
  ev_dose <- dat_raw |> filter(EVID == 1) |>
    transmute(id = ID, time = TIME, evid = 1L, amt = AMT,
              WT = WT, ACE_PROFILE = as.integer(ACE_PROFILE), TAD = NA_real_)
  ev_template <- bind_rows(ev_dose, ev_obs) |> arrange(id, time)

  uniq_ids <- unique(ev_template$id)
  nsim <- 200
  theta_pars <- c(
    tka = unname(pf["tKA"]),
    tcl = unname(pf["tCL"]),
    tv2 = unname(pf["tV2"]),
    tq  = unname(pf["tQ"]),
    tv3 = unname(pf["tV3"]),
    d1  = unname(pf["dCldACEPROFILE1"]),
    d3  = unname(pf["dCldACEPROFILE3"])
  )

  set.seed(42)
  cat("Simulating", nsim, "VPC replicates...\n")
  sim_list <- vector("list", nsim)
  for (k in seq_len(nsim)) {
    eta_draw <- data.frame(id = uniq_ids,
                           eta_cl = rnorm(length(uniq_ids), 0, om_cl_sd))
    ev_k <- ev_template |> left_join(eta_draw, by = "id")
    s <- rxode2::rxSolve(
           sim_mod, params = theta_pars, events = ev_k,
           keep = c("TAD", "ACE_PROFILE"),
           returnType = "data.frame")
    s <- s |> filter(!is.na(TAD)) |>
      mutate(rep = k,
             DV_sim = cp * (1 + rnorm(n(), 0, prop_sd)))
    sim_list[[k]] <- s
    if (k %% 50 == 0) cat("  ", k, "/", nsim, "\n", sep = "")
  }
  sim_all <- bind_rows(sim_list) |>
    mutate(Phenotype = factor(ACE_PROFILE, 1:3,
                              c("Slow", "Intermediate", "Rapid")))
  saveRDS(sim_all, "output/vpc_simulations.rds")

  tad_breaks <- c(-0.01, 0.5, 1.5, 5, 16, 24.5)
  tad_centers <- c(0, 1, 2, 8, 24)
  bin <- function(x) tad_centers[as.integer(cut(x, tad_breaks))]
  sim_all$tad_bin <- bin(sim_all$TAD)
  obs <- fit_df |> filter(EVID == 0) |> mutate(tad_bin = bin(TAD))

  sim_pct <- sim_all |>
    group_by(Phenotype, tad_bin, rep) |>
    summarise(p5  = quantile(DV_sim, .05, na.rm = TRUE),
              p50 = median  (DV_sim,      na.rm = TRUE),
              p95 = quantile(DV_sim, .95, na.rm = TRUE),
              .groups = "drop") |>
    group_by(Phenotype, tad_bin) |>
    summarise(p5_lo  = quantile(p5,  .025), p5_hi  = quantile(p5,  .975),
              p50_med = median(p50), p50_lo = quantile(p50, .025), p50_hi = quantile(p50, .975),
              p95_lo = quantile(p95, .025), p95_hi = quantile(p95, .975),
              .groups = "drop")

  obs_pct <- obs |> group_by(Phenotype, tad_bin) |>
    summarise(o5 = quantile(DV, .05, na.rm = TRUE),
              o50 = median(DV, na.rm = TRUE),
              o95 = quantile(DV, .95, na.rm = TRUE),
              .groups = "drop")

  gg_vpc <- ggplot() +
    geom_ribbon(data = sim_pct, aes(tad_bin, ymin = p5_lo,  ymax = p5_hi),
                fill = "steelblue", alpha = .3) +
    geom_ribbon(data = sim_pct, aes(tad_bin, ymin = p95_lo, ymax = p95_hi),
                fill = "steelblue", alpha = .3) +
    geom_ribbon(data = sim_pct, aes(tad_bin, ymin = p50_lo, ymax = p50_hi),
                fill = "firebrick", alpha = .35) +
    geom_line(data = sim_pct, aes(tad_bin, p50_med), color = "firebrick") +
    geom_line(data = obs_pct, aes(tad_bin, o50),  color = "black", linewidth = .7) +
    geom_line(data = obs_pct, aes(tad_bin, o5),   color = "black", linewidth = .5, linetype = "dashed") +
    geom_line(data = obs_pct, aes(tad_bin, o95),  color = "black", linewidth = .5, linetype = "dashed") +
    geom_point(data = obs, aes(TAD, DV), color = "grey40", alpha = .3, size = .9) +
    facet_wrap(~ Phenotype) +
    scale_y_log10() +
    labs(x = "Time after dose (h)",
         y = expression(INH~plasma~conc.~"("*mu*g/mL*")"),
         title = "Visual predictive check by NAT2 acetylator phenotype",
         subtitle = paste0("Black lines: observed 5/50/95 percentiles. ",
                           "Ribbons: 95% CI of simulated percentiles (n=", nsim, ").")) +
    theme_bw()
  ggsave("output/figures/vpc_by_phenotype.pdf", gg_vpc, width = 10, height = 4)
  cat("VPC saved to output/figures/vpc_by_phenotype.pdf\n")
  TRUE
}, error = function(e) {
  cat("VPC step failed: ", conditionMessage(e), "\n", sep = "")
  cat("GOF panels are still in output/figures/.\n")
  FALSE
})

cat("\nDONE. VPC produced:", vpc_ok, "\n")

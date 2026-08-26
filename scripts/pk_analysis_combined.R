# ============================================================================
#  POCPOI combined-cohort popPK analysis (Campo Grande + Manaus)
#  Run on Sherlock with: ml R/4.3 ; Rscript pk_analysis_combined.R
# ============================================================================
suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(tidyvpc)
})

setwd(Sys.getenv("POCPOI_WORKDIR",
                 "/scratch/users/kesiaeds/POCPOI"))   # adjust on Sherlock
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)

# ------ Load harmonized data ------------------------------------------------
dat <- read_csv("data/pk_nonmem_combined.csv", show_col_types = FALSE)

# nlmixr2 wants numeric IDs; map "CG_n"/"MA_n" -> integer
dat <- dat %>%
  mutate(SUBJID = ID,
         ID = as.integer(factor(ID))) %>%
  filter(WT > 0, !is.na(ACE_PROFILE))

# Sanity check
cat("\n--- Combined dataset summary ---\n")
cat("Subjects:", length(unique(dat$ID)),
    " | obs:",  sum(dat$EVID == 0),
    " | doses:", sum(dat$EVID == 1), "\n")
print(table(site = dat$SITE[!duplicated(dat$ID)],
            phen = dat$ACE_PROFILE[!duplicated(dat$ID)]))

# ============================================================================
#  Two-compartment popPK model (extends Sarkodie/CG model)
#  CL scaled allometrically; ACE_PROFILE shifts CL multiplicatively;
#  proportional residual error.
# ============================================================================
cmt2mod <- function() {
  ini({
    tKA  <-  3.18
    tCL  <- 30.0
    tV2  <- 81.9
    tQ   <-  1.7
    tV3  <- 16.5
    dCldACEPROFILE1 <- -0.55   # slow:   lower CL
    dCldACEPROFILE3 <-  0.30   # rapid:  higher CL
    dVdWT  <- fix(1.0)
    dCldWT <- fix(0.75)
    eta.cl ~ 0.5
    prop.err.p <- 0.30
  })
  model({
    Ka <- tKA
    CL <- tCL * (WT/70)^dCldWT *
            exp(dCldACEPROFILE1 * (ACE_PROFILE == 1)) *
            exp(dCldACEPROFILE3 * (ACE_PROFILE == 3)) *
            exp(eta.cl)
    V2 <- tV2 * (WT/70)^dVdWT
    Q  <- tQ  * (WT/70)^dCldWT
    V3 <- tV3 * (WT/70)^dVdWT
    linCmt() ~ prop(prop.err.p)
  })
}

# ------ Fit FOCEI ----------------------------------------------------------
fit <- nlmixr2(cmt2mod, dat, est = "focei",
               control = foceiControl(print = 5))
saveRDS(fit, "output/fit_combined_focei.rds")

# ------ Parameter table ----------------------------------------------------
ptab <- fit$parFixedDf
write.csv(ptab, "output/popPK_parameters.csv", row.names = TRUE)
print(ptab)

# ------ Empirical Bayes individual CL --------------------------------------
ind <- fit %>% as.data.frame() %>%
  group_by(ID) %>% slice(1) %>%
  select(ID, WT, ACE_PROFILE, eta.cl) %>%
  mutate(CL_i = ptab["tCL", "Estimate"] *
                (WT/70)^0.75 *
                exp(ptab["dCldACEPROFILE1", "Estimate"] * (ACE_PROFILE == 1)) *
                exp(ptab["dCldACEPROFILE3", "Estimate"] * (ACE_PROFILE == 3)) *
                exp(eta.cl))
write.csv(ind, "output/individual_CL.csv", row.names = FALSE)

# ------ NCA-style exposures (AUC0-24, Cmax, C24) for each subject-occasion --
exposures <- dat %>%
  filter(EVID == 0) %>%
  group_by(SUBJID, ID, OCC, DAY, ACE_PROFILE, WT) %>%
  summarise(
    Cmax  = max(DV, na.rm = TRUE),
    C24   = DV[TAD == 24][1],
    AUC24 = {
      tt <- TAD; cc <- DV
      ord <- order(tt); tt <- tt[ord]; cc <- cc[ord]
      sum(diff(tt) * (head(cc, -1) + tail(cc, -1)) / 2)
    },
    .groups = "drop"
  )
write.csv(exposures, "output/exposures.csv", row.names = FALSE)

# ------ Wilcoxon signed-rank: within-subject AUC flat vs guided -------------
paired <- exposures %>%
  select(SUBJID, OCC, AUC24) %>%
  pivot_wider(names_from = OCC, names_prefix = "AUC_OCC",
              values_from = AUC24) %>%
  drop_na()
wt <- wilcox.test(paired$AUC_OCC2, paired$AUC_OCC1, paired = TRUE)
cat("\nWilcoxon signed-rank (Day14 flat vs Day7 guided): p =",
    format.pval(wt$p.value), "\n")

# ------ McNemar on C24 > 0.15 µg/mL (supratherapeutic threshold) ------------
mc_dat <- exposures %>%
  mutate(supra = C24 > 0.15) %>%
  select(SUBJID, OCC, supra) %>%
  pivot_wider(names_from = OCC, names_prefix = "supra_OCC",
              values_from = supra) %>%
  drop_na()
tab <- table(guided = mc_dat$supra_OCC1, flat = mc_dat$supra_OCC2)
print(tab)
mc <- mcnemar.test(tab)
cat("McNemar (C24>0.15): p =", format.pval(mc$p.value), "\n")

# ------ Save tidy exposure summary by phenotype x occasion ------------------
summ <- exposures %>%
  mutate(Phenotype = factor(ACE_PROFILE, 1:3,
                            c("Slow", "Intermediate", "Rapid")),
         Dosing    = factor(OCC, 1:2,
                            c("Genotype-guided", "Flat 900 mg"))) %>%
  group_by(Dosing, Phenotype) %>%
  summarise(across(c(AUC24, Cmax, C24),
                   list(median = ~median(.x, na.rm = TRUE),
                        q1 = ~quantile(.x, .25, na.rm = TRUE),
                        q3 = ~quantile(.x, .75, na.rm = TRUE))),
            n = n(), .groups = "drop")
write.csv(summ, "output/exposure_summary.csv", row.names = FALSE)
print(summ)

# ------ Visual predictive check --------------------------------------------
vpc_obj <- vpcPlot(fit, n = 500, show = list(obs_dv = TRUE))
ggsave("output/figures/vpc.pdf", vpc_obj, width = 7, height = 5)

# ------ Profile plot: median +/- IQR by phenotype/day -----------------------
p <- dat %>% filter(EVID == 0) %>%
  mutate(Phenotype = factor(ACE_PROFILE, 1:3,
                            c("Slow", "Intermediate", "Rapid")),
         Dosing    = factor(OCC, 1:2,
                            c("Genotype-guided", "Flat 900 mg"))) %>%
  group_by(Phenotype, Dosing, TAD) %>%
  summarise(med = median(DV, na.rm = TRUE),
            q1  = quantile(DV, .25, na.rm = TRUE),
            q3  = quantile(DV, .75, na.rm = TRUE),
            .groups = "drop")
gg <- ggplot(p, aes(TAD, med, color = Dosing, fill = Dosing)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = .2, color = NA) +
  geom_line(linewidth = 1) + geom_point() +
  facet_wrap(~Phenotype) +
  scale_y_log10() +
  labs(x = "Time after dose (h)",
       y = expression(INH~plasma~conc.~"("*mu*g/mL*")")) +
  theme_bw()
ggsave("output/figures/profiles_by_phenotype.pdf", gg, width = 9, height = 4)

cat("\nDONE. Artifacts in output/.\n")

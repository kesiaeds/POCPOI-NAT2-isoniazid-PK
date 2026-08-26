# ============================================================================
#  POCPOI primary popPK model: 2-compartment + Michaelis-Menten elimination
#  Reconstructed from briefing/dose_finding_MM_briefing.html (2026-06-14).
#  Produces fit_Km10.rds  ("Km10" = Km initial estimate of 10 ug/mL).
#  Run locally (macOS, R 4.5):  Rscript pk_analysis_MM.R
# ============================================================================
suppressPackageStartupMessages({
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(readr)
})

# Resolve project root regardless of where the script is launched from
here  <- tryCatch(dirname(normalizePath(sub("--file=", "",
            grep("--file=", commandArgs(FALSE), value = TRUE)[1]))),
            error = function(e) getwd())
proj  <- normalizePath(file.path(here, ".."))
datfp <- file.path(proj, "data", "pk_nonmem_combined.csv")
outfp <- file.path(here, "fit_Km10.rds")

# ------ Load harmonized data (same pipeline as the linear model) ------------
dat <- read_csv(datfp, show_col_types = FALSE) %>%
  mutate(SUBJID = ID,
         ID = as.integer(factor(ID))) %>%
  filter(WT > 0, !is.na(ACE_PROFILE))
if ("CENS" %in% names(dat) && !"cens" %in% names(dat)) dat$cens <- dat$CENS

cat("\n--- Dataset ---\n")
cat("Subjects:", length(unique(dat$ID)),
    " | obs:", sum(dat$EVID == 0),
    " | doses:", sum(dat$EVID == 1), "\n")

# ============================================================================
#  2-compartment model, MM elimination from the central compartment.
#  NAT2 phenotype shifts Vmax (enzyme capacity), not Km. Allometric scaling on
#  Vmax/Q (0.75) and V2/V3 (1.0). Single exponential BSV on Vmax. Prop. error.
#  Units: amounts mg, conc ug/mL (= mg/L), time h  ->  Cp = central/V2.
# ============================================================================
mmmod <- function() {
  ini({
    tKA   <-  3.18
    tVmax <- 300.0      # mg/h, 70 kg, intermediate (starting value)
    tKm   <-  10.0      # ug/mL  <-- the "Km10" starting value
    tV2   <-  81.9
    tQ    <-   1.7
    tV3   <-  16.5
    dVmdACEPROFILE1 <- -0.55   # slow:  lower Vmax
    dVmdACEPROFILE3 <-  0.30   # rapid: higher Vmax
    dVdWT  <- fix(1.0)         # allometric exponent, volumes
    dVmdWT <- fix(0.75)        # allometric exponent, Vmax/Q
    eta.vmax ~ 0.1
    prop.err.p <- 0.30
  })
  model({
    Ka   <- tKA
    Vmax <- tVmax * (WT/70)^dVmdWT *
              exp(dVmdACEPROFILE1 * (ACE_PROFILE == 1)) *
              exp(dVmdACEPROFILE3 * (ACE_PROFILE == 3)) *
              exp(eta.vmax)
    Km   <- tKm
    V2   <- tV2 * (WT/70)^dVdWT
    Q    <- tQ  * (WT/70)^dVmdWT
    V3   <- tV3 * (WT/70)^dVdWT
    Cp   <- central / V2
    d/dt(depot)      = -Ka * depot
    d/dt(central)    =  Ka * depot - Vmax * Cp/(Km + Cp) -
                         Q * (Cp - peripheral/V3)
    d/dt(peripheral) =  Q * (Cp - peripheral/V3)
    Cp ~ prop(prop.err.p)
  })
}

# ------ Fit FOCEI -----------------------------------------------------------
fit <- nlmixr2(mmmod, dat, est = "focei",
               control = foceiControl(print = 5))

saveRDS(fit, outfp)
cat("\nSaved:", outfp, "\n")
print(fit$parFixedDf)
cat("\nOFV:", fit$objDf$OBJF, " | AIC:", fit$objDf$AIC, "\n")

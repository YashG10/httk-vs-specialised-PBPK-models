# =============================================================================
# httk_qspr_sweep_combined_fixed.R
# Full PBPK pipeline wrapped as a function and run once per QSPR source.
# Natively integrated with httk package data-loading functions.
# Includes explicit parameterization forcing and integrated observed PK data.
# =============================================================================

library(httk)
library(ggplot2)
library(reshape2)
library(openxlsx)
library(dplyr)

# =============================================================================
# GLOBAL DOSE, OBSERVED DATA & PLOT THEME DEFINITIONS
# =============================================================================
DEFAULT_DOSE_MG_KG <- 1.0

# -----------------------------------------------------------------------------
# Reference hepatic blood flow (well-stirred model upper bound for CL)
# -----------------------------------------------------------------------------
derive_QH_L_H_KG <- function(species, fallback) {
  tryCatch({
    bw <- as.numeric(physiology.data[physiology.data$Parameter == "Average BW", species])
    liver_flow_allo <- as.numeric(
      tissue.data$value[tolower(tissue.data$Tissue) == "liver" &
                          tissue.data$Species == species &
                          grepl("^Flow", tissue.data$variable)]
    )
    if (length(bw) != 1 || length(liver_flow_allo) != 1 ||
        is.na(bw) || is.na(liver_flow_allo)) stop("lookup returned no/ambiguous match")
    # mL/min/kg^0.75 -> L/h/kg BW: *60/1000 (mL/min -> L/h) * BW^-0.25 (deallometrize)
    liver_flow_allo * 0.06 * bw^(-0.25)
  }, error = function(e) {
    cat(sprintf("   [QH] Could not derive %s hepatic flow from httk tables (%s); using literature default %.2f L/h/kg\n",
                species, conditionMessage(e), fallback))
    fallback
  })
}

QH_HUMAN_L_H_KG <- derive_QH_L_H_KG("Human", fallback = 1.24)
QH_RAT_L_H_KG   <- derive_QH_L_H_KG("Rat",   fallback = 4.2)
cat(sprintf("Reference hepatic blood flow: Human = %.3f L/h/kg, Rat = %.3f L/h/kg\n",
            QH_HUMAN_L_H_KG, QH_RAT_L_H_KG))

# Doses scaled to mg/kg assuming a 70 kg body weight.
cpd_obs_dose <- list(
  Paracet    = 14.29,  
  Flutam     = 3.57,   
  
  # --- Integrated Modeled Exposure Scenarios ---
  Aldica     = 17.40,
  Androst    = 1.43,
  Aspirin    = 5.714,
  AZT        = 1.43,
  BaP        = 1.43e-6,
  BPS        = 0.10,
  BPA        = 0.10,
  Caffei     = 4.29,
  Carbamaz   = 11.43,
  Cetiri     = 0.143,
  CPF        = 2.00,
  TCPy       = 2.00,
  Cotinine   = 0.286,
  Acetam     = 0.0214,
  Dexame     = 0.0571,
  Dextrom    = 0.857,
  `24D`      = 5.00,
  Diphenhy   = 0.714,
  Erythro    = 14.29,
  Famoti     = 0.571,
  Flucona    = 2.14,
  Genist     = 0.714,
  Haloper    = 0.50,
  Ibupro     = 5.71,
  Lorata     = 0.143,
  Lovasta    = 0.571,
  MEHP       = 0.693,
  `5OH-MEHP` = 0.693,
  Naprox     = 7.14,
  Nicotine   = 0.0857,
  Nitro_f    = 1.43,
  Omepra     = 0.143,
  Phenyt     = 1.43,
  Prednis    = 0.286,
  Progest    = 0.286,
  PrP        = 2.50,
  Naphtha    = 0.020,
  Pyrene     = 0.030,
  Pyrimet_am = 0.714,
  Raloxi     = 0.857,
  Rifamp     = 6.43,
  Riluzo     = 0.714,
  Simvast    = 0.571,
  Sulfasal   = 28.57,
  Tamox      = 0.286,
  Tizoxanide = 7.14,
  Triamc     = 0.0714,
  Triclop    = 0.50,
  Triclo     = 0.0571,
  Triclocarban = 0.0010,
  UV327      = 0.30,
  Warfarin   = 0.357,
  Zamifen    = 0.30,
  `3-PBA`    = 0.10
)

# Added observed data, raw units converted to uM, with integrated specific references
obs_parent <- list(
  Aspirin = list(
    Human = data.frame(
      time_h      = c(0.167, 0.25, 0.333, 0.5, 0.667, 0.75, 1.0, 1.5, 2.0, 3.0),
      Cplasma_uM  = c(2.5, 4.5, 8.0, 14.0, 9.5, 7.0, 4.5, 2.2, 1.2, 0.45),
      sd_uM       = NA,
      dose_mg_kg  = 5.714,
      linear_pk   = TRUE,
      matrix      = "Plasma — Acetylsalicylic acid (400 mg oral dose; 12 healthy subjects)",
      stringsAsFactors = FALSE)
  ),
  Cetiri = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 24.0),
      Cplasma_uM  = c(0, 40, 100, 225, 250, 200, 170, 140, 125, 108, 93, 70, 22) / 388.89,
      sd_uM       = NA,
      dose_mg_kg  = 0.143,
      linear_pk   = TRUE,
      matrix      = "Plasma — Cetirizine (10 mg Zyrtecset oral dose; 12 healthy volunteers)",
      stringsAsFactors = FALSE)
  ),
  Haloper = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 12.0, 18.0, 24.0, 36.0, 44.0),
      Cplasma_uM  = c(4.4, 16.0, 34.0, 38.0, 34.0, 28.0, 22.0, 19.0, 14.0, 9.5, 10.0, 7.2, 4.2, 3.8) / 375.86,
      sd_uM       = NA,
      dose_mg_kg  = 0.50,
      linear_pk   = TRUE,
      matrix      = "Serum — Haloperidol (0.5 mg/kg oral dose; Subject 2)",
      stringsAsFactors = FALSE)
  ),
  Lorata = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 36.0, 48.0, 72.0, 96.0),
      Cplasma_uM  = c(0.1, 1.5, 2.7, 2.6, 2.2, 1.75, 1.35, 1.1, 0.8, 0.65, 0.4, 0.2, 0.1, 0.05, 0.04, 0.02, 0.0) / 382.88,
      sd_uM       = NA,
      dose_mg_kg  = 0.143,
      linear_pk   = TRUE,
      matrix      = "Plasma — Loratadine (10 mg oral dose)",
      stringsAsFactors = FALSE)
  ),
  MEHP = list(
    Human = data.frame(
      time_h      = c(2.0, 4.0, 6.5, 8.3),
      Cplasma_uM  = c(4.95, 0.57, 0.29, 0.15),
      sd_uM       = NA,
      dose_mg_kg  = 0.693,
      linear_pk   = TRUE,
      matrix      = "Plasma — MEHP (48.5 mg oral dose)",
      stringsAsFactors = FALSE)
  ),
  `5OH-MEHP` = list(
    Human = data.frame(
      time_h      = c(2.0, 4.0, 6.5, 8.3),
      Cplasma_uM  = c(0.20, 0.14, 0.05, 0.03),
      sd_uM       = NA,
      dose_mg_kg  = 0.693,
      linear_pk   = TRUE,
      matrix      = "Plasma — 5OH-MEHP (48.5 mg oral dose)",
      stringsAsFactors = FALSE)
  ),
  Naprox = list(
    Human = data.frame(
      time_h      = c(2.5, 5.0, 10.0, 15.0, 24.0, 28.0, 38.0, 48.0, 72.0),
      Cplasma_uM  = c(55.0, 40.0, 28.0, 20.0, 13.0, 11.0, 8.0, 5.5, 3.0) * 1000 / 230.26,
      sd_uM       = NA,
      dose_mg_kg  = 7.14,
      linear_pk   = TRUE,
      matrix      = "Plasma — Naproxen (500 mg oral dose; Subject A representative)",
      stringsAsFactors = FALSE)
  ),
  Prednis = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 12.0),
      Cplasma_uM  = c(0, 210, 370, 635.2, 625, 490, 390, 315, 250, 180, 120, 100, 50) / 360.44,
      sd_uM       = NA,
      dose_mg_kg  = 0.286,
      linear_pk   = TRUE,
      matrix      = "Serum — Prednisolone (20 mg oral dose; Reference formulation)",
      stringsAsFactors = FALSE)
  ),
  BPA = list(
    Human = data.frame(
      time_h      = c(0,    0.5,    1.0,    1.33,   2.0,   3.0,   4.0,   6.0,   8.0),
      Cplasma_uM  = c(0, 0.00120, 0.00250, 0.00340, 0.00280, 0.00180, 0.00090, 0.00020, 0.00005),
      sd_uM       = c(0, 0.00028, 0.00058, 0.00082, 0.00066, 0.00043, 0.00022, 0.00006,    NA),
      dose_mg_kg  = 0.0714,
      linear_pk   = TRUE,
      matrix      = "Serum — unconjugated BPA; 5 mg oral dose (~0.0714 mg/kg); 5 male volunteers",
      stringsAsFactors = FALSE)
  ),
  Paracet = list(
    Human = data.frame(
      time_h      = c(0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
      # MW = 151.16 g/mol; conversion from ug/mL to uM: (ug/mL) * 1000 / 151.16
      Cplasma_uM  = c(12.0, 14.0, 14.0, 7.5, 5.0, 3.8, 2.8, 2.2) * 1000 / 151.16,
      sd_uM       = NA,
      dose_mg_kg  = 14.29,
      linear_pk   = TRUE,
      matrix      = "Plasma — Paracetamol (1000 mg)",
      stringsAsFactors = FALSE)
  ),
  Riluzo = list(
    Human = data.frame(
      time_h      = c(0.0, 0.33, 0.67, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0, 24.0, 32.0, 36.0, 48.0, 60.0, 72.0),
      Cplasma_uM  = c(0.1879, 0.3757, 0.7173, 0.7344, 0.6063, 0.5764, 0.3928, 0.2989, 0.2476, 0.1964, 0.1623, 0.1281, 0.0939, 0.0683, 0.0598, 0.0512, 0.0384, 0.0342),
      sd_uM       = NA,
      dose_mg_kg  = 0.714,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Nitro_f = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 12.0, 16.0),
      Cplasma_uM  = c(0, 10, 50, 150, 260, 290, 320, 450, 430, 370, 340, 300, 240, 170, 80, 40, 10) / 238.16,
      sd_uM       = NA,
      dose_mg_kg  = 1.43,
      linear_pk   = TRUE,
      matrix      = "Plasma — Nitrofurantoin (reference)",
      stringsAsFactors = FALSE)
  ),
  Flutam = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0),
      Cplasma_uM  = c(0.0, 3.9, 9.6, 18.8, 20.2, 13.2, 8.5, 2.4) / 276.29,
      sd_uM       = NA,
      dose_mg_kg  = 3.57,
      linear_pk   = TRUE,
      matrix      = "Plasma — Flutamide",
      stringsAsFactors = FALSE)
  ),
  Dexame = list(
    Human = data.frame(
      time_h      = c(1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0),
      Cplasma_uM  = c(17.0, 26.0, 24.0, 27.0, 24.0, 16.0, 11.0, 4.5, 1.0) / 392.46,
      sd_uM       = NA,
      dose_mg_kg  = 0.057,
      linear_pk   = TRUE,
      matrix      = "Plasma — Dexamethasone (4 mg)",
      stringsAsFactors = FALSE)
  ),
  BPS = list(
    Human = data.frame(
      time_h      = c(0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 24.0, 48.0),
      Cplasma_uM  = c(0.05, 60.0, 120.0, 130.0, 80.0, 60.0, 45.0, 30.0, 25.0, 20.0, 15.0, 10.0, 9.0, 2.5, 0.4) / 250.27,
      sd_uM       = NA,
      dose_mg_kg  = 0.100,
      linear_pk   = TRUE,
      matrix      = "Plasma — BPS-d8",
      stringsAsFactors = FALSE)
  ),
  `3-PBA` = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 1.5, 4.0, 6.0, 8.0, 10.0, 24.0, 48.0, 72.0),
      Cplasma_uM  = c(0.00085, 0.6, 1.1, 1.45, 0.92, 1.0, 0.95, 0.78, 0.22, 0.016, 0.0035),
      sd_uM       = NA,
      dose_mg_kg  = 0.1,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Naphtha = list(
    Human = data.frame(
      time_h      = c(0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0),
      Cplasma_uM  = c(0.016, 0.0115, 0.0072, 0.0053, 0.0032, 0.0019, 0.00115, 0.0006, 0.00027, 5e-05),
      sd_uM       = NA,
      dose_mg_kg  = 0.02,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  ),
  PhE = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 1.1, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 48.0),
      Cplasma_uM  = c(0.0145, 0.1918, 0.2027, 0.1144, 0.0362, 0.0253, 0.0261, 0.021, 0.0195, 0.0232),
      sd_uM       = NA,
      dose_mg_kg  = 5.0,
      linear_pk   = TRUE,
      matrix      = "Whole Blood",
      stringsAsFactors = FALSE)
  ),
  PrP = list(
    Human = data.frame(
      time_h      = c(0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0),
      Cplasma_uM  = c(0.3, 0.1, 0.028, 0.013, 0.0055, 0.004, 0.0045),
      sd_uM       = NA,
      dose_mg_kg  = 2.5,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  ),
  Progest = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0),
      Cplasma_uM  = c(0.00115, 0.009, 0.0046, 0.00305, 0.0026, 0.00175, 0.0013, 0.00125),
      sd_uM       = NA,
      dose_mg_kg  = 0.286,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  ),
  Pyrene = list(
    Human = data.frame(
      time_h      = c(0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 6.0),
      Cplasma_uM  = c(0.00074, 0.00145, 0.00205, 0.00095, 0.00065, 0.00038, 0.00028, 0.000115),
      sd_uM       = NA,
      dose_mg_kg  = 0.03,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  ),
  Pyrimet_am = list(
    Human = data.frame(
      time_h      = c(0.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0, 72.0, 120.0, 168.0, 336.0, 504.0),
      Cplasma_uM  = c(0.0, 1.2062, 1.6003, 1.347, 1.3067, 1.1861, 1.0856, 0.9047, 0.7438, 0.5549, 0.4222, 0.1367, 0.0402),
      sd_uM       = NA,
      dose_mg_kg  = 0.714,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Raloxi = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 12.0, 24.0, 48.0, 72.0, 96.0),
      Cplasma_uM  = c(0.0, 0.000222, 0.000247, 0.000268, 0.000285, 0.000291, 0.000591, 0.000532, 0.000492, 0.000458, 0.000443, 0.000458, 0.000454, 0.000429, 0.000401, 0.000211, 9.9e-05, 4e-05),
      sd_uM       = NA,
      dose_mg_kg  = 0.857,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Rifamp = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 12.0, 24.0),
      Cplasma_uM  = c(1.944, 6.197, 13.974, 6.683, 3.767, 2.066, 1.094, 0.061),
      sd_uM       = NA,
      dose_mg_kg  = 6.43,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Simvast = list(
    Human = data.frame(
      time_h      = c(0.0, 0.33, 0.67, 1.0, 1.33, 1.67, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 24.0),
      Cplasma_uM  = c(0.0, 0.00108, 0.00526, 0.008, 0.0092, 0.0096, 0.00884, 0.00776, 0.00848, 0.008, 0.00681, 0.00645, 0.00478, 0.00311, 0.00239, 0.00179, 0.00053),
      sd_uM       = NA,
      dose_mg_kg  = 0.571,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Sulfasal = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 9.0, 12.0, 24.0),
      Cplasma_uM  = c(0.0, 3.0121, 8.7854, 23.093, 28.8662, 32.8824, 26.3561, 15.0606, 9.5384, 3.0121),
      sd_uM       = NA,
      dose_mg_kg  = 28.57,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Tamox = list(
    Human = data.frame(
      time_h      = c(1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0),
      Cplasma_uM  = c(0.03768, 0.06595, 0.08613, 0.10229, 0.09286, 0.07537, 0.05922, 0.03095),
      sd_uM       = NA,
      dose_mg_kg  = 0.286,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Tizoxanide = list(
    Human = data.frame(
      time_h      = c(0.0, 0.33, 0.67, 1.0, 1.33, 1.67, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0),
      Cplasma_uM  = c(0.0, 0.6786, 5.6173, 11.9133, 16.0603, 19.0009, 21.6023, 25.2969, 22.8087, 16.8143, 6.3713, 2.6013, 0.9802, 0.377),
      sd_uM       = NA,
      dose_mg_kg  = 7.14,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Triamc = list(
    Human = data.frame(
      time_h      = c(0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0),
      Cplasma_uM  = c(0.0, 0.003797, 0.010817, 0.01565, 0.016571, 0.018412, 0.016571, 0.013809, 0.010587, 0.007825, 0.004258, 0.001956, 0.001105, 0.000552),
      sd_uM       = NA,
      dose_mg_kg  = 0.0714,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Triclo = list(
    Human = data.frame(
      time_h      = c(1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0),
      Cplasma_uM  = c(0.58714, 0.449, 0.25903, 0.16578, 0.12088, 0.08289, 0.04317, 0.01312),
      sd_uM       = NA,
      dose_mg_kg  = 0.0571,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Triclocarban = list(
    Human = data.frame(
      time_h      = c(0.0, 0.67, 2.0, 3.0, 5.5, 9.5, 24.0, 48.0),
      Cplasma_uM  = c(0.285, 0.3, 0.31, 0.53, 0.325, 0.11, 0.048, 0.022),
      sd_uM       = NA,
      dose_mg_kg  = 0.001,
      linear_pk   = TRUE,
      matrix      = "Whole Blood",
      stringsAsFactors = FALSE)
  ),
  Triclop = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0),
      Cplasma_uM  = c(4.8739, 5.8486, 6.0436, 4.484, 3.1193, 2.1445, 1.4817, 0.7798, 0.1872),
      sd_uM       = NA,
      dose_mg_kg  = 0.5,
      linear_pk   = TRUE,
      matrix      = "Whole Blood",
      stringsAsFactors = FALSE)
  ),
  UV327 = list(
    Human = data.frame(
      time_h      = c(0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 24.0, 34.0, 48.0, 72.0),
      Cplasma_uM  = c(0.0, 0.91092, 1.52565, 1.76595, 1.65698, 1.29373, 0.28781, 0.21516, 0.11456, 0.07544),
      sd_uM       = NA,
      dose_mg_kg  = 0.3,
      linear_pk   = TRUE,
      matrix      = "Whole Blood",
      stringsAsFactors = FALSE)
  ),
  Warfarin = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 2.0, 4.0, 7.0, 9.0, 11.0, 24.0, 36.0, 48.0, 60.0, 72.0, 96.0, 120.0, 144.0, 168.0),
      Cplasma_uM  = c(9.0812, 8.2704, 7.4595, 6.8109, 5.8379, 5.6757, 6.0001, 4.3784, 3.7298, 2.5298, 2.1081, 1.6541, 1.0703, 0.6324, 0.4541, 0.2854),
      sd_uM       = NA,
      dose_mg_kg  = 0.357,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Zamifen = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 28.0, 52.0, 74.0, 100.0, 148.0, 192.0),
      Cplasma_uM  = c(2.328, 1.9788, 1.6762, 1.3503, 1.1175, 0.7682, 0.5122, 0.3492, 0.2561, 0.1816, 0.1769, 0.0978, 0.0559, 0.0279, 0.0116, 0.0049),
      sd_uM       = NA,
      dose_mg_kg  = 0.3,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `Aldica` = list(
    Human = data.frame(
      time_h      = c(0.0, 6.0, 12.0, 17.5, 21.5, 44.0, 78.0),
      Cplasma_uM  = c(11.0, 3.4, 2.6, 1.8, 1.1, 0.4, 0.15),
      sd_uM       = NA,
      dose_mg_kg  = 17.4,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `Androst` = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 24.0),
      Cplasma_uM  = c(0.008, 0.052, 0.136, 0.157, 0.164, 0.158, 0.128, 0.11, 0.106, 0.058, 0.041, 0.032, 0.026, 0.01),
      sd_uM       = NA,
      dose_mg_kg  = 1.43,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `AZT` = list(
    Human = data.frame(
      time_h      = c(0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 1.92, 3.08, 4.17),
      Cplasma_uM  = c(0.17, 0.18, 4.04, 1.8, 1.05, 0.79, 0.49, 0.28, 0.12),
      sd_uM       = NA,
      dose_mg_kg  = 1.43,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  ),
  `BaP` = list(
    Human = data.frame(
      time_h      = c(0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 8.0, 24.0, 48.0),
      Cplasma_uM  = c(0.0, 1.6e-09, 9.99e-08, 9.11e-08, 5.19e-08, 1.86e-08, 2.58e-08, 6.34e-09, 5.15e-09, 2.77e-09, 1.19e-09),
      sd_uM       = NA,
      dose_mg_kg  = 1.43e-06,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `Caffei` = list(
    Human = data.frame(
      time_h      = c(0.6, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 24.0),
      Cplasma_uM  = c(38.9, 37.8, 37.6, 38.1, 36.6, 31.9, 29.6, 25.5, 18.8, 13.4, 10.6, 2.57),
      sd_uM       = NA,
      dose_mg_kg  = 4.29,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `Carbamaz` = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 7.0, 23.5, 32.0, 48.0, 57.0, 72.5, 81.0, 98.0, 104.0),
      Cplasma_uM  = c(0.0, 8.25, 32.38, 41.05, 22.43, 20.1, 19.68, 15.45, 12.91, 10.79, 7.19),
      sd_uM       = NA,
      dose_mg_kg  = 11.43,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `CPF` = list(
    Human = data.frame(
      time_h      = c(3.0, 6.5, 10.5),
      Cplasma_uM  = c(0.012, 0.048, 0.007),
      sd_uM       = NA,
      dose_mg_kg  = 2.0,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `TCPy` = list(
    Human = data.frame(
      time_h      = c(2.0, 3.5, 6.5, 10.5, 25.0, 37.0, 47.0),
      Cplasma_uM  = c(0.38, 0.65, 8.4, 8.0, 6.5, 5.5, 5.0),
      sd_uM       = NA,
      dose_mg_kg  = 2.0,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `Cotinine` = list(
    Human = data.frame(
      time_h      = c(0.0, 0.33, 0.67, 1.0, 2.0, 4.0, 6.0, 8.0, 24.0, 30.0, 48.0, 54.0, 72.0, 81.0, 96.0),
      Cplasma_uM  = c(0.0, 0.34, 1.16, 1.69, 1.66, 1.45, 1.18, 1.01, 0.41, 0.31, 0.15, 0.11, 0.062, 0.045, 0.023),
      sd_uM       = NA,
      dose_mg_kg  = 0.286,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `Acetam` = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.15, 1.65, 2.0, 2.6, 3.05, 4.05, 5.0, 6.0, 8.0, 9.85, 24.0),
      Cplasma_uM  = c(0.0, 0.0687, 0.055, 0.0494, 0.0411, 0.0341, 0.0276, 0.0238, 0.0177, 0.013, 0.0072, 0.0043, 2e-04),
      sd_uM       = NA,
      dose_mg_kg  = 0.0214,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Dextrom = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 1.5, 2.5, 3.5, 6.0, 8.0, 10.0, 24.0, 48.0),
      Cplasma_uM  = c(0.0081, 0.0265, 0.0424, 0.0361, 0.031, 0.0203, 0.0184, 0.014, 0.0092, 0.0033),
      sd_uM       = NA,
      dose_mg_kg  = 0.857,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  `24D` = list(
    Human = data.frame(
      time_h      = c(0.5, 1.0, 2.0, 4.0, 7.0, 12.0, 24.0, 48.0, 72.0, 96.0, 120.0, 144.0),
      Cplasma_uM  = c(13.6, 41.6, 83.7, 126.7, 115.4, 88.2, 22.6, 2.94, 0.72, 0.25, 0.095, 0.013),
      sd_uM       = NA,
      dose_mg_kg  = 5.0,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Diphenhy = list(
    Human = data.frame(
      time_h      = c(0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0),
      Cplasma_uM  = c(0.043, 0.045, 0.074, 0.135, 0.137, 0.174, 0.155, 0.133, 0.11, 0.092, 0.074, 0.065, 0.053, 0.045),
      sd_uM       = NA,
      dose_mg_kg  = 0.714,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Erythro = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 9.0, 11.0),
      Cplasma_uM  = c(0.0, 0.57, 5.25, 8.56, 6.95, 5.93, 5.25, 4.91, 3.88, 2.73, 2.34, 1.5, 0.82),
      sd_uM       = NA,
      dose_mg_kg  = 14.29,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  ),
  Famoti = list(
    Human = data.frame(
      time_h      = c(0.5, 1.3, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0),
      Cplasma_uM  = c(0.058, 0.166, 0.201, 0.181, 0.154, 0.099, 0.07, 0.046, 0.029),
      sd_uM       = NA,
      dose_mg_kg  = 0.571,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Flucona = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0),
      Cplasma_uM  = c(0.0, 5.39, 7.31, 8.29, 8.65, 8.78, 8.26, 7.64, 7.08, 6.37, 4.93, 2.84),
      sd_uM       = NA,
      dose_mg_kg  = 2.14,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Genist = list(
    Human = data.frame(
      time_h      = c(0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 12.0, 16.0),
      Cplasma_uM  = c(0.0, 0.003, 0.0067, 0.0093, 0.0107, 0.0126, 0.0152, 0.0181, 0.0414, 0.0155, 0.0141, 0.0144, 0.0118, 0.0178),
      sd_uM       = NA,
      dose_mg_kg  = 0.714,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Ibupro = list(
    Human = data.frame(
      time_h      = c(0.0, 0.17, 0.33, 0.5, 0.67, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0),
      Cplasma_uM  = c(0.0, 33.9, 80.0, 123.6, 147.9, 129.9, 101.8, 84.8, 60.1, 40.7, 18.4, 5.3),
      sd_uM       = NA,
      dose_mg_kg  = 5.71,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Lovasta = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0),
      Cplasma_uM  = c(0.0, 0.0059, 0.0099, 0.0111, 0.0067, 0.0035, 0.003, 0.002, 0.001, 5e-04),
      sd_uM       = NA,
      dose_mg_kg  = 0.571,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Nicotine = list(
    Human = data.frame(
      time_h      = c(0.0, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0),
      Cplasma_uM  = c(0.005, 0.029, 0.044, 0.059, 0.085, 0.064, 0.054, 0.032, 0.022, 0.017, 0.014, 0.009),
      sd_uM       = NA,
      dose_mg_kg  = 0.0857,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Omepra = list(
    Human = data.frame(
      time_h      = c(0.083, 0.167, 0.2, 0.267, 0.333, 0.417, 0.5, 0.667, 0.833, 1.0, 1.5, 2.0),
      Cplasma_uM  = c(0.42, 0.85, 1.15, 1.05, 0.92, 0.42, 0.35, 0.22, 0.16, 0.135, 0.072, 0.058),
      sd_uM       = NA,
      dose_mg_kg  = 0.143,
      linear_pk   = TRUE,
      matrix      = "Plasma",
      stringsAsFactors = FALSE)
  ),
  Phenyt = list(
    Human = data.frame(
      time_h      = c(0.0, 1.0, 1.5, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0, 12.0, 24.0, 32.0, 48.0, 56.0, 78.0),
      Cplasma_uM  = c(0.0, 3.69, 4.8, 6.2, 5.07, 4.12, 4.4, 4.2, 4.16, 4.44, 4.92, 1.33, 0.91, 0.2, 0.1, 0.04),
      sd_uM       = NA,
      dose_mg_kg  = 1.43,
      linear_pk   = TRUE,
      matrix      = "Serum",
      stringsAsFactors = FALSE)
  )
)

# Global plot theme
theme_pbpk <- function(bs = 10) {
  theme_bw(base_size = bs) +
    theme(plot.title    = element_text(face = "bold", size = bs, hjust = 0),
          plot.subtitle = element_text(colour = "grey45", size = bs - 2),
          plot.caption  = element_text(colour = "grey55", size = bs - 3,
                                       hjust = 0, margin = margin(t = 4)),
          plot.margin   = margin(6, 10, 4, 6),
          axis.title    = element_text(face = "bold", size = bs - 1),
          axis.ticks    = element_line(colour = "grey70"),
          panel.border  = element_rect(colour = "grey75", fill = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          legend.position   = "bottom",
          legend.text       = element_text(size = 8),
          legend.key.width  = unit(1.8, "lines"),
          legend.background = element_blank())
}

run_pipeline_for_source <- function(QSPR_SOURCE, QSPR_OUTPUT_SUFFIX = QSPR_SOURCE) {
  
  # =============================================================================
  # Reset httk database natively
  # =============================================================================
  reset_httk()
  
  if (QSPR_SOURCE == "Sipes2017") suppressMessages(load_sipes2017(overwrite = TRUE))
  if (QSPR_SOURCE == "Pradeep2020") suppressMessages(load_pradeep2020(overwrite = TRUE))
  if (QSPR_SOURCE == "Dawson2021") suppressMessages(load_dawson2021(overwrite = TRUE))
  
  PREG_T0_H <- 91
  
  wb <- createWorkbook()
  
  compounds <- list(
    list(abbrev="BBP", cas="85-68-7", name="Benzyl butyl phthalate"),
    list(abbrev="BPA", cas="80-05-7", name="Bisphenol A"),
    list(abbrev="BPB", cas="77-40-7", name="Bisphenol B"),
    list(abbrev="BPS", cas="80-09-1", name="Bisphenol S"),
    list(abbrev="CPF", cas="2921-88-2", name="Chlorpyrifos"),
    list(abbrev="Clothi", cas="210880-92-5", name="Clothianidin"),
    list(abbrev="Cyflut", cas="68359-37-5", name="Cyfluthrin"),
    list(abbrev="Cyprod", cas="121552-61-2", name="Cyprodinil"),
    list(abbrev="DBP", cas="84-74-2", name="Dibutyl phthalate"),
    list(abbrev="DEHP", cas="117-81-7", name="Bis(2-ethylhexyl) phthalate"),
    list(abbrev="DcHP", cas="84-61-7", name="Dicyclohexyl phthalate"),
    list(abbrev="Dexame", cas="50-02-2", name="Dexamethasone"),
    list(abbrev="DiBP", cas="84-69-5", name="Diisobutyl phthalate"),
    list(abbrev="Dichl", cas="15165-67-0", name="Dichlorprop-p"),
    list(abbrev="Dienes", cas="84-17-3", name="Dienestrol"),
    list(abbrev="Dimeth", cas="110488-70-5", name="Dimethomorph"),
    list(abbrev="Flutam", cas="13311-84-7", name="Flutamide"),
    list(abbrev="Imazal", cas="35554-44-0", name="Imazalil"),
    list(abbrev="Linuro", cas="330-55-2", name="Linuron"),
    list(abbrev="Mecop", cas="16484-77-8", name="Mecoprop-p"),
    list(abbrev="Nitro_f", cas="67-20-9", name="Nitrofurantoin"),
    list(abbrev="PCB153", cas="35065-27-1", name="2,2',4,4',5,5'-Hexachlorobiphenyl"),
    list(abbrev="Paracet", cas="103-90-2", name="Acetaminophen"),
    list(abbrev="Permet", cas="52645-53-1", name="Permethrin"),
    list(abbrev="Riluzo", cas="1744-22-5", name="Riluzole"),
    list(abbrev="TCEP", cas="115-96-8", name="Tris(2-chloroethyl) phosphate"),
    list(abbrev="TCIPP", cas="13674-84-5", name="Tris(chloropropyl) phosphate"),
    list(abbrev="TDCIPP", cas="13674-87-8", name="Tris(1,3-dichloro-2-propyl) phosphate"),
    list(abbrev="Tebuco", cas="107534-96-3", name="Tebuconazole"),
    list(abbrev="Thiacl", cas="111988-49-9", name="Thiacloprid"),
    list(abbrev="Vinclo", cas="50471-44-8", name="Vinclozolin"),
    list(abbrev="Azoxy", cas="131860-33-8", name="Azoxystrobin"),
    list(abbrev="Glyph", cas="1071-83-6", name="Glyphosate"),
    list(abbrev="Thiamet", cas="153719-23-4", name="Thiamethoxam"),
    list(abbrev="Chloran", cas="500008-45-7", name="Chlorantraniliprole"),
    list(abbrev="Atraz", cas="1912-24-9", name="Atrazine"),
    list(abbrev="Glufos", cas="51276-47-2", name="Glufosinate"),
    list(abbrev="Prothio", cas="178928-70-6", name="Prothioconazole"),
    list(abbrev="Metola", cas="51218-45-2", name="Metolachlor"),
    list(abbrev="Cadmium", cas="7440-43-9", name="Cadmium"),
    list(abbrev="Mercury", cas="7439-97-6", name="Mercury"),
    list(abbrev="Arsenic", cas="7440-38-2", name="Arsenic"),
    list(abbrev="Lead", cas="7439-92-1", name="Lead"),
    list(abbrev="Ibupro", cas="15687-27-1", name="Ibuprofen"),
    list(abbrev="Aspirin", cas="50-78-2", name="Aspirin"),
    list(abbrev="Naprox", cas="22204-53-1", name="Naproxen"),
    list(abbrev="Diphenhy", cas="58-73-1", name="Diphenhydramine"),
    list(abbrev="Omepra", cas="73590-58-6", name="Omeprazole"),
    list(abbrev="Cetiri", cas="83881-51-0", name="Cetirizine"),
    list(abbrev="Lorata", cas="79794-75-5", name="Loratadine"),
    list(abbrev="Famoti", cas="76824-35-6", name="Famotidine"),
    list(abbrev="Dextrom", cas="125-71-3", name="Dextromethorphan"),
    list(abbrev="DDT", cas="50-29-3", name="DDT"),
    list(abbrev="HxBB", cas="36355-01-8", name="Hexabromobiphenyl"),
    list(abbrev="PFOA", cas="335-67-1", name="PFOA"),
    list(abbrev="PFOS", cas="1763-23-1", name="PFOS"),
    list(abbrev="PFNA", cas="375-95-1", name="PFNA"),
    list(abbrev="PFHxS", cas="355-46-4", name="PFHXS"),
    list(abbrev="Lovasta", cas="75330-75-5", name="Lovastatin"),
    list(abbrev="Benzoph", cas="119-61-9", name="Benzophenone"),
    list(abbrev="Fenari", cas="60168-88-9", name="Fenarimol"),
    list(abbrev="Fipro", cas="120068-37-3", name="Fipronil"),
    list(abbrev="Naphtha", cas="91-20-3", name="Naphthalene"),
    list(abbrev="Triclo", cas="3380-34-5", name="Triclosan"),
    list(abbrev="DEET", cas="134-62-3", name="DEET"),
    list(abbrev="Prochlo", cas="67747-09-5", name="Prochloraz"),
    list(abbrev="DES", cas="56-53-1", name="Diethylstilbestrol"),
    list(abbrev="Promet", cas="1610-18-0", name="Prometon"),
    list(abbrev="Propico", cas="60207-90-1", name="Propiconazole"),
    list(abbrev="Raloxi", cas="84449-90-1", name="Raloxifene"),
    list(abbrev="Difeno", cas="119446-68-3", name="Difenoconazole"),
    list(abbrev="Pyrimet", cas="53112-28-0", name="Pyrimethanil"),
    list(abbrev="Aldica", cas="116-06-3", name="Aldicarb"),
    list(abbrev="Hexazi", cas="51235-04-2", name="Hexazinone"),
    list(abbrev="Spiroxa", cas="118134-30-8", name="Spiroxamine"),
    list(abbrev="Dipheny", cas="122-39-4", name="Diphenylamine"),
    list(abbrev="Clofen", cas="74115-24-5", name="Clofentezine"),
    list(abbrev="Fenbuco", cas="114369-43-6", name="Fenbuconazole"),
    list(abbrev="Fenhexa", cas="126833-17-8", name="Fenhexamid"),
    list(abbrev="Tebufen", cas="112410-23-8", name="Tebufenozide"),
    list(abbrev="Pyridab", cas="96489-71-3", name="Pyridaben"),
    list(abbrev="Fluroxy", cas="69377-81-7", name="Fluroxypyr"),
    list(abbrev="Flusila", cas="85509-19-9", name="Flusilazole"),
    list(abbrev="Fenthio", cas="55-38-9", name="Fenthion"),
    list(abbrev="Simazi", cas="122-34-9", name="Simazine"),
    list(abbrev="Propox", cas="114-26-1", name="Propoxur"),
    list(abbrev="Fenitro", cas="122-14-5", name="Fenitrothion"),
    list(abbrev="Thidia", cas="51707-55-2", name="Thidiazuron"),
    list(abbrev="Pymetro", cas="123312-89-0", name="Pymetrozine"),
    list(abbrev="Triadim", cas="55219-65-3", name="Triadimenol"),
    list(abbrev="Fenami", cas="161326-34-7", name="Fenamidone"),
    list(abbrev="24D", cas="94-75-7", name="2,4-Dichlorophenoxyacetic acid"),
    list(abbrev="Etoxaz", cas="153233-91-1", name="Etoxazole"),
    list(abbrev="Mesotr", cas="104206-82-8", name="Mesotrione"),
    list(abbrev="Hexaco", cas="79983-71-4", name="Hexaconazole"),
    list(abbrev="Lactof", cas="77501-63-4", name="Lactofen"),
    list(abbrev="Thiaben", cas="148-79-8", name="Thiabendazole"),
    list(abbrev="Fluoxa", cas="361377-29-9", name="Fluoxastrobin"),
    list(abbrev="Caffei", cas="58-08-2", name="Caffeine"),
    list(abbrev="MEHP", cas="4376-20-9", name="Mono(2-ethylhexyl) phthalate"),
    list(abbrev="MBP", cas="131-70-4", name="Mono-n-butyl phthalate"),
    list(abbrev="MBzP", cas="2528-16-7", name="Monobenzyl phthalate"),
    list(abbrev="5OH-MEHP", cas="40321-98-0", name="Mono(2-ethyl-5-hydroxyhexyl) phthalate"),
    list(abbrev="5oxo-MEHP", cas="40321-99-1", name="Mono(2-ethyl-5-oxohexyl) phthalate"),
    list(abbrev="MiBP", cas="30833-53-5", name="Monoisobutyl phthalate"),
    list(abbrev="MCHP", cas="7517-36-4", name="Monocyclohexyl phthalate"),
    list(abbrev="BPA-Gluc", cas="63562-33-4", name="BPA glucuronide"),
    list(abbrev="Methamido", cas="10265-92-6", name="Methamidophos"),
    list(abbrev="OH-Flut", cas="52806-53-8", name="Hydroxyflutamide"),
    list(abbrev="Desmethyl-Acetam", cas="194992-44-4", name="N-desmethyl-acetamiprid"),
    list(abbrev="Vinclo-M1", cas="71707-56-7", name="Vinclozolin M1 (butenoic acid metab.)"),
    list(abbrev="Vinclo-M2", cas="66246-88-6", name="Vinclozolin M2 (enanilide metab.)"),
    list(abbrev="6OH-Dexame", cas="2135-17-3", name="6β-hydroxydexamethasone"),
    list(abbrev="APAP-Gluc", cas="120066-54-8", name="Acetaminophen glucuronide"),
    list(abbrev="APAP-Sulf", cas="4410-31-5", name="Acetaminophen sulfate"),
    list(abbrev="BDCIPP", cas="72236-72-7", name="Bis(1,3-dichloro-2-propyl) phosphate"),
    list(abbrev="BCEP", cas="6294-34-4", name="Bis(2-chloroethyl) phosphate"),
    list(abbrev="TCPy", cas="6515-38-4", name="3,5,6-Trichloro-2-pyridinol"),
    list(abbrev="3-PBA", cas="3739-38-6", name="3-Phenoxybenzoic acid"),
    list(abbrev="DCA", cas="95-76-1", name="3,4-Dichloroaniline"),
    list(abbrev="MeP", cas="99-76-3", name="Methylparaben"),
    list(abbrev="PrP", cas="94-13-3", name="Propylparaben"),
    list(abbrev="BuP", cas="94-26-8", name="Butylparaben"),
    list(abbrev="BPF", cas="620-92-8", name="Bisphenol F"),
    list(abbrev="BPAF", cas="1478-61-1", name="Bisphenol AF"),
    list(abbrev="PFHpA", cas="375-85-9", name="Perfluoroheptanoic acid"),
    list(abbrev="GenX", cas="13252-13-6", name="HFPO-DA (GenX)"),
    list(abbrev="PFBS", cas="375-73-5", name="Perfluorobutane sulfonic acid"),
    list(abbrev="Malath", cas="121-75-5", name="Malathion"),
    list(abbrev="Diazin", cas="333-41-5", name="Diazinon"),
    list(abbrev="Triclocarban", cas="101-20-2", name="Triclocarban"),
    list(abbrev="Octocry", cas="6197-30-4", name="Octocrylene"),
    list(abbrev="Oxyben", cas="131-57-7", name="Oxybenzone (BP-3)"),
    list(abbrev="Nicotine", cas="54-11-5", name="Nicotine"),
    list(abbrev="Cotinine", cas="486-56-6", name="Cotinine"),
    list(abbrev="HCB", cas="118-74-1", name="Hexachlorobenzene"),
    list(abbrev="PBDE47", cas="40088-47-9", name="2,2',4,4'-Tetrabromodiphenyl ether"),
    list(abbrev="Tizoxanide", cas="173903-47-4", name="Tizoxanide"),
    list(abbrev="UV327", cas="3864-99-1", name="UV-327"),
    list(abbrev="Triclop", cas="55335-06-3", name="Triclopyr"),
    list(abbrev="BaP", cas="50-32-8", name="Benzo[a]pyrene"),
    list(abbrev="Tamox", cas="10540-29-1", name="Tamoxifen"),
    list(abbrev="Erythro", cas="114-07-8", name="Erythromycin"),
    list(abbrev="Zamifen", cas="127308-82-1", name="Zamifenacin"),
    list(abbrev="Triamc", cas="124-94-7", name="Triamcinolone"),
    list(abbrev="Rifamp", cas="13292-46-1", name="Rifampicin"),
    list(abbrev="Pyrene", cas="129-00-0", name="Pyrene"),
    list(abbrev="Carbamaz", cas="298-46-4", name="Carbamazepine"),
    list(abbrev="AZT", cas="30516-87-1", name="Azidothymidine"),
    list(abbrev="Genist", cas="446-72-0", name="Genistein"),
    list(abbrev="Haloper", cas="52-86-8", name="Haloperidol"),
    list(abbrev="Prednis", cas="53-03-2", name="Prednisone"),
    list(abbrev="Phenyt", cas="57-41-0", name="Phenytoin"),
    list(abbrev="Progest", cas="57-83-0", name="Progesterone"),
    list(abbrev="Pyrimet_am", cas="58-14-0", name="Pyrimethamine"),
    list(abbrev="Sulfasal", cas="599-79-1", name="Sulphasalazine"),
    list(abbrev="Androst", cas="63-05-8", name="Androstenidione"),
    list(abbrev="Simvast", cas="79902-63-9", name="Simvastatin"),
    list(abbrev="Warfarin", cas="81-81-2", name="Warfarin"),
    list(abbrev="Flucona", cas="86386-73-4", name="Fluconazole"),
    list(abbrev="PhE", cas="122-99-6", name="2-Phenoxyethanol")
  )
  
  cat(sprintf("\n== PIPELINE RUNNING: %s ==\n", QSPR_SOURCE))
  flush.console()
  
  safe_numeric <- function(x, fallback = NA_real_) {
    if (is.null(x))   return(fallback)
    if (is.list(x))   x <- unlist(x)
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (length(x) == 0) return(fallback)
    x[1]
  }
  
  calc_auc_manual  <- function(time, conc) {
    n <- length(time); sum((conc[-1] + conc[-n]) / 2 * diff(time))
  }
  calc_cmax        <- function(conc) max(conc, na.rm = TRUE)
  calc_tmax        <- function(time, conc) time[which.max(conc)]
  safe_stat        <- function(expr) tryCatch(expr, error = function(e) NA)
  
  calc_halflife_manual <- function(time, conc) {
    n    <- length(time)
    term <- seq(round(0.7 * n), n)
    valid <- which(is.finite(conc[term]) & conc[term] > 1e-12)
    if (length(valid) < 3) return(NA_real_)
    t_sub <- time[term][valid]
    c_sub <- conc[term][valid]
    fit   <- tryCatch(lm(log(c_sub + 1e-12) ~ t_sub), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    kel <- abs(coef(fit)[2])
    if (!is.finite(kel) || kel < 1e-5) return(NA_real_)
    log(2) / kel
  }
  
  plot_each_compartment <- function(df, species_label, chem_name,
                                    obs_parent_df = NULL,
                                    sim_dose = NULL) {
    skip_cols <- c("time", "Rblood2plasma", "Rfblood2plasma")
    comp_cols <- setdiff(names(df), skip_cols)
    plots <- list()
    for (col in comp_cols) {
      if (all(is.na(df[[col]])) || all(df[[col]] == 0, na.rm = TRUE)) next
      p <- ggplot(df, aes(x = time, y = .data[[col]])) +
        geom_line(colour = "#2E86C1", linewidth = 0.9) +
        labs(title    = col,
             subtitle = paste0(chem_name, " — ", species_label),
             x = "Time (h)", y = paste0(col, " (µM)")) +
        theme_pbpk()
      
      if (col == "Cplasma" && !is.null(obs_parent_df) &&
          !grepl("Pregnancy|Preg|Rat", species_label, ignore.case = TRUE)) {
        sim_d   <- if (!is.null(sim_dose)) sim_dose else 1
        scale_f <- if (!is.null(obs_parent_df$dose_mg_kg) &&
                       obs_parent_df$linear_pk[1] &&
                       obs_parent_df$dose_mg_kg[1] > 0)
          sim_d / obs_parent_df$dose_mg_kg[1] else 1
        obs_sc <- obs_parent_df
        obs_sc$Cplasma_uM <- obs_sc$Cplasma_uM * scale_f
        if (!all(is.na(obs_sc$sd_uM))) obs_sc$sd_uM <- obs_sc$sd_uM * scale_f
        dose_note <- if (abs(scale_f - 1) < 0.01)
          sprintf("Observed @ %.4g mg/kg (matched)", obs_parent_df$dose_mg_kg[1])
        else
          sprintf("Obs %.4g → scaled to %.4g mg/kg",
                  obs_parent_df$dose_mg_kg[1], sim_d)
        p <- p +
          geom_point(data = obs_sc,
                     aes(x = time_h, y = Cplasma_uM),
                     colour = "#C0392B", shape = 17, size = 3, inherit.aes = FALSE) +
          labs(caption = paste0("▲ Observed parent plasma (", dose_note, ")\n",
                                obs_parent_df$matrix[1]))
        if (!all(is.na(obs_sc$sd_uM)))
          p <- p +
          geom_errorbar(data = obs_sc,
                        aes(x = time_h, ymin = Cplasma_uM - sd_uM,
                            ymax = Cplasma_uM + sd_uM),
                        colour = "#C0392B", width = 0.25,
                        linewidth = 0.5, inherit.aes = FALSE)
      }
      plots[[col]] <- p
    }
    return(plots)
  }
  
  mk_hdr <- function(fill)
    createStyle(fontColour = "#FFFFFF", fgFill = fill, halign = "CENTER",
                fontName = "Calibri", fontSize = 10, textDecoration = "Bold")
  hdr_pk    <- mk_hdr("#2E4057"); hdr_human <- mk_hdr("#1A5276")
  hdr_rat   <- mk_hdr("#7B241C"); hdr_preg  <- mk_hdr("#1D6A39")
  hdr_qspr  <- mk_hdr("#6C3483")
  lbl_style <- createStyle(textDecoration = "Bold", fontName = "Calibri",
                           fontSize = 11, fgFill = "#D5D8DC")
  row_even  <- createStyle(fgFill = "#EBF5FB", fontName = "Calibri", fontSize = 9)
  row_odd   <- createStyle(fgFill = "#FDFEFE", fontName = "Calibri", fontSize = 9)
  
  apply_table_style <- function(wb, sheet, df, start_row, start_col, hdr_style) {
    nc <- ncol(df); nr <- nrow(df)
    addStyle(wb, sheet, hdr_style, rows = start_row,
             cols = start_col:(start_col + nc - 1), gridExpand = TRUE)
    if (nr > 0) {
      odd_rows  <- start_row + seq(1, nr, by = 2)
      even_rows <- start_row + seq(2, nr, by = 2)
      if (length(odd_rows) > 0)
        addStyle(wb, sheet, row_odd, rows = odd_rows,
                 cols = start_col:(start_col + nc - 1), gridExpand = TRUE, stack = TRUE)
      if (length(even_rows) > 0)
        addStyle(wb, sheet, row_even, rows = even_rows,
                 cols = start_col:(start_col + nc - 1), gridExpand = TRUE, stack = TRUE)
    }
    setColWidths(wb, sheet, cols = start_col:(start_col + nc - 1), widths = 14)
  }
  
  write_block <- function(wb, sheet, df, start_row, start_col,
                          hdr_style, label = NULL) {
    if (is.null(df) || nrow(df) == 0) return(start_row)
    row <- start_row
    if (!is.null(label)) {
      writeData(wb, sheet, label, startRow = row, startCol = start_col)
      addStyle(wb, sheet, lbl_style, rows = row,
               cols = start_col:(start_col + ncol(df) - 1), gridExpand = TRUE)
      mergeCells(wb, sheet, cols = start_col:(start_col + ncol(df) - 1), rows = row)
      row <- row + 1
    }
    writeData(wb, sheet, df, startRow = row, startCol = start_col)
    apply_table_style(wb, sheet, df, row, start_col, hdr_style)
    return(row + nrow(df) + 2)
  }
  
  embed_plots <- function(wb, plots, sheet_name, w = 6, h = 3.5, rows_per = 22) {
    if (length(plots) == 0) return(invisible(NULL))
    addWorksheet(wb, sheet_name)
    row <- 2
    for (nm in names(plots)) {
      tmp <- tempfile(fileext = ".png")
      ggsave(tmp, plot = plots[[nm]], width = w, height = h, dpi = 150, bg = "white")
      insertImage(wb, sheet_name, tmp, startRow = row, startCol = 2,
                  width = w, height = h, units = "in")
      row <- row + rows_per
    }
  }
  
  all_pk_stats     <- list()
  skipped          <- c()
  curve_store      <- list(parent = list())
  
  cat("── Loading httk chemical database ... ")
  flush.console()
  avail_cas_global <- suppressMessages(get_cheminfo(info = "CAS"))
  cat(sprintf("✓  (%d compounds loaded natively)\n\n", length(avail_cas_global)))
  flush.console()
  
  for (cpd in compounds) {
    if (!cpd$cas %in% avail_cas_global) {
      cat(sprintf("[SKIP] %s (CAS %s not in active httk database)\n\n",
                  cpd$abbrev, cpd$cas))
      flush.console()
      skipped <- c(skipped, cpd$abbrev)
      next
    }
    
    is_qspr <- (QSPR_SOURCE != "Default")
    cat(sprintf("── %s : %s%s\n", cpd$abbrev, cpd$name,
                if (is_qspr) paste0(" [", QSPR_SOURCE, " loaded]") else ""))
    flush.console()
    
    cpd_mw_cached <- safe_numeric(tryCatch(
      suppressMessages(get_physchem_param(param = "MW", chem.cas = cpd$cas)),
      error = function(e) NA), fallback = 300)
    
    cpd_dose <- if (!is.null(cpd_obs_dose[[cpd$abbrev]]))
      cpd_obs_dose[[cpd$abbrev]] else DEFAULT_DOSE_MG_KG
    cat(sprintf("   Dose: %.4g mg/kg%s\n", cpd_dose,
                if (!is.null(cpd_obs_dose[[cpd$abbrev]]))
                  " [matched to observed]" else " [default 1 mg/kg]"))
    flush.console()
    
    h_human <- if (is_qspr) hdr_qspr else hdr_human
    h_rat   <- if (is_qspr) hdr_qspr else hdr_rat
    h_preg  <- if (is_qspr) hdr_qspr else hdr_preg
    
    max_t_h <- 24
    if (!is.null(obs_parent[[cpd$abbrev]]$Human)) {
      max_obs <- max(obs_parent[[cpd$abbrev]]$Human$time_h, na.rm = TRUE)
      if (is.finite(max_obs) && max_obs > max_t_h) max_t_h <- ceiling(max_obs)
    }
    
    # -------------------------------------------------------------------------
    # EXPLICIT PARAMETERIZATION FIX: Guarantees the Schmitt method is forced
    # prior to passing values to the ODE solver.
    # -------------------------------------------------------------------------
    run_sim <- function(cas, species, extra = list()) {
      tryCatch(
        suppressMessages(suppressWarnings({
          param_args <- list(chem.cas = cas, species = species)
          if (isTRUE(extra$default.to.human)) param_args$default.to.human <- TRUE
          
          params <- do.call(parameterize_pbtk, param_args)
          
          base_args <- list(parameters = params,
                            dose       = cpd_dose,
                            times      = seq(0, max_t_h, 0.1) / 24,
                            maxsteps   = 500000)
          
          df <- as.data.frame(do.call(solve_pbtk, base_args))
          df$time <- df$time * 24                          # days -> hours 
          df
        })),
        error = function(e) {
          cat(sprintf("FAIL [%s]: %s\n", species, conditionMessage(e)))
          flush.console()
          if (!isTRUE(extra$default.to.human)) {
            cat(sprintf("   Retrying with default.to.human=TRUE ... "))
            flush.console()
            
            tryCatch(
              suppressMessages(suppressWarnings({
                param_args2 <- list(chem.cas = cas, species = species, default.to.human = TRUE)
                params2 <- do.call(parameterize_pbtk, param_args2)
                
                base_args2 <- list(parameters = params2,
                                   dose       = cpd_dose,
                                   times      = seq(0, max_t_h, 0.1) / 24,
                                   maxsteps   = 500000)
                
                df2 <- as.data.frame(do.call(solve_pbtk, base_args2))
                df2$time <- df2$time * 24
                df2
              })),
              error = function(e2) {
                cat(sprintf("also FAIL: %s\n", conditionMessage(e2)))
                flush.console(); NULL
              })
          } else NULL
        })
    }
    
    cat("   Parent Human     ... "); flush.console()
    out_human <- run_sim(cpd$cas, "Human")
    if (!is.null(out_human)) cat(sprintf("✓ (%d rows)\n", nrow(out_human))) else cat("NULL\n")
    flush.console()
    
    cat("   Parent Rat       ... "); flush.console()
    out_rat <- run_sim(cpd$cas, "Rat", list(default.to.human = TRUE))
    if (!is.null(out_rat)) cat(sprintf("✓ (%d rows)\n", nrow(out_rat))) else cat("NULL\n")
    flush.console()
    
    cat("   Parent Pregnancy ... "); flush.console()
    out_preg <- tryCatch(
      suppressMessages(suppressWarnings({
        preg_params <- parameterize_fetal_pbtk(chem.cas = cpd$cas)
        
        as.data.frame(solve_fetal_pbtk(parameters = preg_params,
                                       dose     = cpd_dose,
                                       times    = seq(PREG_T0_H, 280, 0.1),
                                       maxsteps = 100000))
      })),
      error = function(e) {
        cat(sprintf("FAIL: %s\n", conditionMessage(e)))
        flush.console(); NULL 
      }
    )
    if (!is.null(out_preg)) cat(sprintf("✓ (%d rows)\n", nrow(out_preg))) else cat("NULL\n")
    flush.console()
    
    curve_store$parent[[cpd$abbrev]] <- bind_rows(
      if (!is.null(out_human)) data.frame(time = out_human$time, Cplasma = out_human$Cplasma, Population = "Human"),
      if (!is.null(out_rat))   data.frame(time = out_rat$time,   Cplasma = out_rat$Cplasma,   Population = "Rat"),
      if (!is.null(out_preg))  data.frame(time = out_preg$time - PREG_T0_H, Cplasma = out_preg$Cplasma, Population = "Pregnancy")
    )
    
    CLhep_h <- safe_numeric(safe_stat(suppressMessages(
      calc_hepatic_clearance(chem.cas = cpd$cas, species = "Human",
                             suppress.messages = TRUE))))
    if (is.na(CLhep_h) || CLhep_h == 0)
      CLhep_h <- safe_numeric(safe_stat(suppressMessages(
        calc_total_clearance(chem.cas = cpd$cas, species = "Human"))))
    if (is.na(CLhep_h) || CLhep_h == 0) CLhep_h <- 0.5
    
    CLhep_r <- safe_numeric(safe_stat(suppressMessages(
      calc_hepatic_clearance(chem.cas = cpd$cas, species = "Rat",
                             default.to.human = TRUE,
                             suppress.messages = TRUE))))
    if (is.na(CLhep_r) || CLhep_r == 0)
      CLhep_r <- safe_numeric(safe_stat(suppressMessages(
        calc_total_clearance(chem.cas = cpd$cas, species = "Rat",
                             default.to.human = TRUE))))
    if (is.na(CLhep_r) || CLhep_r == 0) CLhep_r <- 0.8
    
    cat(sprintf("   CLhep H=%.3f R=%.3f\n", CLhep_h, CLhep_r))
    flush.console()
    
    mk_pk <- function(out, sp, CLhep) {
      if (is.null(out) || !"Cplasma" %in% names(out)) return(NULL)
      sp_httk <- if (grepl("Rat", sp, ignore.case = TRUE)) "Rat" else "Human"
      auc <- safe_numeric(safe_stat(suppressMessages(
        calc_tkstats(chem.cas = cpd$cas, species = sp_httk)$AUC)))
      if (is.na(auc) || auc == 0)
        auc <- calc_auc_manual(out$time, out$Cplasma)
      cmax   <- calc_cmax(out$Cplasma)
      tmax   <- calc_tmax(out$time, out$Cplasma)
      t_half <- safe_numeric(safe_stat(suppressMessages(
        calc_half_life(chem.cas = cpd$cas, species = sp_httk))))
      if (is.na(t_half) || t_half <= 0)
        t_half <- calc_halflife_manual(out$time, out$Cplasma)
      if (!is.finite(t_half)) t_half <- NA_real_
      Dose_umol_kg <- cpd_dose * 1000 / cpd_mw_cached
      cl_val <- if (!is.na(auc) && auc > 0) Dose_umol_kg / auc else CLhep
      vd <- if (!is.na(cl_val) && cl_val > 0 &&
                !is.na(t_half) && t_half > 0 && is.finite(t_half))
        cl_val * t_half / log(2)
      else {
        v_httk <- safe_numeric(safe_stat(suppressMessages(
          calc_vdist(chem.cas = cpd$cas, species = sp_httk))))
        if (!is.na(v_httk) && v_httk > 0) v_httk else NA_real_
      }
      qh_lim   <- if (grepl("Rat", sp, ignore.case = TRUE)) QH_RAT_L_H_KG else QH_HUMAN_L_H_KG
      supra_cl <- !is.na(cl_val) && is.finite(cl_val) && cl_val > qh_lim
      cl_flag  <- if (supra_cl)
        sprintf(">%.1f L/h/kg (flow-limited; CL artefact)", qh_lim) else NA_character_
      hl_flag  <- if (is.na(t_half)) {
        if (supra_cl) "NA:SUPRA-CL"
        else if (grepl("Pregnancy|Preg", sp, ignore.case = TRUE)) "NA:PREG-INIT"
        else "NA:WINDOW"
      } else NA_character_
      tmax_pd <- if (grepl("Pregnancy|Preg", sp, ignore.case = TRUE))
        round(tmax - PREG_T0_H, 3) else round(tmax, 3)
      data.frame(Compound         = cpd$abbrev, Full_Name = cpd$name,
                 Species          = sp,
                 Data_Source      = if (is_qspr) QSPR_SOURCE else "httk native",
                 AUC_uM_h         = round(auc, 4),
                 Cmax_uM          = round(cmax, 6),
                 Tmax_h           = round(tmax, 3),
                 Tmax_post_dose_h = tmax_pd,
                 HalfLife_h       = if (is.na(t_half)) NA_real_ else round(t_half, 3),
                 HalfLife_Flag    = hl_flag,
                 CL_L_h_kg        = round(cl_val, 5),
                 CL_Flag          = cl_flag,
                 Vdist_L_kg       = if (is.na(vd)) NA_real_ else round(vd, 4),
                 Fetal_AUC        = NA)
    }
    
    pk_h <- mk_pk(out_human, "Human", CLhep_h)
    pk_r <- mk_pk(out_rat,   "Rat",   CLhep_r)
    
    pk_p <- NULL
    if (!is.null(out_preg) && "Cplasma" %in% names(out_preg)) {
      mat_auc   <- if ("AUC"  %in% names(out_preg)) tail(out_preg$AUC, 1)
      else calc_auc_manual(out_preg$time, out_preg$Cplasma)
      fetal_auc <- if ("fAUC" %in% names(out_preg)) tail(out_preg$fAUC, 1)
      else if ("Cfplasma" %in% names(out_preg))
        calc_auc_manual(out_preg$time, out_preg$Cfplasma)
      else NA
      cmax_p   <- calc_cmax(out_preg$Cplasma)
      tmax_p   <- calc_tmax(out_preg$time, out_preg$Cplasma)
      t_half_p <- calc_halflife_manual(out_preg$time, out_preg$Cplasma)
      if (!is.finite(t_half_p)) t_half_p <- NA_real_
      Dose_umol_kg_p <- cpd_dose * 1000 / cpd_mw_cached
      cl_p <- if (!is.na(mat_auc) && mat_auc > 0) Dose_umol_kg_p / mat_auc
      else CLhep_h * 0.85
      vd_p <- if (!is.na(cl_p) && cl_p > 0 &&
                  !is.na(t_half_p) && t_half_p > 0 && is.finite(t_half_p))
        cl_p * t_half_p / log(2) else NA_real_
      supra_cl_p <- !is.na(cl_p) && is.finite(cl_p) && cl_p > QH_HUMAN_L_H_KG
      cl_flag_p  <- if (supra_cl_p)
        sprintf(">%.1f L/h/kg (flow-limited; CL artefact)", QH_HUMAN_L_H_KG) else NA_character_
      hl_flag_p  <- if (is.na(t_half_p)) {
        if (supra_cl_p) "NA:SUPRA-CL" else "NA:PREG-INIT"
      } else NA_character_
      tmax_pd_p <- round(tmax_p - PREG_T0_H, 3)
      pk_p <- data.frame(Compound         = cpd$abbrev, Full_Name = cpd$name,
                         Species          = "Pregnancy",
                         Data_Source      = if (is_qspr) QSPR_SOURCE else "httk native",
                         AUC_uM_h         = round(mat_auc, 4),
                         Cmax_uM          = round(cmax_p, 6),
                         Tmax_h           = round(tmax_p, 3),
                         Tmax_post_dose_h = tmax_pd_p,
                         HalfLife_h       = if (is.na(t_half_p)) NA_real_ else round(t_half_p, 3),
                         HalfLife_Flag    = hl_flag_p,
                         CL_L_h_kg        = round(cl_p, 5),
                         CL_Flag          = cl_flag_p,
                         Vdist_L_kg       = if (is.na(vd_p)) NA_real_ else round(vd_p, 4),
                         Fetal_AUC        = round(fetal_auc, 4))
    }
    all_pk_stats[[cpd$abbrev]] <- bind_rows(pk_h, pk_r, pk_p)
    
    cat("   Writing raw data  ... "); flush.console()
    sh <- substr(cpd$abbrev, 1, 31); addWorksheet(wb, sh)
    nr <- 1
    if (!is.null(out_human)) nr <- write_block(wb, sh, out_human, nr, 1, h_human,
                                               sprintf("HUMAN PARENT — %s", cpd$name))
    if (!is.null(out_rat))   nr <- write_block(wb, sh, out_rat, nr, 1, h_rat,
                                               sprintf("RAT PARENT — %s", cpd$name))
    if (!is.null(out_preg)) {
      fetal_extra_cols <- setdiff(
        names(out_preg),
        c("time","Cgut","Cliver","Cven","Cart","Clung","Ckidney","Crest",
          "Cplasma","Rblood2plasma","AUC","Qgut","Qliver","Qkidney","Qrest"))
      preg_label <- sprintf(
        "PREGNANCY PARENT — %s%s",
        cpd$name,
        if (length(fetal_extra_cols) > 0)
          paste0("  [fetal cols: ", paste(fetal_extra_cols, collapse=", "), "]")
        else "  [no fetal compartment columns in this simulation]")
      nr <- write_block(wb, sh, out_preg, nr, 1, h_preg, preg_label)
    }
    
    cat("✓\n   Embedding parent plots   ... "); flush.console()
    obs_parent_cpd <- obs_parent[[cpd$abbrev]]
    
    if (!is.null(out_human))
      embed_plots(wb,
                  plot_each_compartment(out_human, "Human", cpd$name,
                                        obs_parent_df = obs_parent_cpd$Human,
                                        sim_dose      = cpd_dose),
                  substr(paste0(cpd$abbrev, "_H_Plots"), 1, 31))
    if (!is.null(out_rat))
      embed_plots(wb,
                  plot_each_compartment(out_rat, "Rat", cpd$name),
                  substr(paste0(cpd$abbrev, "_R_Plots"), 1, 31))
    if (!is.null(out_preg))
      embed_plots(wb,
                  plot_each_compartment(out_preg, "Pregnancy", cpd$name),
                  substr(paste0(cpd$abbrev, "_P_Plots"), 1, 31))
    
    cat(sprintf("✓\n── %s complete ──────────────────────────────────\n\n",
                cpd$abbrev))
    flush.console()
  }
  
  addWorksheet(wb, "Sheet_Map")
  sheet_map_rows <- bind_rows(lapply(compounds, function(cpd_i) {
    data.frame(
      Sheet           = cpd_i$abbrev,
      Compound        = cpd_i$name,
      Parameterisation= if (QSPR_SOURCE != "Default") QSPR_SOURCE else "httk native",
      Banner_row      = NA_integer_,
      Human_Parent_row    = 1,
      Rat_Parent_row      = 490,
      Pregnancy_Parent_row= 980,
      Note = "v28: Parent-only simulation model.",
      stringsAsFactors = FALSE)
  }))
  writeData(wb, "Sheet_Map",
            "Sheet_Map v28 — Parent simulation only.",
            startRow = 1, startCol = 1)
  addStyle(wb, "Sheet_Map",
           createStyle(textDecoration = "Bold", fontSize = 11, fontName = "Calibri",
                       fgFill = "#D6EAF8", wrapText = TRUE),
           rows = 1, cols = 1:ncol(sheet_map_rows), gridExpand = TRUE)
  setRowHeights(wb, "Sheet_Map", rows = 1, heights = 40)
  writeData(wb, "Sheet_Map", sheet_map_rows, startRow = 2, startCol = 1)
  apply_table_style(wb, "Sheet_Map", sheet_map_rows, 2, 1, hdr_pk)
  setColWidths(wb, "Sheet_Map", cols = which(names(sheet_map_rows) == "Note"), widths = 70)
  setColWidths(wb, "Sheet_Map", cols = which(names(sheet_map_rows) == "Compound"), widths = 35)
  
  full_pk <- bind_rows(all_pk_stats)
  addWorksheet(wb, "PK_Summary")
  write_block(wb, "PK_Summary", full_pk, 1, 1, hdr_pk,
              paste0("Parent PK v28 — All Compounds & Populations",
                     "  |  Standardized: 24 h single dose for all compounds",
                     "  |  Tmax_post_dose_h = Tmax-91 for Pregnancy, Tmax for others",
                     "  |  HalfLife_Flag: NA:SUPRA-CL / NA:PREG-INIT / NA:WINDOW",
                     "  |  CL_Flag: flagged when CL > QH (flow-limited artefact)"))
  
  addWorksheet(wb, "Skipped_Compounds")
  skipped_detail <- data.frame(
    Abbrev    = c("PDDP", "MelAlt", "Mancoz", "Cd", "iAs", "UFPs"),
    Full_Name = c("Phenol alkylation products (PDDP)",
                  "Melaleuca alternifolia tea tree oil",
                  "Mancozeb",
                  "Cadmium (metal ion)",
                  "Inorganic Arsenic (As3+/As5+)",
                  "Ultrafine Particles (<100 nm)"),
    CAS       = c("68937-41-7", "68647-73-4", "8018-01-7",
                  "7440-43-9",  "7440-38-2",  "N/A"),
    Skip_Type = c("Genuine — Mixture", "Genuine — Mixture",
                  "Genuine — Mixture/Polymer",
                  "Metal ion — PBTK not applicable",
                  "Metalloid — Specialised As-PBPK required",
                  "Particle mixture — Not a chemical entity"),
    Scientific_Reason = c(
      "MIXTURE — Complex commercial phenol alkylation blend; no single MW/logP.",
      "MIXTURE — Tea tree essential oil (>100 terpene components).",
      "POLYMER — Mn/Zn-EBDC coordination polymer.",
      "METAL ION — Cd2+ has no meaningful logP or Fup.",
      "METALLOID — Inorganic arsenic undergoes methylation via AS3MT enzyme.",
      "PARTICLE MIXTURE — UFPs (<100 nm) are not a single chemical entity."),
    stringsAsFactors = FALSE)
  writeData(wb, "Skipped_Compounds",
            "Skipped Compounds — Mixtures, metal ions, and non-chemical entities (v28)",
            startRow = 1, startCol = 1)
  addStyle(wb, "Skipped_Compounds",
           createStyle(textDecoration = "Bold", fontSize = 13, fontName = "Calibri",
                       fgFill = "#FADBD8", fontColour = "#7B241C"),
           rows = 1, cols = 1:5, gridExpand = TRUE)
  writeData(wb, "Skipped_Compounds", skipped_detail, startRow = 2, startCol = 1)
  apply_table_style(wb, "Skipped_Compounds", skipped_detail, 2, 1, mk_hdr("#7B241C"))
  setColWidths(wb, "Skipped_Compounds",
               cols = which(names(skipped_detail) %in%
                              c("Scientific_Reason", "Full_Name", "Skip_Type")),
               widths = c(70, 35, 25))
  
  all_sh   <- names(wb)
  priority <- c("Sheet_Map", "PK_Summary", "Skipped_Compounds")
  priority <- priority[priority %in% all_sh]
  rest     <- all_sh[!all_sh %in% priority]
  worksheetOrder(wb) <- c(match(priority, all_sh), match(rest, all_sh))
  
  excel_path <- file.path(getwd(), sprintf("PBPK_QSPR_%s_Results.xlsx", QSPR_OUTPUT_SUFFIX))
  cat("\n── Building summary sheets and saving workbook ... ")
  flush.console()
  saveWorkbook(wb, excel_path, overwrite = TRUE)
  if (file.exists(excel_path))
    cat(sprintf("\n✔  Saved → %s\n", normalizePath(excel_path)))
  
  invisible(list(excel_path = excel_path,
                 pk_stats = full_pk,
                 curves = curve_store))
}

# =============================================================================
# DRIVER
# =============================================================================
qspr_sources <- c("Default", "Sipes2017", "Pradeep2020", "Dawson2021")
results_by_source <- list()

for (.src in qspr_sources) {
  cat(sprintf("\n\n#####################################################\n"))
  cat(sprintf("### RUNNING FULL PIPELINE -- QSPR SOURCE: %s\n", .src))
  cat(sprintf("#####################################################\n\n"))
  flush.console()
  res <- tryCatch(
    run_pipeline_for_source(QSPR_SOURCE = .src, QSPR_OUTPUT_SUFFIX = .src),
    error = function(e) {
      cat(sprintf("\n!! PIPELINE FAILED for source %s: %s\n", .src, conditionMessage(e)))
      NULL
    })
  if (!is.null(res)) results_by_source[[.src]] <- res
  gc()
}

cat("\n\nAll QSPR-source pipeline runs complete.\n")

# =============================================================================
# COMBINED COMPARISON WORKBOOK
# =============================================================================
if (length(results_by_source) > 0) {
  
  cat("\nBuilding combined comparison workbook ...\n")
  flush.console()
  
  src_colors <- c(Default = "#4A4A4A", Sipes2017 = "#2E86AB", Pradeep2020 = "#A23B72", Dawson2021 = "#F18F01")
  
  cwb <- createWorkbook()
  
  pk_all <- bind_rows(lapply(names(results_by_source), function(s) {
    df <- results_by_source[[s]]$pk_stats
    if (!is.null(df) && nrow(df) > 0) cbind(Source = s, df) else NULL
  }))
  if (nrow(pk_all) > 0) {
    addWorksheet(cwb, "PK_Summary")
    writeData(cwb, "PK_Summary", pk_all)
    setColWidths(cwb, "PK_Summary", cols = 1:ncol(pk_all), widths = "auto")
    
    if ("Human" %in% pk_all$Species) {
      auc_human <- pk_all %>%
        filter(Species == "Human") %>%
        select(Compound, Full_Name, Source, AUC_uM_h) %>%
        dcast(Compound + Full_Name ~ Source, value.var = "AUC_uM_h")
      
      addWorksheet(cwb, "Human_AUC_Comparison")
      writeData(cwb, "Human_AUC_Comparison", 
                "Human AUC (uM*h) Comparison Across QSPR Parameterizations", 
                startRow = 1, startCol = 1)
      addStyle(cwb, "Human_AUC_Comparison",
               createStyle(textDecoration = "Bold", fontSize = 11, fontName = "Calibri", fgFill = "#D6EAF8"),
               rows = 1, cols = 1:ncol(auc_human), gridExpand = TRUE)
      writeData(cwb, "Human_AUC_Comparison", auc_human, startRow = 2)
      
      addStyle(cwb, "Human_AUC_Comparison",
               createStyle(fontColour = "#FFFFFF", fgFill = "#1A5276", halign = "CENTER", fontName = "Calibri", textDecoration = "Bold"),
               rows = 2, cols = 1:ncol(auc_human), gridExpand = TRUE)
      setColWidths(cwb, "Human_AUC_Comparison", cols = 1:ncol(auc_human), widths = "auto")
    }
  }
  
  addWorksheet(cwb, "Cplasma_Overlay")
  row_cursor <- 2
  writeData(cwb, "Cplasma_Overlay",
            "Parent plasma concentration — QSPR sources overlaid, faceted by population",
            startRow = 1, startCol = 1)
  
  all_abbrevs <- unique(unlist(lapply(results_by_source, function(r) names(r$curves$parent))))
  for (ab in all_abbrevs) {
    curves <- bind_rows(lapply(names(results_by_source), function(s) {
      cd <- results_by_source[[s]]$curves$parent[[ab]]
      if (is.null(cd) || nrow(cd) == 0) return(NULL)
      cbind(cd, Source = s)
    }))
    if (nrow(curves) == 0) next
    
    p <- ggplot(curves, aes(x = time, y = Cplasma, colour = Source)) +
      geom_line(linewidth = 0.8, na.rm = TRUE) +
      scale_colour_manual(values = src_colors) +
      facet_wrap(~ Population, scales = "free", ncol = 3) +
      labs(title = paste0(ab, " — Plasma concentration, QSPR sources"),
           x = "Time (h)", y = "Cplasma (µM)") +
      theme_pbpk()
    
    if (ab %in% names(obs_parent) && !is.null(obs_parent[[ab]]$Human)) {
      obs_df <- obs_parent[[ab]]$Human
      obs_df$Population <- "Human" 
      
      p <- p + geom_point(data = obs_df, aes(x = time_h, y = Cplasma_uM),
                          colour = "black", shape = 16, size = 2, inherit.aes = FALSE) +
        labs(caption = paste0("● Observed parent plasma (", obs_df$matrix[1], ")"))
    }
    
    tmp <- tempfile(fileext = ".png")
    ggsave(tmp, plot = p, width = 12, height = 4.2, dpi = 150, bg = "white")
    insertImage(cwb, "Cplasma_Overlay", tmp, startRow = row_cursor, startCol = 1,
                width = 12, height = 4.2, units = "in")
    row_cursor <- row_cursor + 24
  }
  
  compare_path <- file.path(getwd(), "PBPK_QSPR_Comparison.xlsx")
  saveWorkbook(cwb, compare_path, overwrite = TRUE)
  cat(sprintf("\n✔  Saved combined comparison workbook → %s\n", normalizePath(compare_path)))
} else {
  cat("\nNo successful runs — comparison workbook not built.\n")
}

# =============================================================================
# File Name:    Automated Detection Main
# Author:       Ayooluwa Adeyinka
# Created:      June 25, 2026
# Description:  Automated annotation, spectrographic/cepstral feature
#               generation, and ML training data export. Uses template
#               annotations (.txt) from the working directory to generate
#               automated annotations for other files, exported as
#               tab-delimited .txt.
#
# NOTE: mclapply() (used internally by some `ohun`/`warbleR`/`bioacoustics`
#       functions when running in parallel) only works on Mac/Linux — not
#       Windows.
#
# Dependencies: tuneR, warbleR, ohun, Rraven, readr, reticulate, seewave,
#               ggplot2, parallel, bioacoustics, compiler
#
# NOTE: Several variables referenced below (MaxCores, TemplateSubspaces,
#       nharmonics, parallel, fast, ovlp, fsmooth, wl, wl.freq,
#       MinCallDuration, MaxCallDuration, FreqRange, TemplateDetectionTreshhold)
#       are assumed to be defined in a setup/config script sourced before this
#       file runs.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Python environment (macOS-specific)
# -----------------------------------------------------------------------------
# Mac OS specific python library specification.
# If on another system, point reticulate to the python library on your
# system, or leave commented to use RStudio's default.
# Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/bin/python3.14")


# -----------------------------------------------------------------------------
# 1. Load required libraries
# -----------------------------------------------------------------------------

ProgramPackages = c(
  "tuneR",
  "warbleR",
  "ohun",
  "Rraven",
  "reticulate",
  "readr",
  "seewave",
  "ggplot2",
  "parallel",
  "bioacoustics",
  "compiler"
)
load_packages(ProgramPackages)


# -----------------------------------------------------------------------------
# 2. File paths
# -----------------------------------------------------------------------------

AutoGenFileStoragePath <- "~/OneDrive - UCB-O365/Automated Detection Project/R/AutoGenResults"
WorkingDirectory       <- "~/Documents/Automated-Pika-Detection/Target Files"
ML_TrainingDataStorage <- "~/Documents/Automated-Pika-Detection/ML Training data"  # Currently unused!
PlotStorage            <- "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Plots"

setwd(WorkingDirectory)


# -----------------------------------------------------------------------------
# 3. TRAINING — Load reference annotations and build templates
# -----------------------------------------------------------------------------

# Import all .txt Raven selection tables from the working directory
PikaReference <- imp_raven(
  warbler.format = TRUE,   # converts to warbleR column names
  all.data       = FALSE,  # only keep warbleR-required columns, smaller footprint
  parallel = MaxCores
)

gc()
gc()

# Select representative templates from the acoustic space defined by the
# reference calls
PikaTemplate <- get_templates(
  reference    = PikaReference,
  n.sub.spaces = TemplateSubspaces,
  plot         = TRUE,
  nharmonics   = nharmonics,
  # harmonicity = harmonicity,
  parallel = MaxCores,
  fast     = fast,
  ovlp     = ovlp,
  fsmooth  = fsmooth,
  wl       = wl,
  wl.freq  = wl.freq
  
)


# -----------------------------------------------------------------------------
# 4. ANNOTATING — Run detector across all wav files
# -----------------------------------------------------------------------------

# List all wav files in the working directory
wav_files <- list.files(pattern = "\\.wav$", full.names = TRUE)

# Run detection on each file, collecting results into a list
GenPikaAnn <- lapply(wav_files, function(wav) {
  cat("Processing:", wav, "\n")
  detec <- DetectCalls(
    wavFileName     = wav,
    given_template  = PikaTemplate,
    exportanntable  = TRUE,
    minimum_dur     = MinCallDuration,
    maximum_dur     = MaxCallDuration,
    FrequencyRange  = FreqRange
  )
})

gc()
gc()

# Combine all per-file detection tables into one dataframe
GenPikaAnn <- do.call(rbind, GenPikaAnn)


# -----------------------------------------------------------------------------
# 5. Plot a detection vs. reference spectrogram for visual QC
# -----------------------------------------------------------------------------

# TODO: `detec` is assigned *inside* the lapply() function above, so it's a
# local variable scoped to each function call and does not exist in the
# global environment out here. This will error with "object 'detec' not
# found" unless `detec` is assigned somewhere else first (e.g. pulling one
# result out of GenPikaAnn, or re-running DetectCalls() directly).
PikaWavTest <- readWave("~/Documents/Automated-Pika-Detection/Target Files/Large Training Set.wav")

# Reserve outer margin space above the plot for the title (see note below)
par(oma = c(0, 0, 1, 0))

label_spectro(
  wave = PikaWavTest,
  reference = PikaReference,
  detection = detec,
  smooth = 3,
  # template.correlation = Correlations[[1]],
  flim = c(0, 20),
  tlim = c(73.6, 86.2),
  # threshold = TemplateDetectionTreshhold,
  hop.size = 10,
  cexaxis = 1.2,
  cexlab = 1.2,
  ovlp = 90,
  wl = 1024,
  palette = temp.colors,
  collevels = seq(5, 95, 1),
  dBref = 2 * 10e-5,
  # fastdisp = TRUE
)

# Title drawn into the outer margin reserved by par(oma=...) above, per the
# earlier fix for spectro()/label_spectro() title placement.
title(
  main = "|    Computer detections comapred to reference material",
  cex.main = 1.5,
  line = 0,
  col.main = "black",
  outer = TRUE
)
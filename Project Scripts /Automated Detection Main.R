# File Name:       Automated Detection Main
# Author:          Ayooluwa Adeyinka
# Created:         June 25, 2026
# Description:     Automated annotation, spectrographical/cepstral feature
#                  generation and ML training. Uses template annotations (.txt)
#                  from working directory to generate automated annotations for
#                  other files, exports as tab-delimited .txt.
# NOTE: mclapply is for only mac and linux systems
# General Dependencies: tuneR, warbleR, ohun, Rraven, readr, reticulate, seewave, ggplot2,parallel


# Mac OS specific python library specification.
# If on another system point reticulate to the python library on your system or use defualt with Rstudio
# Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/bin/python3.14")

# Load required libraries
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

LoadPackages(ProgramPackages)

# FILE PATHS
# ==============================================================================
AutoGenFileStoragePath <- "~/OneDrive - UCB-O365/Automated Detection Project/R/AutoGenResults"
WorkingDirectory       <- "~/Documents/Automated-Pika-Detection/Target Files"
ML_TrainingDataStorage <- "~/Documents/Automated-Pika-Detection/ML Training data" #Currently unused!
PlotStorage <- "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R/Plots"

setwd(WorkingDirectory)

# ==============================================================================
# TRAINING — Load reference annotations and build templates
# ==============================================================================

# Import all .txt Raven selection tables from working directory
PikaReference <- imp_raven(
  warbler.format = TRUE,
  # converts to warbleR column names
  all.data       = FALSE,
  # only keep warbleR-required columns smaller footprint
  parallel = MaxCores
)
gc()
gc()

# Select representative templates from acoustic space
PikaTemplate <- get_templates(
  reference   = PikaReference,
  n.sub.spaces = TemplateSubspaces,
  plot        = TRUE,
  nharmonics = nharmonics,
  # harmonicity = harmonicity,
  parallel = parallel,
  fast  = fast,
  ovlp = ovlp,
  fsmooth = fsmooth,
  wl = wl,
  wl.freq = wl.freq
)


# ==============================================================================
# ANNOTATING — Run detector across all wav files
# ==============================================================================

# List all wav files in working directory
wav_files <- list.files(pattern = "\\.wav$", full.names = TRUE)

# Run detection on each file, collect results
GenPikaAnn <- lapply(wav_files, function(wav) {
  cat("Processing:", wav, "\n")
  detec <-  DetectCalls(
    wavFileName    = "Large Training Set.wav",
    given_template = PikaTemplate,
    exportanntable = TRUE,
    minimum_dur = MinCallDuration,
    maximum_dur = MaxCallDuration,
    FrequencyRange = FreqRange
  )
})

gc()
gc()
# Combine all detection tables into one dataframe
GenPikaAnn <- do.call(rbind, GenPikaAnn)


# ==============================================================================
PikaWavTest <- readWave("~/Documents/Automated-Pika-Detection/Target Files/Large Training Set.wav")

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
title(
  main = "|    Computer detections comapred to reference material",
  cex.main = 1.5,
  line = 0,
  col.main = "black",
  outer = TRUE
)
# ==============================================================================

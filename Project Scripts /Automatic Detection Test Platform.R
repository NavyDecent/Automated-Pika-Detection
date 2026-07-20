# =============================================================================
# File Name:    Automated Detection Test Platform
# Author:       Ayooluwa Adeyinka
# Created:      June 25, 2026
# Description:  Testbed for stress-testing/optimizing pika call detection
#               variables (template-based and energy-based detection via
#               the `ohun` package).
#
# NOTE: mclapply() (used internally by some `ohun`/`warbleR` functions when
#       cores > 1) only works on Mac/Linux — not Windows.
#
# Dependencies: tuneR, warbleR, ohun, Rraven, readr, reticulate, seewave,
#               ggplot2, parallel
#
# NOTE: Several variables referenced below (WorkingDirectory, maxcores,
#       Detection_Treshhold, WindowLength, FreqRange, Speed_Fast, Overlap,
#       MaxCores, ML_TRainingDataStorage) are assumed to be defined in a
#       setup/config script sourced before this file runs.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load audio files for spectro analysis
# -----------------------------------------------------------------------------

# Large training clip used to build/test templates and run detection
PikaWavTest <- readWave(
  "~/Documents/Automated-Pika-Detection/Target Files/Large Training Set.wav"
)

# Separate recording used as a second, independent test case
PikaConvoWav <- readWave(
  "~/Documents/Automated-Pika-Detection/Target Files/PikaConvoWav.wav"
)


# -----------------------------------------------------------------------------
# 2. Load reference annotations
# -----------------------------------------------------------------------------

# Import all .txt Raven selection tables from the target folder and convert
# them into warbleR-formatted reference data (ground-truth pika calls)
PikaReference <- imp_raven(
  path = "~/Documents/Automated-Pika-Detection/Target Files",
  warbler.format = TRUE,  # critical - converts to warbleR column names
  all.data = FALSE        # only keep warbleR-relevant columns
)

gc()
gc()


# -----------------------------------------------------------------------------
# 3. (Archived) Energy-based detection attempt
# -----------------------------------------------------------------------------
# Kept for reference — superseded by the template-based approach below.
# Uncomment to re-run energy detection / threshold optimization.
# Will Implement later for mor accurate detections

# Energy_detection <- energy_detector(
#   min.duration = 0.1, max.duration = .5, cores = maxcores,
#   hold.time = 200, threshold = 50, thinning = .9, wl = 1024,
#   smooth = 350, bp = c(2, 18), files = "Large Training Set.wav"
# )
#
# optimization <- optimize_energy_detector(
#   files = "Large Training Set.wav",
#   reference = PikaReference,
#   threshold = seq(40, 60, 1),
#   min.duration = 0.1, max.duration = .5, cores = maxcores,
#   hold.time = 200, wl = 1024, smooth = c(200, 400, 20), bp = c(2, 18)
# )


# -----------------------------------------------------------------------------
# 4. Build detection template
# -----------------------------------------------------------------------------

# Get mean acoustic structure template(s) from the reference calls.
# NOTE: get_templates() errors if a value EQUALS the max rather than
# exceeding it — keep this in mind when filtering imported numeric columns.
template <- get_templates(
  reference = PikaReference,
  path = WorkingDirectory,
  n.sub.spaces = 3
)


# -----------------------------------------------------------------------------
# 5. Template cross-correlation
# -----------------------------------------------------------------------------

# Correlate the template against the large training recording
TestCorrelations1 <- template_correlator(
  templates = template,
  files = c("Large Training Set.wav"),
  path = "~/Documents/Automated-Pika-Detection/Target Files",
  wl = 512,
  wl.freq = 1024
)

# Correlate the template against the second (conversation) recording
TestCorrelations2 <- template_correlator(
  templates = template,
  files = c("PikaConvoWav.wav"),
  path = "~/Documents/Automated-Pika-Detection/Target Files",
  wl = 512,
  wl.freq = 1024
)


# -----------------------------------------------------------------------------
# 6. Run detection — Large Training Set
# -----------------------------------------------------------------------------

detectiontest1 <- template_detector(
  template.correlations = TestCorrelations1,
  threshold = Detection_Treshhold,
  cores = maxcores
)

# Label detections against reference, then consolidate duplicate/overlapping
# detections by keeping the highest-scoring one per event
detectiontest1 <- consensus_detection(
  detection = label_detection(
    reference = PikaReference,
    detection = detectiontest1,
    cores = maxcores
  ),
  by = "scores",
  cores = maxcores
)

# Merge any remaining overlapping detections into single events
detectiontest1 <- merge_overlaps(detectiontest1, pb = TRUE, cores = maxcores)


# -----------------------------------------------------------------------------
# 7. Run detection — Pika Convo recording
# -----------------------------------------------------------------------------

detectiontest2 <- template_detector(
  template.correlations = TestCorrelations2,
  threshold = Detection_Treshhold,
  cores = maxcores
)

detectiontest2 <- consensus_detection(
  detection = label_detection(
    reference = PikaReference,
    detection = detectiontest2,
    cores = maxcores
  ),
  by = "scores",
  cores = maxcores
)

detectiontest2 <- merge_overlaps(detectiontest2, pb = TRUE, cores = maxcores)


# -----------------------------------------------------------------------------
# 8. Diagnose detection performance
# -----------------------------------------------------------------------------

# TODO: `TESTER` is not defined anywhere above 
diagnoses <- diagnose_detection(
  reference = PikaReference,
  detection = TESTER,
  cores = maxcores
)

# Optimize the correlation threshold against the reference annotations
optimization <- optimize_template_detector(
  template.correlations = TestCorrelations1,
  reference = PikaReference,
  threshold = seq(0.3, 0.6, 0.01),
  cores = maxcores,
  min.overlap = 0.6
)


# -----------------------------------------------------------------------------
# 9. (Archived) Export detections/features to Raven-format annotation files
# -----------------------------------------------------------------------------
# Uncomment to re-export. Paths currently point to a OneDrive location that
# may not exist on this machine — update before re-enabling.

# exp_raven(
#   detectiontest1,
#   file.name = "Generated annotations1.txt",
#   khz.to.hz = TRUE,
#   path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
# )
# exp_raven(
#   detectiontest2,
#   file.name = "Generated annotations2.txt",
#   khz.to.hz = TRUE,
#   path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
# )
# Spectrographical_features$start <- PikaReference$start
# Spectrographical_features$end   <- PikaReference$end
# Spectrographical_features$label <- "1" # Pika
# str(Spectrographical_features)
#
# exp_raven(
#   Spectrographical_features,
#   file.name = "training dat.txt",
#   khz.to.hz = TRUE,
#   path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
# )


# -----------------------------------------------------------------------------
# 10. Plot labeled spectrograms
# -----------------------------------------------------------------------------

# Spectrogram for the large training set, with reference + detection +
# correlation score overlays
label_spectro(
  wave = PikaWavTest,
  reference = TESTER,
  detection = detectiontest1,
  template.correlation = TestCorrelations1[[1]],
  flim = c(0, 22),
  tlim = c(0, 40),
  threshold = Detection_Treshhold,
  hop.size = 10,
  ovlp = 50,
  palette = reverse.cm.colors,
  fastdisp = TRUE
)

# Spectrogram for the conversation recording, with detection + correlation
# score overlays (no reference table for this one)
label_spectro(
  wave = PikaConvoWav,
  detection = detectiontest2,
  template.correlation = TestCorrelations2[[1]],
  flim = c(0, 22),
  threshold = Detection_Treshhold,
  hop.size = 10,
  ovlp = 50,
  palette = reverse.cm.colors,
)


# -----------------------------------------------------------------------------
# 11. Compare detection methods
# -----------------------------------------------------------------------------

# Compare spectrographic parameters (SP) vs. frequency-DTW (ffDTW) methods
# across the reference calls, saving spectrogram images for visual review
compare_methods(
  X = PikaReference,
  methods = c("SP", "ffDTW"),
  wl = WindowLength,
  bp = FreqRange,
  fast = Speed_Fast,
  flim = FreqRange,
  it = "jpeg",
  res = 2000,
  length.out = 25,
  path = WorkingDirectory,
  ovlp = Overlap,
  parallel = MaxCores
)


# -----------------------------------------------------------------------------
# 12. Export ML training data
# -----------------------------------------------------------------------------

# Strip identifier columns not needed for model training
Spectrographical_features$sound.files <- NULL
Spectrographical_features$selec <- NULL

# Label all rows as positive (pika) examples
Spectrographical_features$label <- "1"

# Write features out to the ML training data folder
setwd(ML_TRainingDataStorage)
write.table(
  SpectroCoeff,
  file = "SpectroFeatures.txt",
  sep = "\t",
  row.names = FALSE
)

# Return to the main working directory
setwd(WorkingDirectory)


# -----------------------------------------------------------------------------
# 13. Cleanup
# -----------------------------------------------------------------------------

gc()
gc()
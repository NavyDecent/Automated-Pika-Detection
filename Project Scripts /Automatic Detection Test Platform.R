# File Name:       Automated Detection Main
# Author:          Ayooluwa Adeyinka
# Created:         June 25, 2026
# Description:    This is where optimizations to the detection varibales can be stress tested.
.
# NOTE: mclapply is for only mac and linux systems
# General Dependencies: tuneR, warbleR, ohun, Rraven, readr, reticulate, seewave, ggplot2,parallel





PikaWavTest <- readWave(
  "~/Documents/Automated-Pika-Detection/Target Files/Large Training Set.wav"
)

PikaConvoWav <- readWave(
  "~/Documents/Automated-Pika-Detection/Target Files/PikaConvoWav.wav"
)


# Import all .txt selection tables from a folder
PikaReference <- imp_raven(
  path = "~/Documents/Automated-Pika-Detection/Target Files",
  warbler.format = TRUE,
  # critical — converts to warbleR column names
  all.data = FALSE         # only keep warbleR-relevant columns
)
gc()
gc()



# Energy_detection <- energy_detector( min.duration = 0.1, max.duration = .5, cores = maxcores, hold.time = 200, threshold = 50, thinning = .9, wl=1024, smooth = 350, bp=c(2,18), files = "Large Training Set.wav")
#
# optimization <-
#   optimize_energy_detector(files = "Large Training Set.wav",
#     reference = PikaReference,
#     threshold = seq(40, 60,1),
#     min.duration = 0.1, max.duration = .5, cores = maxcores, hold.time = 200, wl=1024, smooth = c(200,400,20), bp=c(2,18),
#   )




# get mean structure template
#keep in mind it errors out if value = max not value > max for imported numbers
template <-
  get_templates(reference = PikaReference,
                path = WorkingDirectory,
                n.sub.spaces = 3)

# get correlations
TestCorrelations1 <-
  template_correlator(
    templates = template,
    files = c("Large Training Set.wav"),
    path = "~/Documents/Automated-Pika-Detection/Target Files",
    wl =512,
    wl.freq=1024
  )


TestCorrelations2 <-
  template_correlator(
    templates = template,
    files = c("PikaConvoWav.wav"),
    path = "~/Documents/Automated-Pika-Detection/Target Files",
    wl = 512,
    wl.freq=1024
  )


# run detection
detectiontest1 <-
  template_detector(
    template.correlations = TestCorrelations1,
    threshold = Detection_Treshhold,
    cores = maxcores
  )

detectiontest1 <- consensus_detection(
  detection = label_detection(
    reference = PikaReference,
    detection = detectiontest1,
    cores = maxcores
  ),
  by = "scores"
  ,
  cores = maxcores
)

detectiontest1 <- merge_overlaps(detectiontest1, pb = TRUE, cores = maxcores)

# run detection
detectiontest2 <-
  template_detector(
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
  by = "scores"
  ,
  cores = maxcores
)

detectiontest2 <- merge_overlaps(detectiontest2, pb = TRUE, cores = maxcores)

#diagnose detection for optimization
diagnoses <- diagnose_detection(reference = PikaReference,
                                detection = TESTER,
                                cores = maxcores)

# run optimization
optimization <-
  optimize_template_detector(
    template.correlations = TestCorrelations1,
    reference = PikaReference,
    threshold = seq(0.3, 0.6, 0.01),
    cores = maxcores,
    min.overlap = 0.6
  )


# exp_raven(
#   detectiontest1,
#   file.name = "Generated annotations1.txt",
#   khz.to.hz = TRUE,
#   path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
#   
# )
# exp_raven(
#   detectiontest2,
#   file.name = "Generated annotations2.txt",
#   khz.to.hz = TRUE,
#   path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
#   
# )
# Spectrographical_features$start <-  PikaReference$start
# Spectrographical_features$end <-  PikaReference$end
# Spectrographical_features$label <- "1"
# str(Spectrographical_features)
# 
# exp_raven(
#   Spectrographical_features,
#   file.name = "training dat.txt",
#   khz.to.hz = TRUE,
#   path = "~/Library/CloudStorage/OneDrive-UCB-O365/Automated Detection Project/R",
# )



# plot spectrogram
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
  fast.dip = TRUE
)

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

compare_methods(
  X = PikaReference,
  methods = c("SP", "ffDTW"),
  wl = WindowLength,
  bp          = FreqRange,
  fast        = Speed_Fast,
  flim       = FreqRange,
  it = "jpeg",
  res = 2000,
  length.out = 25,
  path = WorkingDirectory,
  ovlp = Overlap,
  parallel = MaxCores,
)

# ML training data export
# ==============================================================================
# Remove non useful ML training  data columns
Spectrographical_features$sound.files <- NULL
Spectrographical_features$selec <- NULL
# label for pikas
Spectrographical_features$label <- "1"

#Export ML training data to ML trainig data folder
setwd(ML_TRainingDataStorage)
write.table(
  SpectroCoeff,
  file = "SpectroFeatures.txt",
  sep = "\t",
  row.names = FALSE
)
# reset working directory
setwd(WorkingDirectory)

#   garbage collections

gc()
gc()

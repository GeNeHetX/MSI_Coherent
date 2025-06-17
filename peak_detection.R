library(Cardinal)
library(yaml)
library(feather)


# Read the hyperparameters from the config file
config <- read_yaml("config.yaml")

# Path to the data and results
path1 <- sprintf("%s/%s", config$path_to_data, config$slide1)
path2 <- sprintf("%s/%s", config$path_to_data, config$slide2)

slide_files <- c(sprintf("%s/maldi/mse.imzML", path1),
                 sprintf("%s/maldi/mse.imzML", path2))

msa <- readMSIData(slide_files)

# Read the imzML + ibd object
mse1 <- readMSIData(sprintf("%s/maldi/mse.imzML", path1))
mse2 <- readMSIData(sprintf("%s/maldi/mse.imzML", path2))

# Compute the mean intensity of each m/z value
mse1 <- mse1 |>
  summarizeFeatures(stat = c(Mean = "mean"),
                    BPPARAM = MulticoreParam())
mse2 <- mse2 |>
  summarizeFeatures(stat = c(Mean = "mean"),
                    BPPARAM = MulticoreParam())

# Change the pixel data run factor to the slide name
pixelData(mse1)$run <- factor(config$slide1)
pixelData(mse2)$run <- factor(config$slide2)

# load the pixels data feather file
pixels1 <- read_feather(sprintf("%s/results/mse_pixels.feather", path1))
pixels2 <- read_feather(sprintf("%s/results/mse_pixels.feather", path2))

# Add the lesion and defects columns to the pixelData
pixelData(mse1)$Lesion <- pixels1$Density_Lesion > 0.5
pixelData(mse2)$Lesion <- pixels2$Density_Lesion > 0.5

pixelData(mse1)$Defects <- pixels1$Density_Defects > 0.5
pixelData(mse2)$Defects <- pixels2$Density_Defects > 0.5

# Drop the columns that already exist in the pixelData
pixels1 <- pixels1[, !colnames(pixels1) %in% c("x", "y", "run", "X3DPositionX", "X3DPositionY", "X3DPositionZ", "Density_Lesion", "Density_Defects")]
pixels2 <- pixels2[, !colnames(pixels2) %in% c("x", "y", "run", "X3DPositionX", "X3DPositionY", "X3DPositionZ", "Density_Lesion", "Density_Defects")]

# Add the pixels columns to the pixelData
for (column in colnames(pixels1)) {
  pixelData(mse1)[[column]] <- pixels1[[column]]
}

for (column in colnames(pixels2)) {
  pixelData(mse2)[[column]] <- pixels2[[column]]
}



# Subset the pixels out of the lesion
pid1 <- pixels(mse1, Lesion)
pid2 <- pixels(mse2, Lesion)

# Subset the features in the m/z range of interest
fid1 <- features(mse1, 500 < mz, mz < 3600)
fid2 <- features(mse2, 500 < mz, mz < 3600)

mse_region1 <- mse1[fid1, pid1]
mse_region2 <- mse2[fid2, pid2]

# Bin the m/z values to a resolution of 1 Da
mse_binned1 <- bin(mse_region1, resolution = 0.1, units = "mz")
mse_binned2 <- bin(mse_region2, resolution = 0.1, units = "mz")

# Compute the mean intensity of each m/z value
mse_region1 <- mse_region1 |>
  summarizeFeatures(stat = "mean")

# Normalize the data and reduce the baseline
mse_processed1 <- mse_region1 |>
  normalize(method = config$normalization_method) |>  # Normalize the data
  reduceBaseline(method = config$base_line_method) |>  # Reduce the baseline
  process()

mse_processed2 <- mse_region2 |>
  normalize(method = config$normalization_method) |>  # Normalize the data
  reduceBaseline(method = config$base_line_method) |>  # Reduce the baseline
  process()

mse_peaks <- mse_processed |>
  peakPick(method = config$peak_detection_method,
           SNR = config$signal_to_noise,
           tolerance = config$peak_pick_tolerance,
           units = config$units) |>
  process()

mse_aligned <- mse_peaks |>
  peakAlign()

spectra(mse_region)[1:6, 1:6]
spectra(mse_processed)[1:6, 1:6]


mse_peaks1 <- mse_processed1 |>
  peakPick(method = config$peak_detection_method,
           SNR = config$signal_to_noise,
           tolerance = config$peak_pick_tolerance,
           units = config$units)

mse_peaks2 <- mse_processed2 |>
  peakPick(method = config$peak_detection_method,
           SNR = config$signal_to_noise,
           tolerance = config$peak_pick_tolerance,
           units = config$units)
# Save the peaks data


# # Compute the mean intensity of each m/z value in the processed data
# mse_processed <- mse_processed |>
#   summarizeFeatures(stat = c(Mean_proc = "mean"),
#                     BPPARAM = MulticoreParam())

# Create the directory for the figures if it does not exist
dir.create(sprintf("%s/results", path), showWarnings = FALSE)

# Save the processed data
writeMSIData(mse_processed,
             file = sprintf("%s/results/mse_processed.imzML", path))


msa1 <- convertMSImagingExperiment2Arrays(mse1)
msa2 <- convertMSImagingExperiment2Arrays(mse2)

msa_all <- cbind(msa1, msa2)

# Adjust the figure size
options(repr.plot.width = 20, repr.plot.height = 5)
plot(mse_region1, i = 1000, xlim = c(580, 595))
plot(msa1, mean ~ mz)
plot(mse2, "mean", xlim = c(580, 595))

featureData(mse_region1)

plot(mse1, coord=list(x=16, y=16))

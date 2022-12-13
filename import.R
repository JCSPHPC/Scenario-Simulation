
# List the files and save to an R object
files <- list.files(pattern="testing_resource_.*.csv", full.names = TRUE)

# Extract the scenarios from the file list
scenarios <- sub('.*_(\\d+)_.*', '\\1', files)

# Number of scenarios
nscen <- length(unique(scenarios))

# Loop through the files, saving to a list of scenarios...
sc.data <- lapply(1:nscen, function(i) grep(paste0('.*_', i, '_.*'), files, value=TRUE))

# ... and import the data
data <- lapply(sc.data, function(x) lapply(x, read.csv))

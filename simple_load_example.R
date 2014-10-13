
#Load the data for Dog_1
library(R.matlab)
path <- c("e:/Temp/Dog_1/")

# "Positive" outcome data
filename <- paste0("Dog_1_preictal_segment_000", 1, ".mat") #1st segment of 6
pathname <- file.path(path, filename)    
data <- readMat(pathname)
df_PR <- data.frame(data[[1]][[1]][,])
dim(df_PR)

data_length_sec_PR <- data[[1]][[2]]
sampling_frequency_PR <- data[[1]][[3]]
channels_PR <- data[[1]][[4]]
sequence_PR <- data[[1]][[5]]


# "Negative" outcome data
filename <- paste0("Dog_1_interictal_segment_000", 1, ".mat") #1st segment of 6
pathname <- file.path(path, filename)    
data <- readMat(pathname)
df_IN <- data.frame(data[[1]][[1]][,])
dim(df_IN)

data_length_sec_IN <- data[[1]][[2]]
sampling_frequency_IN <- data[[1]][[3]]
channels_IN <- data[[1]][[4]]
sequence_IN <- data[[1]][[5]]

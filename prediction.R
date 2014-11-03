folder <- "Dog_1"
path <- paste0("j:/Temp/", folder)
fft_comp_n <- 100
N.of.clusters <- 6 # Number of clusters for a parralel execution

## All variables shoulbe be entered above _______________

library(foreach); library(doParallel)
workers <- makeCluster(N.of.clusters);  registerDoParallel(workers)

require(data.table)

MAT <- fread( paste0("fft_MAT_", folder, "_SegmentsSeparate.csv") , sep = ",") )

#Build models for each segment
models <- list(6 elements)
require(caret)
for(seg in 1:6){
    cat("working with segment ", seg)
    Out <- as.factor( MAT[Segment==seg, "Out"]  )    
    mod[1] <- train(MAT[Segment==seg], Out, method="rf")
}

#Make prediction
for(seg in 1:6){
    cat("Predicting with segment ", seg)
    Out <- as.factor( MAT[Segment==seg, "Out"]  )    
    mod[1] <- predict()
}

stopCluster(workers)
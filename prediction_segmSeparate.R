fft_comp_n <- 100 #number of Fourier components to use
N.of.clusters <- 3 # Number of clusters for a parralel execution
data_location <- "j:/Data_Epil/"

model.for.training <- "svmLinear"
#trGrid <- expand.grid(C=c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1, 1e1, 1e2, 1e3, 1e4))
trGrid <- expand.grid(C=c(1e-4, 1e-3, 1e-2, 1e-1, 1))
#trGrid <- expand.grid(C=c(1e-3))
#trGrid <- expand.grid(mtry=c(2, 80))

output.f <- paste0("./Output/FFTsmooth", fft_comp_n,"_segSep_", model.for.training, "/")
progress_file <- "./Output/training_progress.txt"

## All variables shoulbe be entered above _______________
#setup parallel execution
library(Revobase);   setMKLthreads(4)
library(foreach); library(doParallel)
workers <- makeCluster(N.of.clusters);  registerDoParallel(workers)

library(R.matlab)
library(data.table)
library(caret); library(pROC)

base::cat("", file=progress_file, append=F)

for(folder in dir(data_location)){
    #read current file 
    base::cat("reading", folder, "\n", file=progress_file, append=T)
    MAT <- fread( paste0("./MAT/FFTsmooth", fft_comp_n , "_MAT_", folder, "_SegmentsSeparate.csv") , sep = ",") 
        
    library(stringr)
    MAT <- MAT[, lapply(.SD, str_trim)]
    
    # choose data for test
    TestDT <- MAT[Type==c("test"), .SD, .SDcols=-c("Out", "Type", "Segment", "file.name") ]
    TestDT <- TestDT[, lapply(.SD, as.numeric)]
    
    Prediction <- matrix(data=NA, nrow=nrow(TestDT), ncol=7)
    
    #train model and make prediction for each segment separatelly
    for(seg in 1:6){
            
        TrainDT <- MAT[Type %in% c("preictal", "interictal") & Segment == seg,
                       .SD, .SDcols=-c("Out", "Type", "Segment", "file.name") ]
        
        TrainDT <- TrainDT[, lapply(.SD, as.numeric)]
        
            #name of a factor should not start with a number
        Y <- MAT[Type %in% c("preictal", "interictal") & Segment == seg , factor( paste0("P", Out) ) ] 
        
        #TrainDT <- TrainDT[1:50]
        #Y <- Y[1:50]
        
        # http://stats.stackexchange.com/questions/17602/caret-re-sampling-methods
        ctrl <- trainControl(method="boot632", number=25,
                             summaryFunction = twoClassSummary, 
                             classProbs = TRUE, 
                             verboseIter=T)
        
        
        base::cat("training", folder, " segment ", seg, "\n", file=progress_file, append=T)
        Mod <- train(TrainDT, Y, method=model.for.training,
                     tuneGrid=trGrid,
                     preProcess=c("center", "scale"), metric= "ROC", 
                     verbose = TRUE,
                     trControl = ctrl )
        
        capture.output(Mod, file=progress_file, append=T)
        #Save model if needed
        save(Mod, file=paste0(output.f, "Predicted_fft_Model_", fft_comp_n, "_", folder, "_seg", seg, ".model"))
            
        #Make prediction
        base::cat("prediction process...\n", file=progress_file, append=T)
        Prediction[, seg] <- as.numeric( predict(Mod, newdata=TestDT ) ) -1
    }    

    Prediction[, 7] <- as.numeric( rowSums(Prediction[, 1:6]) > 0 )
        
    Predicted_Data <- cbind(MAT[Type==c("test"), file.name], Prediction)
    write.table(Predicted_Data, file=paste0(output.f, "Predicted_fft_SegmSeprate_Full", fft_comp_n, "_", folder, ".csv"), quote=FALSE, 
                row.names=F, col.names = F, append=F, sep=",")    
    write.table(Predicted_Data[, c(1, 8)], file=paste0(output.f, "Predicted_fft_SegmSeprate_Final", fft_comp_n, "_", folder, ".csv"), quote=FALSE, 
                row.names=F, col.names = F, append=F, sep=",")    
}
stopCluster(workers)

system2("C:/cygwin64/bin/cat.exe", 
        args = c("./Output/Prediction_header.txt", 
                 paste0(output.f, grep("Predicted_fft_SegmSeprate_Final*", dir(output.f), value=T) ) ),
        stdout=paste0(output.f, "Together.csv") )
        
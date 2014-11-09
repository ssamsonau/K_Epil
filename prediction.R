folder <- "Dog_2"
path <- "j:/Temp/"
fft_comp_n <- 100
#N.of.clusters <- 2 # Number of clusters for a parralel execution

## All variables shoulbe be entered above _______________

library(Revobase);   setMKLthreads(2)

#library(foreach); library(doParallel)
#workers <- makeCluster(N.of.clusters);  registerDoParallel(workers)

library(data.table)

MAT <- fread( paste0("fft_MAT_", folder, "_SegmentsSeparate.csv") , sep = ",") 
library(stringr)
MAT <- MAT[, lapply(.SD, str_trim)]
setkey(MAT, Type)

TrainDT <- MAT[Type==c("preictal", "interictal") , .SD, .SDcols=-c("Out", "Type", "Segment") ]
TrainDT <- TrainDT[, lapply(.SD, as.numeric)]
    #name of factor should not start with a number
Y <- MAT[Type==c("preictal", "interictal") , factor( paste0("P", Out) ) ] 

#TrainDT <- TrainDT[1:50]
#Y <- Y[1:50]

library(caret); library(pROC)

    # cv does not make sense - for many folds there will be no possitive examples.
    # from experience cv does not give any reasonable results - Sens high, but Spec only 0     
    # bootstrapping ... I need to think... In any case boot632 is a better one
    # boot632 gives reasonable results
    
    # http://stats.stackexchange.com/questions/17602/caret-re-sampling-methods
ctrl <- trainControl(method="boot632", number=10,
                     summaryFunction = twoClassSummary, 
                     classProbs = TRUE, 
                     verboseIter=T)

trGrid <- expand.grid(C=c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1, 1e1, 1e2, 1e3, 1e4))
#trGrid <- expand.grid(mtry=c(2, 80))

Mod <- train(TrainDT, Y, method="svmLinear", 
#Mod <- train(TrainDT, Y, method="rf",
             tuneGrid=trGrid,
             preProcess=c("center", "scale"), metric= "ROC", 
             verbose = TRUE,
             trControl = ctrl )

#Save model if needed

#Make prediction
TestDT <- MAT[Type==c("test") , .SD, .SDcols=-c("Out", "Type", "Segment") ]
TestDT <- TestDT[, lapply(.SD, as.numeric)]

Prediction <- predict(Mod, newdata=TestDT )

files_names <- grep("test", dir(paste0(path, folder)), value=T)

#files_names <- files_names[1:50]

Predicted_Data <- cbind(files_names, Prediction)

write.table(Predicted_Data, file=paste0("./Predicted_", folder, ".txt"), quote=FALSE, 
            row.names=F, col.names = F, append=F, sep=",")

#stopCluster(workers)

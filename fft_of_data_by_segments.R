folder <- "Patient_1"
path <- paste0("e:/Temp/", folder)
fft_comp_n <- 100
#N.of.clusters <- 8 # Number of clusters for a parralel execution

## All variables shoulbe be entered above _______________

library(Revobase);   setMKLthreads(4)

#library(foreach); library(doParallel)
#workers <- makeCluster(N.of.clusters);  registerDoParallel(workers)

library(R.matlab)

ptm <- proc.time()

total_positive <- length(grep("preictal", dir(path))) 
total_negative <- length(grep("interictal", dir(path)) )
total_test <- length(grep("test", dir(path))) 

    # Read first file to determine number of sensors
filename <- paste0(folder, "_", "preictal", "_segment_", "0001", ".mat")
pathname <- file.path(path, filename)    
data_temp <- R.matlab::readMat(pathname)
mat_temp <- as.matrix(data_temp[[1]][[1]][,])  # the data
N.of.sensors <- nrow(mat_temp)

#Data table for input data
MAT <- data.frame( matrix(data = NA, nrow = total_positive+total_negative+total_test,
             ncol = (2*fft_comp_n) *N.of.sensors) )

MAT$Out <- c(rep(1, total_positive), rep(0, total_negative), rep(NA, total_test) )
MAT$Type <- "NA"
MAT$Segment <- NA

library(GeneCycle)

write("", file="./current_progress.txt", append=FALSE)

for(fullCase_cur_number in 1:(total_positive+total_negative+total_test)){
#foreach(fullCase_cur_number = 1:(total_positive+total_negative+total_test) ) %dopar% {
#foreach can not work with sideeffects. One should make a separate code to take this into account    
    ## FORM A FILE NAME TO READ
    if( fullCase_cur_number <= total_positive ){
        file_type <- "preictal"
        file_number <- fullCase_cur_number            
        MAT[fullCase_cur_number, "Out"] <- 1
    }
    else if( fullCase_cur_number <= total_positive+total_negative ){
        file_type <- "interictal"   
        file_number <- fullCase_cur_number-total_positive
        MAT[fullCase_cur_number, "Out"] <- 0
    }
    else{
        file_type <- "test"   
        file_number <- fullCase_cur_number-(total_positive+total_negative)
    }
    
    MAT[fullCase_cur_number, "Type"] <- file_type
        
    file_number_t <- paste0(file_number)
    
    if(file_number <= 9)
        file_number_t <- paste0("000", file_number_t)     
    else if(file_number <= 99)
        file_number_t <- paste0("00", file_number_t)
    else if(file_number <= 999)
        file_number_t <- paste0("0", file_number_t)            
    
    filename <- paste0(folder, "_", file_type, "_segment_", file_number_t, ".mat")

    ## READ FILE
    cat("Working with file", filename, "\n")
    write( paste("Working with file", filename, "\n"), "current_progress.txt", append=TRUE)
    pathname <- file.path(path, filename)    
    data_temp <- R.matlab::readMat(pathname)
    mat_temp <- as.matrix(data_temp[[1]][[1]][,])  # the data
    if(file_type %in% c("preictal", "interictal"))  
        # there is no inforamtion about segment for test data
        MAT[fullCase_cur_number, "Segment"] <- data_temp[[1]][[5]][,] #segment number
    
    ## MAKE FFT transformation for data from a given file
    for(r in 1:N.of.sensors){ # we have N.of.sensors rows in data - one row for every sensor
        pg <- GeneCycle::periodogram(mat_temp[r, ])
        s_ord <- order(pg$spec, decreasing = T)
        
        cur_MAT_index <- 1 # for a given sensor we fill 
                           #in the cur_MAT_index place
        
        for(k in 1:max(s_ord)){ # go through max amplitude values and choose
                                #only those wich are on top
            cur_max_number <- which(s_ord==k)
            
            if(cur_MAT_index > fft_comp_n) break 
                #stop loop when needed number of places filled in

            if(cur_max_number == 1 | cur_max_number == length(pg$spec)) next
                # do not include the point on edge
            
            # check if this is true maximum
            if( pg$spec[cur_max_number] > pg$spec[cur_max_number-1] &
                    pg$spec[cur_max_number] > pg$spec[cur_max_number+1] ){
                # if yes, save the data
                MAT[ fullCase_cur_number, (r-1)*(2*fft_comp_n) + cur_MAT_index] <-
                    pg$spec[cur_max_number] # modes amplitude
                
                MAT[ fullCase_cur_number, (r-1)*(2*fft_comp_n) + fft_comp_n + cur_MAT_index] <-
                    pg$freq[cur_max_number] # correspond. frequencies
            
                #print((r-1)*(2*fft_comp_n) + cur_MAT_index)
                #print((r-1)*(2*fft_comp_n) + fft_comp_n + cur_MAT_index)
                #print(cur_MAT_index)

                cur_MAT_index = cur_MAT_index + 1
            }
        }
    }
}


# Save the data
library(MASS)
write.matrix(MAT, file = paste0("fft_MAT_", folder, "_SegmentsSeparate.csv") , sep = ",")
# next step would be To take learning algorithm and apply to data

write( proc.time() - ptm, "current_progress.txt", append=TRUE)

#stopCluster(workers)

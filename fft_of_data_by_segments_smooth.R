fft_comp_n <- 100
#N.of.clusters <- 2 # Number of clusters for a parralel execution
data_location <- "j:/Data_Epil/"

## All variables shoulbe be entered above _______________
library(Revobase);   setMKLthreads(4)

#library(foreach); library(doParallel)
#workers <- makeCluster(N.of.clusters);  registerDoParallel(workers)

library(R.matlab)
library(data.table)

write("", file="./fft_progress.txt", append=FALSE)

for(folder in dir(data_location)){
#for(folder in c("Patient_2")){

    #folder <- "Dog_1"
    path <- paste0(data_location, folder)
    
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
                              ncol = 3*fft_comp_n *N.of.sensors) )
    
    MAT$Out <- c(rep(1, total_positive), rep(0, total_negative), rep(NA, total_test) )
    MAT$Type <- "NA"
    MAT$Segment <- NA
    MAT$file.name <- "NA"
    
    #library(GeneCycle)
    
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
        MAT[fullCase_cur_number, "file.name"] <- filename
        
        ## READ FILE
        cat("Working with file", filename, "\n")
        write( paste("Working with file", filename, "\n"), "fft_progress.txt", append=TRUE)
        pathname <- file.path(path, filename)    
        data_temp <- R.matlab::readMat(pathname)
        mat_temp <- as.matrix(data_temp[[1]][[1]][,])  # the data
        if(file_type %in% c("preictal", "interictal"))  
            # there is no inforamtion about segment for test data
            MAT[fullCase_cur_number, "Segment"] <- data_temp[[1]][[5]][,] #segment number
        
        ## MAKE FFT transformation for data from a given file
        for(r in 1:N.of.sensors){ # we have N.of.sensors rows in data - one row for every sensor
            #pg <- GeneCycle::periodogram(mat_temp[r, ])
            pg.Mod <- Mod(fft(mat_temp[r, ])) 
            pg.Mod <- pg.Mod[1:(length(pg.Mod)/2)] # we need only unique half
            
            pg.Freq <- (1:length(pg.Mod) ) /length(pg.Mod) #this is not the actual freq scale. But we don't need actual
            
            max_pos <- which(diff(sign(diff( pg.Mod )))==-2)+1
            
            #check
            #pg.Mod.t <- pg.Mod[1:30]
            #plot(pg.Mod.t)
            #max_pos.t <- which(diff(sign(diff( pg.Mod.t )))==-2)+1
            #points((1:29)[max_pos.t], rep(-2, length(max_pos.t)), col="red")
            #
            
            FFT.DT <- data.table(spec = pg.Mod[max_pos], freq = pg.Freq[max_pos])
            setkey(FFT.DT, freq)
            FFT.DT.L <- nrow(FFT.DT)
            
            FFT.DT$group <- (1:FFT.DT.L) %/% (FFT.DT.L/fft_comp_n)
            FFT.DT$group[FFT.DT.L] <- FFT.DT$group[FFT.DT.L-1]
            
#             myFunc <- function(s, f){
#                 #library(caret)
#                 #l.model <- lm(s ~ f)
#                 #r.slope <- l.model$coeff[2]
#                 #r.resid <- summary(l.model)$coeff[1, 2]
#                 m.spec <- mean(spec)
#                 variance <- var(spec)
#                 m.freq <- mean(f)
#                 return(list(m.spec, variance, m.freq))
#             }
            
            FFT.DT[, c("m.spec", "variance", "m.freq") := list(median(spec), var(spec), median(freq)), by=group]
            
                
            save.FFT <- FFT.DT[, lapply(.SD, unique), 
                               by=group, .SDcols = c("m.spec", "variance", "m.freq")]
            
            MAT[ fullCase_cur_number, ((r-1)*(3*fft_comp_n)+1):((r-1)*3*fft_comp_n + 3*fft_comp_n) ] <- 
                unlist(save.FFT[, .SD, .SDcols=-"group"])
            
        }
    }
       
    # Save the data
    library(MASS)
    write.matrix(MAT, file = paste0("./MAT/FFTsmooth", fft_comp_n, "_MAT_", folder, "_SegmentsSeparate.csv") , sep = ",")
    # next step would be To take learning algorithm and apply to data
    
    write( proc.time() - ptm, "fft_progress.txt", append=TRUE)
    cat(proc.time() - ptm)
}

#stopCluster(workers)
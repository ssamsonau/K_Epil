
#Load the data for Dog_1
library(R.matlab)
path <- c("e:/Temp/Dog_1/")

num_positive <- 4
num_negative <- 120 #120

fft_comp_n <- 30

    #Matrix for input data
MAT <- matrix(data = NA, nrow = num_positive+num_negative,
             ncol = (2*fft_comp_n) *16)
    # vector for responses
Out <- c(rep(1, num_positive), rep(0, num_negative) )

library(GeneCycle)

for(i in 1:(num_positive+num_negative)){
    for(s in 1:6){
        
        ###########################separation of "Positive" and "Negative"
        if(i<=num_positive){
            file_type <- "preictal"
            cur_number <- (i-1)*6 + s            
        }
        else{
            file_type <- "interictal"   
            cur_number <- (i-num_positive-1)*6 + s
        }
        
        ###########################Form file name
        cur_number_t <- paste0(cur_number)
        
        if(cur_number <= 9)
            cur_number_t <- paste0("000", cur_number_t)
        else if(cur_number <= 99)
            cur_number_t <- paste0("00", cur_number_t)
        else if(cur_number <= 999)
            cur_number_t <- paste0("0", cur_number_t)            
        
        filename <- paste0("Dog_1_", file_type, "_segment_", cur_number_t, ".mat")

        ###########################Read data
        cat("Working with file", filename, "\n")
        pathname <- file.path(path, filename)    
        data_temp <- readMat(pathname)
        mat_temp <- as.matrix(data_temp[[1]][[1]][,])
        
        
        ###########################FFT        
        for(r in 1:16){ # we have 16 rows in row data - one row for every sensor

            pg <- periodogram(mat_temp[r, ])
            s_ord <- order(pg$spec, decreasing = T)
            
            cur_MAT_index <- 1 # for a given segment we fill 
                               #in the cur_MAT_index place
            
            for(k in 1:max(s_ord)){ # go through max amplitude values and choose
                                    #only those wich are on top
                cur_max_number <- which(s_ord==k)
                
                if(cur_MAT_index > fft_comp_n) break 
                    #stop loop when needed number of places filled in

                if(cur_max_number == 1 | cur_max_number == length(pg$spec)) next
                    # do not include the point on edge
                
                # check if this is true maximum
                if( pg$spec[cur_max_number] - pg$spec[cur_max_number-1] > 0 &
                        pg$spec[cur_max_number] - pg$spec[cur_max_number+1] > 0 ){
                    # if yes, save the data
                    MAT[cur_number, (r-1)*(fft_comp_n) + cur_MAT_index] <-
                        pg$spec[cur_max_number] # modes amplitude
                    
                    MAT[cur_number, (r-1)*(fft_comp_n) + fft_comp_n + cur_MAT_index] <-
                        pg$freq[cur_max_number] # correspond. frequencies
                
                    #print((r-1)*(2*fft_comp_n) + cur_MAT_index)
                    #print((r-1)*(2*fft_comp_n) + fft_comp_n + cur_MAT_index)
                    #print(cur_MAT_index)

                    cur_MAT_index = cur_MAT_index + 1
                }
            }
        }
    }
}


# Save the data
write.matrix(MAT, file = "MAT.csv", sep = ",")

# next step would be To take learning algorithm and apply to data

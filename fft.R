
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
        for(f in 1:16){
            pg <- periodogram(mat_temp[f, ])
            s_ord <- order(pg$spec, decreasing = T)
            MAT[cur_number, (f-1)*(2*fft_comp_n) + seq(1,fft_comp_n)] <-
                pg$spec[s_ord<=fft_comp_n]
            MAT[cur_number, (f-1)*(2*fft_comp_n) + seq(fft_comp_n+1 ,2*fft_comp_n)] <-
                pg$freq[s_ord<=fft_comp_n]
        }
    }
}


# Save the data
write.matrix(MAT, file = "MAT.csv", sep = ",")

# next step would be To take learning algorithm and apply to data

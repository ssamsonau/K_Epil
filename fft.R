
#Load the data for Dog_1
library(R.matlab)
path <- c("e:/Temp/Dog_1/")

num_positive <- 4
num_negative <- 3 #120

fft_comp_n <- 2

    #Matrix for input data
MAT <- matrix(data = NA, nrow = num_positive+num_negative,
             ncol = (2*fft_comp_n) *16)
    # vector for responses
Out <- c(rep(1, num_positive), rep(0, num_negative) )

library(GeneCycle)

# "Positive" outcome data
for(i in 1:num_positive){
    for(s in 1:6){
        cur_number <- (i-1)*6 + s
        if (cur_number <= 9) filename <- 
            paste0("Dog_1_preictal_segment_000", cur_number, ".mat")
        else filename <- paste0("Dog_1_preictal_segment_00", cur_number, ".mat")        
        
        cat("Working with file", filename, "\n")
        pathname <- file.path(path, filename)    
        data_temp <- readMat(pathname)
        mat_temp <- as.matrix(data_temp[[1]][[1]][,])
        
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

# "Negative" outcome data
for(i in (num_positive+1):num_negative){
    for(s in 1:6){
        cur_number <- (i-1)*6 + s
        if (cur_number <= 9) filename <- 
            paste0("Dog_1_interictal_segment_000", cur_number, ".mat")
        else filename <- paste0("Dog_1_interictal_segment_00", cur_number, ".mat")        
        
        cat("Working with file", filename, "\n")
        pathname <- file.path(path, filename)    
        data_temp <- readMat(pathname)
        mat_temp <- as.matrix(data_temp[[1]][[1]][,])
        
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

# American Epilepsy Society Seizure Prediction Challenge

## Use of Fast Fourier transform and SVM for training

[Click-Through Rate Prediction](https://www.kaggle.com/c/seizure-prediction) hosted by [Kaggle](https://www.kaggle.com/)

The goal is to develop algorithm capable of predicting epilepsy seizure. 

Files:

* fft_of_data_by_segments.R - Read data and make Fast Fourier transform for each segment. Extract cetrain number of FFT components with heighest contribution. Save data in corresponding files
* fft_of_data_by_segments_smooth.R - Read data and make Fast Fourier transform for each segment. Smooth data obtained by FFT. Save data in corresponding files
* prediction_segmSeparate.R - based on data delivered by one of fft files make prediction for every segment.
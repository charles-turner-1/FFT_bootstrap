function [smoothed_data] = smooth_for_FFT(input_data,data_freq,max_freq_in_year)
    
    if nargin < 3
        max_freq_in_year = 0.8;   % Default maximum frequency to 0.8 years
    end
   

    times = 1:length(input_data); 
    detrended_input = detrend(input_data);
    end_time = times(end);
    fs = data_freq;
    t = 0:1/fs:(end_time-1)/fs;
    x = detrended_input;
    y = fft(x);


    r = max_freq_in_year*t(end); % range of frequencies we want to preserve

    rectangle = ones(size(y));
    rectangle(1:r+1) = 0;               % preserve low +ve frequencies
    y_half = ifft(y.*rectangle);   % +ve low-pass filtered signal
    rectangle(end-r+1:r+1) = 0;         % preserve low -ve frequencies
    y_rect = ifft(y.*rectangle);   % full low-pass filtered signal

    smoothed_data = real(input_data - y_rect);
    


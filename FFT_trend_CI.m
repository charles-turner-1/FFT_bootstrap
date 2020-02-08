function [trends_and_CIs] = FFT_trend_CI(input_data,data_freq)

warning('off','MATLAB:polyfit:PolyNotUnique')


if size(data_freq(:)) ~=1
    error('Data frequency must be a single value (ie. a scalar)')
end

if data_freq ~= 1 && data_freq~=12
    warning("Data frequency isn't 1 (yearly inputs) or 12 (monthly iputs). Something is probably wrong")
end


input_data = squeeze(input_data); % Make sure we don't have extraneous dimensions
if size(input_data,2) ~= 1
    input_data = input_data'; % Make sure input data is a column vector
end

detrended_var = detrend(input_data);

input_length = length(input_data);

gradient_dist = zeros(1000,input_length/data_freq);

disp('Bootstrapping input data');
for i = input_length/data_freq:-1:1
    X = detrended_var(1:data_freq*i);
    x = [(1:data_freq*i)/data_freq]';
    parfor j = 1:1000
        Y = fft(X); 
        rand_phases = pi  * (2*rand(data_freq*i,1) - 1);
        Z_real = abs(Y).* cos(rand_phases);
        Z_im = abs(Y) .* 1i .* sin(rand_phases);
        Z = Z_real + Z_im;
        shuffled_yvals = real(ifft(Z));
        % Now, for each series of random numbers, ie
        % rand_trend(i,j).random_dist, get the slope of that distribution
        p1 = polyfit(x,shuffled_yvals,1);
        gradient_dist(j,i) = p1(1);
    end
end

% Does the x variable used for polyfit (line 53) need to be random canth
% values in the range of 'real' Canth values in this time frame?
disp('Bootstrap complete, analysing output');

b = sort(gradient_dist,1);
conf_vals = zeros(input_length/data_freq,1);
parfor i = 1:input_length/data_freq
    [~,conf_vals(i)] = normfit(b(:,i)); 
end

confidence_intervals = horzcat(conf_vals, -1*conf_vals);

%

% So conf_vals is the gradient required for the trend to be statistically
% significant at 1 sigma, each entry corresponding to length of period trend
% calculated on in years, ie conf_vals(1) is 1 year etc etc.

trends = NaN(input_length/data_freq,1);
parfor i = 2:input_length/data_freq     % Leave first data point blank since no trend can be discrened from one point
    p1 = polyfit([(1:data_freq*i)./data_freq]',input_data(1:data_freq*i),1);
    trends(i) = p1(1);
end

trends_and_CIs.trends = trends;
trends_and_CIs.one_sigma = confidence_intervals;

warning('on','MATLAB:polyfit:PolyNotUnique')










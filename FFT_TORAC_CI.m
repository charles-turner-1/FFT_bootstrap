function [TORACSs_and_CIs] = FFT_TORAC_CI(input_temp_data,input_carbon_data,input_freq,data_source,output_freq)

warning('off','MATLAB:polyfit:PolyNotUnique')

switch data_source
	case 'NEMO'
		if size(input_freq(:)) ~=1
			error('Data frequency must be a single value (ie. a scalar)')
		end

		if input_freq ~= 1 && input_freq~=12
			warning("Data frequency isn't 1 (yearly inputs) or 12 (monthly inputs). Something is probably wrong")
        end
        
        fprintf("Data source '%s' recognised, continuing to bootstrap\n",data_source);
	case 'RAPID'
		if size(input_freq(:)) ~=1
			error('Data frequency must be a single value (ie. a scalar)')
		end
		
		if input_freq ~= 10 && input_freq ~= 365/10
			warning("Data frequency isn't 10 (given as daily) or 365/10 (given as yearly). Something is probably wrong")
		end

		input_freq(input_freq~=36.5) = 36.5; % Try to sanitize input by turning it into a yearly frequency if it isnt already.
        
        fprintf("Data source '%s' recognised, continuing to bootstrap\n",data_source);

	otherwise
		fprintf("Data Source '%s' not recognised. Check spelling or add to source switch\n",data_source); 
end

if ~isscalar(output_freq)
	error('Output Frequency must be scalar');
end

input_temp_data = squeeze(input_temp_data); % Make sure we don't have extraneous dimensions
if size(input_temp_data,2) ~= 1
    input_temp_data = input_temp_data'; % Make sure input data is a column vector
end

input_carbon_data = squeeze(input_carbon_data); % Make sure we don't have extraneous dimensions
if size(input_carbon_data,2) ~= 1
    input_carbon_data = input_carbon_data'; % Make sure input data is a column vector
end

detrended_var = detrend(input_temp_data);

input_length = length(input_temp_data);
output_length = floor(output_freq*input_length/input_freq);
gradient_dist = zeros(1000,output_length);

disp('Bootstrapping input data');
for i = output_length:-1:1
    X = detrended_var(1:floor(input_freq*i/output_freq));
    x = input_carbon_data(1:length(X));
    X = repmat(X,[1 1000]);
    Y = fft(X);                      
    rand_phases = pi  * (2*rand(size(X)) - 1); 
    Z_real = abs(Y).* cos(rand_phases);        
    Z_im = abs(Y) .* 1i .* sin(rand_phases);
    shuffled_yvals = real(ifft(Z_real + Z_im));
    x = repmat(x,[1 1000]);
    k = multpolyfit(x,shuffled_yvals,1);
    gradient_dist(:,i) = k(2,:);
end
% Does the x variable used for polyfit (line 53) need to be random canth
% values in the range of 'real' Canth values in this time frame?
disp('Bootstrap complete, analysing output');

conf_vals = std(gradient_dist)';
confidence_intervals = horzcat(conf_vals, -1*conf_vals);

%

% So conf_vals is the gradient required for the trend to be statistically
% significant at 1 sigma, each entry corresponding to length of period trend
% calculated on in years, ie conf_vals(1) is 1 year etc etc.


TORAC = NaN(output_length,1);
nloops = output_length;

if input_freq <= 1
    parfor i = 2:nloops    % Leave first data point blank since no trend can be discerned from one point
        warning('off','all');
        p1 = polyfit(input_carbon_data(1:input_freq*i/output_freq),input_temp_data(1:input_freq*i/output_freq),1);
        TORAC(i) = p1(1);
        warning('on','all');
    end
else
    parfor i = 1:nloops
        warning('off','all');
        p1 = polyfit(input_carbon_data(1:input_freq*i/output_freq),input_temp_data(1:input_freq*i/output_freq),1);
        TORAC(i) = p1(1);
        warning('on','all');
    end
end

TORACSs_and_CIs.TORAC = TORAC;
TORACSs_and_CIs.one_sigma = confidence_intervals;

end


function [c] = multpolyfit(x,y,n)
m = size(x,2);
c = zeros(n+1,m);
for k = 1:m
 M = repmat(x(:,k),1,n+1);
 M = bsxfun(@power,M,0:n);
 c(:,k) = M\y(:,k);
end
end









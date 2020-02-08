function [output] = determine_trend_emergence(input_struct,n_sigma,criteria,input_data_dates)
% This function will only take outputs from FFT_trend_CI.m, which are in a
% structure. You also will need to supply it with input data dates, so it
% can give a year back. If you don't give it dates, it will just give back
% indices, and years will need to be matched up manually.

% Criteria is a scalar value ranging from 1 to 4. Select a value
% to determine what criteria this function will use to determine trend
% emergence. More criteria can be added if necessary

% 1: Trend/TORAC is greater than n_sigma. This is the most basic possible 
% test of trend emergence. It simply finds the last time at which the trend
% is not significant, and defines trend emergence as the next time after 
% this time

% All further trend emergence criteria will also require that the trend is
% statistically significant at this level, and simply add extra tests to
% determine trend emergence. 

% 2: Trend is within 12.5% of its final value (at end of timeseries of
% trends/TORAC). 12.5% is an arbitrary choice, but seems pretty reasonable.
% It's also 1/8 which is kinda mathematically pleasing.

% 3: Trend is increasing/decreasing at a rate that is less than 1% of its
% current value.

% 4: We require the conditions of both 2 and 3.

% Current criteria don't seem very sophisticated (2,3,4), and I get the feeling
% like they are making trend emergence artificially great. I need to think
% of a more sophisticated way of determining whether a trend has actually
% emerged.

% 5: Perform a linear fit from each year to the end of the run. Then
% require that the gradient of trends/TORAC is changing at less than 10% of
% the final value

warning('off','MATLAB:polyfit:PolyNotUnique')


if find(round(input_data_dates) ~= input_data_dates)
    error('Input data dates must be years only, given as a list, ie. 1980:2099');
end

if length(n_sigma) ~= 1
    error('n_sigma must be a scalar');
end

if length(input_data_dates) ~= length(input_struct.TORAC)
   error("Dates don't match up with length of timeseries"); 
end

% For the sake of this function, rename TORAC trends so that it recognises
% it as a trend to emerge, if necessary.
try
    input_struct.trends = input_struct.TORAC;
    input_struct = rmfield(input_struct,TORAC);
catch
    % If the field was already named trends then theres no need to do
    % anything.
end

if nargin < 4
    input_data_dates = 1:length(input_struct.trends);
end

trends = input_struct.trends;
one_sigma = input_struct.one_sigma;


trend_minus_CI = abs(trends) - one_sigma(:,1).*n_sigma;

% This chunk belongs to all criteria, as we always require a statistically
% significant trend/TORAC
trend_minus_CI(trend_minus_CI< 0) = NaN;
non_sig_vals = find(isnan(trend_minus_CI));
if non_sig_vals
    last_non_sig = non_sig_vals(end);
    emergence_index = last_non_sig + 2;
else
    warning('No insignificant values found')
    output = -1;
end

if criteria > 1
    emergence_index_min = emergence_index;
    clear emergence_index
end

% Following block: we require that we have a trend/TORAC within 5% of its
% final (estimated) value.
if criteria == 2
    final_trend_val = trends(end);
    frac_diff = (trends - repmat(final_trend_val, size(trends)))./final_trend_val;
    frac_diff(abs(frac_diff) > 0.125) = NaN;
    
    non_sig_vals = find(isnan(frac_diff));
    last_non_sig = non_sig_vals(end);
    emergence_index = last_non_sig + 2;
    
    if emergence_index < emergence_index_min
        emergence_index = emergence_index_min; % If it's stable before significant, take when significant, not when stable
    end
end

% Following block: we require that we have a trend/TORAC that is stable to
% within 1% of its current value

if criteria == 3
    trend_grad = gradient(trends);
    deriv_frac = trend_grad./trends;
    deriv_frac(abs(deriv_frac)>0.01) = NaN;
    
    non_sig_vals = find(isnan(deriv_frac));
    last_non_sig = non_sig_vals(end);
    emergence_index = last_non_sig + 2;
    
    if emergence_index < emergence_index_min
        emergence_index = emergence_index_min; % If it's stable before significant, take when significant, not when stable
    end
end

if criteria == 4
    trend_grad = gradient(trends);
    deriv_frac = trend_grad./trends;
    deriv_frac(abs(deriv_frac)>0.01) = NaN;

    final_trend_val = trends(end);
    frac_diff = (trends - repmat(final_trend_val, size(trends)))./final_trend_val;
    frac_diff(abs(frac_diff) > 0.125) = NaN;

    non_sig_vals = find(isnan(deriv_frac) | isnan(frac_diff));
    last_non_sig = non_sig_vals(end);
    emergence_index = last_non_sig + 2;
    
    if emergence_index < emergence_index_min
        emergence_index = emergence_index_min; % If it's stable before significant, take when significant, not when stable
    end
end

if criteria == 5
    grads = NaN(length(trends)-1,1);
    for i = 3:length(grads) %Determine the drift of the trend. Trends only 
        % exist from 2nd value, so start fitting at 3rd time step
        p1 = polyfit([2:i]',trends(2:i),1);
        grads(i) = p1(1);
    end
    drift = abs(grads/trends(end));
    drift(drift > 0.1) = NaN;
        
    non_sig_vals = find(isnan(drift));
    last_non_sig = non_sig_vals(end);
    emergence_index = last_non_sig + 2;
    
    if emergence_index < emergence_index_min
        emergence_index = emergence_index_min; % If it's stable before significant, take when significant, not when stable
    end
    
end


try 
    output = input_data_dates(emergence_index);
catch
    if emergence_index == length(trends) - 1
        warning('No significant values found')
        output = Inf;
    else
        output = NaN;
    end
        
end

warning('on','MATLAB:polyfit:PolyNotUnique')

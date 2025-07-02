function [lambda] = sgs_compute_lambda(data, data_time, mo, varargin)

% Description:
% Computes the decorrelation timescale and lambda based on autolag of timeseries
% (1) Create x(t) and x(t+tau), 
% (2) Find tc = (auto)correlation < 1/e, 
% (3) lambda = 1/tc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% data = hourly time series of water nontidal residuals (can have NaNs)
% data_time = timeseries associated with the nontidal residual (datetime)
% mo = month for which we are computing lambda for (number 1 through 12)
% varagin{1} - 1 if last year of data_time is not complete and we don't want
%   it to be included in the lambda. For example, Jan data can be used for
%   the December lag calculations but not Jan decorrelation

% Outputs:
% lamda = 1/decorrelation timescale

% Provenance: 
% Affiliation: Physical Sciences Laboratory, NOAA and Cooperative Institue for Research in
% Environmental Sciences, CU Boulder
% Written: Jun 2025 by Paige Hovenga (contact paige.hovenga@noaa.gov)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare variables
lag = 0; %counter
lags = []; myautocorr = []; check = []; ii = [];

while isempty(ii);

    lag = lag + 1;
    lags = [lags lag];

    %loop through each year and construct data1 and data2 (lagged hourly)
    data1 = [];
    data2 = []; %pull 1hr lagged + 1hr extra from next month

    %find data1 and data2 which is lagged
    if varargin{1} == 1; %exclude the last year of data for incomplete years
        idx = find(month(data_time) == mo & year(data_time) < max(year(data_time)));
    else
        idx = find(month(data_time) == mo);
    end

    data1 = data(idx);
    data2 = data(idx+lag);

    % Remove NaNs found in Data1 then Data2
    remove = find(isnan(data1));
    data1(remove) = [];
    data2(remove) = [];
    remove = find(isnan(data2));
    data1(remove) = [];
    data2(remove) = [];

    %compute the correlation (this is essentially lagged auto)
    myautocorr(lag) = corr(data1, data2);
    check(lag) = mean(data1.*data2)/mean(data1.^2); %autocorrelation check

    %find the first instance for when autocorrelation goes to 1/e and linearly interpolate tau_c
    ii = min(find(myautocorr <= 1/exp(1)));
end

%if decorrelated after one time unit (happen for daily means)
if ii == 1;
    tau_c = interp1([1, myautocorr], [0, lags], 1/exp(1)); 
else %interpolate (linear but could be more accurate)
    tau_c = interp1(myautocorr(ii-1:ii), lags(ii-1:ii), 1/exp(1));
end
lambda = 1/tau_c;

% figure; plot(lags, myautocorr);
% hold on;
% yline(1/exp(1));
% xline(tau_c);
% ylabel('AutoCorr');
% xlabel('Tau (hours)');
% grid on;
% title(['month = ' mon_text]);

end
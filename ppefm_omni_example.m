%% ppefm_omni_example.m
% Clean example script: download OMNI 5-minute GSM By/Bz + Vsw and run PPEFM
%
% Requirements:
%   1) ppefm_time_domain.m on MATLAB path
%   2) Internet access (OMNIWeb)
%
% This script demonstrates how to call ppefm_time_domain.m and plot results.
%
% Data source:
%   https://omniweb.gsfc.nasa.gov/form/omni_min.html

clear;
clc;
close all;

%% User settings
start_time_utc = datetime(2024,5,11,0,0,0,'TimeZone','UTC');
end_time_utc   = datetime(2024,5,11,2,0,0,'TimeZone','UTC');
longitude_deg  = -105; % example longitude

% Optional: add path to ppefm_time_domain.m
% addpath('/Users/manojnair/m');

%% Build OMNIWeb URL (5-min data: By GSM, Bz GSM, Vsw)
start_str = datestr(start_time_utc, 'yyyymmddHH');
end_str   = datestr(end_time_utc, 'yyyymmddHH');

base_url = 'https://omniweb.gsfc.nasa.gov/cgi/nx1.cgi';
query = [ ...
    'activity=retrieve', ...
    '&res=5min', ...
    '&spacecraft=omni_5min', ...
    '&vars=17', ... % By GSM
    '&vars=18', ... % Bz GSM
    '&vars=21', ... % Vsw
    '&start_date=', start_str, ...
    '&end_date=', end_str, ...
    '&view=0' ...
];

full_url = [base_url, '?', query];
fprintf('Requesting OMNI data: %s\n', full_url);

%% Download and parse OMNI response
raw_text = webread(full_url);
lines = regexp(raw_text, '\r?\n', 'split');

rows = {};
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if ~isempty(line) && ~isempty(regexp(line, '^\d{4}\s+\d{1,3}\s+\d{1,2}\s+\d{1,2}', 'once'))
        rows{end+1} = line; 
    end
end

if isempty(rows)
    error('No data rows parsed from OMNI response.');
end

n = numel(rows);
data = nan(n, 7);
for i = 1:n
    cols = regexp(rows{i}, '\s+', 'split');
    if numel(cols) < 7
        continue;
    end
    data(i,1:7) = str2double(cols(1:7));
end

year = data(:,1);
doy  = data(:,2);
hour = data(:,3);
minu = data(:,4);
By_GSM = data(:,5);
Bz_GSM = data(:,6);
Vsw = data(:,7);

% Replace OMNI missing data flags with NaN
By_GSM(By_GSM == 9999.99) = NaN;
Bz_GSM(Bz_GSM == 9999.99) = NaN;
Vsw(Vsw == 99999.9) = NaN;

% Build time vector
time_utc = datetime(year,1,1,0,0,0,'TimeZone','UTC') + days(doy - 1) + hours(hour) + minutes(minu);

% Interpolate missing values (simple linear)
if any(isnan(Vsw)), Vsw = fillmissing(Vsw, 'linear'); end
if any(isnan(By_GSM)), By_GSM = fillmissing(By_GSM, 'linear'); end
if any(isnan(Bz_GSM)), Bz_GSM = fillmissing(Bz_GSM, 'linear'); end

%% Run PPEFM
[EEF, IEF_Ey, IEF_Ez] = ppefm_time_domain( ...
    Vsw, By_GSM, Bz_GSM, ...
    'cadence_minutes', 5, ...
    'longitude_deg', longitude_deg, ...
    'start_time_utc', time_utc(1), ...
    'apply_lt', true, ...
    'apply_delay', true, ...
    'delay_minutes', 17, ...
    'verbose', true);

%% Plot results
figure('Color','w','Position',[100 100 900 800]);

subplot(4,1,1);
plot(time_utc, Bz_GSM, 'b-', 'LineWidth', 1.2); hold on;
plot(time_utc, By_GSM, 'r--', 'LineWidth', 1.2);
ylabel('IMF [nT]');
legend('Bz GSM','By GSM','Location','best');
title('OMNI GSM IMF and PPEFM Output');
grid on;
yline(0, 'k:', 'LineWidth', 0.5);

subplot(4,1,2);
plot(time_utc, Vsw, 'k-', 'LineWidth', 1.2);
ylabel('V_{sw} [km/s]');
grid on;

subplot(4,1,3);
plot(time_utc, IEF_Ey, 'b-', 'LineWidth', 1.2); hold on;
plot(time_utc, IEF_Ez, 'r--', 'LineWidth', 1.0);
ylabel('IEF [mV/m]');
legend('IEF_{Ey}','IEF_{Ez}','Location','best');
grid on;
yline(0, 'k:', 'LineWidth', 0.5);

subplot(4,1,4);
plot(time_utc, EEF, 'Color', [0.2 0.6 0.2], 'LineWidth', 1.6);
xlabel('Time [UTC]');
ylabel('EEF [mV/m]');
title('Equatorial Ionospheric Electric Field (PPEFM)');
grid on;
yline(0, 'k:', 'LineWidth', 0.5);

addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

% Load conditions and filter for 'Desync'
condition_mat = load('D:/src/conditions.mat');
condition_data = condition_mat.data;
filtered_cond_rows = condition_data(strcmp(condition_data.condition, "Desync"), :);
file_list = filtered_cond_rows.orig_file;

% Run analysis on each file
oscillatory_data = cell(1, numel(file_list));

for i = 1:numel(file_list)
    filename = fullfile(DIR, file_list{i});
    load(filename, 'reref');

    data = reref(3); 
    data.label = reref(3).label; 

    % FOOOF aperiodic
    cfg               = [];
    cfg.foilim        = [1 30];
    cfg.pad           = 4;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    fractal = ft_freqanalysis(cfg, data);

    % Original power
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, data);

    % Oscillatory power
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';
    oscillatory = ft_math(cfg, fractal, original);

    oscillatory_data{i} = oscillatory;
end

% Grand average across sessions
cfg = [];
cfg.parameter = 'powspctrm';
avg_oscillatory = ft_freqgrandaverage(cfg, oscillatory_data{:});

% Extract frequency and power
freq_oscillatory = avg_oscillatory.freq;
pow  = avg_oscillatory.powspctrm;
chan_labels = avg_oscillatory.label;
nchan = size(pow, 1);

% Alpha band
alpha_band = [7 13];
alpha_idx = find(freq_oscillatory >= alpha_band(1) & freq_oscillatory <= alpha_band(2));
alpha_freqs = freq_oscillatory(alpha_idx);

alpha_peak_freqs = zeros(nchan, 1);
alpha_peak_vals = zeros(nchan, 1);
alpha_mean_pow = zeros(nchan, 1);

for c = 1:nchan
    alpha_power = pow(c, alpha_idx);
    [peak_val, peak_idx] = max(alpha_power);
    alpha_peak_freqs(c) = alpha_freqs(peak_idx);
    alpha_peak_vals(c) = peak_val;
    alpha_mean_pow(c) = mean(alpha_power);
end

[~, sorted_idx] = sort(alpha_mean_pow, 'descend');
max_chan_idx = sorted_idx(1);
max_alpha_channel = chan_labels{max_chan_idx};
fprintf('\nTop alpha channel: %s\n', max_alpha_channel);

% Define alpha band from that top channel
MAX_ALPHA_RANGE = 3;
THETA_WIDTH     = 3;

alpha_power = avg_oscillatory.powspctrm(ismember(avg_oscillatory.label, max_alpha_channel), :);
[peaks, locs] = findpeaks(alpha_power, freq_oscillatory, 'MinPeakProminence', 1, 'SortStr', 'descend');

peaks = peaks(1);
locs = locs(1);
peak_tenth = 0.1 * peaks;

alpha_range = find(alpha_power > peak_tenth & ...
    freq_oscillatory >= locs - MAX_ALPHA_RANGE & ...
    freq_oscillatory <= locs + MAX_ALPHA_RANGE);

lower_alpha = freq_oscillatory(alpha_range(1));
upper_alpha = freq_oscillatory(alpha_range(end));
alpha_band_def = [lower_alpha, upper_alpha];

% Theta band
lower_theta = lower_alpha - THETA_WIDTH;
theta_band_def = [lower_theta, lower_alpha];
fprintf('Theta band: %.2f to %.2f Hz\n', theta_band_def(1), theta_band_def(2));

% Theta power
theta_idx = find(freq_oscillatory >= theta_band_def(1) & freq_oscillatory <= theta_band_def(2));
theta_freqs = freq_oscillatory(theta_idx);

theta_peak_freqs = zeros(nchan, 1);
theta_peak_vals = zeros(nchan, 1);
theta_mean_pow = zeros(nchan, 1);

for c = 1:nchan
    theta_power = pow(c, theta_idx);
    [peak_val, peak_idx] = max(theta_power);
    theta_peak_freqs(c) = theta_freqs(peak_idx);
    theta_peak_vals(c) = peak_val;
    theta_mean_pow(c) = mean(theta_power);
end

% Determine Top 5 indices for alpha and theta
[~, sorted_alpha_idx] = sort(alpha_mean_pow, 'descend');
[~, sorted_theta_idx] = sort(theta_mean_pow, 'descend');

% Initialize logical columns
Top5Alpha = false(nchan,1);
Top5Theta = false(nchan,1);
Top5Alpha(sorted_alpha_idx(1:5)) = true;
Top5Theta(sorted_theta_idx(1:5)) = true;

% Create final results table
results = table(chan_labels(:), ...
    alpha_mean_pow, alpha_peak_vals, alpha_peak_freqs, ...
    repmat(lower_alpha, nchan, 1), repmat(upper_alpha, nchan, 1), ...
    theta_mean_pow, theta_peak_vals, theta_peak_freqs, ...
    repmat(theta_band_def(1), nchan, 1), repmat(theta_band_def(2), nchan, 1), ...
    Top5Alpha, Top5Theta, ...
    'VariableNames', {'Channel', ...
    'AlphaMeanPower', 'AlphaPeakPower', 'AlphaPeakFreq', ...
    'AlphaLower', 'AlphaUpper', ...
    'ThetaMeanPower', 'ThetaPeakPower', 'ThetaPeakFreq', ...
    'ThetaLower', 'ThetaUpper', ...
    'Top5Alpha', 'Top5Theta'});

% Save to single Excel file
writetable(results, 'D:/post_desync_grandavg.xlsx');
fprintf('Results table saved to D:/post_desync_grandavg.xlsx\n');

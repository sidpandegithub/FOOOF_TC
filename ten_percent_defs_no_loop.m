addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

FILES = {
    '20190717_AMLEO_S1.mat';
    '20190723_AMLEO_S2.mat';
    '20190726_AMLEO_S3.mat'
};

oscillatory_data = cell(1, numel(FILES));

for i = 1:numel(FILES)
    filename = fullfile(DIR, FILES{i});
    load(filename, 'reref');

    data = reref(1); 
    data.label = reref(1).label; 
   
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

% Find top alpha channels, peak search between 7-13Hz
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

% Define alpha peak range from top alpha channel
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

figure;
plot(freq_oscillatory, alpha_power, 'k', 'LineWidth', 1.5);
hold on;
plot(freq_oscillatory(alpha_range), alpha_power(alpha_range), 'r', 'LineWidth', 1.5);
title(sprintf('Alpha band peak range - %s', max_alpha_channel));
xlabel('Frequency (Hz)');
ylabel('Power');
grid on;

lower_alpha = freq_oscillatory(alpha_range(1));
upper_alpha = freq_oscillatory(alpha_range(end));
alpha_band_def = [lower_alpha, upper_alpha];
% Define theta band based on alpha range
lower_theta = lower_alpha - THETA_WIDTH;
theta_band_def = [lower_theta, lower_alpha];
fprintf('Theta band: %.2f to %.2f Hz\n', theta_band_def(1), theta_band_def(2));

% Find top theta channels based on this theta band
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

[~, sorted_theta_idx] = sort(theta_mean_pow, 'descend');
fprintf('\nTop 5 channels by theta power (%.2f–%.2f Hz):\n', theta_band_def(1), theta_band_def(2));
for i = 1:5
    idx = sorted_theta_idx(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz\n', ...
        i, chan_labels{idx}, theta_mean_pow(idx), theta_peak_vals(idx), theta_peak_freqs(idx));
end
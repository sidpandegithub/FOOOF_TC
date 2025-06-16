addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\12-reref';
FILE= '20210420_JATSO_S1_reref.mat';
filename = fullfile(DIR, FILE);
load(filename, 'reref');

data = reref(1); 

% FOOOF aperiodic component
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

% Oscillatory power = Original / Aperiodic
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2./x1';
oscillatory = ft_math(cfg, fractal, original);

% Mean across channels
mean_original   = mean(original.powspctrm, 1);
mean_fractal    = mean(fractal.powspctrm, 1);
mean_oscillatory = mean(oscillatory.powspctrm, 1);

freq_oscillatory = oscillatory.freq;
pow = oscillatory.powspctrm;
freq_original = original.freq;
chan_labels = data.label;  % FIXED: Use channel labels from the data
nchan = size(pow, 1);

% Parameters
MAX_ALPHA_RANGE = 3;
THETA_WIDTH     = 3;

% --- Alpha band stats
alpha_band = [7 13];
alpha_idx = find(freq_oscillatory >= alpha_band(1) & freq_oscillatory <= alpha_band(2));
alpha_freqs = freq_oscillatory(alpha_idx);

alpha_peak_freqs = zeros(nchan, 1);
alpha_peak_vals  = zeros(nchan, 1);
alpha_mean_pow   = zeros(nchan, 1);

for c = 1:nchan
    alpha_power = pow(c, alpha_idx);
    [peak_val, peak_idx] = max(alpha_power);
    alpha_peak_freqs(c) = alpha_freqs(peak_idx);
    alpha_peak_vals(c)  = peak_val;
    alpha_mean_pow(c)   = mean(alpha_power);
end

% --- Top alpha channel
[~, sorted_alpha_idx] = sort(alpha_mean_pow, 'descend');
max_chan_idx = sorted_alpha_idx(1);
max_alpha_channel = chan_labels{max_chan_idx};

alpha_power = pow(ismember(chan_labels, max_alpha_channel), :);
[peaks, locs] = findpeaks(alpha_power, freq_oscillatory, ...
    'MinPeakProminence', 1, 'SortStr', 'descend');

% Base name fallback
if ~exist('base_name', 'var')
    base_name = FILE;
end

if isempty(peaks)
    warning('No strong alpha peak in %s. Skipping.', base_name);
end

peak_val = peaks(1);
peak_freq = locs(1);
peak_tenth = 0.1 * peak_val;

alpha_range = find(alpha_power > peak_tenth & ...
    freq_oscillatory >= peak_freq - MAX_ALPHA_RANGE & ...
    freq_oscillatory <= peak_freq + MAX_ALPHA_RANGE);

lower_alpha = freq_oscillatory(alpha_range(1));
upper_alpha = freq_oscillatory(alpha_range(end));
alpha_band_def = [lower_alpha, upper_alpha];

% --- Theta band from alpha edge
lower_theta = lower_alpha - THETA_WIDTH;
theta_band_def = [lower_theta, lower_alpha];

% --- Theta stats
theta_idx = find(freq_oscillatory >= theta_band_def(1) & freq_oscillatory <= theta_band_def(2));
theta_freqs = freq_oscillatory(theta_idx);

theta_peak_freqs = zeros(nchan, 1);
theta_peak_vals  = zeros(nchan, 1);
theta_mean_pow   = zeros(nchan, 1);

for c = 1:nchan
    theta_power = pow(c, theta_idx);
    [peak_val, peak_idx] = max(theta_power);
    theta_peak_freqs(c) = theta_freqs(peak_idx);
    theta_peak_vals(c)  = peak_val;
    theta_mean_pow(c)   = mean(theta_power);
end

[~, sorted_theta_idx] = sort(theta_mean_pow, 'descend');
Top5Alpha = false(nchan,1);
Top5Theta = false(nchan,1);
Top5Alpha(sorted_alpha_idx(1:5)) = true;
Top5Theta(sorted_theta_idx(1:5)) = true;

% Plot
figure;
semilogy(freq_original, mean_original, 'k-', 'LineWidth', 2); hold on;
semilogy(freq_original, mean_fractal, 'r--', 'LineWidth', 2);
semilogy(freq_oscillatory, mean_oscillatory, 'b-', 'LineWidth', 2);
legend({'Original', 'Fractal (Aperiodic)', 'Oscillatory'}, 'Location', 'northeast');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Power Spectra: Original vs Fractal vs Oscillatory');
grid on;

% Layout
topo_data = [];
topo_data.label = chan_labels;
topo_data.dimord = 'chan_time';
topo_data.time = 0;

cfg = [];
cfg.layout = 'easycapM1.mat';  
layout = ft_prepare_layout(cfg, topo_data);

% Alpha Topoplot
topo_data.avg = alpha_mean_pow(:)';  % make row vector
cfg = [];
cfg.layout = layout;
cfg.marker = 'on';
cfg.comment = 'no';
cfg.colorbar = 'yes';
cfg.zlim = [0 max(abs(topo_data.avg(:)))];
cfg.colormap = jet;
figure;
ft_topoplotER(cfg, topo_data);

% Theta Topoplot
topo_data.avg = theta_mean_pow(:)';  % make row vector
cfg = [];
cfg.layout = layout;
cfg.marker = 'on';
cfg.comment = 'no';
cfg.colorbar = 'yes';
cfg.zlim = [0 max(abs(topo_data.avg(:)))];
cfg.colormap = jet;

figure;
ft_topoplotER(cfg, topo_data);


% --- Top 5 Alpha Channels

% Find channel with max alpha power
[~, max_chan_idx] = max(alpha_mean_pow);
max_alpha_channel = chan_labels{max_chan_idx};

fprintf('Channel with max alpha power: %s (%.2f µV^2/Hz at %.2f Hz)\n', ...
    max_alpha_channel, alpha_peak_vals(max_chan_idx), alpha_peak_freqs(max_chan_idx));

% Sort and display top N alpha channels
[~, sorted_idx] = sort(alpha_mean_pow, 'descend');
topN = 5;
fprintf('\nTop %d channels by alpha power:\n', topN);
for i = 1:topN
    idx = sorted_idx(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz\n', ...
        i, chan_labels{idx}, alpha_mean_pow(idx), alpha_peak_vals(idx), alpha_peak_freqs(idx));
end

% Plot top 5 alpha channels
figure;
hold on;
colors = lines(topN);
for i = 1:topN
    idx = sorted_idx(i);
    semilogy(freq_oscillatory, pow(idx, :), ...
        'LineWidth', 2, 'Color', colors(i, :));
end
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Full Spectrum for Top 5 Channels by Alpha Power (7–13 Hz)');
legend(chan_labels(sorted_idx(1:topN)), 'Location', 'northeast');
grid on;
xticks(1:1:35);



% --- Top 5 Theta Channels

% Find channel with max theta power
[~, max_chan_idx_theta] = max(theta_mean_pow);
max_theta_channel = chan_labels{max_chan_idx_theta};

fprintf('\nChannel with max theta power: %s (%.2f µV^2/Hz at %.2f Hz)\n', ...
    max_theta_channel, theta_peak_vals(max_chan_idx_theta), theta_peak_freqs(max_chan_idx_theta));

% Sort and display top N theta channels
[~, sorted_idx_theta] = sort(theta_mean_pow, 'descend');
topN = 5;
fprintf('\nTop %d channels by theta power:\n', topN);
for i = 1:topN
    idx = sorted_idx_theta(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz\n', ...
        i, chan_labels{idx}, theta_mean_pow(idx), theta_peak_vals(idx), theta_peak_freqs(idx));
end

% Plot top 5 theta channels
figure;
hold on;
colors = lines(topN);
for i = 1:topN
    idx = sorted_idx_theta(i);
    semilogy(freq_oscillatory, pow(idx, :), ...
        'LineWidth', 2, 'Color', colors(i, :));
end
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Full Spectrum for Top 5 Channels by Theta Power');
legend(chan_labels(sorted_idx_theta(1:topN)), 'Location', 'northeast');
grid on;
xticks(1:1:35);

addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

FILES = {
    '20210420_JATSO_S1.mat';
    '20210426_JATSO_S2.mat';
    '20210503_JATSO_S3.mat'
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

% Extract info
freq_oscillatory = avg_oscillatory.freq;
pow = avg_oscillatory.powspctrm;
chan_labels = avg_oscillatory.label;
nchan = size(pow, 1);

%% --- Alpha stats
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

[~, sorted_alpha_idx] = sort(alpha_mean_pow, 'descend');
max_chan_idx = sorted_alpha_idx(1);
max_alpha_channel = chan_labels{max_chan_idx};

fprintf('\nTop alpha channel: %s\n', max_alpha_channel);

% Define alpha peak range
MAX_ALPHA_RANGE = 3;
THETA_WIDTH     = 3;

alpha_power = pow(max_chan_idx, :);
[peaks, locs] = findpeaks(alpha_power, freq_oscillatory, ...
    'MinPeakProminence', 1, 'SortStr', 'descend');
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
lower_theta = lower_alpha - THETA_WIDTH;
theta_band_def = [lower_theta, lower_alpha];
fprintf('Theta band: %.2f to %.2f Hz\n', theta_band_def(1), theta_band_def(2));

%% --- Theta stats
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

%% --- Label Top 5
Top5Alpha = false(nchan,1);
Top5Theta = false(nchan,1);
Top5Alpha(sorted_alpha_idx(1:5)) = true;
Top5Theta(sorted_theta_idx(1:5)) = true;

%% --- Plot spectrum (average across channels)
mean_oscillatory = mean(pow, 1);
figure;
semilogy(freq_oscillatory, mean_oscillatory, 'b-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Power Spectrum (Oscillatory, Grand Average)');
grid on;

%% --- Topoplot layout setup
topo_data = [];
topo_data.label = chan_labels;
topo_data.dimord = 'chan_time';
topo_data.time = 0;

cfg = [];
cfg.layout = 'easycapM1.mat';  
layout = ft_prepare_layout(cfg, topo_data);

%% --- Alpha Topoplot
topo_data.avg = alpha_mean_pow(:)'; 
cfg = [];
cfg.layout = layout;
cfg.marker = 'on';
cfg.comment = 'no';
cfg.colorbar = 'yes';
cfg.zlim = [0 max(abs(topo_data.avg(:)))];
cfg.colormap = jet;
figure;
ft_topoplotER(cfg, topo_data);
title('Topography - Alpha Mean Power');

%% --- Theta Topoplot
topo_data.avg = theta_mean_pow(:)'; 
cfg = [];
cfg.layout = layout;
cfg.marker = 'on';
cfg.comment = 'no';
cfg.colorbar = 'yes';
cfg.zlim = [0 max(abs(topo_data.avg(:)))];
cfg.colormap = jet;
figure;
ft_topoplotER(cfg, topo_data);
title('Topography - Theta Mean Power');

%% --- Print and plot top 5 Alpha
fprintf('\nTop 5 channels by alpha power:\n');
for i = 1:5
    idx = sorted_alpha_idx(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz\n', ...
        i, chan_labels{idx}, alpha_mean_pow(idx), alpha_peak_vals(idx), alpha_peak_freqs(idx));
end

figure;
hold on;
colors = lines(5);
for i = 1:5
    idx = sorted_alpha_idx(i);
    semilogy(freq_oscillatory, pow(idx, :), ...
        'LineWidth', 2, 'Color', colors(i, :));
end
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Top 5 Channels by Alpha Power (7–13 Hz)');
legend(chan_labels(sorted_alpha_idx(1:5)), 'Location', 'northeast');
grid on;
xticks(1:1:35);

%% --- Print and plot top 5 Theta
fprintf('\nTop 5 channels by theta power:\n');
for i = 1:5
    idx = sorted_theta_idx(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz\n', ...
        i, chan_labels{idx}, theta_mean_pow(idx), theta_peak_vals(idx), theta_peak_freqs(idx));
end

figure;
hold on;
colors = lines(5);
for i = 1:5
    idx = sorted_theta_idx(i);
    semilogy(freq_oscillatory, pow(idx, :), ...
        'LineWidth', 2, 'Color', colors(i, :));
end
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Top 5 Channels by Theta Power');
legend(chan_labels(sorted_theta_idx(1:5)), 'Location', 'northeast');
grid on;
xticks(1:1:35);

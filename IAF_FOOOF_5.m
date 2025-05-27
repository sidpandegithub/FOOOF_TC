addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

FILES = {
    '20190715_ANSCH_S1.mat';
    '20190717_ANSCH_S2.mat';
    '20190723_ANSCH_S3.mat'
};

oscillatory_data = cell(1, numel(FILES));

for i = 1:numel(FILES)
    filename = fullfile(DIR, FILES{i});
    load(filename, 'reref');
    data = reref(1);  % "PRE RS-EEG"
    
    % FOOOF aperiodic component
    cfg               = [];
    cfg.foilim        = [1 35];
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

    oscillatory_data{i} = oscillatory;
end

% Grand average across all sessions
cfg = [];
cfg.parameter = 'powspctrm';
avg_oscillatory = ft_freqgrandaverage(cfg, oscillatory_data{:});

% Extract frequency and power
freq = avg_oscillatory.freq;
pow  = avg_oscillatory.powspctrm;

% Plot average spectrum across channels
% Get frequency and average power across channels
freq = original.freq;

% Mean across channels
mean_original   = mean(original.powspctrm, 1);
mean_fractal    = mean(fractal.powspctrm, 1);
mean_oscillatory = mean(oscillatory.powspctrm, 1);

% Plot
figure;
plot(freq, mean_original, 'k-', 'LineWidth', 2); hold on;
plot(freq, mean_fractal, 'r--', 'LineWidth', 2);
plot(freq, mean_oscillatory, 'b-', 'LineWidth', 2);
legend({'Original', 'Fractal (Aperiodic)', 'Oscillatory'}, 'Location', 'northeast');

xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Power Spectra: Original vs Fractal vs Oscillatory');
xlim([1 35]);
grid on;


% Define alpha band
alpha_band = [7 13];
alpha_idx = find(freq >= alpha_band(1) & freq <= alpha_band(2));
alpha_freqs = freq(alpha_idx);

% Initialize arrays
nchan = size(pow, 1);
alpha_peak_freqs = zeros(nchan, 1);
alpha_peak_vals = zeros(nchan, 1);
alpha_mean_pow = zeros(nchan, 1);

% Loop over channels
for c = 1:nchan
    alpha_power = pow(c, alpha_idx);
    
    % Find peak within alpha band
    [peak_val, peak_idx] = max(alpha_power);
    alpha_peak_freqs(c) = alpha_freqs(peak_idx);
    alpha_peak_vals(c) = peak_val;

    % Mean alpha power
    alpha_mean_pow(c) = mean(alpha_power);
end

% Get channel names
chan_labels = avg_oscillatory.label;

% Find channel(s) with most alpha power
[~, max_chan_idx] = max(alpha_mean_pow);
max_alpha_channel = chan_labels{max_chan_idx};

fprintf('Channel with max alpha power: %s (%.2f µV^2/Hz at %.2f Hz)\n', ...
    max_alpha_channel, alpha_peak_vals(max_chan_idx), alpha_peak_freqs(max_chan_idx));

% Sort and display top N
[~, sorted_idx] = sort(alpha_mean_pow, 'descend');
topN = 5;
fprintf('\nTop %d channels by alpha power:\n', topN);
for i = 1:topN
    idx = sorted_idx(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz\n', ...
        i, chan_labels{idx}, alpha_mean_pow(idx), alpha_peak_vals(idx), alpha_peak_freqs(idx));
end

%% Bar Plot of Mean Alpha Power
figure;
bar(alpha_mean_pow(sorted_idx));
set(gca, 'XTick', 1:nchan, 'XTickLabel', chan_labels(sorted_idx));
xtickangle(45);
ylabel('Mean Alpha Power (7–13 Hz)');
title('Alpha Power by Channel (Sorted, Post-FOOOF)');
grid on;
%% 

% Build data structure for topoplot
topo_data = [];
topo_data.label = chan_labels;
topo_data.dimord = 'chan_time';
topo_data.time = 0;
topo_data.avg = alpha_mean_pow(:)';  % This is correct structure for chan_time

% Load layout
cfg = [];
cfg.layout = 'easycapM1.mat';  % Explicitly use your layout file
layout = ft_prepare_layout(cfg, topo_data);

% Plot topoplot
cfg = [];
cfg.layout = layout;
cfg.marker = 'on';
cfg.comment = 'no';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';  % or set manually like [min max]
cfg.colormap = parula;

figure;
ft_topoplotER(cfg, topo_data);
title('Topographical Distribution of Mean Alpha Power (7–13 Hz)');

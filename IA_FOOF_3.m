addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

FILES = {
    '20190716_LIGOH_S1.mat';
    '20190719_LIGOH_S2.mat';
    '20190724_LIGOH_S3.mat'
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

% Extract data
frequencies = avg_oscillatory.freq;
power_matrix = avg_oscillatory.powspctrm;  % channels x freqs
channel_labels = avg_oscillatory.label;

figure;
plot(frequencies,power_matrix);

% Define alpha band
alpha_idx = find(frequencies >= 8 & frequencies <= 13);

channels_alpha = channel_labels(alpha_idx);

% Mean alpha power per channel
mean_alpha_power = mean(power_matrix(:, alpha_idx), 2);


% Sort channels by alpha power
[sorted_power, sorted_idx] = sort(mean_alpha_power, 'descend');
sorted_channels = channel_labels(sorted_idx);

% Display top N channels
top_N = 10;
fprintf('\nTop %d channels with highest alpha power (8–13 Hz):\n', top_N);
for i = 1:min(top_N, numel(sorted_channels))
    fprintf('%2d. %s: %.3f\n', i, sorted_channels{i}, sorted_power(i));
end

% Optional: Plot all channels' alpha power
figure;
bar(mean_alpha_power(sorted_idx));
set(gca, 'XTickLabel', sorted_channels, 'XTick', 1:numel(sorted_channels));
xtickangle(45);
ylabel('Mean Alpha Power (8–13 Hz)');
title('Alpha Power by Channel (Post-FOOOF)');
grid on;

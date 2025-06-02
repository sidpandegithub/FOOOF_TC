clear all;
addpath('C:\fieldtrip-20250106');
% addpath('C:\fieldtrip-20240113');
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
% 
%     % Keep only the first 64 EEG channels in the labels
%     reref(1).label = NEWEST_LABELS(1:64);  
% 
    % Loop through each trial and slice out the first 64 channels
 % for j = 1:numel(data(1).trial)
 %        data(1).trial{j} = data(1).trial{j}(1:64, :);  
 % end
% 
    data = reref(1); 
    data.label=reref(1).label; 
   
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

    oscillatory_data{i} = oscillatory;
end


% Grand average across all sessions
cfg = [];
cfg.parameter = 'powspctrm';
avg_oscillatory = ft_freqgrandaverage(cfg, oscillatory_data{:});

% Extract frequency and power
freq_oscillatory = avg_oscillatory.freq;
pow  = avg_oscillatory.powspctrm;


% Get frequency and average power across channels
freq_original = original.freq;

% Mean across channels
mean_original   = mean(original.powspctrm, 1);
mean_fractal    = mean(fractal.powspctrm, 1);
mean_oscillatory = mean(oscillatory.powspctrm, 1);

% Define alpha band
alpha_band = [7 13];
alpha_idx = find(freq_oscillatory >= alpha_band(1) & freq_oscillatory <= alpha_band(2));
alpha_freqs = freq_oscillatory(alpha_idx);
alpha_pow_all_channels = pow(:, alpha_idx);

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

% Define ±2 Hz around peak alpha frequency of the max channel
peak_freq = alpha_peak_freqs(max_chan_idx);
individual_alpha_band = [peak_freq - 2, peak_freq + 2];

% Display IAF band numerically
fprintf('Individual Alpha Band for %s: %.2f Hz to %.2f Hz\n', ...
    max_alpha_channel, individual_alpha_band(1), individual_alpha_band(2));

% Find indices of frequencies in this personalized alpha band
indiv_idx = find(freq_oscillatory >= individual_alpha_band(1) & ...
                 freq_oscillatory <= individual_alpha_band(2));
indiv_freqs = freq_oscillatory(indiv_idx);

indiv_power = pow(max_chan_idx, indiv_idx);

% Plot
figure;
plot(indiv_freqs, indiv_power, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title(sprintf('Individual Alpha Band: %s (%.2f ± 2 Hz)', ...
    max_alpha_channel, peak_freq));
grid on;

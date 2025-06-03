clear all; close all;
addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';
image_DIR = 'E:\BDI\Twincoil\harddrive\12_IAF\post';
FILES = {
    '20190722_TAWEI_S1';
    '20190725_TAWEI_S2.mat';
    '20190730_TAWEI_S3.mat'
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

    all_data = [data.trial{:}];  % 68 x (2048*186)
    all_time = [data.time{:}];   % 1 x (2048*186)

    data_continuous = [];
    data_continuous.label = data.label;
    data_continuous.fsample = data.fsample;
    data_continuous.trial = {all_data};  
    data_continuous.time = {all_time};
    data_continuous.sampleinfo = [1, size(all_data,2)];


    % chunk into 2-second segments  
    cfg               = [];
    cfg.length        = 2;
    cfg.overlap       = 0.5;
    data_epoched       = ft_redefinetrial(cfg, data_continuous);
   
    % FOOOF aperiodic component
    cfg               = [];
    cfg.foilim        = [1 30];
    cfg.pad           = 4;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    fractal = ft_freqanalysis(cfg, data_epoched);

    % Original power
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, data_epoched);

    % Oscillatory power = Original / Aperiodic
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';   
    %cfg.operation     = 'x2-x1';
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
chan_labels = avg_oscillatory.label;


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
% alpha_peak_freqs = zeros(nchan, 1);
% alpha_peak_vals = zeros(nchan, 1);
% alpha_mean_pow = zeros(nchan, 1);
% Initialize arrays
alpha_peak_freqs = nan(nchan, 1);
alpha_peak_vals  = nan(nchan, 1);
alpha_mean_pow   = nan(nchan, 1);
alpha_peak_width = nan(nchan, 1);

% 初始化结构体数组
alpha_results = struct('channel', [], 'peak_freq', [], 'peak_power', [], 'mean_power', [], 'peak_width', [], 'freq_band', []);

for c = 1:nchan
    alpha_power = pow(c, alpha_idx); 
    f_alpha     = alpha_freqs;

    % prominence > 0.05 µV²/Hz）
    [peaks, locs, widths, proms] = findpeaks(alpha_power, f_alpha, 'MinPeakProminence', 0.05);

    if isempty(peaks)
        alpha_results(c).channel = chan_labels{c};
        alpha_results(c).peak_freq = NaN;
        alpha_results(c).peak_power = NaN;
        alpha_results(c).mean_power = mean(alpha_power);
        alpha_results(c).peak_width = NaN;
        alpha_results(c).freq_band = [NaN NaN];
        continue;
    end

    [~, max_idx] = max(proms);
    peak_freq = locs(max_idx);
    peak_val  = peaks(max_idx);
    peak_prom = proms(max_idx);


    peak_width = widths(max_idx);
    if isempty(peak_width) || peak_width > 6
        peak_width = 3;  % 宽度设为3Hz（±1.5Hz）
    end

    freq_band = [peak_freq - peak_width/2, peak_freq + peak_width/2];


    alpha_results(c).channel = chan_labels{c};
    alpha_results(c).peak_freq = peak_freq;
    alpha_results(c).peak_power = peak_val;
    alpha_results(c).mean_power = mean(alpha_power);
    alpha_results(c).peak_width = peak_width;
    alpha_results(c).freq_band = freq_band;
end


[~, max_chan_idx] = max([alpha_results.mean_power]);
max_alpha_channel = alpha_results(max_chan_idx).channel;

fprintf('Channel with max alpha power: %s (%.2f µV^2/Hz at %.2f Hz)\n', ...
    max_alpha_channel, alpha_results(max_chan_idx).peak_power, alpha_results(max_chan_idx).peak_freq);

topN = 5;  
[~, sorted_idx] = sort([alpha_results.mean_power], 'descend');

for i = 1:topN
    idx = sorted_idx(i);
    fprintf('%2d. %s | Mean: %.2f | Peak: %.2f at %.2f Hz | Width: %.2f Hz | Band: %.2f-%.2f Hz\n', ...
        i, alpha_results(idx).channel, alpha_results(idx).mean_power, ...
        alpha_results(idx).peak_power, alpha_results(idx).peak_freq, ...
        alpha_results(idx).peak_width, alpha_results(idx).freq_band(1), alpha_results(idx).freq_band(2));
end

% Define Individual Theta Band: Alpha peak -6 to -2 Hz
theta_band = [peak_freq - 6, peak_freq - 2];

% Display
fprintf('Individual Theta Band for %s: %.2f Hz to %.2f Hz\n', ...
    max_alpha_channel, theta_band(1), theta_band(2));

% Find indices of theta band frequencies
theta_idx = find(freq_oscillatory >= theta_band(1) & ...
                 freq_oscillatory <= theta_band(2));
theta_freqs = freq_oscillatory(theta_idx);
theta_power = pow(max_chan_idx, theta_idx);

figure;
plot(theta_freqs, theta_power, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title(sprintf('Individual Theta Band: %s (%.2f–%.2f Hz)', ...
    max_alpha_channel, theta_band(1), theta_band(2)));
grid on;

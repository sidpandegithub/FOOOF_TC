clear all;
addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

FILE = '20190715_ANSCH_S1.mat';
filename = fullfile(DIR, FILE);
load(filename, 'reref');

data = reref(1); 


% LABELS = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4',...
%     'O1','O2','F7','F8','T7','T8','P7','P8',...
%     'Fz','Cz','Pz','Iz','FC1','FC2','CP1','CP2',...
%     'FC5','FC6','CP5','CP6','TP9','TP10','AFz','FCz',...
%     'F1','F2','C1','C2','P1','P2','AF3','AF4',...
%     'FC3','FC4','CP3','CP4','PO3','PO4','F5','F6',...
%     'C5','C6','P5','P6','AF7','AF8','FT7','FT8',...
%     'TP7','TP8','PO7','PO8','Fpz','CPz','POz','Oz',...
%     'ECG','HEOG','VEOG','Marker'}';

OLD_LABELS  = {'Fp1','P9','Fp2','F7','F3','Fz','F4','F8',...
    'FC5','FC1','FC2','FC6','PO9','T7','C3','Cz',...
    'C4','T8','PO10','CP5','CP1','CP2','CP6','P7',...
    'P3','Pz','P4','P8','POz','O1','Oz','O2',...
    'AF7','AF3','AF4','AF8','F5','F1','F2','F6',...
    'FC3','FCz','FC4','C5','C1','C2','C6','CP3',...
    'CPz','CP4','P5','P1','P2','P6','FT9','PO3',...
    'PO4','FT10','FT7','FT8','TP7','TP8','PO7','PO8',...
    'ECG','HEOG','VEOG','Marker'}';

% NEWEST_LABELS  = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4',...
%     'O1','O2','F7','F8','T7','T8','P7','P8',...
%     'Fz','Cz','Pz','Iz','FC1','FC2','CP1','CP2',...
%     'FC5','FC6','CP5','CP6','TP9','TP10','AFz','FCz',...
%     'F1','F2','C1','C2','P1','P2','AF3','AF4',...
%     'FC3','FC4','CP3','CP4','PO3','PO4','F5','F6',...
%     'C5','C6','P5','P6','AF7','AF8','FT7','FT8',...
%     'TP7','TP8','PO7','PO8','Fpz','CPz','POz','Oz',...
%     'ECG','HEOG','VEOG','Marker'}';

%if changing to one of the labels above
data.label = OLD_LABELS(1:64);

%if using default labels 
%data.label=reref.label 

%make 1:68 if using default labelling 
 for j = 1:numel(data(1).trial)
        data(1).trial{j} = data(1).trial{j}(1:64, :);  
 end
    
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

% Extract frequency and power
freq_oscillatory = oscillatory.freq;
pow = oscillatory.powspctrm;

% Plot average spectrum across channels
freq_original = original.freq;

% Mean across channels
mean_original   = mean(original.powspctrm, 1);
mean_fractal    = mean(fractal.powspctrm, 1);
mean_oscillatory = mean(oscillatory.powspctrm, 1);

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
chan_labels = oscillatory.label;

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

% Plot top 5 alpha channels using semilogy
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

% Bar Plot of Mean Alpha Power
figure;
bar(alpha_mean_pow(sorted_idx));
set(gca, 'XTick', 1:nchan, 'XTickLabel', chan_labels(sorted_idx));
xtickangle(45);
ylabel('Mean Alpha Power (7–13 Hz)');
title('Alpha Power by Channel (Sorted, Post-FOOOF)');
grid on;

% Topographical plot
topo_data = [];
topo_data.label = chan_labels;
topo_data.dimord = 'chan_time';
topo_data.time = 0;
topo_data.avg = alpha_mean_pow(:)';  

% Load layout
cfg = [];
cfg.layout = 'easycapM1.mat';  
layout = ft_prepare_layout(cfg, topo_data);

% Plot topoplot
cfg = [];
cfg.layout = layout;
cfg.marker = 'on';
cfg.comment = 'no';
cfg.colorbar = 'yes';
cfg.zlim = [0 max(abs(topo_data.avg(:)))];  
cfg.colormap = jet;

figure;
ft_topoplotER(cfg, topo_data);
title('Topographical Distribution of Mean Alpha Power (7–13 Hz)');

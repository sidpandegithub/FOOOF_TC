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


% Remove the file extension
image_DIR = 'D:/src/IAF_avgs';
[~, name, ~] = fileparts(FILES{1});
parts = split(name, '_');
core_name = parts{2}; 


% Build data structure for topoplot
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
cfg.zlim = [0 max(abs(topo_data.avg(:)))] ;  
cfg.colormap = jet;

figure;
ft_topoplotER(cfg, topo_data);
title('Topographical Distribution of Mean Alpha Power (7–13 Hz)');
exportgraphics(gcf, fullfile(image_DIR, [core_name '_topo.jpg']), 'ContentType', 'vector');




figure ;
semilogy(freq_original, mean_original, 'k-', 'LineWidth', 2); hold on;
semilogy(freq_original, mean_fractal, 'r--', 'LineWidth', 2);
semilogy(freq_oscillatory, mean_oscillatory, 'b-', 'LineWidth', 2);
legend({'Original', 'Fractal (Aperiodic)', 'Oscillatory'}, 'Location', 'southwest');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title('Power Spectra');
grid on;
exportgraphics(gcf, fullfile(image_DIR, [core_name '_FOOOF_plot.jpg']), 'ContentType', 'vector');

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
title('Top 5 Channels by Alpha Power');
legend(chan_labels(sorted_idx(1:topN)), 'Location', 'northeast');
grid on;
xticks(1:1:35);
exportgraphics(gcf, fullfile(image_DIR, [core_name '_top_channels_plot.jpg']), 'ContentType', 'vector');

figure;
bar(alpha_mean_pow(sorted_idx));
set(gca, 'XTick', 1:nchan, 'XTickLabel', chan_labels(sorted_idx));
xtickangle(45);
ylabel('Mean Alpha Power (7–13 Hz)');
title('Alpha Power by Channel');
grid on;
exportgraphics(gcf, fullfile(image_DIR, [core_name '_power_bar_plot.jpg']), 'ContentType', 'vector');


%% width of alpha 

% Get the power spectrum for the channel with max alpha power
max_chan_power = pow(max_chan_idx, alpha_idx);

% Get the peak power and frequency
peak_val = alpha_peak_vals(max_chan_idx);
peak_freq = alpha_peak_freqs(max_chan_idx);

% Find index of the peak frequency within the alpha band
[~, peak_idx_within_alpha] = max(max_chan_power);

% Calculate 10% of peak (prominence threshold)
threshold = 0.10 * peak_val;

% Search to the left of the peak for crossing
left_idx = peak_idx_within_alpha;
while left_idx > 1 && max_chan_power(left_idx) > threshold
    left_idx = left_idx - 1;
end

% Search to the right of the peak for crossing
right_idx = peak_idx_within_alpha;
while right_idx < length(max_chan_power) && max_chan_power(right_idx) > threshold
    right_idx = right_idx + 1;
end

% Get frequencies at boundaries
left_freq = alpha_freqs(left_idx);
right_freq = alpha_freqs(right_idx);

% Compute width
bandwidth = right_freq - left_freq;

% Apply 6 Hz max width rule
if bandwidth > 6
    bandwidth = 3;  % Apply 3 Hz symmetric width
    band_low = peak_freq - 1.5;
    band_high = peak_freq + 1.5;
else
    band_low = left_freq;
    band_high = right_freq;
end

fprintf('\nDefined alpha peak band around %s: %.2f–%.2f Hz (width = %.2f Hz)\n', ...
    max_alpha_channel, band_low, band_high, bandwidth);

figure;
hold on;
plot(alpha_freqs, max_chan_power, 'b-', 'LineWidth', 2);
plot(peak_freq, peak_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
y_limits = ylim;
plot([band_low band_low], y_limits, 'r--', 'LineWidth', 1.5);
plot([band_high band_high], y_limits, 'r--', 'LineWidth', 1.5);

xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');
title(sprintf('Alpha Peak Band for Channel %s', max_alpha_channel));
legend('Power Spectrum', 'Alpha Peak', 'Band Edges');
grid on;
hold off;

%% 


% Extract the oscillatory power for the selected channel
chan_power = pow(max_chan_idx, :);
log_power = log10(chan_power);

% Plot log-power spectrum for that channel
figure;
plot(freq_oscillatory, log_power, 'b-', 'LineWidth', 2); hold on;
yline(0, 'k--', 'LineWidth', 1.5); % log10(1) = 0 reference line
xlabel('Frequency (Hz)');
ylabel('Log Power (log_{10}(\muV^2/Hz))');
title(sprintf('Log Power Spectrum - %s', chan_labels{max_chan_idx}));
grid on;

% Find where the curve crosses 0 (i.e., log10(1))
above_thresh = log_power > 0;
crossings = diff(above_thresh);

% Indexes of rising and falling edges
start_idx = find(crossings == 1, 1, 'first') + 1;
end_idx = find(crossings == -1 & (1:numel(crossings))' > start_idx, 1, 'first') + 1;

% Check if both were found and are within bounds
if ~isempty(start_idx) && ~isempty(end_idx) && ...
   start_idx <= length(freq_oscillatory) && end_idx <= length(freq_oscillatory)

    f_start = freq_oscillatory(start_idx);
    f_end = freq_oscillatory(end_idx);
    alpha_peak_width = f_end - f_start;

    % Plot markers
    xline(f_start, 'g--', 'LineWidth', 1.5, 'Label', 'Start');
    xline(f_end, 'r--', 'LineWidth', 1.5, 'Label', 'End');

    fprintf('\nAlpha peak width estimated from log crossing:\n');
    fprintf('Start: %.2f Hz, End: %.2f Hz, Width: %.2f Hz\n', f_start, f_end, alpha_peak_width);
else
    fprintf('\nCould not identify valid start/end of log crossing in alpha band.\n');
end


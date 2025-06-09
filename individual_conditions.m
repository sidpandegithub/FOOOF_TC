addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

% Load condition data
condition_mat = load('D:/src/conditions.mat');
condition_data = condition_mat.data;

% Choose condition: 'Desync', 'Sync', or 'Sham'
TARGET_CONDITION = "Desync"; 
filtered_cond_rows = condition_data(strcmp(condition_data.condition, TARGET_CONDITION), :);
file_list = filtered_cond_rows.orig_file;

% Parameters
MAX_ALPHA_RANGE = 3;
THETA_WIDTH     = 3;

for i = 1:numel(file_list)
    file_name = file_list{i};
    [~, base_name, ~] = fileparts(file_name);
    fprintf('\nProcessing %s (%d of %d)\n', file_name, i, numel(file_list));
    
    full_path = fullfile(DIR, file_name);
    if ~isfile(full_path)
        warning('File not found: %s', full_path);
        continue;
    end
    
    load(full_path, 'reref');
    data = reref(3); 
    data.label = reref(3).label;

    % --- FOOOF aperiodic
    cfg               = [];
    cfg.foilim        = [1 30];
    cfg.pad           = 4;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    fractal = ft_freqanalysis(cfg, data);

    % --- Original power
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, data);

    % --- Oscillatory power
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';
    oscillatory = ft_math(cfg, fractal, original);

    % --- Extract power spectrum
    freq_oscillatory = oscillatory.freq;
    pow              = oscillatory.powspctrm;
    chan_labels      = oscillatory.label;
    nchan            = size(pow, 1);

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

    if isempty(peaks)
        warning('No strong alpha peak in %s. Skipping.', base_name);
        continue;
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

    % --- Save results as Excel
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

    outfile = sprintf('D:/individual_results_post_desync/%s_%s_results.xlsx', TARGET_CONDITION, base_name);
    writetable(results, outfile);
    fprintf('Saved results to: %s\n', outfile);
end

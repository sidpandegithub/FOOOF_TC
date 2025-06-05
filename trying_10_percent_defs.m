addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';
OUTFILE = 'D:/output_figures/osc_power_summary.csv';

% Get all .mat files
all_files = dir(fullfile(DIR, '*.mat'));
file_names = {all_files.name};

% Group files by participant
participant_map = containers.Map;

for i = 1:length(file_names)
    parts = split(file_names{i}, '_');
    if numel(parts) >= 3
        participant_id = parts{2};
        if isKey(participant_map, participant_id)
        tmp = participant_map(participant_id);  % get current list
        tmp{end+1} = file_names{i};             % append to it
        participant_map(participant_id) = tmp;  % store it back
else
    participant_map(participant_id) = {file_names{i}};
end

    end
end

% Initialize results storage
all_results = {
    'Participant', 'AlphaPeakHz', 'AlphaBandHz', ...
    'ThetaPeakHz', 'ThetaBandHz', 'Top5AlphaChans', 'Top5ThetaChans'
};

% Loop over each participant
participant_ids = keys(participant_map);

for p = 1:length(participant_ids)
    participant_id = participant_ids{p};
    files = participant_map(participant_id);
    
    if length(files) < 3
        fprintf('Skipping %s: fewer than 3 sessions\n', participant_id);
        continue
    end

    fprintf('\nProcessing participant: %s\n', participant_id);
    
    % Sort files to ensure S1, S2, S3 order
    [~, sort_idx] = sort(files);
    files = files(sort_idx);
    
    oscillatory_data = cell(1, 3);
    
    for i = 1:3
        load(fullfile(DIR, files{i}), 'reref');
        data = reref(1); 
        data.label = reref(1).label;

        % FOOOF aperiodic
        cfg = [];
        cfg.foilim = [1 30];
        cfg.pad = 4;
        cfg.tapsmofrq = 2;
        cfg.method = 'mtmfft';
        cfg.output = 'fooof_aperiodic';
        fractal = ft_freqanalysis(cfg, data);

        % Original power
        cfg.output = 'pow';
        original = ft_freqanalysis(cfg, data);

        % Oscillatory = original / aperiodic
        cfg = [];
        cfg.parameter = 'powspctrm';
        cfg.operation = 'x2./x1';
        oscillatory = ft_math(cfg, fractal, original);

        oscillatory_data{i} = oscillatory;
    end

    % Grand average
    cfg = [];
    cfg.parameter = 'powspctrm';
    avg_oscillatory = ft_freqgrandaverage(cfg, oscillatory_data{:});
    freq = avg_oscillatory.freq;
    pow = avg_oscillatory.powspctrm;
    labels = avg_oscillatory.label;
    nchan = size(pow, 1);

    % Alpha band analysis
    alpha_idx = find(freq >= 7 & freq <= 13);
    alpha_freqs = freq(alpha_idx);
    alpha_peak_freqs = zeros(nchan,1);
    alpha_peak_vals = zeros(nchan,1);
    alpha_mean_pow = zeros(nchan,1);

    for c = 1:nchan
        alpha_power = pow(c, alpha_idx);
        [peak_val, peak_idx] = max(alpha_power);
        alpha_peak_freqs(c) = alpha_freqs(peak_idx);
        alpha_peak_vals(c) = peak_val;
        alpha_mean_pow(c) = mean(alpha_power);
    end

    [~, sorted_idx] = sort(alpha_mean_pow, 'descend');
    max_chan_idx = sorted_idx(1);
    max_alpha_channel = labels{max_chan_idx};
    peak_alpha_freq = alpha_peak_freqs(max_chan_idx);
    top5_alpha_chans = labels(sorted_idx(1:5));

    % Get alpha range
    MAX_ALPHA_RANGE = 3;
    THETA_WIDTH = 3;
    alpha_power = pow(max_chan_idx, :);
    [peaks, locs] = findpeaks(alpha_power, freq, 'MinPeakProminence', 1, 'SortStr', 'descend');

    if isempty(peaks)
        fprintf('No alpha peak found for %s — skipping.\n', participant_id);
        continue
    end

    locs = locs(1);
    peak_tenth = 0.1 * peaks(1);
    alpha_range = find(alpha_power > peak_tenth & ...
        freq >= locs - MAX_ALPHA_RANGE & ...
        freq <= locs + MAX_ALPHA_RANGE);

    lower_alpha = freq(alpha_range(1));
    upper_alpha = freq(alpha_range(end));
    alpha_band_def = [lower_alpha, upper_alpha];

    % Define theta band
    lower_theta = lower_alpha - THETA_WIDTH;
    theta_band_def = [lower_theta, lower_alpha];
    theta_idx = find(freq >= theta_band_def(1) & freq <= theta_band_def(2));
    theta_freqs = freq(theta_idx);

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
    max_theta_chan_idx = sorted_theta_idx(1);
    peak_theta_freq = theta_peak_freqs(max_theta_chan_idx);
    top5_theta_chans = labels(sorted_theta_idx(1:5));

    % Add row to results
    row = {
        participant_id, ...
        peak_alpha_freq, ...
        sprintf('%.2f–%.2f', alpha_band_def(1), alpha_band_def(2)), ...
        peak_theta_freq, ...
        sprintf('%.2f–%.2f', theta_band_def(1), theta_band_def(2)), ...
        strjoin(top5_alpha_chans, ', '), ...
        strjoin(top5_theta_chans, ', ')
    };

    all_results(end+1, :) = row;
end

% Save to CSV
OUTFILE = fullfile(DIR, 'post_results-10_percent_def.xlsx');
results_tbl = cell2table(all_results(2:end,:), 'VariableNames', all_results(1,:));
writetable(results_tbl, OUTFILE);
fprintf('\nSaved all results to: %s\n', OUTFILE);


addpath('C:/fieldtrip-20240113');
ft_defaults;

DIR = 'D:/src/11-reref';
FILES = {
    '20190716_LIGOH_S1.mat';
    '20190719_LIGOH_S2.mat';
    '20190724_LIGOH_S3.mat'
};

% Step 1: Load and compute power spectra for all sessions
all_pow = cell(1, length(FILES));
for i = 1:length(FILES)
    filename = fullfile(DIR, FILES{i});
    load(filename, 'reref');
    data = reref(1); % Baseline resting-state EEG

    cfg = [];
    cfg.foilim = [1 35];
    cfg.method = 'mtmfft';
    cfg.pad = 4;
    cfg.tapsmofrq = 2;
    cfg.output = 'pow';
    pow = ft_freqanalysis(cfg, data);

    all_pow{i} = pow;
end

% Step 2: Grand average power spectrum across sessions (per channel)
cfg = [];
cfg.parameter = 'powspctrm';
avg_pow = ft_freqgrandaverage(cfg, all_pow{:});

% Step 3: Compute aperiodic (1/f) component on the average spectrum
cfg = [];
cfg.foilim = [1 35];
cfg.pad = 4;
cfg.tapsmofrq = 2;
cfg.method = 'mtmfft';
cfg.output = 'fooof_aperiodic';
aperiodic = ft_freqanalysis(cfg, avg_pow);

% Step 4: Calculate residual oscillatory spectrum by dividing original avg_pow by aperiodic
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2./x1';
residual = ft_math(cfg, aperiodic, avg_pow);

freq = residual.freq;
alpha_band = [7 13];
alpha_idx = find(freq >= alpha_band(1) & freq <= alpha_band(2));

% For each channel, find alpha peak and alpha power AUC
results = struct();
for ch = 1:length(residual.label)
    chan_name = residual.label{ch};
    psd = residual.powspctrm(ch, :);
    alpha_psd = psd(alpha_idx);
    alpha_freqs = freq(alpha_idx);

    % Use findpeaks to detect peaks in alpha band
    [pks, locs, widths, prominences] = findpeaks(alpha_psd, alpha_freqs, ...
        'MinPeakProminence', 0.05);

    if isempty(pks)
        % No clear alpha peak found
        results.(chan_name).IAF = NaN;
        results.(chan_name).band = [NaN NaN];
        results.(chan_name).AUC = NaN;
    else
        % Take peak with highest prominence
        [~, max_i] = max(prominences);
        peak_freq = locs(max_i);
        peak_width = widths(max_i);

        % Define integration window for AUC (max 6 Hz width, centered on peak)
        left = max(alpha_band(1), peak_freq - peak_width/2);
        right = min(alpha_band(2), peak_freq + peak_width/2);
        if (right - left) > 6
            left = peak_freq - 3;
            right = peak_freq + 3;
            % Clamp to alpha band edges
            if left < alpha_band(1), left = alpha_band(1); end
            if right > alpha_band(2), right = alpha_band(2); end
        end

        % Compute area under curve for oscillatory power in that window
        range_idx = find(freq >= left & freq <= right);
        AUC = trapz(freq(range_idx), psd(range_idx));

        % Store results
        results.(chan_name).IAF = peak_freq;
        results.(chan_name).band = [left right];
        results.(chan_name).AUC = AUC;
    end
end

% Optional: Display top 10 channels by alpha AUC
all_channels = fieldnames(results);
auc_values = nan(length(all_channels),1);
for k = 1:length(all_channels)
    auc_values(k) = results.(all_channels{k}).AUC;
end
[sorted_auc, sort_idx] = sort(auc_values, 'descend');

fprintf('Top 10 channels by alpha power (AUC in residual spectrum):\n');
for i = 1:min(10,length(all_channels))
    ch_name = all_channels{sort_idx(i)};
    fprintf('%2d. %s: IAF = %.2f Hz, AUC = %.3f\n', i, ch_name, ...
        results.(ch_name).IAF, results.(ch_name).AUC);
end

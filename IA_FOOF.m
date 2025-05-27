addpath('C:\fieldtrip-20240113');
DIR = 'D:\src\11-reref';

FILES = {
    '20190716_LIGOH_S1.mat',
    '20190719_LIGOH_S2.mat',
    '20190724_LIGOH_S3.mat'
};

oscillatory_data = cell(1, numel(FILES));

for i = 1:numel(FILES)
    filename = fullfile(DIR, FILES{i});
    load(filename, 'reref');
    data = reref(1);  % Baseline RS-EEG
    
    % Compute aperiodic component using FOOOF settings
    cfg               = [];
    cfg.foilim        = [1 35];
    cfg.pad           = 4;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    fractal = ft_freqanalysis(cfg, data);

    % Compute original power
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, data);

    % Get oscillatory component by removing aperiodic
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';  % Divide original by fractal
    oscillatory = ft_math(cfg, fractal, original);

    % Keep only Oz channel
    cfg = [];
    cfg.channel = 'Oz';
    oscillatory_Oz = ft_selectdata(cfg, oscillatory);

    oscillatory_data{i} = oscillatory_Oz;
end

% Average across the three sessions
cfg = [];
cfg.parameter = 'powspctrm';
avg_oscillatory = ft_freqgrandaverage(cfg, oscillatory_data{:});

% Find IAF in the 8–13 Hz range
frequencies = avg_oscillatory.freq;
power = avg_oscillatory.powspctrm(1,:);  % Only Oz channel

alpha_idx = find(frequencies >= 8 & frequencies <= 13);
[~, max_idx] = max(power(alpha_idx));
IAF = frequencies(alpha_idx(max_idx));

fprintf('Individual Alpha Frequency (Oz, post-FOOOF): %.2f Hz\n', IAF);


figure;
plot(frequencies, power, 'LineWidth', 2);
hold on;
plot(IAF, power(frequencies == IAF), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlim([5 15]);
xlabel('Frequency (Hz)');
ylabel('Power (Oscillatory)');
title(sprintf('Oscillatory Power Spectrum at Oz (IAF = %.2f Hz)', IAF));
grid on;

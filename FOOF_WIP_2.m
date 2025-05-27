addpath('C:\fieldtrip-20240113');
% addpath('C:\fieldtrip-20250106');
DIR = 'D:\src\11-reref';
save_DIR = 'D:\src\FOOOF_division_smallerfreq';


fileList = dir(fullfile(DIR, '*.mat'));

for i = 1:length(fileList)
    FILE = fileList(i).name;
    filename = fullfile(DIR, FILE);
    load(filename);
    
    try
        data = reref(2);
    catch
        fprintf('Skipping %s: reref(2) not available.\n', FILE);
        continue
    end

    % FOOOF (fractal) and original power spectra
    cfg               = [];
    cfg.foilim        = [1 35];
    cfg.pad           = 4;
    cfg.tapsmofrq     = 2;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    fractal = ft_freqanalysis(cfg, data);
    cfg.output        = 'pow';
    original = ft_freqanalysis(cfg, data);

    % Oscillatory power (division)
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';
    oscillatory_alt = ft_math(cfg, fractal, original);

    % Channel selection
    frontal_channels  = ismember(original.label, {'F4','F2','FC2'});
    parietal_channels = ismember(original.label, {'P4','P2','PO4'});

    % Plotting
    fig = figure('Visible', 'off');
    subplot(1,2,1);
    semilogy(original.freq, squeeze(mean(original.powspctrm(frontal_channels,:),1)), 'k', 'linewidth', 1.5); hold on
    semilogy(fractal.freq, squeeze(mean(fractal.powspctrm(frontal_channels,:),1)), 'r', 'linewidth', 1.5);
    semilogy(oscillatory_alt.freq, squeeze(mean(oscillatory_alt.powspctrm(frontal_channels,:),1)), 'g', 'linewidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('log(Power)'); title('Frontal Channels');
    legend({'original','fractal','oscillatory\_alt'}, 'Location','southwest');
    grid on;
    xticks(0:2:50);  
    axis tight;
    
    subplot(1,2,2);
    semilogy(original.freq, squeeze(mean(original.powspctrm(parietal_channels,:),1)), 'k', 'linewidth', 1.5); hold on
    semilogy(fractal.freq, squeeze(mean(fractal.powspctrm(parietal_channels,:),1)), 'r', 'linewidth', 1.5);
    semilogy(oscillatory_alt.freq, squeeze(mean(oscillatory_alt.powspctrm(parietal_channels,:),1)), 'g', 'linewidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('log(Power)'); title('Parietal Channels');
    legend({'original','fractal','oscillatory\_alt'},'Location', 'southwest');
    grid on;
    xticks(0:2:50); 
    axis tight;

    % Save figure
    [~, name, ~] = fileparts(FILE);
    saveas(fig, fullfile(save_DIR, [name '.png']));
    close(fig); 
    fprintf('Processed and saved: %s\n', FILE);
end


%% 

figure;
subplot(1,2,1);
%subplot(2,1,2)
semilogy(original.freq, squeeze(mean(original.powspctrm(frontal_channels,:),1)), 'k', 'linewidth', 1.5);
hold on
semilogy(fractal.freq, squeeze(mean(fractal.powspctrm(frontal_channels,:),1)), 'r', 'linewidth', 1.5);
semilogy(oscillatory.freq, squeeze(mean(oscillatory.powspctrm(frontal_channels,:),1)), 'g', 'linewidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('log(Power)');
title('Frontal Channels');
legend({'original','fractal','oscillatory\_alt'});
grid on; axis tight;

subplot(1,2,2);
%subplot(2,1,2)
semilogy(original.freq, squeeze(mean(original.powspctrm(parietal_channels,:),1)), 'k', 'linewidth', 1.5);
hold on
semilogy(fractal.freq, squeeze(mean(fractal.powspctrm(parietal_channels,:),1)), 'r', 'linewidth', 1.5);
semilogy(oscillatory.freq, squeeze(mean(oscillatory.powspctrm(parietal_channels,:),1)), 'g', 'linewidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('log(Power)');
title('Parietal Channels');
legend({'original','fractal','oscillatory'});
grid on; axis tight;



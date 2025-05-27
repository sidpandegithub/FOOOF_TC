addpath ('C:\fieldtrip-20240113');
% addpath('C:\fieldtrip-20250106');
DIR = 'D:\src\11-reref';
FILE = '20190716_LIGOH_S1.mat';
filename = fullfile(DIR,FILE);

load(filename);
data = reref(2);

%% just plotting raw PSD
cfg = [];
cfg.foilim = [1 200];       
cfg.method = 'mtmfft';       
cfg.output = 'pow';          
cfg.tapsmofrq = 2;           
cfg.pad = 4;                 
psd_raw = ft_freqanalysis(cfg, data);  

figure;
plot(psd_raw.freq, log10(squeeze(psd_raw.powspctrm)));
xlabel('Frequency (Hz)');
ylabel('log_{10} Power');
title('Raw Power Spectral Density');
xlim([1 100]);  

%% FOOOF tutorial from fieldtrip

cfg               = [];
cfg.foilim        = [1 200];
cfg.pad           = 4;
cfg.tapsmofrq     = 2;
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, data);

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
oscillatory = ft_math(cfg, fractal, original);

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
oscillatory_alt = ft_math(cfg, fractal, original);

% display the spectra on a log-log scale
figure();
subplot(1,2,1); hold on;
plot(log(original.freq), log(original.powspctrm),'k');
plot(log(fractal.freq), log(fractal.powspctrm));
plot(log(fractal.freq), log(oscillatory.powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum-fractal'},'location','southwest');
F = mean(fractal.powspctrm(:));
O = mean(oscillatory.powspctrm(:));
if F~=0 && O==0
  title('pure fractal signal');
elseif F==0 && O~=0
  title('pure oscillatory signal');
elseif F~=0 && O~=0
  title('mixed signal');
end
subplot(1,2,2); hold on;
plot(log(original.freq), log(original.powspctrm),'k');
plot(log(fractal.freq), log(fractal.powspctrm));
plot(log(oscillatory_alt.freq), log(oscillatory_alt.powspctrm));
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
title('oscillatory = spectrum / fractal');







%% average across all channels for a frequency bin
DIR = 'D:\src\11-reref';
FILE = '20190715_ANSCH_S1.mat';
filename = fullfile(DIR, FILE);
load(filename);
data = reref(2);

cfg = [];
cfg.foilim = [1 200];
cfg.pad = 4;
cfg.tapsmofrq = 2;
cfg.method = 'mtmfft';
cfg.output = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, data);
cfg.output = 'pow';
original = ft_freqanalysis(cfg, data);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2 - x1';
oscillatory = ft_math(cfg, fractal, original);

cfg = [];
cfg.parameter     = 'powspctrm';
cfg.operation = 'x2 ./ x1';
oscillatory_alt = ft_math(cfg, fractal, original);



avg_original = mean(original.powspctrm, 1);
avg_fractal = mean(fractal.powspctrm, 1);
avg_oscillatory = mean(oscillatory.powspctrm, 1);
avg_oscillatory_alt = mean(oscillatory_alt.powspctrm, 1);


figure();
hold on;
plot(original.freq, log(avg_original), 'k', 'LineWidth', 1.5);
plot(fractal.freq, log(avg_fractal), 'b', 'LineWidth', 1.5);
plot(oscillatory.freq, log(avg_oscillatory), 'r', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('log-Power');
title('Average PSD: Original, Fractal (Aperiodic), and Oscillatory');
legend({'Original', 'Fractal', 'Oscillatory'}, 'Location', 'SouthWest');
xlim([1 100]);
grid on;


figure();
hold on;
plot(log(original.freq), log(avg_original), 'k');
plot(log(fractal.freq), log(avg_fractal), 'b');
plot(log(oscillatory_alt.freq), log(avg_oscillatory_alt), 'm');
xlabel('log-freq'); ylabel('log-power'); grid on;
legend({'Original','Fractal','Oscillatory = Spectrum / Fractal'},'location','southwest');
title('Oscillatory = Spectrum / Fractal');


%% 

cfg               = [];
cfg.foilim        = [1 40];
cfg.pad           = 4;
cfg.tapsmofrq     = 2;
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'pow';
original = ft_freqanalysis(cfg, data);

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
oscillatory = ft_math(cfg, fractal, original);

cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
oscillatory_alt = ft_math(cfg, fractal, original);
frontal_channels = ismember(original.label, {'F4','F2','FC2'});
parietal_channels = ismember(original.label, {'P4','P2','PO4'});

% Plot power spectra for the frontal ROI

frontal_channels = ismember(original.label, {'F4','F2','FC2'});
parietal_channels = ismember(original.label, {'P4','P2','PO4'});
%% log

figure;

subplot(1,2,1);
plot(original.freq, log(squeeze(mean(original.powspctrm(frontal_channels,:), 1))), 'k', 'LineWidth', 1.5); hold on;
plot(fractal.freq, log(squeeze(mean(fractal.powspctrm(frontal_channels,:), 1))), 'r', 'LineWidth', 1.5);
plot(oscillatory_alt.freq, log(squeeze(mean(oscillatory_alt.powspctrm(frontal_channels,:), 1))), 'g', 'LineWidth', 1.5);
plot(oscillatory.freq, log(squeeze(mean(oscillatory.powspctrm(frontal_channels,:), 1))), 'm', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('log(Power)');
title('Frontal Channels');
legend({'Original', 'Fractal', 'Alt Oscillatory', 'Oscillatory'}, 'Location', 'southwest');
grid on;
xticks(0:2:50);
axis tight;

subplot(1,2,2);
plot(original.freq, log(squeeze(mean(original.powspctrm(parietal_channels,:), 1))), 'k', 'LineWidth', 1.5); hold on;
plot(fractal.freq, log(squeeze(mean(fractal.powspctrm(parietal_channels,:), 1))), 'r', 'LineWidth', 1.5);
plot(oscillatory_alt.freq, log(squeeze(mean(oscillatory_alt.powspctrm(parietal_channels,:), 1))), 'g', 'LineWidth', 1.5);
plot(oscillatory.freq, log(squeeze(mean(oscillatory.powspctrm(parietal_channels,:), 1))), 'm', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('log(Power)');
title('Parietal Channels');
legend({'Original', 'Fractal', 'Alt Oscillatory', 'Oscillatory'}, 'Location', 'southwest');
grid on;
xticks(0:2:50);
axis tight;
%% semilogy


figure;

% Frontal Channels subplot
subplot(1,2,1);
semilogy(original.freq, squeeze(mean(original.powspctrm(frontal_channels,:),1)), 'k', 'LineWidth', 1.5); hold on;
semilogy(fractal.freq, squeeze(mean(fractal.powspctrm(frontal_channels,:),1)), 'r', 'LineWidth', 1.5);
semilogy(oscillatory_alt.freq, squeeze(mean(oscillatory_alt.powspctrm(frontal_channels,:),1)), 'g', 'LineWidth', 1.5);
semilogy(oscillatory.freq, squeeze(mean(oscillatory.powspctrm(frontal_channels,:),1)), 'm', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (log scale)');
title('Frontal Channels');
legend({'original', 'fractal', 'oscillatory\_alt', 'oscillatory'}, 'Location', 'southwest');
grid on;
xticks(0:2:50);
axis tight;

% Parietal Channels subplot
subplot(1,2,2);
semilogy(original.freq, squeeze(mean(original.powspctrm(parietal_channels,:),1)), 'k', 'LineWidth', 1.5); hold on;
semilogy(fractal.freq, squeeze(mean(fractal.powspctrm(parietal_channels,:),1)), 'r', 'LineWidth', 1.5);
semilogy(oscillatory_alt.freq, squeeze(mean(oscillatory_alt.powspctrm(parietal_channels,:),1)), 'g', 'LineWidth', 1.5);
semilogy(oscillatory.freq, squeeze(mean(oscillatory.powspctrm(parietal_channels,:),1)), 'm', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Power (log scale)');
title('Parietal Channels');
legend({'original', 'fractal', 'oscillatory\_alt', 'oscillatory'}, 'Location', 'southwest');
grid on;
xticks(0:2:50);
axis tight;


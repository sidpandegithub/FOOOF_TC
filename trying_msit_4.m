% data_folder = 'D:\src\12-reref';  % Folder with individual .mat files
% data_DIR    = 'D:\src';           % Folder with twincoil_data.mat
% save_folder = 'D:\src\preMSIT';   % Output folder

data_folder = 'E:\src\12-reref';  % Folder with individual .mat files
data_DIR    = 'E:\src';           % Folder with twincoil_data.mat
save_folder = 'E:\src\preMSIT';   % Output folder

% Load twincoil data
load(fullfile(data_DIR, 'twincoil_data.mat'));

% Extract core names and assign participant IDs
core_names = cellfun(@(x) extractBetween(x, '_', '_S'), {twincoil_data.name}, 'UniformOutput', false);
core_names = cellfun(@(x) x{1}, core_names, 'UniformOutput', false);
[~, ~, participant_ids] = unique(core_names, 'stable');

for p = 1:length(twincoil_data)
    twincoil_data(p).participant_id = participant_ids(p);
end

% Get unique participant IDs
unique_ids = unique(participant_ids);

% Loop over each participant
for pid = unique_ids'
    % Find all entries for this participant
    idxs = find(participant_ids == pid);
    data_for_subj = [];

    for j = 1:length(idxs)
        entry = twincoil_data(idxs(j));
        filename = [entry.name '.mat'];
        filepath = fullfile(data_folder, filename);

        if ~isfile(filepath)
            warning('Missing file: %s', filepath);
            continue;
        end

        % Load reref
        S = load(filepath, 'reref');
        reref = S.reref;

        % Extract desired block
        target_block_idx = 4;
        block_data = reref(target_block_idx);

        % Attach session and name info
        tokens = regexp(filename, '^\d+_([^_]+)_S(\d+)\.mat$', 'tokens');
        if isempty(tokens)
            warning('empty: %s', filename);
            continue;
        end

        subject = tokens{1}{1};
        session_num = str2double(tokens{1}{2});
        block_data.session = session_num;
        block_data.name = subject;

        data_for_subj = [data_for_subj, block_data];
    end

    % Save to one file per subject
    if ~isempty(data_for_subj)
        save_path = fullfile(save_folder, [subject '.mat']);
        save(save_path, 'data_for_subj');
        fprintf('Saved all sessions for %s to %s\n', subject, save_path);
    end
end



%%

clear all;
savefolder = 'E:\src\preMSIT'; 
filename = 'VIGUP.mat';          
fullpath = fullfile(savefolder, filename);
[~, subject_id, ~] = fileparts(filename);  

% Make folders
fig_save_folder = fullfile(savefolder, 'figures');
if ~exist(fig_save_folder, 'dir')
    mkdir(fig_save_folder);
end

% Load data
load(fullpath); 
reref = data_for_subj;

% Select HIT trials from each session
for i = 1:3
    hit_trials = find(contains({reref(i).event.type}, 'HIT'));  
    cfg = [];
    cfg.trials = hit_trials;
    filtered_reref(i) = ft_selectdata(cfg, reref(i));
end

reref = filtered_reref;
cfg = [];  % Reset cfg
combined_data = ft_appenddata(cfg, reref(1), reref(2), reref(3));
reref_all = combined_data;

% Baseline correction 
BASE = [-0.2 0.0];
cfg = [];
cfg.channel = 'all';
cfg.demean = 'yes';
cfg.baselinewindow = BASE;
temp_msit = ft_preprocessing(cfg, reref_all);
all_trials = temp_msit;

% ERP Analysis 
cfg = [];
cfg.covariance = 'yes';
tlock = ft_timelockanalysis(cfg, all_trials);

% ERP Component Settings 
N2_COI = 'Cz';
P3_COI = 'Pz';
N2_TOI = [0.25 0.35];
P3_TOI = [0.30 0.80];

n2_chan = strcmp(tlock.label, N2_COI);
p3_chan = strcmp(tlock.label, P3_COI);
n2_time = tlock.time >= N2_TOI(1) & tlock.time <= N2_TOI(2);
p3_time = tlock.time >= P3_TOI(1) & tlock.time <= P3_TOI(2);

% ERP Waveforms and Peaks
n2_wave = tlock.avg(n2_chan, n2_time);
p3_wave = tlock.avg(p3_chan, p3_time);
n2_times = tlock.time(n2_time);
p3_times = tlock.time(p3_time);

[n2_peak_amp, n2_peak_idx] = min(n2_wave);
[p3_peak_amp, p3_peak_idx] = max(p3_wave);
n2_peak_latency = n2_times(n2_peak_idx);
p3_peak_latency = p3_times(p3_peak_idx);

% Mean Around Peaks
n2_window = [n2_peak_latency - 0.020, n2_peak_latency + 0.020];
p3_window = [p3_peak_latency - 0.050, p3_peak_latency + 0.050];
n2_peak_range = tlock.time >= n2_window(1) & tlock.time <= n2_window(2);
p3_peak_range = tlock.time >= p3_window(1) & tlock.time <= p3_window(2);
n2_mean_around_peak = mean(tlock.avg(n2_chan, n2_peak_range), 2);
p3_mean_around_peak = mean(tlock.avg(p3_chan, p3_peak_range), 2);

%Print to Console
fprintf('Subject: %s\n', subject_id);
fprintf('N200 peak at %s: %.3f µV (%.3f s)\n', N2_COI, n2_peak_amp, n2_peak_latency);
fprintf('P300 peak at %s: %.3f µV (%.3f s)\n', P3_COI, p3_peak_amp, p3_peak_latency);
fprintf('N200 mean ±20ms: %.3f µV\n', n2_mean_around_peak);
fprintf('P300 mean ±50ms: %.3f µV\n', p3_mean_around_peak);

% Plot and Save Figure 
figure('Visible', 'off');

subplot(2,1,1)
plot(tlock.time, tlock.avg(n2_chan, :), 'b', 'LineWidth', 1.5); hold on;
xline(N2_TOI(1), 'r--'); xline(N2_TOI(2), 'r--');
plot(n2_peak_latency, n2_peak_amp, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
title(['N200 at ', N2_COI]); grid on;

subplot(2,1,2)
plot(tlock.time, tlock.avg(p3_chan, :), 'g', 'LineWidth', 1.5); hold on;
xline(P3_TOI(1), 'r--'); xline(P3_TOI(2), 'r--');
plot(p3_peak_latency, p3_peak_amp, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
title(['P300 at ', P3_COI]); grid on;

saveas(gcf, fullfile(fig_save_folder, [subject_id '_ERP.png']));
close;

% Save Results to Excel
results_table = table({subject_id}, n2_peak_amp, n2_peak_latency, n2_mean_around_peak, ...
                      p3_peak_amp, p3_peak_latency, p3_mean_around_peak, ...
                      'VariableNames', ...
                      {'Subject', 'N2_Peak_Amp', 'N2_Latency', 'N2_MeanAroundPeak', ...
                       'P3_Peak_Amp', 'P3_Latency', 'P3_MeanAroundPeak'});

writetable(results_table, fullfile(savefolder, [subject_id '_ERP_summary.xlsx']));
fprintf('Saved: %s_ERP.png and %s_ERP_summary.xlsx\n', subject_id, subject_id);

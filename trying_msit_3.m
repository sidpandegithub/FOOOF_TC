clear;

data_folder = 'D:\src\12-reref';  % Folder with individual .mat files
data_DIR    = 'D:\src';           % Folder with twincoil_data.mat
save_folder = 'D:\src\preMSIT';   % Output folder

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
savefolder = 'D:\src\preMSIT';  
filename = 'IPAN.mat';          
fullpath = fullfile(savefolder, filename);

load(fullpath); 
reref = data_for_subj;

for i = 1:3
    hit_trials = find(contains({reref(i).event.type}, 'HIT'));  
    cfg = [];
    cfg.trials = hit_trials;
    filtered_reref(i) = ft_selectdata(cfg, reref(i));
end

reref = filtered_reref;

combined_data = ft_appenddata(cfg, reref(1), reref(2), reref(3));
 
reref_all = combined_data;


BASE    = [-0.2 0.0]; % or [-0.5 0.0];
% m = 1:length(msit_index);
cfg                 = [];
cfg.channel         = 'all';
cfg.demean          = 'yes';
cfg.baselinewindow  = BASE;
% temp_msit           = ft_preprocessing(cfg,reref(msit_index(m)));
temp_msit           = ft_preprocessing(cfg,reref_all);
all_trials = temp_msit;
cfg = [];
cfg.covariance = 'yes';
tlock = ft_timelockanalysis(cfg, all_trials);


% Define channel and time-of-interest
N2_COI = 'Cz';
P3_COI = 'Pz';
N2_TOI = [0.25 0.35];
P3_TOI = [0.30 0.80];

% Find channel indices and time windows
n2_chan = strcmp(tlock.label, N2_COI);
p3_chan = strcmp(tlock.label, P3_COI);

n2_time = tlock.time >= N2_TOI(1) & tlock.time <= N2_TOI(2);
p3_time = tlock.time >= P3_TOI(1) & tlock.time <= P3_TOI(2);

% Extract mean amplitudes
n2_data = mean(tlock.avg(n2_chan, n2_time), 2);
p3_data = mean(tlock.avg(p3_chan, p3_time), 2);

% Plotting
figure;

% N200 subplot 
subplot(2,1,1)
plot(tlock.time, tlock.avg(n2_chan, :), 'b', 'LineWidth', 1.5); hold on;
xline(N2_TOI(1), 'r--');
xline(N2_TOI(2), 'r--');
xticks(min(tlock.time):0.1:max(tlock.time));  % More x-axis increments
xlabel('Time (s)');
ylabel('Amplitude (µV)');
title(['N200 at ', N2_COI]);
legend('ERP', 'N200 window');
grid on;

% P300 subplot 
subplot(2,1,2)
plot(tlock.time, tlock.avg(p3_chan, :), 'g', 'LineWidth', 1.5); hold on;
xline(P3_TOI(1), 'r--');
xline(P3_TOI(2), 'r--');
xticks(min(tlock.time):0.1:max(tlock.time));  % More x-axis increments
xlabel('Time (s)');
ylabel('Amplitude (µV)');
title(['P300 at ', P3_COI]);
legend('ERP', 'P300 window');
grid on;


% Extract ERP waveforms and time vectors within each TOI
n2_wave  = tlock.avg(n2_chan, n2_time);
p3_wave  = tlock.avg(p3_chan, p3_time);
n2_times = tlock.time(n2_time);
p3_times = tlock.time(p3_time);

% Find N200 negative peak
[n2_peak_amp, n2_peak_idx] = min(n2_wave);
n2_peak_latency = n2_times(n2_peak_idx);

% Find P300 positive peak
[p3_peak_amp, p3_peak_idx] = max(p3_wave);
p3_peak_latency = p3_times(p3_peak_idx);

% Print peak results
fprintf('N200 peak amplitude at %s: %.3f µV (latency: %.3f s)\n', N2_COI, n2_peak_amp, n2_peak_latency);
fprintf('P300 peak amplitude at %s: %.3f µV (latency: %.3f s)\n', P3_COI, p3_peak_amp, p3_peak_latency);

% Mark peaks on the plots
subplot(2,1,1)
plot(n2_peak_latency, n2_peak_amp, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
subplot(2,1,2)
plot(p3_peak_latency, p3_peak_amp, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');


% N200 ±20 ms window
n2_window = [n2_peak_latency - 0.020, n2_peak_latency + 0.020];
n2_peak_range = tlock.time >= n2_window(1) & tlock.time <= n2_window(2);
n2_mean_around_peak = mean(tlock.avg(n2_chan, n2_peak_range), 2);

% P300 ±50 ms window
p3_window = [p3_peak_latency - 0.050, p3_peak_latency + 0.050];
p3_peak_range = tlock.time >= p3_window(1) & tlock.time <= p3_window(2);
p3_mean_around_peak = mean(tlock.avg(p3_chan, p3_peak_range), 2);

% Print results
fprintf('N200 average amplitude (±20ms around peak): %.3f µV\n', n2_mean_around_peak);
fprintf('P300 average amplitude (±50ms around peak): %.3f µV\n', p3_mean_around_peak);

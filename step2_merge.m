function step2_merge(patient)
% STEP2_MERGE  Merge overlapping UFOs from ufo_inspect_log.mat and recreate figures.
%
% This script:
%   1) Loads the existing log from Results/<patient>/ufo_inspect_log.mat
%   2) Groups entries by micro_channel + macro_channel.
%   3) Within each group, merges any entries whose [new_start_us, new_end_us] overlap.
%   4) Overwrites the original ufo_inspect_log.mat with the merged set.
%   5) Re-exports all PNG figures to Results/<pat>/UFO_figures/.

if nargin < 1 || isempty(patient)
    error('step2_merge: patient ID required, e.g. step2_merge(''pat001'')');
end

[base_dir, matmef_path] = ufo_config();
if exist(matmef_path,'dir'), addpath(matmef_path); end

res_dir  = fullfile(base_dir, 'Results', patient);
log_path = fullfile(res_dir, 'ufo_inspect_log.mat');

if ~exist(log_path,'file')
    error('step2_merge: Log file not found at %s', log_path);
end

% Output folders
res_fig_dir = fullfile(res_dir, 'UFO_figures');
out_disc = fullfile(res_fig_dir, 'discard');
out_t1   = fullfile(res_fig_dir, 'type1');
out_t2   = fullfile(res_fig_dir, 'type2');

%% ------------------------------------------------------------------------
% 1) Merge overlaps
%% ------------------------------------------------------------------------
% Load Results
fprintf('Loading log from %s...\n', log_path);
S = load(log_path,'Results');
Results = S.Results;
n_orig = numel(Results);
fprintf('  Found %d entries.\n', n_orig);

fprintf('\nChecking for overlapping UFOs...\n');

% Group by micro+macro
[groups, ~, group_idx] = unique(strcat({Results.micro_channel}, '_', {Results.macro_channel}));
merged_results = struct('patient',{},'micro_channel',{},'macro_channel',{}, ...
    'ufo_idx',{},'uuid',{},'orig_start_us',{},'orig_end_us',{}, ...
    'orig_freq_hz',{},'freq_class',{},'decision',{}, ...
    'new_start_us',{},'new_end_us',{},'type',{});

any_merged = false;
for g = 1:numel(groups)
    R_group = Results(group_idx == g);
    
    % Sort by start time
    [~, ord] = sort([R_group.new_start_us]);
    R_group = R_group(ord);
    
    if isempty(R_group), continue; end
    
    current = R_group(1);
    for i = 2:numel(R_group)
        next = R_group(i);
        
        % Check for overlap
        if next.new_start_us <= current.new_end_us
            % Overlap! Merge them.
            fprintf('  [merge] %s: merging idx %d and %d\n', ...
                current.micro_channel, current.ufo_idx, next.ufo_idx);
            
            % Union of windows
            current.new_start_us = min(current.new_start_us, next.new_start_us);
            current.new_end_us   = max(current.new_end_us, next.new_end_us);
            
            % For other fields, take the "most classified" one or the first one
            if current.type == 0 && next.type ~= 0
                current.type = next.type;
            end
            if strcmp(current.decision,'skip') && ~strcmp(next.decision,'skip')
                current.decision = next.decision;
            end
            any_merged = true;
        else
            % No overlap, push current and start new
            merged_results(end+1) = current;
            current = next;
        end
    end
    merged_results(end+1) = current;
end

if ~any_merged
    fprintf('  No overlaps found. Nothing to merge.\n');
    return;
end

n_merged = numel(merged_results);
fprintf('  Reduction: %d -> %d entries.\n', n_orig, n_merged);

% Overwrite log
Results = merged_results;
save(log_path,'Results');
fprintf('  Updated log saved to: %s\n', log_path);

% Re-create/clean figure folders (ONLY if we merged something)
fprintf('\nCleaning and preparing figure folders in Results...\n');
if exist(out_disc,'dir'), rmdir(out_disc,'s'); end
if exist(out_t1,'dir'),   rmdir(out_t1,'s');   end
if exist(out_t2,'dir'),   rmdir(out_t2,'s');   end

if ~exist(res_fig_dir,'dir'), mkdir(res_fig_dir); end
mkdir(out_disc);
mkdir(out_t1);
mkdir(out_t2);

%% ------------------------------------------------------------------------
% 2) Re-export all figures
%% ------------------------------------------------------------------------
fprintf('\nRe-exporting all figures for merged set (to Results only)...\n');

data_dir = fullfile(base_dir, 'Data', patient);
mef_path = fullfile(data_dir, 'sub.mefd');
pw       = 'bemena';

if ~exist(mef_path,'dir')
    warning('MEF directory not found at: %s.', mef_path);
    return;
end

for k = 1:numel(Results)
    R = Results(k);
    [trace, t_ms] = read_ufo_window_simple(mef_path, pw, R.micro_channel, R.new_start_us, R.new_end_us);
    
    if ~isempty(trace)
        hFig = figure('Visible','off','position',[0 387 1024 150]);
        plot(t_ms, trace, 'k-'); hold on;
        dur_ms = (R.new_end_us - R.new_start_us)/1e3;
        plot(t_ms(t_ms>=0 & t_ms<=dur_ms), trace(t_ms>=0 & t_ms<=dur_ms), 'r-', 'LineWidth', 1.5);
        title(sprintf('MERGED: %s | %s | idx %d', R.patient, R.micro_channel, R.ufo_idx));
        
        id_str = sprintf('%s_%s_%s_idx%04d', R.patient, R.micro_channel, R.macro_channel, R.ufo_idx);
        if strcmp(R.decision,'discard'), sd='discard'; elseif R.type==1, sd='type1'; elseif R.type==2, sd='type2'; else, sd='.'; end
        exportgraphics(hFig, fullfile(res_fig_dir, sd, [id_str '.png']), 'Resolution',150);
        close(hFig);
    end
end

fprintf('\nDone!\n');

end

function [trace, t_ms] = read_ufo_window_simple(mef_path, pw, ch_name, start_us, end_us)
trace = []; t_ms = []; pad_us = 20e3;
win_start = int64(start_us - pad_us); win_end = int64(end_us + pad_us);
try
    [~, tr] = readMef3(mef_path, pw, ch_name, 'time', win_start, win_end);
    if isempty(tr), return; end
    tr = double(tr(:)); dt_us = double(win_end - win_start) / max(1, numel(tr));
    fs = 1e6/dt_us; t_ms = (0:numel(tr)-1)'/fs*1e3 + double(win_start-start_us)/1e3;
    trace = tr;
catch, end
end


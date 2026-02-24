function ufo_edit_master(patient)
% UFO_EDIT_MASTER  One-stop UFO editing + logging + figure export by ID.
%
% Usage:
%   ufo_edit_master('pat079')
%
% Workflow:
%   - Loads prep_basics.mat (Dat, UFOdat) and optional existing UFOdat_edited/UFOlog.
%   - Creates/maintains UFOlog with stable integer ID per UFO.
%   - Auto-merges overlapping UFOs per channel (union window, keep longer).
%   - Lets you:
%       * list UFOs
%       * edit single UFO by ID:
%           - classify (Type 1/2)
%           - trim (clicks)
%           - keep/discard
%       * regenerate that UFO's PNG live (Desktop/UFOfigs_by_id).
%
%   - All information is stored in:
%       Data/<pat>/UFO_edited/UFOdat_edited.mat
%       Data/<pat>/UFO_edited/UFOlog.mat
%
%   - PNGs are named: patXXX_ufoID####_<status>_<type>.png
%       status: kept / disc
%       type: t0 / t1L / t1H / t2L / t2H
%

if nargin < 1 || isempty(patient)
    error('ufo_edit_master: patient ID required, e.g. ufo_edit_master(''pat079'')');
end

[base_dir, matmef_path] = ufo_config();
res_dir   = fullfile(base_dir, 'Results', patient);
data_dir  = fullfile(base_dir, 'Data',    patient, 'UFO_edited');
if ~exist(data_dir,'dir'), mkdir(data_dir); end

prep_file   = fullfile(res_dir,  'prep_basics.mat');
edited_file = fullfile(data_dir, 'UFOdat_edited.mat');
log_file    = fullfile(data_dir, 'UFOlog.mat');

if ~exist(prep_file,'file')
    error('ufo_edit_master: prep_basics.mat not found at %s', prep_file);
end

Sprep = load(prep_file,'Dat','UFOdat');
Dat   = Sprep.Dat;
UFO0  = Sprep.UFOdat;   % original UFOs (cell array)
clear Sprep;

% Ensure matmef on path for readMef3
if exist(matmef_path,'dir'), addpath(matmef_path); end

%% ------------------------------------------------------------------------
% 1) Load or initialise UFOdat_edited + UFOlog
%% ------------------------------------------------------------------------

if exist(edited_file,'file')
    Sed = load(edited_file,'UFOdat_edited');
    UFOed = Sed.UFOdat_edited;
else
    % Initialise from original UFO0: copy [start,end,freq] and add type=0
    UFOed = cell(size(UFO0));
    for ch = 1:numel(UFO0)
        U = UFO0{ch};
        if isempty(U), continue; end
        if size(U,2) == 3
            type_col = zeros(size(U,1),1,'int8');
            UFOed{ch} = [U, type_col];
        elseif size(U,2) >= 4
            UFOed{ch} = U(:,1:4);
        end
    end
end

% Ensure 4 columns
for ch = 1:numel(UFOed)
    U = UFOed{ch};
    if isempty(U), continue; end
    if size(U,2) == 3
        type_col = zeros(size(U,1),1,'int8');
        UFOed{ch} = [U, type_col];
    elseif size(U,2) > 4
        UFOed{ch} = U(:,1:4);
    end
end

% Load or create UFOlog
Log = [];
if exist(log_file,'file')
    Slog = load(log_file,'UFOlog');
    Log  = Slog.UFOlog;
end

if isempty(Log)
    % Create fresh log from UFOed and original UFO0
    fprintf('ufo_edit_master: Initialising UFOlog from current edited set...\n');
    Log = struct('id',{},'chan',{},'orig_start_us',{},'orig_end_us',{}, ...
                 'orig_freq_hz',{},'cur_start_us',{},'cur_end_us',{}, ...
                 'cur_freq_hz',{},'type',{},'status',{});
    next_id = 1;
    for ch = 1:numel(UFOed)
        Ue = UFOed{ch};
        if isempty(Ue), continue; end
        for i = 1:size(Ue,1)
            row = Ue(i,:);
            % Best-effort match to original for orig_* (may be equal)
            orig_row = row(1,1:3);
            Log(end+1) = struct( ...
                'id',           next_id, ...
                'chan',         ch, ...
                'orig_start_us',double(orig_row(1)), ...
                'orig_end_us',  double(orig_row(2)), ...
                'orig_freq_hz', double(orig_row(3)), ...
                'cur_start_us', double(row(1)), ...
                'cur_end_us',   double(row(2)), ...
                'cur_freq_hz',  double(row(3)), ...
                'type',         int8(row(4)), ...
                'status',       'kept'); %#ok<AGROW>
            next_id = next_id + 1;
        end
    end
else
    % Ensure consistency (e.g. types) from UFOed into Log.cur_*
    fprintf('ufo_edit_master: Loaded existing UFOlog with %d entries\n', numel(Log));
    for k = 1:numel(Log)
        ch = Log(k).chan;
        if ch > numel(UFOed) || isempty(UFOed{ch}), continue; end
        % Find closest by start time
        U = UFOed{ch};
        diffs = abs(double(U(:,1)) - Log(k).cur_start_us);
        [dmin, ei] = min(diffs);
        if dmin < 5e3 % 5 ms tolerance
            row = U(ei,:);
            Log(k).cur_start_us = double(row(1));
            Log(k).cur_end_us   = double(row(2));
            Log(k).cur_freq_hz  = double(row(3));
            Log(k).type         = int8(row(4));
            if strcmp(Log(k).status,'kept') && row(4)==0
                % still unclassified but kept
            end
        end
    end
    if isempty(Log)
        warning('UFOlog loaded but empty; reinitialising.');
    end
end

%% ------------------------------------------------------------------------
% 2) Auto-merge overlaps (once per run)
%% ------------------------------------------------------------------------

[UFOed, Log] = auto_merge_overlaps(UFOed, Log);

%% ------------------------------------------------------------------------
% 3) Main interactive loop
%% ------------------------------------------------------------------------

out_root = fullfile(home_dir,'Desktop','UFOfigs_by_id');
if ~exist(out_root,'dir'), mkdir(out_root); end
% Ensure all subfolders exist
mkdir(fullfile(out_root,'type1'));
mkdir(fullfile(out_root,'type2'));
mkdir(fullfile(out_root,'discarded'));

done = false;
valid_chans = 1:numel(Dat.allmi);
while ~done
    fprintf('\n=== UFO EDIT MASTER (%s) ===\n', patient);
    fprintf('Total UFOs (log): %d\n', numel(Log));
    % Count only valid channels (in Dat.allmi)
    valid_mask = ismember([Log.chan], valid_chans);
    n_kept_valid  = sum(strcmp({Log.status},'kept') & valid_mask(:));
    n_disc_valid  = sum(strcmp({Log.status},'discarded') & valid_mask(:));
    n_invalid = sum(~valid_mask);
    fprintf('  Kept (valid): %d | Discarded (valid): %d', n_kept_valid, n_disc_valid);
    if n_invalid > 0
        fprintf(' | Invalid channels: %d', n_invalid);
    end
    fprintf('\n');
    fprintf('\nCommands:\n');
    fprintf('  S  = step through all UFOs (sequential classify)\n');
    fprintf('  L  = list first N UFOs\n');
    fprintf('  G  = goto UFO by ID (edit)\n');
    fprintf('  E  = export all figs by ID\n');
    fprintf('  Q  = quit (save & exit)\n');
    cmd = lower(strtrim(input('Choice [S/L/G/E/Q]: ','s')));
    if isempty(cmd), cmd = 's'; end

    switch cmd
        case 's'
            % Sequential pass over all kept UFOs
            [UFOed, Log] = sequential_pass(patient, UFOed, Log, Dat, out_root);

        case 'l'
            valid_chans = 1:numel(Dat.allmi);
            N = min(20, numel(Log));
            fprintf(' id   ch  status   type   freq(Hz)   t_start(s)\n');
            shown = 0;
            for k = 1:numel(Log)
                if shown >= N, break; end
                Lk = Log(k);
                if ~ismember(Lk.chan, valid_chans)
                    continue; % skip invalid channels
                end
                fprintf('%4d  %2d  %-8s  %2d    %7.1f   %.6f\n', ...
                    Lk.id, Lk.chan, Lk.status, Lk.type, ...
                    Lk.cur_freq_hz, Lk.cur_start_us/1e6);
                shown = shown + 1;
            end

        case 'g'
            id = str2double(input('  Enter UFO ID: ','s'));
            if isnan(id) || ~ismember(id,[Log.id])
                fprintf('  → Invalid ID.\n'); continue;
            end
            k = find([Log.id]==id,1);
            if ~ismember(Log(k).chan, valid_chans)
                fprintf('  → WARNING: Channel %d is not in Dat.allmi (only %d valid channels).\n', ...
                    Log(k).chan, numel(Dat.allmi));
                fprintf('  → Plot may be empty, but you can still classify/discard this UFO.\n');
            end
            [UFOed, Log] = edit_single_ufo_by_id(patient, id, UFOed, Log, Dat, out_root);

        case 'e'
            export_all_figs_by_id(patient, Dat, UFOed, Log, out_root);

        case 'q'
            done = true;

        otherwise
            fprintf('  → Unknown command.\n');
    end

    % Save state after each top-level command
    UFOdat_edited = UFOed; %#ok<NASGU>
    UFOlog        = Log;   %#ok<NASGU>
    save(edited_file,'UFOdat_edited','-v7.3');
    save(log_file,   'UFOlog',       '-v7.3');
    fprintf('  ✓ State saved.\n');
end

fprintf('ufo_edit_master: done.\n');

end  % main


%% ===== Helper: auto-merge overlaps =====================================
function [UFOed, Log] = auto_merge_overlaps(UFOed, Log)
fprintf('\n[auto-merge] Checking for overlapping UFOs per channel...\n');
for ch = 1:numel(UFOed)
    U = UFOed{ch};
    if isempty(U) || size(U,1) < 2, continue; end
    % Sort by start
    [~, ord] = sort(double(U(:,1)));
    U = U(ord,:);
    changed = false;

    j = 1;
    while j < size(U,1)
        s1 = double(U(j,1)); e1 = double(U(j,2));
        s2 = double(U(j+1,1)); e2 = double(U(j+1,2));
        if e1 > s2  % overlap
            % Union window
            new_start = min(s1,s2);
            new_end   = max(e1,e2);
            dur1 = e1-s1; dur2 = e2-s2;
            if dur1 >= dur2
                new_freq = double(U(j,3));
                new_type = U(j,4);
            else
                new_freq = double(U(j+1,3));
                new_type = U(j+1,4);
            end
            if new_type<0 || new_type>4, new_type = 0; end
            % Update row j, delete j+1
            U(j,1) = new_start;
            U(j,2) = new_end;
            U(j,3) = new_freq;
            U(j,4) = new_type;
            U(j+1,:) = [];
            changed = true;
            % In log: mark later-starting ID as discarded_overlap
            ids = find([Log.chan]==ch);
            % approximate matching by start times
            for k = ids
                if abs(Log(k).cur_start_us - s2) < 5e3 % 5 ms
                    Log(k).status = 'discarded';
                end
            end
        else
            j = j+1;
        end
    end
    if changed
        UFOed{ch} = U;
    end
end
fprintf('[auto-merge] Done.\n');
end

%% ===== Helper: sequential pass over all UFOs ============================
function [UFOed, Log] = sequential_pass(patient, UFOed, Log, Dat, out_root)
% SEQUENTIAL_PASS  Step through all "kept" UFOs in ID order.
% Skips channels not in Dat.allmi (invalid channels).

ids = [Log.id];
[~, order] = sort(ids);
valid_chans = 1:numel(Dat.allmi);  % Only process channels that exist in Dat.allmi
skipped = 0;
for idx = 1:numel(order)
    k = order(idx);
    if ~strcmp(Log(k).status,'kept')
        continue; % skip discarded
    end
    if ~ismember(Log(k).chan, valid_chans)
        skipped = skipped + 1;
        continue; % skip invalid channels (not in Dat.allmi)
    end
    fprintf('\n[seq] %d/%d  (ID %d, ch %d)\n', idx-skipped, numel(order)-skipped, Log(k).id, Log(k).chan);
    [UFOed, Log, stop_all] = edit_single_ufo_by_id(patient, Log(k).id, UFOed, Log, Dat, out_root);
    if stop_all
        fprintf('[seq] User requested quit. Stopping sequential pass.\n');
        break;
    end
end
if skipped > 0
    fprintf('[seq] Skipped %d UFOs from invalid channels (not in Dat.allmi)\n', skipped);
end
end


%% ===== Helper: edit single UFO by ID ===================================
function [UFOed, Log, stop_all] = edit_single_ufo_by_id(patient, id, UFOed, Log, Dat, out_root)
k  = find([Log.id]==id,1);
Lk = Log(k);
ch = Lk.chan;

fprintf('\n[edit] UFO ID %d, channel %d\n', id, ch);
fprintf('  status = %s, type = %d, freq = %.1f Hz\n', ...
    Lk.status, Lk.type, Lk.cur_freq_hz);
fprintf('  current window: [%.6f %.6f] s\n', ...
    Lk.cur_start_us/1e6, Lk.cur_end_us/1e6);

% Ensure this UFO exists (or recreate if discarded)
U = UFOed{ch};
if strcmp(Lk.status,'discarded')
    % Re-insert into edited list
    row = [Lk.cur_start_us Lk.cur_end_us Lk.cur_freq_hz double(Lk.type)];
    U = [U; int64(row)];
    [~, ord] = sort(double(U(:,1)));
    U = U(ord,:);
    UFOed{ch} = U;
    Lk.status = 'kept';
    Log(k)    = Lk;
end

% Find row index in U again
diffs = abs(double(U(:,1)) - Lk.cur_start_us);
[~, ei] = min(diffs);

% Plot
[trace, t_ms] = read_ufo_window(Dat, ch, Lk.cur_start_us, Lk.cur_end_us);
fig = figure('Name',sprintf('UFO ID %d', id),'Position',[100 100 900 300]);
plot(t_ms, trace,'k-'); hold on;
yl = ylim;
e_start_ms = 0;
e_end_ms   = (Lk.cur_end_us - Lk.cur_start_us)/1e3;
fill([e_start_ms e_end_ms e_end_ms e_start_ms], ...
     [yl(1) yl(1) yl(2) yl(2)], 'r','FaceAlpha',0.15,'EdgeColor','none');
xlabel('Time (ms, rel to UFO start)'); ylabel('Amplitude (\muV)');
title(sprintf('ID %d | Ch %d | freq=%.1f Hz', id, ch, Lk.cur_freq_hz));
grid on;

% Simple edit loop
done = false;
stop_all = false;
while ~done
    fprintf('\nActions for ID %d:\n', id);
    fprintf('  1 = Type 1\n');
    fprintf('  2 = Type 2\n');
    fprintf('  t = trim\n');
    fprintf('  k = keep (no change)\n');
    fprintf('  d = discard\n');
    fprintf('  q = quit editing this UFO\n');
    a = lower(strtrim(input('Choice [1/2/t/k/d/q]: ','s')));
    if isempty(a), a = 'k'; end

    switch a
        case '1'
            % Type 1: auto low/high based on 2kHz
            if Lk.cur_freq_hz < 2000, Lk.type = int8(1);
            else,                    Lk.type = int8(2);
            end
            % Update matrix immediately
            U(ei,4) = Lk.type;
            UFOed{ch} = U;
            fprintf('  → Set type to %d (Type 1)\n', Lk.type);
        case '2'
            if Lk.cur_freq_hz < 2000, Lk.type = int8(3);
            else,                    Lk.type = int8(4);
            end
            % Update matrix immediately
            U(ei,4) = Lk.type;
            UFOed{ch} = U;
            fprintf('  → Set type to %d (Type 2)\n', Lk.type);
        case 't'
            fprintf('  → Trim: click new start and end on the plot.\n');
            [x,~] = ginput(2);
            if numel(x) < 2
                fprintf('  → Need 2 clicks; cancelling trim.\n');
            else
                x = sort(x);
                new_start_us = Lk.cur_start_us + x(1)*1e3;
                new_end_us   = Lk.cur_start_us + x(2)*1e3;
                if new_end_us > new_start_us
                    Lk.cur_start_us = new_start_us;
                    Lk.cur_end_us   = new_end_us;
                    U(ei,1) = int64(new_start_us);
                    U(ei,2) = int64(new_end_us);
                    UFOed{ch} = U;
                    fprintf('  → Trimmed to [%.6f %.6f] s\n', ...
                        new_start_us/1e6, new_end_us/1e6);
                else
                    fprintf('  → Invalid trim; ignoring.\n');
                end
            end
        case 'k'
            fprintf('  → Keeping with current settings.\n');
            done = true;
        case 'd'
            fprintf('  → Discarding UFO ID %d\n', id);
            Lk.status = 'discarded';
            U(ei,:)   = [];
            UFOed{ch} = U;
            done = true;
        case 'q'
            % Quit this UFO and signal caller to stop the whole sequential pass
            fprintf('  → Quit requested. Will stop sequential pass after this UFO.\n');
            stop_all = true;
            done = true;
        otherwise
            fprintf('  → Unknown.\n');
    end
end

close(fig);

% Update log with latest cur_* and type/status (if still kept)
if strcmp(Lk.status,'kept') && ~isempty(U)
    % Refresh from edited matrix
    diffs = abs(double(U(:,1)) - Lk.cur_start_us);
    [~, ei2] = min(diffs);
    row = U(ei2,:);
    Lk.cur_start_us = double(row(1));
    Lk.cur_end_us   = double(row(2));
    Lk.cur_freq_hz  = double(row(3));
    Lk.type         = int8(row(4));
end
Log(k) = Lk;

% Export single fig for this ID
export_one_fig_by_id(Dat, patient, UFOed, Log(k), out_root);

end


%% ===== Helper: export all figs =========================================
function export_all_figs_by_id(patient, Dat, UFOed, Log, out_root)
if exist(out_root,'dir')
    fprintf('[export] Deleting old folder %s\n', out_root);
    rmdir(out_root,'s');
end
mkdir(out_root);
mkdir(fullfile(out_root,'type1'));
mkdir(fullfile(out_root,'type2'));
mkdir(fullfile(out_root,'discarded'));
for k = 1:numel(Log)
    export_one_fig_by_id(Dat, patient, UFOed, Log(k), out_root);
end
fprintf('[export] All figs written under %s\n', out_root);
end


%% ===== Helper: export one fig ==========================================
function export_one_fig_by_id(Dat, patient, UFOed, Lk, out_root)
ch = Lk.chan;
U  = UFOed{ch};
% If the edited matrix for this channel is empty (e.g. fully discarded),
% fall back to the log's current window for plotting/export.
if isempty(U)
    row = int64([Lk.cur_start_us, Lk.cur_end_us, Lk.cur_freq_hz, Lk.type]);
else
% find closest row by start time
diffs = abs(double(U(:,1)) - Lk.cur_start_us);
[~, ei] = min(diffs);
row = U(ei,:);
end

[trace, t_ms] = read_ufo_window(Dat, ch, double(row(1)), double(row(2)));

fig = figure('Visible','off','Position',[100 100 900 300]);
plot(t_ms, trace,'k-'); hold on;
yl = ylim;
e_start_ms = 0;
e_end_ms   = (double(row(2)) - double(row(1)))/1e3;
fill([e_start_ms e_end_ms e_end_ms e_start_ms], ...
     [yl(1) yl(1) yl(2) yl(2)], 'r','FaceAlpha',0.15,'EdgeColor','none');
xlabel('Time (ms, rel to UFO start)'); ylabel('Amplitude (\muV)');
title(sprintf('ID %d | Ch %d | freq=%.1f Hz | status=%s | type=%d', ...
    Lk.id, ch, double(row(3)), Lk.status, Lk.type), 'Interpreter','none');
grid on;

% Build filename and choose subfolder
is_kept = strcmp(Lk.status,'kept');
status_tag = ternary(is_kept,'kept','disc');
switch Lk.type
    case 1, type_tag = 't1L';
    case 2, type_tag = 't1H';
    case 3, type_tag = 't2L';
    case 4, type_tag = 't2H';
    otherwise, type_tag = 't0';
end
fname = sprintf('%s_ufoID%04d_%s_%s.png', patient, Lk.id, status_tag, type_tag);
if ~is_kept
    subdir_path = fullfile(out_root, 'discarded');
else
    % kept: route by type into type1/type2; untyped kept go to type1 by default
    if Lk.type==1 || Lk.type==2 || Lk.type==0
        subdir_path = fullfile(out_root, 'type1');
    else % types 3/4
        subdir_path = fullfile(out_root, 'type2');
    end
end
if exist(subdir_path,'dir') ~= 7
    mkdir(subdir_path);
end
out_path = fullfile(subdir_path, fname);
exportgraphics(fig, out_path, 'Resolution',150);
close(fig);
end


%% ===== Helper: read window from MEF ====================================
function [trace, t_ms] = read_ufo_window(Dat, ch, start_us, end_us)
trace = zeros(100,1); t_ms = linspace(-20,20,100)';
ch_name = local_get_micro_name(Dat,ch);
try
    if isfield(Dat,'meta_path') && ~isempty(Dat.meta_path)
        mef_path = Dat.meta_path;
    else
        mef_path = fullfile(fileparts(Dat.basedir),'sub.mefd');
    end
    if ~exist(mef_path,'dir')
        warning('read_ufo_window: MEF path %s not found for ch %d', mef_path, ch);
        return;
    end
    if isfield(Dat,'password') && ~isempty(Dat.password)
        pw = Dat.password;
    else
        pw = 'bemena';
    end
    pad_us   = 20e3;
    win_start = int64(start_us - pad_us);
    win_end   = int64(end_us   + pad_us);
    if ch > numel(Dat.allmi) || (ch <= numel(Dat.allmi) && isempty(Dat.allmi{ch}))
        warning('read_ufo_window: Channel %d not in Dat.allmi (only %d channels). Using fallback name "%s"', ...
            ch, numel(Dat.allmi), ch_name);
    end
    [~, tr] = readMef3(mef_path, pw, ch_name, 'time', win_start, win_end);
    if isempty(tr)
        warning('read_ufo_window: readMef3 returned empty for ch %d ("%s")', ch, ch_name);
        return;
    end
    tr = double(tr(:));
    dt_us = double(win_end - win_start)/max(1,numel(tr));
    fs    = 1e6/dt_us;
    t_ms  = (0:numel(tr)-1)'/fs*1e3 + double(win_start - start_us)/1e3;
    trace = tr;
catch ME
    warning('read_ufo_window: Error reading MEF for ch %d ("%s"): %s', ch, ch_name, ME.message);
end
end


%% ===== Tiny helpers =====================================================
function name = local_get_micro_name(Dat, idx)
if ~isfield(Dat,'allmi') || idx<1 || idx>numel(Dat.allmi)
    name = sprintf('mic%d',idx); return;
end
mi = Dat.allmi{idx};
if iscell(mi), name = char(mi{1}); else, name = char(mi); end
end

function s = patient_from_Dat(~) %#ok<INUSD>
s = ''; % not needed if you pass patient separately
end

function y = ternary(cond,a,b)
if cond, y = a; else, y = b; end
end
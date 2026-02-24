function export_ufo_figs(patient)
% EXPORT_UFO_FIGS  Export UFO figures by outcome/type for a single patient.
%
% Usage:
%   export_ufo_figs('pat079')
%
% Creates (or recreates) the folder:
%   ~/Desktop/UFOfigs/
% with subfolders:
%   - discarded/     : UFOs present in original UFOdat but NOT in UFOdat_edited
%   - type1/         : kept UFOs with type 1 (<2kHz or >2kHz)
%   - type2/         : kept UFOs with type 2 (<2kHz or >2kHz)
%
% Figures are generated with 'Visible','off' and saved as PNGs.
%

if nargin < 1 || isempty(patient)
    error('export_ufo_figs: patient ID required');
end

% Paths (from ufo_config)
[base_dir, matmef_path] = ufo_config();
basedir    = base_dir;
resultsDir = fullfile(basedir, 'Results', patient);
dataDir    = fullfile(basedir, 'Data',    patient, 'UFO_edited');

prep_file   = fullfile(resultsDir, 'prep_basics.mat');
edited_file = fullfile(dataDir,   'UFOdat_edited.mat');

if exist(prep_file,'file') ~= 2
    error('export_ufo_figs: prep_basics.mat not found at %s', prep_file);
end
if exist(edited_file,'file') ~= 2
    error('export_ufo_figs: UFOdat_edited.mat not found at %s', edited_file);
end

Sprep = load(prep_file, 'Dat', 'UFOdat');
if ~isfield(Sprep, 'Dat') || ~isfield(Sprep, 'UFOdat')
    error('export_ufo_figs: prep_basics.mat must contain Dat and UFOdat');
end
Dat      = Sprep.Dat;
UFO_orig = Sprep.UFOdat;      % original UFOs (before editing)

Sedit = load(edited_file, 'UFOdat_edited');
if ~isfield(Sedit, 'UFOdat_edited')
    error('export_ufo_figs: UFOdat_edited missing from %s', edited_file);
end
UFO_edit = Sedit.UFOdat_edited;

% Ensure cell arrays
if ~iscell(UFO_orig), error('export_ufo_figs: UFOdat in prep_basics must be a cell array'); end
if ~iscell(UFO_edit), error('export_ufo_figs: UFOdat_edited must be a cell array'); end

% Ensure 4 columns in edited UFOs (start_us, end_us, freq, type)
UFO_edit = ensure_4_columns_local(UFO_edit);

% Output folders
root_out = fullfile(getenv('HOME'), 'Desktop', 'UFOfigs');
if exist(root_out, 'dir')
    fprintf('export_ufo_figs: deleting existing folder %s\n', root_out);
    rmdir(root_out, 's');
end
mkdir(root_out);
out_discard = fullfile(root_out, 'discarded'); mkdir(out_discard);
out_type1   = fullfile(root_out, 'type1');     mkdir(out_type1);
out_type2   = fullfile(root_out, 'type2');     mkdir(out_type2);

% Add matmef for readMef3
if exist(matmef_path,'dir')
    addpath(matmef_path);
else
    warning('export_ufo_figs: matmef path %s not found; traces may not load.', matmef_path);
end

% Resolve MEF path
if isfield(Dat, 'meta_path') && ~isempty(Dat.meta_path)
    mef_path = Dat.meta_path;
else
    mef_path = fullfile(basedir, 'Data', patient, 'sub.mefd');
end
if exist(mef_path, 'dir') ~= 7
    warning('export_ufo_figs: MEF directory %s not found. Figures will be blank waveforms.', mef_path);
end

if ~isfield(Dat, 'password') || isempty(Dat.password)
    pw = 'bemena';
else
    pw = Dat.password;
end

% Helper: channel name
get_ch_name = @(idx) local_get_micro_name(Dat, idx);

% --- Build index of kept events for type1/type2 and compute discarded ---

fprintf('export_ufo_figs: building kept/discarded sets for %s...\n', patient);

% For quick matching: per-channel edited events [start,end,freq]
U_edit_SEF = cellfun(@(u) double(u(:,1:3)), UFO_edit, 'UniformOutput', false);

tol_us  = 2000;   % 2 ms tolerance for matching original->edited
tol_hz  = 1e-3;

% Iterate channels
for ch = 1:numel(UFO_orig)
    Uo = UFO_orig{ch};
    if isempty(Uo), continue; end

    Ue = [];
    type_e = [];
    if ch <= numel(UFO_edit)
        Ue = UFO_edit{ch};
    end
    if ~isempty(Ue)
        type_e = Ue(:,4);
    end

    % Prepare matching arrays
    Se_orig = double(Uo(:,1:3));  % [start,end,freq]
    Se_edit = [];
    if ch <= numel(U_edit_SEF)
        Se_edit = U_edit_SEF{ch};
    end

    n_orig = size(Se_orig,1);

    % Determine which original events are "kept" (overlap any edited UFO)
    % NOTE: We treat an original UFO as "kept" if its time window overlaps
    %       any edited UFO on the same channel. This corresponds to
    %       UFOs you *did not* discard with 'd' in the editor. Trimmed or
    %       re-bounded events will still overlap their originals and thus
    %       will NOT be counted as discarded.
    kept_mask = false(n_orig,1);
    for i = 1:n_orig
        s = Se_orig(i,1);
        e = Se_orig(i,2);
        if isempty(Se_edit)
            kept_mask(i) = false;
            continue;
        end
        % Overlap in time if min(e, e_edit) - max(s, s_edit) > 0
        overlap_us = min(e, Se_edit(:,2)) - max(s, Se_edit(:,1));
        kept_mask(i) = any(overlap_us > 0);
    end

    % Discarded events: originals with no match in edited
    disc_idx = find(~kept_mask);

    % Export discarded UFOs for this channel
    for idx = disc_idx(:)'
        row = Uo(idx,:);
        fig_name = sprintf('pat%s_ch%02d_disc_idx%03d.png', patient(end-2:end), ch, idx);
        out_path = fullfile(out_discard, fig_name);
        export_single_ufo_fig(Dat, mef_path, pw, ch, row, out_path);
    end

    % Export kept UFOs by type from edited file (use edited times & trims)
    if isempty(Ue), continue; end
    for j = 1:size(Ue,1)
        row_e = Ue(j,:);
        tval  = row_e(4);
        if tval == 1 || tval == 2   % Type 1
            fig_name = sprintf('pat%s_ch%02d_type1_idx%03d.png', patient(end-2:end), ch, j);
            out_path = fullfile(out_type1, fig_name);
            export_single_ufo_fig(Dat, mef_path, pw, ch, row_e, out_path);
        elseif tval == 3 || tval == 4   % Type 2
            fig_name = sprintf('pat%s_ch%02d_type2_idx%03d.png', patient(end-2:end), ch, j);
            out_path = fullfile(out_type2, fig_name);
            export_single_ufo_fig(Dat, mef_path, pw, ch, row_e, out_path);
        end
    end
end

fprintf('export_ufo_figs: done. Figures saved under %s\n', root_out);

end  % main

% -------------------------------------------------------------------------
function export_single_ufo_fig(Dat, mef_path, pw, micro_idx, row, out_path)
% Plot a single UFO window and save to out_path (figure is not shown).

start_us = double(row(1));
end_us   = double(row(2));
freq_hz  = double(row(3)); %#ok<NASGU> % not used in plotting, kept for context

% Window: 20 ms before, 20 ms after
pad_us = 20e3;
win_start = int64(start_us - pad_us);
win_end   = int64(end_us   + pad_us);

ch_name = local_get_micro_name(Dat, micro_idx);

trace = [];
t_ms  = [];

if exist(mef_path,'dir') == 7
    try
        [~, trace] = readMef3(mef_path, pw, ch_name, 'time', win_start, win_end);
        trace = double(trace(:));
        if ~isempty(trace)
            dt_us = double(win_end - win_start) / max(1, numel(trace));
            fs    = 1e6 / dt_us;
            t_ms  = (0:numel(trace)-1)'/fs*1e3 + double(win_start - start_us)/1e3;
        end
    catch ME
        warning('export_ufo_figs: error reading MEF for %s: %s', ch_name, ME.message);
    end
end

fig = figure('Visible','off', 'Position',[100 100 800 300]);
if ~isempty(trace)
    plot(t_ms, trace, 'k-', 'LineWidth', 1); hold on;
    tr_min = min(trace); tr_max = max(trace);
    if isfinite(tr_min) && isfinite(tr_max) && tr_max > tr_min
        pad = 0.1 * (tr_max - tr_min);
        ylim([tr_min - pad, tr_max + pad]);
    end
else
    plot([0 1],[0 0],'k-'); hold on;
end

% Mark UFO window (relative to first event start)
e_start_ms = 0;
e_end_ms   = (end_us - start_us)/1e3;
yl = ylim;
fill([e_start_ms e_end_ms e_end_ms e_start_ms], ...
     [yl(1) yl(1) yl(2) yl(2)], ...
     'r', 'FaceAlpha', 0.15, 'EdgeColor','none');

xlabel('Time relative to UFO start (ms)');
ylabel('Amplitude (\muV)');
title(sprintf('Ch %d (%s) | UFO [%.3f–%.3f s]', ...
    micro_idx, ch_name, start_us/1e6, end_us/1e6), 'Interpreter','none');
grid on;

[out_dir,~,~] = fileparts(out_path);
if exist(out_dir,'dir') ~= 7
    mkdir(out_dir);
end

try
    exportgraphics(fig, out_path, 'Resolution',150);
catch
    saveas(fig, out_path);
end
close(fig);

end

% -------------------------------------------------------------------------
function UFOdat_out = ensure_4_columns_local(UFOdat_in)
% Local copy of ensure_4_columns logic for UFOdat_edited
UFOdat_out = UFOdat_in;
for ii = 1:numel(UFOdat_out)
    if ~isempty(UFOdat_out{ii})
        n_cols = size(UFOdat_out{ii}, 2);
        if n_cols == 3
            if isinteger(UFOdat_out{ii})
                type_col = zeros(size(UFOdat_out{ii}, 1), 1, 'like', UFOdat_out{ii});
            else
                type_col = zeros(size(UFOdat_out{ii}, 1), 1, 'int8');
            end
            UFOdat_out{ii} = [UFOdat_out{ii}, type_col];
        elseif n_cols > 4
            UFOdat_out{ii} = UFOdat_out{ii}(:,1:4);
        end
    end
end
end

% -------------------------------------------------------------------------
function name = local_get_micro_name(Dat, idx)
if ~isfield(Dat, 'allmi') || idx < 1 || idx > numel(Dat.allmi)
    name = sprintf('mic%d', idx);
    return;
end
mi = Dat.allmi{idx};
if iscell(mi)
    name = char(mi{1});
else
    name = char(mi);
end
end



function step1_inspect(patient)
% STEP1_INSPECT  Interactive UFO inspector.
% Displays UFOs one at a time, lets you discard/keep/trim and classify.
% Saves results to ufo_inspect_log.mat and PNG figures to Results/<pat>/UFO_figures/.

[base_dir, matmef_path] = ufo_config();
if exist(matmef_path,'dir'), addpath(matmef_path); end

if nargin < 1 || isempty(patient)
    error('step1_inspect: patient ID required, e.g. step1_inspect(''pat001'')');
end

% Basic patient lists (extend as needed)
patIDs = {'pat001','pat003','pat004','pat061','pat072','pat079','pat038','pat397','pat448','pat452','pat453'};
p = find(strcmp(patIDs, patient),1);
if isempty(p)
    error('step1_inspect: Patient %s not recognized.', patient);
end

% Micro/macro mappings
mics = cell(numel(patIDs),1);
macs = cell(numel(patIDs),1);

% pat001: AMG micros map to AMG1, others to HH1
mics{1} = {'CSC1', 'CSC2', 'CSC3', 'CSC4', 'CSC5', 'CSC6', 'CSC7', 'CSC8', ...
    'CSC10', 'CSC1_AMG', 'CSC1_HH', 'CSC2_AMG', 'CSC2_HH', 'CSC3_AMG', ...
    'CSC3_HH', 'CSC4_AMG', 'CSC4_HH', 'CSC5_AMG', 'CSC5_HH', 'CSC6_AMG', ...
    'CSC6_HH', 'CSC7_AMG', 'CSC7_HH', 'CSC8_AMG', 'CSC8_HH'};
macs{1} = cell(size(mics{1}));
for i = 1:numel(mics{1})
    if contains(mics{1}{i}, 'AMG')
        macs{1}{i} = 'AMG1';
    else
        macs{1}{i} = 'HH1';
    end
end

% pat003: CSC64-71 -> SMAAL1, CSC72-79 -> SMAL1
mics{2} = arrayfun(@(x) sprintf('CSC%d', x), 64:79, 'UniformOutput', false);
macs{2} = cell(size(mics{2}));
for i = 1:numel(mics{2})
    ch_num = str2double(mics{2}{i}(4:end));
    if ch_num >= 64 && ch_num <= 71
        macs{2}{i} = 'SMAAL1';
    elseif ch_num >= 72 && ch_num <= 79
        macs{2}{i} = 'SMAL1';
    else
        error('Unexpected channel number for pat003: %s', mics{2}{i});
    end
end

% pat004: CSC64-71 -> AMG1, CSC72-79 -> HH1, CSC80-87 -> HB1, CSC96-105 -> TP1
mics{3} = arrayfun(@(x) sprintf('CSC%d', x), [64:71, 72:79, 80:87, 96:105], 'UniformOutput', false);
macs{3} = cell(size(mics{3}));
for i = 1:numel(mics{3})
    ch_num = str2double(mics{3}{i}(4:end));
    if ch_num >= 64 && ch_num <= 71
        macs{3}{i} = 'AMG1';
    elseif ch_num >= 72 && ch_num <= 79
        macs{3}{i} = 'HH1';
    elseif ch_num >= 80 && ch_num <= 87
        macs{3}{i} = 'HB1';
    elseif ch_num >= 96 && ch_num <= 105
        macs{3}{i} = 'TP1';
    else
        error('Unexpected channel number for pat004: %s', mics{3}{i});
    end
end

% pat061, pat072, pat079: Bm*a micros map to multiple macros (B-group)
% Each micro appears twice (once for each macro it maps to)
% NOTE: We keep macro names here as B1/B2/B3 for logging and stats;
% step3_statistics translates them to B'1/B'2/B'3 when reading MEF.
micro_names = {'Bm1a', 'Bm2a', 'Bm3a', 'Bm4a', 'Bm5a', 'Bm6a'};
macro_pairs = {{'B1','B2'}, {'B1','B2'}, {'B1','B2'}, {'B2','B3'}, {'B2','B3'}, {'B2','B3'}};

% pat061 (index 4)
mics{4} = {};
macs{4} = {};
for i = 1:numel(micro_names)
    for j = 1:numel(macro_pairs{i})
        mics{4}{end+1} = micro_names{i};
        macs{4}{end+1} = macro_pairs{i}{j};
    end
end

% pat072 (index 5) - same mapping
mics{5} = mics{4};
macs{5} = macs{4};

% pat079 (index 6) - same mapping
mics{6} = mics{4};
macs{6} = macs{4};

% pat038 (index 7): Mayo clinic - structure similar to pat453 but with LShaftMicro/RShaftMicro
% LBundle_01-08 -> LMacro_01
mics{7} = arrayfun(@(x) sprintf('LBundle_%02d', x), 1:8, 'UniformOutput', false);
macs{7} = repmat({'LMacro_01'}, 1, 8);
% RBundle_01-08 (missing 06) -> RMacro_01
mics{7} = [mics{7}, arrayfun(@(x) sprintf('RBundle_%02d', x), [1:5 7:8], 'UniformOutput', false)];
macs{7} = [macs{7}, repmat({'RMacro_01'}, 1, 7)];
% LShaftMicro_01-6 -> LMacro_01 and LMacro_02
for i = 1:6
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_01';
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_02';
end
% LShaftMicro_07-12 -> LMacro_02 and LMacro_03
for i = 7:12
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_02';
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_03';
end
% LShaftMicro_13-18 -> LMacro_03 and LMacro_04
for i = 13:18
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_03';
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_04';
end
% LShaftMicro_19-24 -> LMacro_04 only (if needed, or map to LMacro_04 and beyond)
for i = 19:24
    mics{7}{end+1} = sprintf('LShaftMicro_%02d', i);
    macs{7}{end+1} = 'LMacro_04';
end
% RShaftMicro_01-6 -> RMacro_01 and RMacro_02
for i = 1:6
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_01';
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_02';
end
% RShaftMicro_07-12 -> RMacro_02 and RMacro_03
for i = 7:12
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_02';
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_03';
end
% RShaftMicro_13-18 -> RMacro_03 and RMacro_04
for i = 13:18
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_03';
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_04';
end
% RShaftMicro_19-24 -> RMacro_04 only
for i = 19:24
    mics{7}{end+1} = sprintf('RShaftMicro_%02d', i);
    macs{7}{end+1} = 'RMacro_04';
end

% pat397 (index 8): Mayo clinic - similar structure to pat453: LBundle/RBundle and LShaft/RShaft
% RBundle_01-07 -> RMacro_01
mics{8} = arrayfun(@(x) sprintf('RBundle_%02d', x), 1:7, 'UniformOutput', false);
macs{8} = repmat({'RMacro_01'}, 1, 7);
% LBundle_01-08 -> LMacro_01
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LBundle_%02d', x), 1:8, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_01'}, 1, 8)];
% RShaft_01-3 -> RMacro_01 and RMacro_02
mics{8} = [mics{8}, arrayfun(@(x) sprintf('RShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'RMacro_01'}, 1, 3)];
mics{8} = [mics{8}, arrayfun(@(x) sprintf('RShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'RMacro_02'}, 1, 3)];
% RShaft_05-08 -> RMacro_02 and RMacro_03 (note: RShaft_04 missing in JSON)
mics{8} = [mics{8}, arrayfun(@(x) sprintf('RShaft_%02d', x), 5:8, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'RMacro_02'}, 1, 4)];
mics{8} = [mics{8}, arrayfun(@(x) sprintf('RShaft_%02d', x), 5:8, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'RMacro_03'}, 1, 4)];
% LShaft_01-3 -> LMacro_01 and LMacro_02
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_01'}, 1, 3)];
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_02'}, 1, 3)];
% LShaft_04-6 -> LMacro_02 and LMacro_03
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LShaft_%02d', x), 4:6, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_02'}, 1, 3)];
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LShaft_%02d', x), 4:6, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_03'}, 1, 3)];
% LShaft_07-09 -> LMacro_03 and LMacro_04
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LShaft_%02d', x), 7:9, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_03'}, 1, 3)];
mics{8} = [mics{8}, arrayfun(@(x) sprintf('LShaft_%02d', x), 7:9, 'UniformOutput', false)];
macs{8} = [macs{8}, repmat({'LMacro_04'}, 1, 3)];

% pat448 (index 9): Each ADShaft micro maps to multiple macros (similar to pat061/pat072/pat079)
% ADShaft_01-04 are between PDMacro_01 and PDMacro_02 -> map to both
% ADShaft_05-08 are between PDMacro_02 and PDMacro_03 -> map to both
% ADShaft_09-12 are between PDMacro_03 and PDMacro_04 -> map to both
% ADShaft_13-16 are positioned behind PDMacro_04 -> map to PDMacro_04 only
mics{9} = {};
macs{9} = {};
% ADShaft_01-04 -> PDMacro_01 and PDMacro_02
for i = 1:4
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_01';
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_02';
end
% ADShaft_05-08 -> PDMacro_02 and PDMacro_03
for i = 5:8
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_02';
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_03';
end
% ADShaft_09-12 -> PDMacro_03 and PDMacro_04
for i = 9:12
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_03';
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_04';
end
% ADShaft_13-16 -> PDMacro_04 only
for i = 13:16
    mics{9}{end+1} = sprintf('ADShaft_%02d', i);
    macs{9}{end+1} = 'PDMacro_04';
end

% pat452 (index 10): Mayo clinic - similar structure to pat453: LBundle/RBundle
% RBundle_01-08 -> RMacro_01 (note: some may be missing in JSON)
mics{10} = arrayfun(@(x) sprintf('RBundle_%02d', x), 1:8, 'UniformOutput', false);
macs{10} = repmat({'RMacro_01'}, 1, 8);
% LBundle_01-08 -> LMacro_01
mics{10} = [mics{10}, arrayfun(@(x) sprintf('LBundle_%02d', x), 1:8, 'UniformOutput', false)];
macs{10} = [macs{10}, repmat({'LMacro_01'}, 1, 8)];

% pat453 (index 11): Extended mappings
% RBundle_01-08 -> RMacro_01
mics{11} = arrayfun(@(x) sprintf('RBundle_%02d', x), 1:8, 'UniformOutput', false);
macs{11} = repmat({'RMacro_01'}, 1, 8);
% LBundle_01-08 -> LMacro_01
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LBundle_%02d', x), 1:8, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_01'}, 1, 8)];
% RShaft_01-3 -> RMacro_01 and RMacro_02
mics{11} = [mics{11}, arrayfun(@(x) sprintf('RShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'RMacro_01'}, 1, 3)];
mics{11} = [mics{11}, arrayfun(@(x) sprintf('RShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'RMacro_02'}, 1, 3)];
% RShaft_04-6 -> RMacro_02 and RMacro_03
mics{11} = [mics{11}, arrayfun(@(x) sprintf('RShaft_%02d', x), 4:6, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'RMacro_02'}, 1, 3)];
mics{11} = [mics{11}, arrayfun(@(x) sprintf('RShaft_%02d', x), 4:6, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'RMacro_03'}, 1, 3)];
% RShaft_07-9 -> RMacro_03 and RMacro_04
mics{11} = [mics{11}, arrayfun(@(x) sprintf('RShaft_%02d', x), 7:9, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'RMacro_03'}, 1, 3)];
mics{11} = [mics{11}, arrayfun(@(x) sprintf('RShaft_%02d', x), 7:9, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'RMacro_04'}, 1, 3)];
% LShaft_01-3 -> LMacro_01 and LMacro_02
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_01'}, 1, 3)];
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LShaft_%02d', x), 1:3, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_02'}, 1, 3)];
% LShaft_04-6 -> LMacro_02 and LMacro_03
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LShaft_%02d', x), 4:6, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_02'}, 1, 3)];
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LShaft_%02d', x), 4:6, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_03'}, 1, 3)];
% LShaft_07-9 -> LMacro_03 and LMacro_04
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LShaft_%02d', x), 7:9, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_03'}, 1, 3)];
mics{11} = [mics{11}, arrayfun(@(x) sprintf('LShaft_%02d', x), 7:9, 'UniformOutput', false)];
macs{11} = [macs{11}, repmat({'LMacro_04'}, 1, 3)];

if isempty(mics{p})
    error('step1_inspect: No micro/macro mapping defined for %s.', patient);
end

mic_list = mics{p};
mac_list = macs{p};

% Get unique micro channels for searching (avoid duplicates from multi-macro mappings)
[unique_mics, unique_idx] = unique(mic_list, 'stable');
% Create a mapping from unique micro to all its associated macros
mic_to_macros = cell(numel(unique_mics), 1);
for i = 1:numel(unique_mics)
    mic_to_macros{i} = mac_list(strcmp(mic_list, unique_mics{i}));
end

% Paths (base_dir from ufo_config)
data_dir = fullfile(base_dir, 'Data', patient);
res_dir  = fullfile(base_dir, 'Results', patient);
if ~exist(res_dir,'dir'), mkdir(res_dir); end

ufo_json_path = fullfile(data_dir, 'UFO.json');

% Auto-detect MEF directory (may be named sub.mefd, sub-*.mefd, or *.mefd)
mef_path = fullfile(data_dir, 'sub.mefd');
if ~exist(mef_path, 'dir')
    % Try to find any .mefd directory
    dir_list = dir(data_dir);
    mefd_dirs = dir_list([dir_list.isdir] & endsWith({dir_list.name}, '.mefd'));
    if ~isempty(mefd_dirs)
        mef_path = fullfile(data_dir, mefd_dirs(1).name);
        fprintf('[step1_inspect] Found MEF directory: %s\n', mefd_dirs(1).name);
    else
        error('step1_inspect: No .mefd directory found in %s', data_dir);
    end
end

pw            = 'bemena';
log_path      = fullfile(res_dir, 'ufo_inspect_log.mat');

% Performance warning for pat448
if strcmp(patient, 'pat448')
    fprintf('\n[WARNING] pat448 has 203 channels in MEF. readMef3 may be very slow.\n');
    fprintf('          Each channel read may take 1-5 minutes. Please be patient.\n');
    fprintf('          If it hangs for >10 minutes, interrupt (Ctrl+C) and check MEF file integrity.\n\n');
end

% Output folders for figures (only in Results)
res_fig_dir = fullfile(res_dir, 'UFO_figures');
out_disc = fullfile(res_fig_dir, 'discard');
out_t1   = fullfile(res_fig_dir, 'type1');
out_t2   = fullfile(res_fig_dir, 'type2');
if ~exist(res_fig_dir,'dir'), mkdir(res_fig_dir); end
if ~exist(out_disc,'dir'),    mkdir(out_disc);    end
if ~exist(out_t1,'dir'),      mkdir(out_t1);      end
if ~exist(out_t2,'dir'),      mkdir(out_t2);      end

% Initialise or load results log
Results = struct('patient',{},'micro_channel',{},'macro_channel',{},'ufo_idx',{}, ...
    'uuid',{},'orig_start_us',{},'orig_end_us',{},'orig_freq_hz',{}, ...
    'freq_class',{},'decision',{},'new_start_us',{},'new_end_us',{},'type',{});

if exist(log_path,'file')
    Slog = load(log_path,'Results');
    Results = Slog.Results;
    fprintf('Loaded existing log for %s with %d entries.\n', patient, numel(Results));
end

% Choose mode
mode = lower(strtrim(input('\nMode: (a)ll/resume or (s)ingle UFO edit? [a/s, default a]: ','s')));
if isempty(mode), mode = 'a'; end

if strcmp(mode,'s')
    edit_single_ufo_reedit(patient, mics{p}, macs{p}, ufo_json_path, mef_path, pw, log_path, res_fig_dir);
    return;
end

% --- Main Inspector Loop (All/Resume) ---
% Check if UFO.json exists
if ~exist(ufo_json_path, 'file')
    error('step1_inspect: UFO.json not found at %s\nPlease ensure UFO detection has been run for this patient.', ufo_json_path);
end

txt = fileread(ufo_json_path);
if isempty(txt) || isempty(strtrim(txt))
    error('step1_inspect: UFO.json is empty at %s\nPlease ensure UFO detection has been run for this patient.', ufo_json_path);
end

raw = jsondecode(txt);
if isempty(raw)
    error('step1_inspect: UFO.json contains no UFOs for %s\nPlease ensure UFO detection has been run for this patient.', patient);
end

fprintf('\n=== %s UFO INSPECTOR ===\n', patient);

% Diagnostic: show what channels are actually in the JSON file
if numel(raw) > 0
    all_json_channels = {};
    for i = 1:min(100, numel(raw))  % Check first 100 UFOs to get channel names
        chs = raw(i).channels;
        if iscell(chs)
            all_json_channels = [all_json_channels, chs];
        else
            all_json_channels{end+1} = chs;
        end
    end
    unique_json_channels = unique(all_json_channels);
    fprintf('Sample channels found in JSON (first 100 UFOs): %s\n', strjoin(unique_json_channels(1:min(20, numel(unique_json_channels))), ', '));
    if numel(unique_json_channels) > 20
        fprintf('  ... and %d more unique channels\n', numel(unique_json_channels) - 20);
    end
end

% Calculate how many UFO entries remain to inspect for this patient
total_pending_patient = 0;
for m = 1:numel(unique_mics)
    micro = unique_mics{m};
    is_array = numel(raw);
    keep_idx = false(is_array,1);
    for i = 1:is_array
        chs = raw(i).channels;
        if iscell(chs), keep_idx(i) = any(strcmp(chs, micro));
        else,           keep_idx(i) = strcmp(chs, micro); end
    end
    idx_all = find(keep_idx);
    % Remove any that are already in the Results log
    if ~isempty(Results)
        mask_pat = strcmp({Results.patient}, patient);
        mask_mic = strcmp({Results.micro_channel}, micro);
        done_idx = [Results(mask_pat & mask_mic).ufo_idx];
        idx_all = setdiff(idx_all, done_idx);
    end
    total_pending_patient = total_pending_patient + numel(idx_all);
end

% Check if there are any UFOs to inspect
if total_pending_patient == 0
    fprintf('\nNo UFOs found to inspect for %s.\n', patient);
    fprintf('This could mean:\n');
    fprintf('  1. All UFOs have already been inspected (check %s)\n', log_path);
    fprintf('  2. UFO.json does not contain UFOs for the mapped micro channels\n');
    fprintf('  3. UFO detection has not been run for this patient\n');
    fprintf('\nTotal UFOs in JSON file: %d\n', numel(raw));
    fprintf('Mapped micro channels (unique): %s\n', strjoin(unique_mics, ', '));
    return;
end

fprintf('Total UFOs remaining to inspect: %d\n\n', total_pending_patient);

for m = 1:numel(unique_mics)
    micro = unique_mics{m};
    % Get all macros associated with this micro
    associated_macros = mic_to_macros{m};

    is_array = numel(raw);
    keep_idx = false(is_array,1);
    for i = 1:is_array
        chs = raw(i).channels;
        if iscell(chs), keep_idx(i) = any(strcmp(chs, micro));
        else,           keep_idx(i) = strcmp(chs, micro); end
    end
    idx_list = find(keep_idx);

    % Skip already done
    if ~isempty(Results)
        mask_pat = strcmp({Results.patient}, patient);
        mask_mic = strcmp({Results.micro_channel}, micro);
        done_idx = [Results(mask_pat & mask_mic).ufo_idx];
        idx_list = setdiff(idx_list, done_idx);
    end

    n_ufo = numel(idx_list);
    if n_ufo == 0, continue; end
    macro_str = strjoin(associated_macros, '/');
    fprintf('  micro %s (macros: %s): %d remaining UFOs\n', micro, macro_str, n_ufo);

    for k = 1:n_ufo
        u_idx = idx_list(k);
        u     = raw(u_idx);
        start_us = double(u.uutc_left);
        end_us   = double(u.uutc_right);
        freq_hz  = double(u.frequency);
        dur_ms   = (end_us - start_us) / 1e3;
        start_s  = start_us / 1e6;

        fprintf('    UFO %d/%d (JSON idx %d): start=%.6f s, dur=%.3f ms, freq=%.1f Hz\n', ...
            k, n_ufo, u_idx, start_s, dur_ms, freq_hz);

        % Diagnostic output for slow reads (especially pat448)
        if strcmp(patient, 'pat448')
            fprintf('    [Reading MEF data for channel %s...] ', micro);
        end
        [trace, t_ms] = read_ufo_window_simple(mef_path, pw, micro, start_us, end_us);
        if strcmp(patient, 'pat448')
            fprintf('Done.\n');
        end
        if ~isempty(trace)
            % Calculate how many UFOs remain to inspect (including current one)
            if ~isempty(Results)
                current_processed = sum(strcmp({Results.patient}, patient));
            else
                current_processed = 0;
            end
            remaining_incl_current = total_pending_patient - current_processed;
            
            figure(1); clf; set(gcf,'position',[0 387 1024 150])
            plot(t_ms, trace, 'k-'); hold on;
            highlight_mask = t_ms >= 0 & t_ms <= dur_ms;
            plot(t_ms(highlight_mask), trace(highlight_mask), 'r-', 'LineWidth', 1.5);
            xlabel('Time (ms)'); ylabel('uV');
            title(sprintf('%s | %s->%s | Micro: %d/%d | Remaining (incl this): %d | %.1f Hz', ...
                patient, micro, macro_str, k, n_ufo, remaining_incl_current, freq_hz));
            xlim([min(t_ms) max(t_ms)]);
        else
            fprintf('      (No MEF trace available.)\n');
        end

        % Decision
        decision = '';
        while isempty(decision)
            resp = lower(strtrim(input('      Choice [d=disc, k=keep, t=trim, s=skip, q=quit]: ','s')));
            if any(strcmp(resp, {'d','k','t','s','q'})), decision = resp; end
        end

        if strcmp(decision,'q'), return; end

        new_start_us = start_us; new_end_us = end_us; final_decision = 'skip';
        if strcmp(decision,'t')
            fprintf('      Click twice: START then END.\n');
            if ~isempty(trace)
                [x_click, ~] = ginput(2);
                if numel(x_click)==2
                    off_s = min(x_click); off_e = max(x_click);
                    off_s = max(0, min(dur_ms, off_s)); off_e = max(0, min(dur_ms, off_e));
                    new_start_us = start_us + off_s*1e3; new_end_us = start_us + off_e*1e3;
                end
            end
            final_decision = 'keep_trimmed';
        elseif strcmp(decision,'k'), final_decision = 'keep';
        elseif strcmp(decision,'d'), final_decision = 'discard';
        end

        ufo_type = 0;
        if ~strcmp(final_decision,'discard') && ~strcmp(final_decision,'skip')
            t_resp = str2double(strtrim(input('      Type [1/2/0]: ','s')));
            if ismember(t_resp,[0 1 2]), ufo_type = t_resp; end
        end

        freq_class = 'low'; if freq_hz >= 2000, freq_class = 'high'; end

        % Save PNG (use first macro for filename)
        primary_macro = associated_macros{1};
        if ~isempty(trace)
            id_str = sprintf('%s_%s_%s_idx%04d', patient, micro, primary_macro, u_idx);
            if strcmp(final_decision,'discard'), subdir = out_disc;
            elseif ufo_type == 1, subdir = out_t1;
            elseif ufo_type == 2, subdir = out_t2;
            else, subdir = res_fig_dir; end
            exportgraphics(gcf, fullfile(subdir, [id_str '.png']), 'Resolution',150);
        end

        % Log: create one entry per macro this micro maps to
        for mac_idx = 1:numel(associated_macros)
            macro = associated_macros{mac_idx};
            R = struct('patient',patient,'micro_channel',micro,'macro_channel',macro, ...
                'ufo_idx',u_idx,'uuid',char(u.uuid),'orig_start_us',start_us, ...
                'orig_end_us',end_us,'orig_freq_hz',freq_hz,'freq_class',freq_class, ...
                'decision',final_decision,'new_start_us',new_start_us,'new_end_us',new_end_us, ...
                'type',ufo_type);
            Results(end+1) = R; %#ok<AGROW>
        end
        save(log_path,'Results');
    end
end
end

function edit_single_ufo_reedit(patient, mic_list, mac_list, ufo_json_path, mef_path, pw, log_path, res_fig_dir)
    S = load(log_path,'Results'); Results = S.Results;
    micro = strtrim(input('Micro channel: ','s'));
    macro = strtrim(input('Macro channel: ','s'));
    u_idx = str2double(strtrim(input('UFO idx in JSON: ','s')));
    
    k = find(strcmp({Results.patient},patient) & strcmp({Results.micro_channel},micro) & ...
             strcmp({Results.macro_channel},macro) & [Results.ufo_idx] == u_idx, 1);
    if isempty(k), fprintf('Not found in log.\n'); return; end
    
    raw = jsondecode(fileread(ufo_json_path)); u = raw(u_idx);
    start_us = double(u.uutc_left); end_us = double(u.uutc_right);
    freq_hz = double(u.frequency); dur_ms = (end_us-start_us)/1e3;
    
    [trace, t_ms] = read_ufo_window_simple(mef_path, pw, micro, start_us, end_us);
    if ~isempty(trace)
        figure(1); clf; set(gcf,'position',[0 387 1024 150])
        plot(t_ms, trace, 'k-'); hold on;
        plot(t_ms(t_ms>=0 & t_ms<=dur_ms), trace(t_ms>=0 & t_ms<=dur_ms), 'r-', 'LineWidth', 1.5);
        xlim([min(t_ms) max(t_ms)]); title('RE-EDIT');
    end
    
    resp = lower(strtrim(input('Choice [d/k/t/s/q]: ','s')));
    if strcmp(resp,'q'), return; end
    
    new_start_us = start_us; new_end_us = end_us; final_decision = 'skip';
    if strcmp(resp,'t')
        [xc,~] = ginput(2); if numel(xc)==2
            new_start_us = start_us + max(0,min(dur_ms,min(xc)))*1e3;
            new_end_us   = start_us + max(0,min(dur_ms,max(xc)))*1e3;
        end
        final_decision = 'keep_trimmed';
    elseif strcmp(resp,'k'), final_decision = 'keep';
    elseif strcmp(resp,'d'), final_decision = 'discard';
    end
    
    ufo_type = 0;
    if ~strcmp(final_decision,'discard') && ~strcmp(final_decision,'skip')
        t_resp = str2double(strtrim(input('Type [1/2/0]: ','s')));
        if ismember(t_resp,[0 1 2]), ufo_type = t_resp; end
    end
    
    % Update Log
    Results(k).decision = final_decision; Results(k).new_start_us = new_start_us;
    Results(k).new_end_us = new_end_us; Results(k).type = ufo_type;
    save(log_path,'Results');
    
    % Re-save PNG (cleanup old ones first)
    id_str = sprintf('%s_%s_%s_idx%04d', patient, micro, macro, u_idx);
    subdirs = {'discard','type1','type2','.'};
    for i=1:numel(subdirs), f=fullfile(res_fig_dir, subdirs{i}, [id_str '.png']); if exist(f,'file'), delete(f); end; end
    
    if ~isempty(trace)
        if strcmp(final_decision,'discard'), sd='discard'; elseif ufo_type==1, sd='type1'; elseif ufo_type==2, sd='type2'; else, sd='.'; end
        exportgraphics(gcf, fullfile(res_fig_dir, sd, [id_str '.png']), 'Resolution',150);
    end
    fprintf('Updated.\n');
end

function [trace, t_ms] = read_ufo_window_simple(mef_path, pw, ch_name, start_us, end_us)
% READ_UFO_WINDOW_SIMPLE  Read a time window from MEF for a specific channel.
% For pat448, uses Python/pymef for faster reading (readMef3 is very slow).
% For other patients, uses standard readMef3.
trace = []; t_ms = []; pad_us = 20e3;
win_start = int64(start_us - pad_us); win_end = int64(end_us + pad_us);

% Quick check: verify channel directory exists (helps diagnose issues)
ch_dir = fullfile(mef_path, [ch_name '.timd']);
if ~exist(ch_dir, 'dir')
    warning('read_ufo_window_simple: Channel directory not found: %s', ch_dir);
    return;
end

% Determine patient from mef_path
[data_dir, ~] = fileparts(mef_path);
[~, patient] = fileparts(data_dir);

try
    % For pat448, use Python for much faster reading
    if strcmp(patient, 'pat448')
        [meta, tr] = readMef3_python(mef_path, pw, ch_name, 'time', win_start, win_end);
    else
        % For other patients, use standard readMef3
        [~, tr] = readMef3(mef_path, pw, ch_name, 'time', win_start, win_end);
        meta = [];
    end
    
    if isempty(tr), return; end
    tr = double(tr(:));
    
    % Calculate time vector
    if ~isempty(meta) && isfield(meta, 'sampling_frequency') && meta.sampling_frequency > 0
        % Use sampling frequency from Python
        fs = meta.sampling_frequency;
        dt_us = 1e6 / fs;
    else
        % Fallback: estimate from data length
        dt_us = double(win_end - win_start) / max(1, numel(tr));
        fs = 1e6/dt_us;
    end
    
    t_ms = (0:numel(tr)-1)'/fs*1e3 + double(win_start-start_us)/1e3;
    trace = tr;
catch ME
    % If there's an error, at least report it
    warning('read_ufo_window_simple: Error reading %s: %s', ch_name, ME.message);
end
end


function edit_ufos_interactive(patient, results_dir, mef_path, password)
% EDIT_UFOS_INTERACTIVE  Interactive UFO review and editing tool
%
% Usage:
%   edit_ufos_interactive('pat001')
%   edit_ufos_interactive('pat001', results_dir, mef_path, password)
%
% For each UFO, displays the micro channel signal and prompts for action:
%   - '1' or 'type1': Classify as Type 1 (genuine) - automatically uses frequency to determine <2kHz or >2kHz
%   - '2' or 'type2': Classify as Type 2 (artifact) - automatically uses frequency to determine <2kHz or >2kHz
%   - 't' or 'trim': Trim UFO (select new start/end times with cursor)
%   - 'k' or 'keep': Keep UFO and move to next
%   - 'd' or 'discard': Remove UFO entirely
%   - 'n' or 'next': Move to next UFO without action
%   - 'b' or 'back': Go back to previous UFO (undo last action)
%   - 'q' or 'quit': Save and exit
%
% You can do multiple actions on the same UFO:
%   Example: Press '1' to classify, then 't' to trim, then 'k' to keep&next
%
% UFO Type Classification (4 categories, automatically determined by frequency):
%   - Type 1, < 2kHz: Genuine UFO, low frequency (if freq < 2000 Hz)
%   - Type 1, > 2kHz: Genuine UFO, high frequency (if freq >= 2000 Hz)
%   - Type 2, < 2kHz: Artifact, low frequency (if freq < 2000 Hz)
%   - Type 2, > 2kHz: Artifact, high frequency (if freq >= 2000 Hz)
%
% Type encoding in 4th column:
%   0 = unclassified
%   1 = Type 1, < 2kHz
%   2 = Type 1, > 2kHz
%   3 = Type 2, < 2kHz
%   4 = Type 2, > 2kHz
%
% Edited UFOs are saved to: BASE_DIR/Data/<patient>/UFO_edited/
%
% Set CERD_EXCLUDE_TYPE2=true to automatically exclude Type 2 UFOs from analysis

if nargin < 1 || isempty(patient)
    error('edit_ufos_interactive: patient ID required');
end

[base_dir, matmef_path] = ufo_config();
if exist(matmef_path,'dir'), addpath(matmef_path); end

% Check for patients with known readMef3 issues
skip_readmef3_patients = getin('CERD_SKIP_READMEF3_PATIENTS', {});
if ischar(skip_readmef3_patients) || isstring(skip_readmef3_patients)
    skip_readmef3_patients = cellstr(skip_readmef3_patients);
end
if any(strcmpi(skip_readmef3_patients, patient)) || any(strcmpi(skip_readmef3_patients, 'ALL'))
    error('edit_ufos_interactive: Patient %s is in CERD_SKIP_READMEF3_PATIENTS list. readMef3 is known to hang for this patient. Use FAST META mode or skip this patient.', patient);
end

% Minimum frequency filter (default 0 = disabled, set > 0 to enable)
% UFOs below this frequency will be automatically discarded
% Set CERD_MIN_UFO_FREQ_HZ = 0 to disable the filter
min_freq_hz = getin('CERD_MIN_UFO_FREQ_HZ', 0);
if ~isscalar(min_freq_hz) || ~isnumeric(min_freq_hz) || min_freq_hz < 0
    warning('edit_ufos_interactive: Invalid CERD_MIN_UFO_FREQ_HZ, using default 0 (disabled)');
    min_freq_hz = 0;
end
if min_freq_hz > 0
    fprintf('Minimum UFO frequency filter: %.1f Hz (UFOs below this will be automatically discarded)\n', min_freq_hz);
    fprintf('  (Set CERD_MIN_UFO_FREQ_HZ = 0 to disable this filter)\n');
else
    fprintf('Minimum UFO frequency filter: DISABLED (all UFOs will be available for review)\n');
    fprintf('  (Set CERD_MIN_UFO_FREQ_HZ > 0 to enable filtering)\n');
end

% Type 2 exclusion filter (default: false, set to true to exclude Type 2 artifacts)
% Type 2 UFOs are thought to be artifacts and can be automatically excluded
exclude_type2 = getin('CERD_EXCLUDE_TYPE2', false);
if exclude_type2
    fprintf('Type 2 exclusion: ENABLED (Type 2 UFOs will be automatically excluded from analysis)\n');
else
    fprintf('Type 2 exclusion: DISABLED (all classified UFOs will be included)\n');
    fprintf('  (Set CERD_EXCLUDE_TYPE2=true to exclude Type 2 artifacts)\n');
end

% Load Dat and UFOdat from prep_basics.mat
Dat = [];
UFOdat = [];

% Determine results_dir
if nargin < 2 || isempty(results_dir)
    results_dir = fullfile(base_dir, 'Results', patient);
end

prep_file = fullfile(results_dir, 'prep_basics.mat');
if exist(prep_file, 'file')
    try
        Prep = load(prep_file);
        if isfield(Prep, 'Dat')
            Dat = Prep.Dat;
        end
        if isfield(Prep, 'UFOdat')
            UFOdat = Prep.UFOdat;
        end
        fprintf('Loaded Dat and UFOdat from %s\n', prep_file);
    catch ME
        error('Failed to load prep_basics.mat: %s', ME.message);
    end
else
    % prep_basics.mat doesn't exist - we need to create it
    % Since prep_basics and extract_UFOs are local functions in cerd_main,
    % the simplest way is to run cerd_main which will create it automatically.
    % However, the user wants to edit UFOs BEFORE running analysis.
    % Solution: run cerd_main in a minimal way that just does prep.
    
    fprintf('prep_basics.mat not found. Creating it now...\n');
    fprintf('(This requires running prep steps from cerd_main)\n\n');
    
    % Ensure results directory exists
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    
    % The cleanest solution: run cerd_main with the patient, which will
    % call prep_basics and extract_UFOs, then we can save prep_basics.mat
    % But we need to do this without running the full analysis.
    % We can do this by running cerd_main and catching early, or by
    % creating a wrapper that calls the prep functions.
    
    % Since prep_basics is local to cerd_main, we'll use a workaround:
    % Run cerd_main with a special setup that stops after prep.
    % Actually, the simplest: just run cerd_main('pat448') which will
    % do prep, then we extract Dat and UFOdat from the base workspace.
    
    fprintf('Running prep steps via cerd_main...\n');
    fprintf('(This will load patient metadata and extract UFOs from UFO.json)\n\n');
    
    try
        % Run cerd_main which will call prep_basics and extract_UFOs
        % We'll set it up to stop after prep by using a mode that just does prep
        % Actually, cerd_main always does prep first, so we can run it and
        % extract the variables before it does heavy computation.
        
        % Check if cerd_main exists
        if ~exist('cerd_main', 'file')
            error('cerd_main.m not found in path.');
        end
        
        % We'll run cerd_main in generate mode with Nnull=0 to just do prep
        % But that requires modifying parameters. Simpler: just run it normally
        % and it will create prep_basics.mat in 'g' mode (line 186-192 of cerd_main)
        
        % Actually, looking at cerd_main, it creates prep_basics.mat in 'g' mode
        % at lines 186-192. So we can just run:
        %   cerd_main(patient, 'g', 1)
        % But that generates nulls which takes time.
        
        % Better: create prep_basics.mat by calling the prep functions directly
        % Since they're local, we need a workaround. Let's use eval to call
        % them from within cerd_main's file scope.
        
        % Simplest immediate solution: provide clear error with instructions
        error(['prep_basics.mat not found for %s.\n\n' ...
            'QUICK FIX - Choose one:\n\n' ...
            'Option 1 (Recommended): Run the helper script:\n' ...
            '  create_prep_for_editing(''%s'')\n\n' ...
            'Option 2: Run cerd_main in generate mode:\n' ...
            '  cerd_main(''%s'', ''g'', 1)\n' ...
            '  (This creates prep_basics.mat and generates 1 null - you can cancel after prep)\n\n' ...
            'Option 3: Run analysis once (also creates prep_basics.mat):\n' ...
            '  cerd_main(''%s'', ''a'')\n\n' ...
            'After any of these, run: edit_ufos_interactive(''%s'')\n'], ...
            patient, patient, patient, patient, patient);
    catch ME
        if contains(ME.message, 'prep_basics.mat not found')
            rethrow(ME);
        else
            error('Error creating prep_basics.mat: %s', ME.message);
        end
    end
end

% Ensure required fields
if isempty(Dat) || isempty(UFOdat)
    error('Dat or UFOdat not found in prep_basics.mat');
end

if ~isfield(Dat, 'allmi') || isempty(Dat.allmi)
    error('Dat.allmi not found');
end

if ~isfield(Dat, 'meta_path') || isempty(Dat.meta_path)
    if nargin >= 3 && ~isempty(mef_path)
        Dat.meta_path = mef_path;
    else
        error('Dat.meta_path not set and mef_path not provided');
    end
end

% If provided mef_path parameter exists, use it (overrides Dat.meta_path)
if nargin >= 3 && ~isempty(mef_path) && exist(mef_path, 'dir')
    Dat.meta_path = mef_path;
    fprintf('Using provided MEF path: %s\n', Dat.meta_path);
end

% If MEF path doesn't exist, try the default path (same logic as prep_basics)
if ~exist(Dat.meta_path, 'dir')
    default_path = fullfile(base_dir, 'Data', patient, 'sub.mefd');
    if exist(default_path, 'dir')
        fprintf('MEF path from prep_basics not found, using default: %s\n', default_path);
        Dat.meta_path = default_path;
    else
        fprintf('Warning: MEF path not found at %s\n', Dat.meta_path);
        fprintf('  Also checked default path: %s\n', default_path);
        fprintf('  You can provide the correct path: edit_ufos_interactive(''%s'', [], ''/path/to/mef'', [])\n', patient);
    end
end

if ~isfield(Dat, 'password') || isempty(Dat.password)
    if nargin >= 4 && ~isempty(password)
        Dat.password = password;
    else
        Dat.password = 'bemena';  % Default
        fprintf('Using default password\n');
    end
end

% Setup output directory
output_dir = fullfile(base_dir, 'Data', patient, 'UFO_edited');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

% Check for existing edited UFO file
output_file = fullfile(output_dir, 'UFOdat_edited.mat');
load_existing_file = false;
if exist(output_file, 'file')
    fprintf('\n=== Found existing edited UFO file ===\n');
    try
        Existing = load(output_file);
        if isfield(Existing, 'UFOdat_edited')
            % Count UFOs in existing file
            existing_count = 0;
            for ii = 1:numel(Existing.UFOdat_edited)
                if ~isempty(Existing.UFOdat_edited{ii}) && size(Existing.UFOdat_edited{ii}, 1) > 0
                    existing_count = existing_count + size(Existing.UFOdat_edited{ii}, 1);
                end
            end
            fprintf('Existing file contains: %d channels, %d total UFOs\n', ...
                numel(Existing.UFOdat_edited), existing_count);
            fprintf('\nOptions:\n');
            fprintf('  [f] Start fresh from original UFOs (will overwrite existing file)\n');
            fprintf('  [l] Load existing edited file and continue from there\n');
            resp = input('Choose option [f/l]: ', 's');
            if strcmpi(strtrim(resp), 'l')
                load_existing_file = true;
                fprintf('Loading existing edited file...\n');
            else
                fprintf('Starting fresh from original UFOs (existing file will be overwritten when you save)\n');
                load_existing_file = false;
            end
        end
    catch ME
        fprintf('Warning: Could not load existing file: %s\n', ME.message);
        fprintf('Starting fresh\n');
        load_existing_file = false;
    end
end

% Check for existing checkpoint
checkpoint_file = fullfile(output_dir, 'UFO_editing_checkpoint.mat');
resume_editing = false;
current_idx = 1;
history = struct('micro_idx', {}, 'event_idx', {}, 'action', {}, 'old_value', {});
prev_indices = [];  % stack of previously visited UFO indices for 'b' navigation

if exist(checkpoint_file, 'file')
    fprintf('\n=== Found existing checkpoint ===\n');
    try
        Checkpoint = load(checkpoint_file);
        if isfield(Checkpoint, 'UFOdat_edited') && isfield(Checkpoint, 'current_idx')
            fprintf('Checkpoint found: %d UFOs processed\n', Checkpoint.current_idx - 1);
            fprintf('Checkpoint has %d channels in UFOdat_edited\n', numel(Checkpoint.UFOdat_edited));
            resp = input('Resume from checkpoint? [y/n]: ', 's');
            if strcmpi(strtrim(resp), 'y')
                UFOdat_edited = Checkpoint.UFOdat_edited;
                current_idx = Checkpoint.current_idx;
                if isfield(Checkpoint, 'history')
                    history = Checkpoint.history;
                end
                if isfield(Checkpoint, 'prev_indices')
                    prev_indices = Checkpoint.prev_indices;
                end
                if isfield(Checkpoint, 'ufo_list')
                    ufo_list = Checkpoint.ufo_list;
                else
                    % Rebuild ufo_list from current UFOdat_edited state
                    ufo_list = [];
                end
                resume_editing = true;
                fprintf('Resuming from UFO %d\n', current_idx);
            else
                fprintf('Starting fresh (checkpoint will be overwritten)\n');
                % Clear any loaded checkpoint data
                clear Checkpoint;
            end
        end
    catch ME
        fprintf('Warning: Could not load checkpoint: %s\n', ME.message);
        fprintf('Starting fresh\n');
    end
end

% Debug: Check if resuming
if resume_editing
    fprintf('[DEBUG] Resuming mode: current_idx=%d\n', current_idx);
    if exist('ufo_list', 'var')
        fprintf('[DEBUG] ufo_list exists, size: %s\n', mat2str(size(ufo_list)));
    else
        fprintf('[DEBUG] ufo_list does not exist yet\n');
    end
end

% Create a copy of UFOdat for editing (if not resuming)
if ~resume_editing
    % Always load existing edited file when it exists; NEVER overwrite it here.
    if exist(output_file,'file') == 2
        Existing = load(output_file);
        if isfield(Existing, 'UFOdat_edited')
            UFOdat_edited = ensure_4_columns(Existing.UFOdat_edited);
            fprintf('Loaded existing edited file with %d channels (types preserved)\n', numel(UFOdat_edited));
        else
            error('edit_ufos_interactive: Existing edited file %s does not contain UFOdat_edited. Refusing to overwrite automatically.', output_file);
        end
    else
        % No existing edited file: initialise from original UFOdat once
        fprintf('No existing UFOdat_edited.mat found. Initialising from original UFOdat.\n');
    UFOdat_edited = cell(size(UFOdat));
for ii = 1:numel(UFOdat)
    if ~isempty(UFOdat{ii}) && size(UFOdat{ii}, 1) > 0
                % Add 4th column (type, default 0 = unclassified)
                if isinteger(UFOdat{ii})
                    type_col = zeros(size(UFOdat{ii}, 1), 1, 'like', UFOdat{ii});
                else
                    type_col = zeros(size(UFOdat{ii}, 1), 1, 'int8');
                end
                UFOdat_edited{ii} = [UFOdat{ii}, type_col];
            else
                UFOdat_edited{ii} = zeros(0,4,'int64');
            end
        end
    end
else
    % When resuming, check if UFOdat_edited matches Dat.allmi structure.
    % In some configurations, UFOdat_edited may have more channels than
    % Dat.allmi. That is safe: later, when building ufo_list, channels with
    % ii > numel(Dat.allmi) are skipped with a warning. We should NOT throw
    % away existing manual edits in this case.
    if numel(UFOdat_edited) ~= numel(Dat.allmi)
        fprintf('\n=== WARNING: Structure mismatch detected ===\n');
        fprintf('UFOdat_edited has %d channels but Dat.allmi has %d channels.\n', ...
            numel(UFOdat_edited), numel(Dat.allmi));
        fprintf('Proceeding anyway: channels beyond Dat.allmi will be ignored when building ufo_list.\n');
        fprintf('Existing classifications and trims are preserved.\n\n');
        % Do NOT modify UFOdat_edited, current_idx, or resume_editing here.
        % Resuming from checkpoint remains valid; extra channels are just
        % ignored later.
    end
end

% Apply minimum frequency filter to UFOdat_edited (if enabled)
% Remove UFOs below min_freq_hz from UFOdat_edited (they are automatically discarded)
total_before_filter = sum(cellfun(@(u) size(u,1), UFOdat_edited));
if min_freq_hz > 0
    filtered_count = 0;
    for ii = 1:numel(UFOdat_edited)
        if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 1) > 0
            % Filter out rows where frequency (column 3) < min_freq_hz
            keep_mask = UFOdat_edited{ii}(:, 3) >= min_freq_hz;
            filtered_count = filtered_count + sum(~keep_mask);
            UFOdat_edited{ii} = UFOdat_edited{ii}(keep_mask, :);
        end
    end
    total_after_filter = sum(cellfun(@(u) size(u,1), UFOdat_edited));
    if filtered_count > 0
        fprintf('\n=== Frequency Filter Applied ===\n');
        fprintf('Automatically discarded %d UFOs below minimum frequency (%.1f Hz)\n', filtered_count, min_freq_hz);
        fprintf('UFOs remaining: %d (down from %d)\n\n', total_after_filter, total_before_filter);
    else
        fprintf('\n=== Frequency Filter Applied ===\n');
        fprintf('No UFOs were below the minimum frequency threshold (%.1f Hz)\n', min_freq_hz);
        fprintf('All %d UFOs are available for review\n\n', total_before_filter);
    end
else
    fprintf('\n=== Frequency Filter Disabled ===\n');
    fprintf('All %d UFOs are available for review (no frequency filtering applied)\n\n', total_before_filter);
end

% Ensure all UFOdat_edited entries have 4 columns before Type 2 filtering
UFOdat_edited = ensure_4_columns(UFOdat_edited);

% Apply Type 2 exclusion filter (if enabled)
% Remove Type 2 UFOs (artifacts) from UFOdat_edited
% Type encoding: 1=Type1<2kHz, 2=Type1>2kHz, 3=Type2<2kHz, 4=Type2>2kHz
if exclude_type2
    total_before_type2_filter = sum(cellfun(@(u) size(u,1), UFOdat_edited));
    filtered_type2_count = 0;
    for ii = 1:numel(UFOdat_edited)
        if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 1) > 0 && size(UFOdat_edited{ii}, 2) >= 4
            % Filter out rows where type (column 4) == 3 or 4 (Type 2)
            keep_mask = (UFOdat_edited{ii}(:, 4) ~= 3) & (UFOdat_edited{ii}(:, 4) ~= 4);
            filtered_type2_count = filtered_type2_count + sum(~keep_mask);
            UFOdat_edited{ii} = UFOdat_edited{ii}(keep_mask, :);
        end
    end
    total_after_type2_filter = sum(cellfun(@(u) size(u,1), UFOdat_edited));
    if filtered_type2_count > 0
        fprintf('\n=== Type 2 Exclusion Filter Applied ===\n');
        fprintf('Automatically excluded %d Type 2 UFOs (artifacts)\n', filtered_type2_count);
        fprintf('UFOs remaining: %d (down from %d)\n\n', total_after_type2_filter, total_before_type2_filter);
    else
        fprintf('\n=== Type 2 Exclusion Filter Applied ===\n');
        fprintf('No Type 2 UFOs found to exclude\n');
        fprintf('All %d UFOs are available for review\n\n', total_before_type2_filter);
    end
end

% Build list of all UFOs with indices (always rebuild to match current state)
fprintf('[DEBUG] About to rebuild ufo_list...\n');
fprintf('[DEBUG] resume_editing = %d\n', resume_editing);
fprintf('[DEBUG] Rebuilding ufo_list from UFOdat_edited...\n');
fprintf('[DEBUG] UFOdat_edited has %d cells\n', numel(UFOdat_edited));
fprintf('[DEBUG] Dat.allmi has %d elements\n', numel(Dat.allmi));

try
    ufo_list = [];
    for ii = 1:numel(UFOdat_edited)
        if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 1) > 0
            % Validate that this micro index exists in Dat.allmi
            if ii > numel(Dat.allmi)
                fprintf('[DEBUG] Warning: UFOdat_edited{%d} exists but Dat.allmi only has %d elements, skipping...\n', ...
                    ii, numel(Dat.allmi));
                continue;
            end
            if isempty(Dat.allmi{ii})
                fprintf('[DEBUG] Warning: Dat.allmi{%d} is empty, skipping UFOs for this channel...\n', ii);
                continue;
            end
            for kk = 1:size(UFOdat_edited{ii}, 1)
                ufo_list(end+1, :) = [ii, kk];  % [micro_idx, event_idx]
            end
        end
    end
    fprintf('[DEBUG] Rebuilt ufo_list: %d entries\n', size(ufo_list, 1));
catch ME
    fprintf('[DEBUG] Error rebuilding ufo_list: %s\n', ME.message);
    fprintf('[DEBUG] Stack trace:\n');
    for k = 1:numel(ME.stack)
        fprintf('  %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
    rethrow(ME);
end

if isempty(ufo_list)
    fprintf('No UFOs found to review.\n');
    return;
end

total_ufos = size(ufo_list, 1);
fprintf('[DEBUG] total_ufos = %d, current_idx = %d\n', total_ufos, current_idx);

% If resuming, ensure current_idx is valid for the rebuilt list
if resume_editing
    if current_idx > total_ufos
        % Clamp to last valid UFO instead of restarting from the beginning
        fprintf('Warning: Saved index %d exceeds total UFOs (%d). Clamping to last UFO.\n', ...
            current_idx, total_ufos);
        current_idx = total_ufos;
    elseif current_idx < 1
    current_idx = 1;
    else
        % Validate that the current_idx entry is valid
        if current_idx <= size(ufo_list, 1)
            test_micro_idx = ufo_list(current_idx, 1);
            if test_micro_idx > numel(Dat.allmi) || test_micro_idx > numel(UFOdat_edited)
                fprintf('Warning: Checkpoint index %d points to invalid micro_idx %d. Clamping to last UFO.\n', ...
                    current_idx, test_micro_idx);
                current_idx = total_ufos;
            end
        end
    end
end

fprintf('\n=== Interactive UFO Editor ===\n');
fprintf('Total UFOs to review: %d\n', total_ufos);
fprintf('Current index: %d\n', current_idx);
fprintf('Commands: k=keep, d=discard, t=trim, b=back, q=quit (saves checkpoint)\n\n');

% Process each UFO
while current_idx <= total_ufos
    fprintf('Processing UFO %d/%d...\n', current_idx, total_ufos);
    micro_idx = ufo_list(current_idx, 1);
    event_idx = ufo_list(current_idx, 2);
    
    % Validate micro_idx is within bounds
    if micro_idx < 1 || micro_idx > numel(UFOdat_edited)
        fprintf('UFO %d/%d: Invalid micro_idx %d (max: %d), skipping...\n', ...
            current_idx, total_ufos, micro_idx, numel(UFOdat_edited));
        current_idx = current_idx + 1;
        continue;
    end
    
    % Validate micro_idx matches Dat.allmi
    if micro_idx > numel(Dat.allmi)
        fprintf('UFO %d/%d: micro_idx %d exceeds Dat.allmi size (%d), skipping...\n', ...
            current_idx, total_ufos, micro_idx, numel(Dat.allmi));
        current_idx = current_idx + 1;
        continue;
    end
    
    % Check if this UFO still exists (might have been removed)
    if isempty(UFOdat_edited{micro_idx}) || event_idx > size(UFOdat_edited{micro_idx}, 1) || size(UFOdat_edited{micro_idx}, 1) == 0
        % This micro channel has no more UFOs, skip to next
        fprintf('UFO %d/%d: Channel %d has no UFOs or event %d doesn''t exist, skipping...\n', ...
            current_idx, total_ufos, micro_idx, event_idx);
        current_idx = current_idx + 1;
        continue;
    end
    
    U = UFOdat_edited{micro_idx};
    if event_idx > size(U, 1)
        fprintf('UFO %d/%d: Event %d doesn''t exist in channel %d (max: %d), skipping...\n', ...
            current_idx, total_ufos, event_idx, micro_idx, size(U, 1));
        current_idx = current_idx + 1;
        continue;
    end
    
    % Get UFO info
    ufo_start_us = double(U(event_idx, 1));
    ufo_end_us = double(U(event_idx, 2));
    ufo_freq_hz = U(event_idx, 3);
    
    % Validate micro channel exists
    if micro_idx > numel(Dat.allmi) || isempty(Dat.allmi{micro_idx})
        fprintf('UFO %d/%d: Dat.allmi{%d} is empty or invalid, skipping...\n', ...
            current_idx, total_ufos, micro_idx);
        current_idx = current_idx + 1;
        continue;
    end
    
    micro_ch = Dat.allmi{micro_idx};
    if iscell(micro_ch)
        micro_ch_name = micro_ch{1};
    else
        micro_ch_name = micro_ch;
    end
    
    % Validate channel name
    if isempty(micro_ch_name) || ~ischar(micro_ch_name)
        fprintf('UFO %d/%d: Invalid channel name, skipping...\n', current_idx, total_ufos);
        current_idx = current_idx + 1;
        continue;
    end
    
    % Validate time range
    if ~isfinite(ufo_start_us) || ~isfinite(ufo_end_us) || ufo_end_us <= ufo_start_us
        fprintf('UFO %d/%d: Invalid time range (start=%.1f, end=%.1f us), skipping...\n', ...
            current_idx, total_ufos, ufo_start_us, ufo_end_us);
        current_idx = current_idx + 1;
        continue;
    end
    
    % Read and plot the UFO
    window_before_us = 20e3;  % 20 ms before
    window_after_us = 20e3;   % 20 ms after
    
    win_start = int64(ufo_start_us - window_before_us);
    win_end = int64(ufo_end_us + window_after_us);
    
    % Validate window
    if win_end <= win_start
        fprintf('UFO %d/%d: Invalid read window, skipping...\n', current_idx, total_ufos);
            current_idx = current_idx + 1;
            continue;
        end
    
    fprintf('  Reading micro channel %s from MEF...\n', micro_ch_name);
    fprintf('  Time window: %.3f to %.3f seconds (%.1f ms window)\n', ...
        double(win_start)/1e6, double(win_end)/1e6, double(win_end - win_start)/1e3);
    fprintf('  MEF path: %s\n', Dat.meta_path);
    
    % Check if MEF path exists
    mef_available = exist(Dat.meta_path, 'dir');
    if ~mef_available
        fprintf('  WARNING: MEF directory not found: %s\n', Dat.meta_path);
        fprintf('  Plot cannot be displayed, but you can still keep/discard this UFO\n');
    end
    
    % Try to read with a timeout-like approach: use a smaller test read first
    % This helps identify if the channel exists and is readable
    micro = [];
    read_success = false;
    
    % Only try to read if MEF is available
    if mef_available
        % Check for slow read threshold (configurable, default 30 seconds)
        slow_read_threshold = getin('CERD_SLOW_READ_THRESHOLD_SEC', 30);
        skip_slow_reads = getin('CERD_SKIP_SLOW_READS', false);
        
        test_win_size = min(1e6, double(win_end - win_start) / 10);  % 1 second or 10% of window
        test_start = win_start;
        test_end = int64(double(win_start) + test_win_size);
        
        try
            fprintf('  Testing read with small window first...\n');
            fprintf('  (If this hangs, it may be a MATLAB/readMef3 issue for this patient)\n');
            fprintf('  To skip this patient, set: CERD_SKIP_READMEF3_PATIENTS = {''%s''};\n', patient);
            fprintf('  Or use FAST META mode if available for this patient.\n');
            tic;
            [~, micro_test] = readMef3(Dat.meta_path, Dat.password, micro_ch_name, 'time', test_start, test_end);
            test_time = toc;
            fprintf('  Test read completed in %.2f seconds, got %d samples\n', test_time, numel(micro_test));
            
            if test_time > slow_read_threshold
                fprintf('  WARNING: Test read took %.1f seconds (threshold: %.1f s)\n', test_time, slow_read_threshold);
                if skip_slow_reads
                    fprintf('  Skipping slow read (CERD_SKIP_SLOW_READS=true)...\n');
                    read_success = false;  % Don't skip UFO, just don't show plot
                else
                    fprintf('  Continuing anyway (set CERD_SKIP_SLOW_READS=true to skip slow reads)\n');
                end
            end
            
            % If test read succeeded and we're not skipping slow reads, proceed with full read
            if skip_slow_reads
                fprintf('  Skipping full read (CERD_SKIP_SLOW_READS=true)\n');
                read_success = false;
            elseif isempty(micro_test)
                fprintf('  Test read returned empty data - plot unavailable\n');
                read_success = false;
            else
                % Test read succeeded, do the full read
                fprintf('  Reading full window...\n');
                fprintf('  (This may take a while if readMef3 is slow for this patient)\n');
                try
                    tic;
                    [~, micro] = readMef3(Dat.meta_path, Dat.password, micro_ch_name, 'time', win_start, win_end);
                    read_time = toc;
                    fprintf('  Full read completed in %.2f seconds, got %d samples\n', read_time, numel(micro));
                    
                    if read_time > slow_read_threshold
                        fprintf('  WARNING: Full read took %.1f seconds (threshold: %.1f s)\n', read_time, slow_read_threshold);
                        fprintf('  Consider setting CERD_SKIP_SLOW_READS=true to skip slow reads\n');
                    end
                    
                    if ~isempty(micro) && ~any(isnan(micro), 'all')
                        read_success = true;
                    else
                        fprintf('  Full read returned empty or invalid data\n');
                        read_success = false;
                    end
                catch ME2
                    fprintf('  Error during full read: %s\n', ME2.message);
                    read_success = false;
                end
            end
            
    catch ME
            fprintf('  Error reading micro data: %s\n', ME.message);
            fprintf('  Channel: %s\n', micro_ch_name);
            fprintf('  Error type: %s\n', ME.identifier);
            if ~isempty(ME.stack)
                fprintf('  Error location: %s at line %d\n', ME.stack(1).name, ME.stack(1).line);
            end
            fprintf('  Plot unavailable, but you can still keep/discard this UFO\n');
            read_success = false;
        end
        
        % Check if we have valid data
        if read_success && exist('micro', 'var') && ~isempty(micro) && ~any(isnan(micro), 'all')
            % Data read successfully, proceed with plotting
            fprintf('  Data read successfully: %d samples\n', numel(micro));
        else
            read_success = false;
            if exist('micro', 'var') && isempty(micro)
                fprintf('  Full read returned empty data - plot unavailable\n');
            elseif exist('micro', 'var') && any(isnan(micro), 'all')
                fprintf('  Full read contains NaN values - plot unavailable\n');
            else
                fprintf('  Cannot display plot (data unavailable), but you can still keep/discard this UFO\n');
    end
        end
    end
    
    % Only plot if data was successfully read
    if read_success && ~isempty(micro) && ~any(isnan(micro), 'all')
    % Convert to double and estimate sampling rate
    micro = double(micro(:));
    dt_us = double(win_end - win_start) / numel(micro);
    fs = 1e6 / dt_us;
        
        fprintf('  Sampling rate: %.1f Hz, %d samples in window\n', fs, numel(micro));
    
    % Create time array (relative to UFO start, in ms)
        % Time array goes from 0 to window_duration, then we center at UFO start
    t_ms = (0:numel(micro)-1)' / fs * 1000;
        
        % Calculate where UFO start/end fall within the read window
        ufo_start_offset_ms = double(ufo_start_us - win_start) / 1e3;  % Should be ~20ms
        ufo_end_offset_ms = double(ufo_end_us - win_start) / 1e3;
        ufo_duration_ms = (ufo_end_us - ufo_start_us) / 1e3;
        
        fprintf('  UFO start offset in window: %.2f ms (should be ~20 ms)\n', ufo_start_offset_ms);
        fprintf('  UFO end offset in window: %.2f ms\n', ufo_end_offset_ms);
        fprintf('  UFO duration: %.2f ms\n', ufo_duration_ms);
        
        % Center time array at UFO start (so UFO starts at t=0)
        t_ms = t_ms - ufo_start_offset_ms;
        
        % UFO boundaries in centered time (UFO starts at 0)
    ufo_start_rel_ms = 0;
        ufo_end_rel_ms = ufo_duration_ms;
    
        % Create mask for UFO period (should be from 0 to ufo_duration_ms)
    ufo_mask = t_ms >= ufo_start_rel_ms & t_ms <= ufo_end_rel_ms;
        
        fprintf('  Time axis range: %.2f to %.2f ms (centered at UFO start)\n', min(t_ms), max(t_ms));
        fprintf('  UFO period: %.2f to %.2f ms\n', ufo_start_rel_ms, ufo_end_rel_ms);
        fprintf('  Samples in UFO period: %d / %d\n', nnz(ufo_mask), numel(micro));
        
        % Additional diagnostics: signal properties within UFO period
        if any(ufo_mask)
            ufo_signal = micro(ufo_mask);
            ufo_amplitude_range = [min(ufo_signal), max(ufo_signal)];
            ufo_amplitude_span = diff(ufo_amplitude_range);
            ufo_rms = sqrt(mean(ufo_signal.^2));
            ufo_std = std(ufo_signal);
            
            % Compare to baseline (before UFO)
            baseline_mask = t_ms < -5;  % Use signal before -5ms as baseline
            if any(baseline_mask)
                baseline_signal = micro(baseline_mask);
                baseline_std = std(baseline_signal);
                baseline_rms = sqrt(mean(baseline_signal.^2));
                snr_estimate = ufo_std / max(baseline_std, eps);
                fprintf('  UFO signal stats: range=%.3f µV, RMS=%.3f µV, std=%.3f µV\n', ...
                    ufo_amplitude_span, ufo_rms, ufo_std);
                fprintf('  Baseline (pre-UFO) std=%.3f µV, RMS=%.3f µV\n', baseline_std, baseline_rms);
                fprintf('  Estimated SNR (std ratio): %.2f\n', snr_estimate);
            end
            
            % Check if this looks like a typical UFO (high frequency, brief burst)
            % Typical UFOs: 2-8 kHz, 1-10 ms duration, high amplitude
            if ufo_freq_hz < 500 || ufo_freq_hz > 10000
                fprintf('  WARNING: Frequency %.1f Hz is outside typical UFO range (500-10000 Hz)\n', ufo_freq_hz);
            end
            if ufo_duration_ms > 20
                fprintf('  WARNING: Duration %.2f ms is longer than typical UFOs (usually <10 ms)\n', ufo_duration_ms);
            end
        end
        
        % Verify UFO mask is valid
        if ~any(ufo_mask)
            warning('UFO period mask is empty - no samples match UFO time range!');
            warning('This might indicate a timestamp mismatch. Check UFO timestamps.');
        end
    
    % Display plot
    fig = figure(1);
    clf(fig);
    set(fig, 'Position', [100 100 1200 400], 'Name', sprintf('UFO %d/%d', current_idx, total_ufos), ...
        'Visible', 'on');
    
    % Plot entire signal in black first
    plot(t_ms, micro, 'k-', 'LineWidth', 1);
    hold on;
    
        % Overlay UFO period in red (thicker line)
        if any(ufo_mask)
            plot(t_ms(ufo_mask), micro(ufo_mask), 'r-', 'LineWidth', 3);
        else
            warning('No samples in UFO period - check timestamps!');
        end
    
    % Add vertical lines at UFO boundaries
        xline(ufo_start_rel_ms, 'r--', 'LineWidth', 2, 'Label', 'UFO start');
        xline(ufo_end_rel_ms, 'r--', 'LineWidth', 2, 'Label', 'UFO end');
        
        % Add shaded region for UFO period
        if any(ufo_mask)
            y_range = [min(micro), max(micro)];
            fill([ufo_start_rel_ms, ufo_end_rel_ms, ufo_end_rel_ms, ufo_start_rel_ms], ...
                [y_range(1), y_range(1), y_range(2), y_range(2)], ...
                'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
    
    % Formatting
    xlabel('Time (ms, relative to UFO start)');
        ylabel('Amplitude (μV)');
    title(sprintf('UFO %d/%d: Micro %d (%s) | Event %d | Freq=%.1f Hz | Duration=%.2f ms', ...
            current_idx, total_ufos, micro_idx, micro_ch_name, event_idx, ufo_freq_hz, ufo_duration_ms));
    grid on;
    
        % Set reasonable x-axis limits (show ~50ms around UFO)
        xlim([-30, max(30, ufo_duration_ms + 10)]);
        
        % Add text annotation at UFO center
        if any(ufo_mask)
            ufo_center_ms = ufo_duration_ms / 2;
            text(ufo_center_ms, max(micro)*0.9, 'UFO', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
        end
        
        % Add info text
        info_str = sprintf('Window: %.1f ms before, %.1f ms after', window_before_us/1e3, window_after_us/1e3);
        text(0.02, 0.98, info_str, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', [0.5 0.5 0.5]);
    
    % Force display update and bring figure to front
    drawnow;
    figure(fig);  % Bring to front
    else
        % Calculate duration for display even without plot
        ufo_duration_ms = (ufo_end_us - ufo_start_us) / 1e3;
        fprintf('  UFO duration: %.2f ms\n', ufo_duration_ms);
        
        % Create a placeholder figure with message
        fprintf('  Creating placeholder figure...\n');
        fig = figure(1);
        clf(fig);
        set(fig, 'Position', [100 100 1200 400], 'Name', sprintf('UFO %d/%d', current_idx, total_ufos), ...
            'Visible', 'on');
        
        % Display message explaining why plot is unavailable
        axes('Position', [0.1 0.1 0.8 0.8]);
        axis off;
        
        msg_text = {
            sprintf('UFO %d/%d: Micro %d (%s), Event %d', current_idx, total_ufos, micro_idx, micro_ch_name, event_idx);
            '';
            'Plot Unavailable';
            '';
            sprintf('MEF directory not found: %s', Dat.meta_path);
            '';
            sprintf('UFO Properties:');
            sprintf('  Frequency: %.1f Hz', ufo_freq_hz);
            sprintf('  Duration: %.2f ms', ufo_duration_ms);
            sprintf('  Start time: %.3f seconds', double(ufo_start_us)/1e6);
            sprintf('  End time: %.3f seconds', double(ufo_end_us)/1e6);
            '';
            'You can still keep or discard this UFO';
            'Trim is not available without data access';
        };
        
        text(0.5, 0.5, msg_text, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 12, 'FontName', 'Courier');
        
        title(sprintf('UFO %d/%d: Data Unavailable', current_idx, total_ufos), 'FontSize', 14, 'FontWeight', 'bold');
        
        % Force display update and bring figure to front
        drawnow;
        figure(fig);  % Bring to front
    end
    
    % Prompt for action
    fprintf('\nUFO %d/%d: Micro %d (%s), Event %d\n', current_idx, total_ufos, micro_idx, micro_ch_name, event_idx);
    fprintf('  Frequency: %.1f Hz\n', ufo_freq_hz);
    % Show current type if classified
    if size(U, 2) >= 4 && event_idx <= size(U, 1)
        ufo_type = U(event_idx, 4);
        type_str = '';
        if ufo_type == 1
            type_str = 'Type 1, < 2kHz (genuine, low freq)';
        elseif ufo_type == 2
            type_str = 'Type 1, > 2kHz (genuine, high freq)';
        elseif ufo_type == 3
            type_str = 'Type 2, < 2kHz (artifact, low freq)';
        elseif ufo_type == 4
            type_str = 'Type 2, > 2kHz (artifact, high freq)';
        else
            type_str = 'Unclassified';
        end
        fprintf('  Current classification: %s\n', type_str);
    else
        fprintf('  Current classification: Unclassified\n');
    end
    % Loop to allow multiple actions on the same UFO (classify, then trim, etc.)
    ufo_done = false;
    while ~ufo_done && current_idx <= total_ufos
        % Reload U from UFOdat_edited in case it was modified (e.g., by trim)
        U = UFOdat_edited{micro_idx};
        if isempty(U) || event_idx > size(U, 1)
            % UFO was removed, break out
            break;
        end
        
        % Re-display current state
        fprintf('\nUFO %d/%d: Micro %d (%s), Event %d\n', current_idx, total_ufos, micro_idx, micro_ch_name, event_idx);
        fprintf('  Frequency: %.1f Hz\n', U(event_idx, 3));
        % Show current type if classified
        if size(U, 2) >= 4 && event_idx <= size(U, 1)
            ufo_type = U(event_idx, 4);
            type_str = '';
            if ufo_type == 1
                type_str = 'Type 1, < 2kHz (genuine, low freq)';
            elseif ufo_type == 2
                type_str = 'Type 1, > 2kHz (genuine, high freq)';
            elseif ufo_type == 3
                type_str = 'Type 2, < 2kHz (artifact, low freq)';
            elseif ufo_type == 4
                type_str = 'Type 2, > 2kHz (artifact, high freq)';
            else
                type_str = 'Unclassified';
            end
            fprintf('  Current classification: %s\n', type_str);
        else
            fprintf('  Current classification: Unclassified\n');
        end
        
        % Determine frequency threshold for automatic low/high classification
        freq_threshold_hz = 2000;
        is_low_freq = ufo_freq_hz < freq_threshold_hz;
        freq_label = ternary(is_low_freq, '< 2kHz', '> 2kHz');
        
        if read_success && ~isempty(micro)
            action = input(sprintf('Action [1=Type1, 2=Type2, t=trim, k=keep&next, d=discard, n=next, b=back, g=goto, q=quit] (freq=%.1f Hz %s): ', ufo_freq_hz, freq_label), 's');
        else
            fprintf('  (Plot unavailable - trim not available)\n');
            action = input(sprintf('Action [1=Type1, 2=Type2, k=keep&next, d=discard, n=next, b=back, g=goto, q=quit] (freq=%.1f Hz %s): ', ufo_freq_hz, freq_label), 's');
        end
    action = lower(strtrim(action));
    
    if isempty(action)
            action = 'n';  % Default to next
    end
    
    switch action
        case {'k', 'keep'}
            % Keep as is (ensure 4 columns), move to next
            fprintf('  → Keeping UFO as is, moving to next\n');
            % Ensure 4 columns
            if size(U, 2) == 3
                if isinteger(U)
                    type_col = zeros(size(U, 1), 1, 'like', U);
                else
                    type_col = zeros(size(U, 1), 1, 'int8');
                end
                U = [U, type_col];
            end
            UFOdat_edited{micro_idx} = U;
            % Record current index for 'b' back navigation
            prev_indices(end+1) = current_idx;
            ufo_done = true;  % Move to next UFO
            current_idx = current_idx + 1;
            % Save checkpoint periodically (every 10 UFOs)
            if mod(current_idx - 1, 10) == 0
                try
                    save(checkpoint_file, 'UFOdat_edited', 'current_idx', 'history', 'ufo_list', '-v7.3');
                catch ME
                    warning('Could not save checkpoint: %s', ME.message);
                end
            end
            
        case {'n', 'next'}
            % Move to next UFO without any action
            fprintf('  → Moving to next UFO\n');
            % Record current index for 'b' back navigation
            prev_indices(end+1) = current_idx;
            ufo_done = true;  % Move to next UFO
            current_idx = current_idx + 1;
            
        case {'1', 'type1'}
            % Classify as Type 1 (genuine) - automatically determine low/high based on frequency
            % Type encoding: 1=Type1<2kHz, 2=Type1>2kHz
            if is_low_freq
                type_value = 1;  % Type 1, < 2kHz
                type_desc = 'Type 1, < 2kHz (genuine, low frequency)';
            else
                type_value = 2;  % Type 1, > 2kHz
                type_desc = 'Type 1, > 2kHz (genuine, high frequency)';
            end
            fprintf('  → Classified as %s\n', type_desc);
            fprintf('  → (Press t=trim to adjust boundaries, k=keep&next to move on, or continue with more actions)\n');
            % Ensure 4 columns
            if size(U, 2) == 3
                if isinteger(U)
                    type_col = zeros(size(U, 1), 1, 'like', U);
                else
                    type_col = zeros(size(U, 1), 1, 'int8');
                end
                U = [U, type_col];
            end
            if size(U, 2) >= 4 && event_idx <= size(U, 1)
                U(event_idx, 4) = type_value;
                UFOdat_edited{micro_idx} = U;
            else
                warning('Could not set type: invalid index or structure');
            end
            % Stay on same UFO to allow more actions (don't set ufo_done)
            
        case {'2', 'type2'}
            % Classify as Type 2 (artifact) - automatically determine low/high based on frequency
            % Type encoding: 3=Type2<2kHz, 4=Type2>2kHz
            if is_low_freq
                type_value = 3;  % Type 2, < 2kHz
                type_desc = 'Type 2, < 2kHz (artifact, low frequency)';
            else
                type_value = 4;  % Type 2, > 2kHz
                type_desc = 'Type 2, > 2kHz (artifact, high frequency)';
            end
            fprintf('  → Classified as %s\n', type_desc);
            if exclude_type2
                fprintf('  → WARNING: Type 2 exclusion is enabled. This UFO will be filtered out in analysis.\n');
            end
            fprintf('  → (Press t=trim to adjust boundaries, k=keep&next to move on, or continue with more actions)\n');
            % Ensure 4 columns
            if size(U, 2) == 3
                if isinteger(U)
                    type_col = zeros(size(U, 1), 1, 'like', U);
                else
                    type_col = zeros(size(U, 1), 1, 'int8');
                end
                U = [U, type_col];
            end
            if size(U, 2) >= 4 && event_idx <= size(U, 1)
                U(event_idx, 4) = type_value;
                UFOdat_edited{micro_idx} = U;
            else
                warning('Could not set type: invalid index or structure');
            end
            % Stay on same UFO to allow more actions (don't set ufo_done)
            
        case {'d', 'discard'}
            % Remove this UFO
            fprintf('  → Discarding UFO\n');
            
            % Save to history for undo
            history(end+1).micro_idx = micro_idx;
            history(end).event_idx = event_idx;
            history(end).action = 'discard';
            history(end).old_value = U(event_idx, :);
            
            % Remove from UFOdat_edited
            if size(U, 1) == 1
                % Last UFO in this channel, make empty
                UFOdat_edited{micro_idx} = zeros(0, 4, 'like', U);
            else
                U(event_idx, :) = [];
                UFOdat_edited{micro_idx} = U;
            end
            
            % Update ufo_list - remove this entry and adjust indices
            ufo_list(current_idx, :) = [];
            total_ufos = total_ufos - 1;
            
            % Adjust remaining event indices for this micro channel
            mask = ufo_list(:, 1) == micro_idx & ufo_list(:, 2) > event_idx;
            ufo_list(mask, 2) = ufo_list(mask, 2) - 1;
            
            % Save checkpoint after discard
            try
                save(checkpoint_file, 'UFOdat_edited', 'current_idx', 'history', 'ufo_list', 'prev_indices', '-v7.3');
            catch ME
                warning('Could not save checkpoint: %s', ME.message);
            end
            
            % Move to next UFO (discard removes this one)
            ufo_done = true;
            % Don't increment current_idx since we removed an entry
            
        case {'t', 'trim'}
            % Trim UFO - user selects new boundaries
            if ~read_success || isempty(micro)
                fprintf('  → Cannot trim: plot unavailable (MEF data not accessible)\n');
                fprintf('  → Use k=keep or d=discard instead\n');
                continue;
            end
            fprintf('  → Trimming UFO\n');
            fprintf('  Click two points on the plot to select new start and end times\n');
            fprintf('  (First click = new start, Second click = new end)\n');
            
            try
                [x_click, ~] = ginput(2);
                if numel(x_click) < 2
                    fprintf('  → Need 2 clicks. Cancelling trim.\n');
                    % Stay on same UFO (don't set ufo_done)
                    continue;
                end
                
                % Sort clicks (first = start, second = end)
                x_click = sort(x_click);
                new_start_rel_ms = x_click(1);
                new_end_rel_ms = x_click(2);
                
                % Convert back to absolute microseconds (preserve original data type)
                orig_type = class(U);
                new_start_us = ufo_start_us + new_start_rel_ms * 1e3;
                new_end_us = ufo_start_us + new_end_rel_ms * 1e3;
                
                % Cast to original type
                if strcmp(orig_type, 'int64')
                    new_start_us = int64(new_start_us);
                    new_end_us = int64(new_end_us);
                elseif strcmp(orig_type, 'double')
                    new_start_us = double(new_start_us);
                    new_end_us = double(new_end_us);
                end
                
                % Validate
                if new_end_us <= new_start_us
                    fprintf('  → Invalid selection (end <= start). Cancelling trim.\n');
                    % Stay on same UFO (don't set ufo_done)
                    continue;
                end
                
                if new_start_us < win_start || new_end_us > win_end
                    fprintf('  → Selection outside data window. Cancelling trim.\n');
                    % Stay on same UFO (don't set ufo_done)
                    continue;
                end
                
                % Save to history
                history(end+1).micro_idx = micro_idx;
                history(end).event_idx = event_idx;
                history(end).action = 'trim';
                history(end).old_value = U(event_idx, :);
                
                % Update UFO (preserve frequency in column 3)
                U(event_idx, 1) = new_start_us;
                U(event_idx, 2) = new_end_us;
                % Column 3 (frequency) remains unchanged
                UFOdat_edited{micro_idx} = U;
                
                fprintf('  → Trimmed: %.2f ms -> %.2f ms (duration: %.2f ms)\n', ...
                    new_start_rel_ms, new_end_rel_ms, (new_end_us - new_start_us)/1e3);
                fprintf('  → (Press k=keep&next or n=next to move on, or continue with more actions)\n');
                
                % Save checkpoint after trim
                try
                    save(checkpoint_file, 'UFOdat_edited', 'current_idx', 'history', 'ufo_list', '-v7.3');
                catch ME
                    warning('Could not save checkpoint: %s', ME.message);
                end
                
                % Re-read data for updated plot (will happen on next iteration)
                % Stay on same UFO to allow more actions (don't set ufo_done)
                
            catch ME
                fprintf('  → Error during trim: %s. Cancelling trim.\n', ME.message);
                % Stay on same UFO (don't set ufo_done)
            end
            
        case {'b', 'back'}
            % Back: either undo last structural action (discard/trim),
            % or navigate to previously viewed UFO using prev_indices.
            if ~isempty(history)
                % Undo last discard/trim from history
            last = history(end);
            history(end) = [];  % Remove from history
            
            fprintf('  → Undoing last action (%s)\n', last.action);
            
            switch last.action
                case 'discard'
                    % Restore the discarded UFO
                    U_restore = UFOdat_edited{last.micro_idx};
                    if isempty(U_restore)
                        U_restore = last.old_value;
                    else
                        % Insert at original position (sorted by time)
                        U_restore = [U_restore; last.old_value];
                        [~, ord] = sort(U_restore(:, 1));
                        U_restore = U_restore(ord, :);
                    end
                    UFOdat_edited{last.micro_idx} = U_restore;
                    
                    % Rebuild ufo_list
                    ufo_list = [];
                    for ii = 1:numel(UFOdat_edited)
                        if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 1) > 0
                            for kk = 1:size(UFOdat_edited{ii}, 1)
                                ufo_list(end+1, :) = [ii, kk];
                            end
                        end
                    end
                    total_ufos = size(ufo_list, 1);
                    
                    % Find position of restored UFO
                    mask = ufo_list(:, 1) == last.micro_idx & ufo_list(:, 2) == last.event_idx;
                    if any(mask)
                        current_idx = find(mask, 1);
                    end
                    
                case 'trim'
                    % Restore original boundaries
                    U_restore = UFOdat_edited{last.micro_idx};
                    U_restore(last.event_idx, :) = last.old_value;
                    UFOdat_edited{last.micro_idx} = U_restore;
                    
                    % Stay at current UFO to review again
                    % (current_idx already points to this UFO)
                end
            else
                % No undoable history; use navigation stack
                if isempty(prev_indices)
                    fprintf('  → No previous UFO to go back to\n');
                else
                    prev_idx = prev_indices(end);
                    prev_indices(end) = [];
                    % Clamp to valid range
                    prev_idx = max(1, min(prev_idx, total_ufos));
                    fprintf('  → Going back to previous UFO (index %d of %d)\n', prev_idx, total_ufos);
                    current_idx = prev_idx;
                end
            end
            % After handling 'back', force re-draw of the (possibly new) current_idx
            % by exiting the inner action loop and letting the outer loop
            % reload micro_idx/event_idx and replot.
            ufo_done = true;
            
        case {'q', 'quit'}
            % Save and exit
            fprintf('\n=== Saving edited UFOs ===\n');
            % Save final checkpoint before exiting
            try
                save(checkpoint_file, 'UFOdat_edited', 'current_idx', 'history', 'ufo_list', '-v7.3');
            catch ME
                warning('Could not save checkpoint: %s', ME.message);
            end
            ufo_done = true;  % Exit inner loop
            % Mark that we're quitting early (not completing all UFOs)
            user_quit_early = true;  % Set flag to distinguish from natural completion
            current_idx = total_ufos + 1;  % Exit outer loop
            
        case {'g', 'goto'}
            % Jump directly to a specific [micro_idx, event_idx] in the current ufo_list
            fprintf('  → Go to specific UFO\n');
            mi = str2double(input('    Enter micro index: ', 's'));
            ei = str2double(input('    Enter event index within that micro: ', 's'));
            if isnan(mi) || isnan(ei)
                fprintf('  → Invalid input, staying on current UFO\n');
            else
                tgt = find(ufo_list(:,1) == mi & ufo_list(:,2) == ei, 1);
                if isempty(tgt)
                    fprintf('  → No matching [micro=%d, event=%d] in current list (may have been discarded or filtered)\n', mi, ei);
                else
                    current_idx = tgt;
                    fprintf('  → Jumping to UFO %d/%d (micro %d, event %d)\n', current_idx, total_ufos, mi, ei);
                    ufo_done = true;  % break inner loop and redraw at new index
                end
            end
            
        otherwise
            fprintf('  → Unknown action "%s". Use 1=Type1, 2=Type2, t=trim, k=keep&next, d=discard, n=next, b=back, g=goto, q=quit\n', action);
            fprintf('  → You can do multiple actions: e.g., "1" to classify, then "t" to trim, then "k" to keep&next\n');
            % Stay on same UFO (don't set ufo_done)
        end  % End of switch action
        
        % Save checkpoint periodically (every 10 actions, not just every 10 UFOs)
        if mod(current_idx, 10) == 0
            try
                save(checkpoint_file, 'UFOdat_edited', 'current_idx', 'history', 'ufo_list', '-v7.3');
            catch ME
                warning('Could not save checkpoint: %s', ME.message);
    end
end
        
    end  % End of while ~ufo_done (allows multiple actions on same UFO)
end  % End of while current_idx <= total_ufos (main loop)

% Save edited UFOdat with validation
% (output_file already defined earlier)

% Validate structure before saving
fprintf('\n=== Validating edited UFO structure ===\n');
if ~iscell(UFOdat_edited)
    error('UFOdat_edited must be a cell array');
end

% Diagnostic: Check type distribution before save
type_counts_before = zeros(1, 5);
has_types_before = false;
for ii = 1:numel(UFOdat_edited)
    if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 2) >= 4
        has_types_before = true;
        types = UFOdat_edited{ii}(:, 4);
        for t = 0:4
            type_counts_before(t+1) = type_counts_before(t+1) + sum(types == t);
        end
    end
end
if has_types_before
    fprintf('[save] Type distribution before save: Unclassified=%d, Type1_low=%d, Type1_high=%d, Type2_low=%d, Type2_high=%d\n', ...
        type_counts_before(1), type_counts_before(2), type_counts_before(3), type_counts_before(4), type_counts_before(5));
end

total_ufos_check = 0;
for ii = 1:numel(UFOdat_edited)
    if ~isempty(UFOdat_edited{ii})
        n_cols = size(UFOdat_edited{ii}, 2);
        if n_cols ~= 3 && n_cols ~= 4
            error('UFOdat_edited{%d} must have 3 or 4 columns [start_us, end_us, freq_hz, type]', ii);
        end
        % Ensure 4 columns for consistency (preserve existing types if present)
        if n_cols == 3
            if isinteger(UFOdat_edited{ii})
                type_col = zeros(size(UFOdat_edited{ii}, 1), 1, 'like', UFOdat_edited{ii});
            else
                type_col = zeros(size(UFOdat_edited{ii}, 1), 1, 'int8');
            end
            UFOdat_edited{ii} = [UFOdat_edited{ii}, type_col];
            fprintf('[save] Channel %d: Added 4th column (all unclassified)\n', ii);
        elseif n_cols >= 4
            % Already has 4 columns - preserve existing types
            fprintf('[save] Channel %d: Already has %d columns, preserving types\n', ii, n_cols);
        end
        % Check for valid timestamps (end > start)
        for kk = 1:size(UFOdat_edited{ii}, 1)
            if UFOdat_edited{ii}(kk, 2) <= UFOdat_edited{ii}(kk, 1)
                warning('UFOdat_edited{%d}(%d,:) has invalid timestamps (end <= start)', ii, kk);
            end
        end
        total_ufos_check = total_ufos_check + size(UFOdat_edited{ii}, 1);
    end
end

fprintf('Validation passed: %d micro channels, %d total UFOs\n', numel(UFOdat_edited), total_ufos_check);

% Verify types are preserved after validation
type_counts_after = zeros(1, 5);
has_types_after = false;
for ii = 1:numel(UFOdat_edited)
    if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 2) >= 4
        has_types_after = true;
        types = UFOdat_edited{ii}(:, 4);
        for t = 0:4
            type_counts_after(t+1) = type_counts_after(t+1) + sum(types == t);
        end
    end
end
if has_types_after
    fprintf('[save] Type distribution after validation: Unclassified=%d, Type1_low=%d, Type1_high=%d, Type2_low=%d, Type2_high=%d\n', ...
        type_counts_after(1), type_counts_after(2), type_counts_after(3), type_counts_after(4), type_counts_after(5));
    if has_types_before && sum(type_counts_before(2:5)) > 0 && sum(type_counts_after(2:5)) == 0
        warning('[save] WARNING: Types were lost during validation! This should not happen.');
    end
end

% Final safety check: Ensure ALL channels have 4 columns before saving
% This is critical - we must never save with only 3 columns
fprintf('[save] Final safety check: Ensuring all channels have 4 columns...\n');
for ii = 1:numel(UFOdat_edited)
    if ~isempty(UFOdat_edited{ii})
        n_cols = size(UFOdat_edited{ii}, 2);
        if n_cols == 3
            warning('[save] Channel %d still has 3 columns! Adding 4th column now.', ii);
            if isinteger(UFOdat_edited{ii})
                type_col = zeros(size(UFOdat_edited{ii}, 1), 1, 'like', UFOdat_edited{ii});
            else
                type_col = zeros(size(UFOdat_edited{ii}, 1), 1, 'int8');
            end
            UFOdat_edited{ii} = [UFOdat_edited{ii}, type_col];
        elseif n_cols < 3
            error('[save] Channel %d has only %d columns - this should not happen!', ii, n_cols);
        end
    end
end

% Verify final structure one more time
final_cols = cellfun(@(u) size(u,2), UFOdat_edited(~cellfun(@isempty, UFOdat_edited)));
if any(final_cols == 3)
    error('[save] CRITICAL: Some channels still have 3 columns after validation! This is a bug.');
elseif any(final_cols < 4)
    error('[save] CRITICAL: Some channels have fewer than 4 columns!');
else
    fprintf('[save] ✓ All channels verified to have 4 columns\n');
end

% Save with metadata
save(output_file, 'UFOdat_edited', 'Dat', '-v7.3');
fprintf('Saved edited UFOs to: %s\n', output_file);
fprintf('  File contains: UFOdat_edited (cell array), Dat (struct with metadata)\n');
fprintf('  Use load_edited_ufos(''%s'') to load into pipeline\n', patient);

% Final verification: Load back and check
try
    Verify = load(output_file, 'UFOdat_edited');
    verify_cols = cellfun(@(u) size(u,2), Verify.UFOdat_edited(~cellfun(@isempty, Verify.UFOdat_edited)));
    if any(verify_cols == 3)
        warning('[save] WARNING: File was saved but verification shows some channels have 3 columns!');
        warning('[save] This suggests a save/load issue. Please check the file manually.');
    else
        fprintf('[save] ✓ Verification: File saved correctly with 4 columns\n');
    end
catch ME
    warning('[save] Could not verify saved file: %s', ME.message);
end

% Check if we completed all UFOs
% Initialize user_quit_early flag if not set (for natural completion)
if ~exist('user_quit_early', 'var')
    user_quit_early = false;
end

if current_idx > total_ufos && ~user_quit_early
    % Completed all UFOs naturally (didn't press 'q')
    if exist(checkpoint_file, 'file')
        resp = input('\nAll UFOs reviewed. Delete checkpoint? [y/n]: ', 's');
        if strcmpi(strtrim(resp), 'y')
            delete(checkpoint_file);
            fprintf('Checkpoint deleted\n');
        else
            fprintf('Checkpoint kept for safety\n');
        end
    end
else
    % Quit early or didn't complete all UFOs
    fprintf('\nProgress saved. Resume anytime by running this script again.\n');
    fprintf('Checkpoint: %s\n', checkpoint_file);
    if exist('user_quit_early', 'var') && user_quit_early
        fprintf('You quit early - %d/%d UFOs reviewed\n', current_idx - 1, total_ufos);
    else
        fprintf('Progress: %d/%d UFOs reviewed\n', current_idx - 1, total_ufos);
    end
end

% Print summary
% Count from the actual structures being used (not from saved file)
original_count = 0;
for ii = 1:numel(UFOdat)
    if ~isempty(UFOdat{ii}) && size(UFOdat{ii}, 1) > 0
        original_count = original_count + size(UFOdat{ii}, 1);
    end
end

edited_count = 0;
type1_low_count = 0;   % Type 1, < 2kHz (type = 1)
type1_high_count = 0;  % Type 1, > 2kHz (type = 2)
type2_low_count = 0;   % Type 2, < 2kHz (type = 3)
type2_high_count = 0;  % Type 2, > 2kHz (type = 4)
unclassified_count = 0;
for ii = 1:numel(UFOdat_edited)
    if ~isempty(UFOdat_edited{ii}) && size(UFOdat_edited{ii}, 1) > 0
        edited_count = edited_count + size(UFOdat_edited{ii}, 1);
        % Count types if 4th column exists
        if size(UFOdat_edited{ii}, 2) >= 4
            type1_low_count = type1_low_count + sum(UFOdat_edited{ii}(:, 4) == 1);
            type1_high_count = type1_high_count + sum(UFOdat_edited{ii}(:, 4) == 2);
            type2_low_count = type2_low_count + sum(UFOdat_edited{ii}(:, 4) == 3);
            type2_high_count = type2_high_count + sum(UFOdat_edited{ii}(:, 4) == 4);
            unclassified_count = unclassified_count + sum(UFOdat_edited{ii}(:, 4) == 0);
        else
            unclassified_count = edited_count;  % All unclassified if no type column
        end
    end
end

fprintf('\n=== Summary ===\n');
fprintf('Original UFOs: %d\n', original_count);
fprintf('Edited UFOs: %d\n', edited_count);
fprintf('Removed: %d\n', original_count - edited_count);
if edited_count > 0 && (type1_low_count > 0 || type1_high_count > 0 || type2_low_count > 0 || type2_high_count > 0 || unclassified_count > 0)
    fprintf('\nType Classification:\n');
    fprintf('  Type 1, < 2kHz (genuine, low freq): %d\n', type1_low_count);
    fprintf('  Type 1, > 2kHz (genuine, high freq): %d\n', type1_high_count);
    fprintf('  Type 2, < 2kHz (artifact, low freq): %d\n', type2_low_count);
    fprintf('  Type 2, > 2kHz (artifact, high freq): %d\n', type2_high_count);
    if unclassified_count > 0
        fprintf('  Unclassified: %d\n', unclassified_count);
    end
    total_type2 = type2_low_count + type2_high_count;
    if exclude_type2 && total_type2 > 0
        fprintf('  (Type 2 UFOs (%d total) will be excluded from analysis)\n', total_type2);
    end
end
if current_idx <= total_ufos
    fprintf('\nProgress: %d/%d UFOs reviewed\n', current_idx - 1, total_ufos);
end

close all;
fprintf('\nDone!\n');

end

% Local helper function to get variable from base workspace with default
function v = getin(name, default)
    if evalin('base', sprintf('exist(''%s'',''var'')', name))
        v = evalin('base', name);
    else
        v = default;
    end
end

% Helper function for ternary operator
function y = ternary(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end

% Helper function to ensure UFOdat_edited has 4 columns [start_us, end_us, freq_hz, type]
% Type encoding: 0=unclassified, 1=Type1<2kHz, 2=Type1>2kHz, 3=Type2<2kHz, 4=Type2>2kHz
function UFOdat_out = ensure_4_columns(UFOdat_in)
    UFOdat_out = UFOdat_in;
    for ii = 1:numel(UFOdat_out)
        if ~isempty(UFOdat_out{ii})
            n_cols = size(UFOdat_out{ii}, 2);
            if n_cols == 3
                % Add 4th column (type, default 0 = unclassified)
                if isinteger(UFOdat_out{ii})
                    % For integer types, use same type for type column
                    type_col = zeros(size(UFOdat_out{ii}, 1), 1, 'like', UFOdat_out{ii});
                else
                    % For double/float, use int8 for type (0, 1, 2, 3, or 4)
                    type_col = zeros(size(UFOdat_out{ii}, 1), 1, 'int8');
                end
                UFOdat_out{ii} = [UFOdat_out{ii}, type_col];
            elseif n_cols == 4
                % Already has 4 columns, nothing to do
                continue;
            elseif n_cols == 5
                % Has 5 columns - likely from old format, keep only first 4
                warning('UFOdat_edited{%d} has 5 columns, keeping only first 4 columns', ii);
                UFOdat_out{ii} = UFOdat_out{ii}(:, 1:4);
            else
                error('UFOdat_edited{%d} has %d columns (expected 3, 4, or 5)', ii, n_cols);
            end
        end
    end
end


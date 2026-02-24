function summary_results = step3_statistics_fixed(patient, choices, N, equalize, min_duration_ms)
% STEP3_STATISTICS_FIXED  Calculate echo scores using fixed frequency-domain Green's envelope.
%
% This implements the main analysis with fixed parameters:
%   - Frequency band: [60, 100] Hz
%   - Macro natural frequency: f0 = 75 Hz
%   - Damping: gamma = 40 (1/s)
%   - Baseline window: pre_ms = 200 ms
%   - Post window: post_ms = 200 ms
%
% Statistical test: One p-value per micro-macro pair using matched non-UFO time nulls.
%
% Degenerate score (canonical definition): empirical statistic cannot be computed or
% yields NaN/Inf because of missing data / zero variance / empty window / empty null
% distribution. Degenerate pairs are excluded from p-value and FDR and counted as skipped.
%
% Inputs:
%   patient - Patient ID (required)
%   choices - Selection options [1-7] (optional, will prompt if not provided)
%   N       - Number of null/bootstrap iterations (optional, will prompt if not provided)
%   equalize - If true, equalize sample sizes across all selections (optional, default false)
%   min_duration_ms - Minimum UFO duration in milliseconds (optional, default: no filter)
%
% Outputs:
%   summary_results - Cell array of results for summary table

if nargin < 1 || isempty(patient)
    error('step3_statistics_fixed: patient ID required, e.g. step3_statistics_fixed(''pat001'')');
end

if nargin < 4 || isempty(equalize)
    equalize = false;
end

if nargin < 5 || isempty(min_duration_ms)
    min_duration_ms = 0;  % No filtering by default
end

% Fixed analysis parameters
ANALYSIS_CONFIG = struct();
ANALYSIS_CONFIG.band_hz = [60, 100];  % Frequency band (Hz)
ANALYSIS_CONFIG.f0_hz = 75;           % Macro natural frequency (Hz)
ANALYSIS_CONFIG.gamma = 40;           % Damping (1/s)
ANALYSIS_CONFIG.pre_ms = 200;         % Baseline window (ms)
ANALYSIS_CONFIG.post_ms = 200;        % Post-UFO window (ms)
ANALYSIS_CONFIG.nperseg_ms = 100;     % Welch segment length (ms) - defensible for 60-100 Hz
ANALYSIS_CONFIG.overlap = 0.5;        % Welch overlap fraction
ANALYSIS_CONFIG.nfft_min = 512;       % Minimum FFT length
ANALYSIS_CONFIG.one_sided = true;     % One-sided test (greater) for post > pre
% Some MATLAB MEF3 library builds have compatibility issues on certain platforms.
% Set force_python_mef = true to use the Python fallback for all patients.
ANALYSIS_CONFIG.force_python_mef = false;

summary_results = {};

% Paths (from ufo_config)
[base_dir, matmef_path] = ufo_config();
res_dir  = fullfile(base_dir, 'Results', patient);
log_path = fullfile(res_dir, 'ufo_inspect_log.mat');
data_dir = fullfile(base_dir, 'Data', patient);
ufo_json_path = fullfile(data_dir, 'UFO.json');

if ~exist(log_path,'file')
    error('step3_statistics_fixed: Log file not found at %s. Please run step1/2 first.', log_path);
end

% MEF defaults - auto-detect MEF directory
mef_path = fullfile(data_dir, 'sub.mefd');
if ~exist(mef_path, 'dir')
    % Try to find any .mefd directory in the data folder
    dir_list = dir(data_dir);
    mefd_dirs = dir_list([dir_list.isdir] & endsWith({dir_list.name}, '.mefd'));
    if ~isempty(mefd_dirs)
        mef_path = fullfile(data_dir, mefd_dirs(1).name);
    end
end
pw = 'bemena';

% Add matmef for readMef3
if exist(matmef_path,'dir'), addpath(matmef_path); end

% Load data
Slog = load(log_path,'Results'); Results = Slog.Results;

% --- PATIENT-SPECIFIC MACRO MAPPING FIXES ---
% pat001: ensure AMG micros map to AMG1, others to HH1
if strcmp(patient, 'pat001')
    fixed_count = 0;
    for i = 1:numel(Results)
        old_mac = Results(i).macro_channel;
        if contains(Results(i).micro_channel, 'AMG'), new_mac = 'AMG1';
        else, new_mac = 'HH1'; end
        if ~strcmp(old_mac, new_mac)
            Results(i).macro_channel = new_mac;
            fixed_count = fixed_count + 1;
        end
    end
    if fixed_count > 0
        fprintf('  Fixed %d macro channel mappings for pat001.\n', fixed_count);
        save(log_path, 'Results');
    end
end

% pat061 / pat072 / pat079: ensure macros are stored as B1/B2/B3 in Results
if any(strcmp(patient, {'pat061','pat072','pat079'}))
    fixed_count = 0;
    for i = 1:numel(Results)
        mac = Results(i).macro_channel;
        if strcmp(mac, 'B''1')
            Results(i).macro_channel = 'B1';
            fixed_count = fixed_count + 1;
        elseif strcmp(mac, 'B''2')
            Results(i).macro_channel = 'B2';
            fixed_count = fixed_count + 1;
        elseif strcmp(mac, 'B''3')
            Results(i).macro_channel = 'B3';
            fixed_count = fixed_count + 1;
        end
    end
    if fixed_count > 0
        fprintf('  Normalised %d B'' macros back to B1/B2/B3 for %s.\n', fixed_count, patient);
        save(log_path, 'Results');
    end
end

raw_json = jsondecode(fileread(ufo_json_path));

% Fallback bounds (UFO span only) - compute early so we can use when Python/MEF fails (e.g. pat448)
ufo_starts = [raw_json.uutc_left];
ufo_ends = [raw_json.uutc_right];
global_t_min_fallback = min(ufo_starts);
global_t_max_fallback = max(ufo_ends);

% Build cache of recording bounds per macro channel (MEF name -> [t_min, t_max])
% Get unique macro channels from Results
macro_channels = unique({Results.macro_channel});
macro_channels = macro_channels(~cellfun(@isempty, macro_channels));
if isempty(macro_channels)
    error('No macro channels found in Results for %s', patient);
end

% Build cache: containers.Map mapping mac_mef -> [t_min, t_max]
bounds_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
force_python_mef = ANALYSIS_CONFIG.force_python_mef;
fprintf('  Building recording bounds cache for %d macro channel(s)...\n', numel(macro_channels));

for mac_idx = 1:numel(macro_channels)
    mac_log = macro_channels{mac_idx};
    mac_mef = mac_log;
    if any(strcmp(patient, {'pat061','pat072','pat079'})) && startsWith(mac_log,'B')
        mac_mef = sprintf('B''%s', mac_log(2:end));
    end
    
    % Skip if already cached
    if bounds_cache.isKey(mac_mef)
        continue;
    end
    
    % Try to read MEF bounds for this macro channel
    try
        if force_python_mef
            % Use Python: bounds-only request (0, 0) returns session start/end
            [meta, ~] = readMef3_python(mef_path, pw, mac_mef, 'time', 0, 0);
            if isfield(meta, 'session_start_time') && isfield(meta, 'session_end_time')
                t_min = meta.session_start_time;
                t_max = meta.session_end_time;
            else
                error('Could not get session bounds from Python');
            end
        else
            m_info = readMef3(mef_path, pw, mac_mef);
            if isfield(m_info, 'uutc_start')
                t_min = m_info.uutc_start;
                t_max = m_info.uutc_end;
            elseif isfield(m_info, 'earliest_start_time')
                t_min = m_info.earliest_start_time;
                t_max = m_info.latest_end_time;
            else
                error('Could not determine recording bounds');
            end
        end
        % Ensure bounds are in microseconds (readMef3 returns microseconds, Python might return seconds)
        t_min = double(t_min);
        t_max = double(t_max);
        if t_min < 1e12  % Likely in seconds (Unix timestamp), convert to microseconds
            t_min = t_min * 1e6;
            t_max = t_max * 1e6;
        end
        bounds_cache(mac_mef) = [t_min, t_max];
        fprintf('    %s (MEF: %s): %.1f - %.1f seconds\n', mac_log, mac_mef, t_min/1e6, t_max/1e6);
    catch ME
        fprintf('    WARNING: Could not get bounds for %s (MEF: %s): %s\n', mac_log, mac_mef, ME.message);
        % If Python/MEF failed: cache fallback so analysis can run
        if force_python_mef
            force_python_mef = false;  % Disable Python for data reads to avoid repeated failures
            bounds_cache(mac_mef) = [global_t_min_fallback, global_t_max_fallback];
            fprintf('    %s (MEF: %s): using fallback bounds (UFO span)\n', mac_log, mac_mef);
        end
    end
end

% Global UFO exclusion windows: all UFOs from JSON with margins matching scored window
% Exclusion uses same window definition as scoring: [ufo_left - pre_us, ufo_right + post_us]
pre_us = round(ANALYSIS_CONFIG.pre_ms * 1e3);
post_us = round(ANALYSIS_CONFIG.post_ms * 1e3);
global_ufo_starts = [raw_json.uutc_left] - pre_us;
global_ufo_ends   = [raw_json.uutc_right] + post_us;

% Merge overlapping UFO intervals for faster overlap testing
% Sort by start time
[sorted_starts, sort_idx] = sort(global_ufo_starts);
sorted_ends = global_ufo_ends(sort_idx);
% Merge overlapping intervals
merged_starts = [];
merged_ends = [];
if ~isempty(sorted_starts)
    merged_starts(1) = sorted_starts(1);
    merged_ends(1) = sorted_ends(1);
    for i = 2:numel(sorted_starts)
        if sorted_starts(i) <= merged_ends(end)
            % Overlaps with current merged interval - extend it
            merged_ends(end) = max(merged_ends(end), sorted_ends(i));
        else
            % New non-overlapping interval
            merged_starts(end+1) = sorted_starts(i);
            merged_ends(end+1) = sorted_ends(i);
        end
    end
end

% Spectral config with fixed parameters
cfg = struct();
cfg.pre_ms      = ANALYSIS_CONFIG.pre_ms;
cfg.post_ms     = ANALYSIS_CONFIG.post_ms;
cfg.band_hz     = ANALYSIS_CONFIG.band_hz;
cfg.nperseg_ms  = ANALYSIS_CONFIG.nperseg_ms;
cfg.overlap     = ANALYSIS_CONFIG.overlap;
cfg.nfft_min    = ANALYSIS_CONFIG.nfft_min;
cfg.f0_hz       = ANALYSIS_CONFIG.f0_hz;
cfg.gamma       = ANALYSIS_CONFIG.gamma;
cfg.detrend     = true;
cfg.one_sided   = ANALYSIS_CONFIG.one_sided;

% Stats directory
stats_root = fullfile(res_dir, 'stats_fixed');
if ~exist(stats_root,'dir'), mkdir(stats_root); end

%% ------------------------------------------------------------------------
% 1) Selection and Check for Existing Results
%% ------------------------------------------------------------------------
% Prompt for choices and N if not provided
if nargin < 2 || isempty(choices)
    fprintf('  1) Type 1 only, all frequencies\n');
    fprintf('  2) Type 2 only, all frequencies\n');
    fprintf('  3) Both Type 1 & 2, all frequencies\n');
    fprintf('  4) Type 1 only, only >2kHz\n');
    fprintf('  5) Type 2 only, only >2kHz\n');
    fprintf('  6) Both Type 1 & 2, only >2kHz\n');
    fprintf('  7) ALL UFOs (original detections, ignore post-hoc labels)\n');
    choices = input('\nSelect option [1-7] (can be a vector, e.g. [4 5 7]): ');
    if isempty(choices), return; end
end

if nargin < 3 || isempty(N)
    N = input('\nEnter number of null/bootstrap iterations (e.g. 1000): ');
    if isempty(N) || N < 1, return; end
end

%% ------------------------------------------------------------------------
% 2) Process Each Selection
%% ------------------------------------------------------------------------
% If equalizing, first pass: find minimum n_sel across all selections
equalize_N = [];
if equalize
    min_n_sel = Inf;
    for choice = choices
        keep_mask = ~strcmp({Results.decision}, 'discard') & ~strcmp({Results.decision}, 'skip');
        if choice == 7
            R_set = Results;
        else
            R_set = Results(keep_mask);
        end
        switch choice
            case 1, R_sel = R_set([R_set.type] == 1);
            case 2, R_sel = R_set([R_set.type] == 2);
            case 3, R_sel = R_set([R_set.type] == 1 | [R_set.type] == 2);
            case 4, R_sel = R_set([R_set.type] == 1 & strcmp({R_set.freq_class}, 'high'));
            case 5, R_sel = R_set([R_set.type] == 2 & strcmp({R_set.freq_class}, 'high'));
            case 6, R_sel = R_set(([R_set.type] == 1 | [R_set.type] == 2) & strcmp({R_set.freq_class}, 'high'));
            case 7, R_sel = R_set;
        end
        n_sel = numel(R_sel);
        if n_sel > 0 && n_sel < min_n_sel
            min_n_sel = n_sel;
        end
    end
    if isfinite(min_n_sel)
        equalize_N = min_n_sel;
    else
        fprintf('  Warning: No valid selections found. Disabling equalization.\n');
        equalize = false;
    end
end

fprintf('\n>>> Processing %d selection(s) for %s...\n', numel(choices), patient);
choice_idx = 0;
for choice = choices
    choice_idx = choice_idx + 1;
    % Validate choice
    if choice < 1 || choice > 7
        fprintf('\n>>> Invalid choice %d (must be 1-7). Skipping.\n', choice);
        continue;
    end
    fprintf('\n--- Selection %d/%d (choice %d) ---\n', choice_idx, numel(choices), choice);
    
    keep_mask = ~strcmp({Results.decision}, 'discard') & ~strcmp({Results.decision}, 'skip');
    if choice == 7
        R_set = Results;
    else
        R_set = Results(keep_mask);
    end
    switch choice
        case 1, R_sel = R_set([R_set.type] == 1); desc = 'Type1_AllFreq';
        case 2, R_sel = R_set([R_set.type] == 2); desc = 'Type2_AllFreq';
        case 3, R_sel = R_set([R_set.type] == 1 | [R_set.type] == 2); desc = 'Type12_AllFreq';
        case 4, R_sel = R_set([R_set.type] == 1 & strcmp({R_set.freq_class}, 'high')); desc = 'Type1_HighFreq';
        case 5, R_sel = R_set([R_set.type] == 2 & strcmp({R_set.freq_class}, 'high')); desc = 'Type2_HighFreq';
        case 6, R_sel = R_set(([R_set.type] == 1 | [R_set.type] == 2) & strcmp({R_set.freq_class}, 'high')); desc = 'Type12_HighFreq';
        case 7
            R_sel = R_set;
            desc = 'AllUFOs_OrigJSON';
            for ii = 1:numel(R_sel)
                if isfield(R_sel, 'orig_start_us') && ~isempty(R_sel(ii).orig_start_us)
                    R_sel(ii).new_start_us = R_sel(ii).orig_start_us;
                end
            end
    end
    
    % Filter by minimum duration if specified
    if min_duration_ms > 0
        durations_ms = zeros(numel(R_sel), 1);
        for ii = 1:numel(R_sel)
            if isfield(R_sel, 'new_start_us') && isfield(R_sel, 'new_end_us') && ...
               ~isempty(R_sel(ii).new_start_us) && ~isempty(R_sel(ii).new_end_us)
                durations_ms(ii) = (R_sel(ii).new_end_us - R_sel(ii).new_start_us) / 1e3;
            else
                durations_ms(ii) = 0;
            end
        end
        R_sel = R_sel(durations_ms >= min_duration_ms);
        if ~isempty(R_sel) && min_duration_ms > 0
            desc = sprintf('%s_minDur%.1fms', desc, min_duration_ms);
        end
    end
    
    n_sel_orig = numel(R_sel);
    if n_sel_orig == 0
        desc_final = desc;
        if equalize && ~isempty(equalize_N)
            desc_final = sprintf('%s_Equalized', desc);
        end
        summary_results{end+1} = {patient, desc_final, 0, NaN, NaN};
        continue;
    end
    
    % Determine final description
    desc_final = desc;
    if equalize
        desc_final = [desc '_Equalized'];
    end
    
    % Equalize: randomly subsample to minimum size if enabled
    idx_keep = [];
    if equalize && n_sel_orig > equalize_N
        rng('shuffle');
        idx_keep = randperm(n_sel_orig, equalize_N);
        R_sel = R_sel(idx_keep);
        n_sel = equalize_N;
    else
        n_sel = n_sel_orig;
    end

    % Use desc_final for directory
    stats_dir = fullfile(stats_root, desc_final);
    if ~exist(stats_dir,'dir'), mkdir(stats_dir); end
    save_path = fullfile(stats_dir, 'results.mat');
    
    % Diagnostic: show cache status
    if exist(save_path, 'file')
        fprintf('  Cache file exists: %s\n', save_path);
    else
        fprintf('  No cache file found. Will calculate from scratch.\n');
    end

    %% ------------------------------------------------------------------------
    % 3) Calculate Empirical Scores
    %% ------------------------------------------------------------------------
    empirical_scores = NaN(n_sel, 1);
    
    if exist(save_path, 'file')
        % Try to load existing scores
        warning('off', 'MATLAB:load:variableNotFound');
        S_loaded = load(save_path, 'empirical_scores', 'desc_final', 'n_sel');
        warning('on', 'MATLAB:load:variableNotFound');
        if isfield(S_loaded, 'empirical_scores')
            loaded_scores = S_loaded.empirical_scores;
            loaded_desc = '';
            if isfield(S_loaded, 'desc_final')
                loaded_desc = S_loaded.desc_final;
            end
            loaded_n_sel = [];
            if isfield(S_loaded, 'n_sel')
                loaded_n_sel = S_loaded.n_sel;
            end
            if numel(loaded_scores) == n_sel && strcmp(loaded_desc, desc_final)
                empirical_scores = loaded_scores;
            else
                if numel(loaded_scores) ~= n_sel
                    fprintf('  Cache mismatch: loaded %d scores, need %d. Recalculating...\n', numel(loaded_scores), n_sel);
                elseif ~strcmp(loaded_desc, desc_final)
                    fprintf('  Cache mismatch: loaded desc="%s", need "%s". Recalculating...\n', loaded_desc, desc_final);
                end
            end
        end
    end
    
    if any(isnan(empirical_scores))
        n_missing = sum(isnan(empirical_scores));
        fprintf('  Calculating empirical scores (%d total, %d missing)...\n', n_sel, n_missing);
        
        % Pre-calculate constants to avoid repeated calculations
        pre_us_const = round(cfg.pre_ms * 1e3);
        post_us_const = round(cfg.post_ms * 1e3);
        read_before = max(cfg.pre_ms + 2, 10);
        read_after = max(cfg.post_ms + 2, 10);
        read_before_us = round(read_before * 1e3);
        read_after_us = round(read_after * 1e3);
        min_required_window_us = round((cfg.pre_ms + cfg.post_ms + 4) * 1e3);
        
        last_progress = 0;
        for i = 1:n_sel
            if ~isnan(empirical_scores(i)), continue; end
            % Progress reporting every 10%
            progress_pct = floor(100 * (i - 1) / n_sel);
            if progress_pct >= last_progress + 10
                fprintf('    Progress: %d%% (%d/%d)\n', progress_pct, i-1, n_sel);
                last_progress = progress_pct;
            end
            R = R_sel(i);
            try
                % Validate required fields
                if ~isfield(R, 'new_start_us') || isempty(R.new_start_us)
                    continue;
                end
                if ~isfield(R, 'macro_channel') || isempty(R.macro_channel)
                    continue;
                end
                
                mac_log = R.macro_channel;
                mac_mef = mac_log;
                if any(strcmp(patient, {'pat061','pat072','pat079'})) && startsWith(mac_log,'B')
                    mac_mef = sprintf('B''%s', mac_log(2:end));
                end
                
                % Get recording bounds for this macro channel (cached lookup is fast)
                if ~bounds_cache.isKey(mac_mef)
                    if i <= 3
                        fprintf('  WARNING: UFO %d - no bounds cache for macro %s\n', i, mac_mef);
                    end
                    continue;  % Skip if no bounds available
                end
                bounds = bounds_cache(mac_mef);
                read_t_min = double(bounds(1));
                read_t_max = double(bounds(2));
                
                % Quick bounds check: ensure UFO time is within bounds with margin
                ufo_time_us = double(R.new_start_us);
                if ufo_time_us < (read_t_min + pre_us_const) || ufo_time_us > (read_t_max - post_us_const)
                    if i <= 3
                        fprintf('  WARNING: UFO %d time (%.1f s) outside MEF bounds [%.1f, %.1f s] (margin: %.1f ms pre, %.1f ms post)\n', ...
                            i, ufo_time_us/1e6, read_t_min/1e6, read_t_max/1e6, cfg.pre_ms, cfg.post_ms);
                    end
                    continue;
                end
                
                % Calculate and clamp read window (ensure double precision)
                read_start = double(max(read_t_min, ufo_time_us - read_before_us));
                read_end = double(min(read_t_max, ufo_time_us + read_after_us));
                
                % Read MEF data
                if force_python_mef
                    [meta, tr] = readMef3_python(mef_path, pw, mac_mef, 'time', read_start, read_end);
                    if isempty(tr)
                        [~, tr] = readMef3(mef_path, pw, mac_mef, 'time', read_start, read_end);
                    end
                    if ~isempty(tr) && ~isempty(meta) && isfield(meta, 'sampling_frequency')
                        fs = double(meta.sampling_frequency);
                    elseif ~isempty(tr)
                        read_window_us = double(read_end - read_start);
                        if read_window_us > 0 && isfinite(read_window_us)
                            read_window_sec = read_window_us / 1e6;
                            if read_window_sec > 0 && isfinite(read_window_sec)
                                fs = double(numel(tr)) / read_window_sec;
                            else
                                fs = NaN;
                            end
                        else
                            fs = NaN;
                        end
                    else
                        fs = NaN;
                    end
                else
                    [~, tr] = readMef3(mef_path, pw, mac_mef, 'time', read_start, read_end);
                    read_window_us = double(read_end - read_start);
                    if read_window_us > 0 && isfinite(read_window_us)
                        read_window_sec = read_window_us / 1e6;
                        if read_window_sec > 0 && isfinite(read_window_sec)
                            fs = double(numel(tr)) / read_window_sec;
                        else
                            fs = NaN;
                        end
                    else
                        fs = NaN;
                    end
                end
                
                % Extract window if data is valid
                if ~isempty(tr) && isfinite(fs) && fs > 0 && fs <= 1e6
                    % Calculate sample index (ufo_time_us is already validated to be in [read_start, read_end])
                    time_diff_us = ufo_time_us - read_start;
                    time_diff_sec = time_diff_us / 1e6;
                    ufo_start_idx = round(time_diff_sec * fs) + 1;
                    
                    % Validate and extract window
                    if ufo_start_idx >= 1 && ufo_start_idx <= numel(tr)
                        pre_samples = round(cfg.pre_ms * fs / 1000);
                        post_samples = round(cfg.post_ms * fs / 1000);
                        
                        % Exact-length extraction: skip if out of bounds
                        idx_start = ufo_start_idx - pre_samples + 1;
                        idx_end = ufo_start_idx + post_samples;
                        
                        if idx_start >= 1 && idx_end <= numel(tr)
                            idx_range = idx_start:idx_end;
                            empirical_scores(i) = greens_envelope_score(tr(idx_range), fs, cfg);
                        elseif i <= 3
                            fprintf('  WARNING: UFO %d window out of bounds (idx_start=%d, idx_end=%d, tr_length=%d, ufo_idx=%d)\n', ...
                                i, idx_start, idx_end, numel(tr), ufo_start_idx);
                        end
                    elseif i <= 3
                        fprintf('  WARNING: UFO %d invalid start index (idx=%d, tr_length=%d, time_diff_sec=%.6f, fs=%.1f)\n', ...
                            i, ufo_start_idx, numel(tr), time_diff_sec, fs);
                    end
                elseif i <= 3
                    if isempty(tr)
                        fprintf('  WARNING: UFO %d returned empty data from MEF read\n', i);
                    elseif ~isfinite(fs) || fs <= 0 || fs > 1e6
                        fprintf('  WARNING: UFO %d invalid fs (fs=%.1f, tr_length=%d, read_window_us=%.1f)\n', ...
                            i, fs, numel(tr), read_end - read_start);
                    end
                end
            catch ME
                if i <= 3 || mod(i, 50) == 0
                    fprintf('  ERROR processing UFO %d: %s\n', i, ME.message);
                end
            end
        end
        fprintf('    Empirical scores complete (%d/%d calculated)\n', sum(isfinite(empirical_scores)), n_sel);
        % IMPORTANT: preserve any existing cached null/pair results.
        % A previous run may have already computed and saved pair_pvals, etc.
        try
            if exist(save_path, 'file')
                S_existing = load(save_path);
                S_existing.empirical_scores = empirical_scores;
                S_existing.desc_final = desc_final;
                S_existing.n_sel = n_sel;
                save(save_path, '-struct', 'S_existing', '-v7.3');
            else
                save(save_path, 'empirical_scores', 'desc_final', 'n_sel', '-v7.3');
            end
        catch ME_save_emp
            fprintf('  WARNING: Failed to save empirical scores cache: %s\n', ME_save_emp.message);
        end
    else
        fprintf('  Using cached empirical scores (%d total)\n', n_sel);
    end
    
    emp_avg = mean(empirical_scores, 'omitnan');
    
    if isnan(emp_avg)
        fprintf('  WARNING: Empirical average is NaN for %s. Skipping nulls for this selection.\n', desc);
        summary_results{end+1} = {patient, desc_final, n_sel, NaN, NaN, NaN};
        continue;
    end

    %% ------------------------------------------------------------------------
    % 4) Statistical Test: One p-value per micro-macro pair
    %% ------------------------------------------------------------------------
    % Group UFOs by micro-macro pair (needed to check completeness)
    pair_keys = cell(n_sel, 1);
    for i = 1:n_sel
        pair_keys{i} = sprintf('%s_%s', R_sel(i).micro_channel, R_sel(i).macro_channel);
    end
    unique_pairs = unique(pair_keys);
    N_pairs_attempted = numel(unique_pairs);
    
    % Check if complete results already exist (all pairs have a result, N_eff consistent)
    % Note: null_scores_* are kept for extend (1k->10k); do not treat as incomplete just because they exist.
    results_complete = false;
    has_checkpoint_vars = false;  % used only for diagnostic message
    if exist(save_path, 'file')
        warning('off', 'MATLAB:load:variableNotFound');
        S_check = load(save_path, 'pair_pvals', 'pair_null_means', 'pair_N_eff', 'N', 'desc_final', 'n_total_pairs');
        warning('on', 'MATLAB:load:variableNotFound');
        if isfield(S_check, 'pair_pvals') && isfield(S_check, 'pair_null_means') && ...
               isfield(S_check, 'N') && isfield(S_check, 'desc_final') && ...
               S_check.N == N && strcmp(S_check.desc_final, desc_final)
                % Only complete if every pair has an entry (resume uses partial)
                all_have_result = true;
                for pp = 1:numel(unique_pairs)
                    pn = matlab.lang.makeValidName(unique_pairs{pp});
                    if ~isfield(S_check.pair_pvals, pn)
                        all_have_result = false;
                        break;
                    end
                end
                % Require stored null run to match requested N (don't use old N=1000 cache for N=10000)
                N_eff_consistent = false;
                if all_have_result && isfield(S_check, 'pair_N_eff')
                    for pp = 1:numel(unique_pairs)
                        pn = matlab.lang.makeValidName(unique_pairs{pp});
                        if isfield(S_check.pair_pvals, pn) && isfinite(S_check.pair_pvals.(pn)) && ...
                           isfield(S_check.pair_N_eff, pn) && S_check.pair_N_eff.(pn) >= N * 0.5
                            N_eff_consistent = true;
                            break;
                        end
                    end
                    % If no tested pair had N_eff loaded, allow (e.g. old file without pair_N_eff)
                    if ~N_eff_consistent && all_have_result
                        tested_any = false;
                        for pp = 1:numel(unique_pairs)
                            pn = matlab.lang.makeValidName(unique_pairs{pp});
                            if isfield(S_check.pair_pvals, pn) && isfinite(S_check.pair_pvals.(pn))
                                tested_any = true;
                                break;
                            end
                        end
                        if tested_any
                            N_eff_consistent = false;  % have tested pairs but none with N_eff >= 0.5*N
                        else
                            N_eff_consistent = true;  % no tested pairs, allow
                        end
                    end
                elseif all_have_result
                    N_eff_consistent = true;  % no pair_N_eff in file, allow (legacy)
                end
                if all_have_result && N_eff_consistent
                    results_complete = true;
                    n_cached_pairs = 0;
                    if isfield(S_check, 'n_total_pairs')
                        n_cached_pairs = S_check.n_total_pairs;
                    end
                    fprintf('  Found complete cached results (N=%d, %d pairs). Loading...\n', N, n_cached_pairs);
                elseif all_have_result && ~N_eff_consistent
                    fprintf('  Cache N_eff inconsistent with requested N=%d (e.g. from older run). Will recalculate nulls.\n', N);
                end
            end
        if ~results_complete
            if ~isfield(S_check, 'pair_pvals')
                fprintf('  Cache missing pair_pvals. Will recalculate nulls.\n');
            elseif ~isfield(S_check, 'N') || S_check.N ~= N
                fprintf('  Cache N mismatch (cached: %d, requested: %d). Will recalculate nulls.\n', ...
                    isfield(S_check, 'N') && S_check.N, N);
            elseif ~isfield(S_check, 'desc_final') || ~strcmp(S_check.desc_final, desc_final)
                % Changing equalize or min_duration_ms changes desc_final and the UFO set, so per-pair K_eff
                % and the null distribution (mean of K_eff scores) change; nulls cannot be reused from a base run.
                fprintf('  Cache desc mismatch (cached: "%s", requested: "%s"). Will recalculate nulls.\n', ...
                    isfield(S_check, 'desc_final') && S_check.desc_final, desc_final);
            else
                fprintf('  Incomplete cache (not all pairs have results). Will resume or recalculate.\n');
            end
        end
    else
        fprintf('  No cache file found at %s. Will calculate from scratch.\n', save_path);
    end
    
    % Skip counters for pair accounting (attempted vs testable vs tested)
    % Degenerate score = empirical statistic cannot be computed or yields NaN/Inf because of
    % missing data / zero variance / empty window / empty null distribution / etc.
    % Degenerate pairs are excluded from p-value and FDR and counted under "Skipped: Degenerate empirical score".
    skip_no_usable_events = 0;
    skip_out_of_bounds = 0;
    skip_degenerate = 0;
    skip_other = 0;
    
    % Calculate mean score per pair (will be recomputed per pair with K_eff)
    pair_scores = struct();
    pair_ufo_counts = struct();
    pair_ufo_counts_raw = struct();
    pair_ufo_counts_eff = struct();
    
    % For each pair, generate null distribution
    pair_pvals = struct();
    pair_null_means = struct();
    pair_null_stds = struct();
    pair_N_eff = struct();
    
    % Load cached results if complete, or partial results for resume-after-interrupt
    if results_complete
        warning('off', 'MATLAB:load:variableNotFound');
        S_loaded = load(save_path, 'pair_scores', 'pair_pvals', 'pair_null_means', 'pair_null_stds', ...
            'pair_ufo_counts', 'pair_ufo_counts_raw', 'pair_ufo_counts_eff', 'pair_N_eff', ...
            'n_significant', 'n_total_pairs', 'fdr_corrected', 'n_fdr_significant');
        warning('on', 'MATLAB:load:variableNotFound');
        if isfield(S_loaded, 'pair_pvals')
            pair_pvals = S_loaded.pair_pvals;
            pair_null_means = S_loaded.pair_null_means;
            pair_null_stds = S_loaded.pair_null_stds;
            pair_scores = S_loaded.pair_scores;
            pair_ufo_counts = S_loaded.pair_ufo_counts;
            pair_ufo_counts_raw = S_loaded.pair_ufo_counts_raw;
            pair_ufo_counts_eff = S_loaded.pair_ufo_counts_eff;
            pair_N_eff = S_loaded.pair_N_eff;
            fprintf('  Using cached null results. Skipping null sampling.\n');
            % Skip to results display section
            skip_null_sampling = true;
        else
            skip_null_sampling = false;
        end
    else
        skip_null_sampling = false;
        % Load any partial pair results so we can resume after interrupt (pick up at next pair)
        if exist(save_path, 'file')
            warning('off', 'MATLAB:load:variableNotFound');
            S_partial = load(save_path, 'pair_pvals', 'pair_scores', 'pair_null_means', 'pair_null_stds', ...
                'pair_ufo_counts', 'pair_ufo_counts_raw', 'pair_ufo_counts_eff', 'pair_N_eff', ...
                'N', 'desc_final');
            warning('on', 'MATLAB:load:variableNotFound');
            if isfield(S_partial, 'pair_pvals') && isfield(S_partial, 'N') && isfield(S_partial, 'desc_final') && ...
               S_partial.N == N && strcmp(S_partial.desc_final, desc_final)
                pair_pvals = S_partial.pair_pvals;
                if isfield(S_partial, 'pair_scores'), pair_scores = S_partial.pair_scores; end
                if isfield(S_partial, 'pair_null_means'), pair_null_means = S_partial.pair_null_means; end
                if isfield(S_partial, 'pair_null_stds'), pair_null_stds = S_partial.pair_null_stds; end
                if isfield(S_partial, 'pair_ufo_counts'), pair_ufo_counts = S_partial.pair_ufo_counts; end
                if isfield(S_partial, 'pair_ufo_counts_raw'), pair_ufo_counts_raw = S_partial.pair_ufo_counts_raw; end
                if isfield(S_partial, 'pair_ufo_counts_eff'), pair_ufo_counts_eff = S_partial.pair_ufo_counts_eff; end
                if isfield(S_partial, 'pair_N_eff'), pair_N_eff = S_partial.pair_N_eff; end
                % Do not load skip_* — count only this run's skips so pair accounting is correct
                n_resumed = sum(cellfun(@(pk) isfield(pair_pvals, matlab.lang.makeValidName(pk)) && isfinite(pair_pvals.(matlab.lang.makeValidName(pk))), unique_pairs));
                if n_resumed > 0
                    fprintf('  Resuming: %d pair(s) already have results; will compute remaining pairs.\n', n_resumed);
                end
            end
        end
    end
    
    fprintf('  Enumerating candidate pairs: %d\n', N_pairs_attempted);
    
    % Only run null sampling if we don't have complete cached results
    if ~skip_null_sampling
        for p = 1:numel(unique_pairs)
        pair_key = unique_pairs{p};
        pair_name = matlab.lang.makeValidName(pair_key);
        
        % Get UFOs for this pair
        pair_mask = strcmp(pair_keys, pair_key);
        pair_ufos = R_sel(pair_mask);
        K_raw = sum(pair_mask);
        
        if K_raw == 0
            skip_no_usable_events = skip_no_usable_events + 1;
            save_partial_pair_results(save_path, pair_pvals, pair_scores, pair_null_means, pair_null_stds, ...
                pair_ufo_counts, pair_ufo_counts_raw, pair_ufo_counts_eff, pair_N_eff, ...
                skip_no_usable_events, skip_out_of_bounds, skip_degenerate, skip_other, N, desc_final);
            continue;
        end
        
        % Skip only if this pair has a result computed with the requested N (not from older/smaller N)
        if isfield(pair_pvals, pair_name)
            n_eff_this = [];
            if isfield(pair_N_eff, pair_name)
                n_eff_this = pair_N_eff.(pair_name);
            end
            if ~isempty(n_eff_this) && n_eff_this >= N * 0.5
                fprintf('    Pair %d/%d: %s — using saved result (skip).\n', p, numel(unique_pairs), pair_key);
                continue;
            end
            % Has entry but N_eff too small (e.g. 1k when asking 10k) — will extend below
        end
        
        fprintf('    Pair %d/%d: %s (K_raw=%d)...\n', p, numel(unique_pairs), pair_key, K_raw);
        
        % Get valid (finite) scores for this pair
        pair_scores_vec = empirical_scores(pair_mask);
        valid_pair_scores = pair_scores_vec(isfinite(pair_scores_vec));
        K_eff = numel(valid_pair_scores);
        
        % Guard: need at least 3 valid events for meaningful statistics
        if K_eff < 3
            skip_no_usable_events = skip_no_usable_events + 1;
            fprintf('  WARNING: Pair %s has too few valid events (K_raw=%d, K_eff=%d). Skipping nulls.\n', ...
                pair_key, K_raw, K_eff);
            pair_pvals.(pair_name) = NaN;
            pair_null_means.(pair_name) = NaN;
            pair_null_stds.(pair_name) = NaN;
            pair_ufo_counts.(pair_name) = K_raw;  % Store raw count
            pair_ufo_counts_raw.(pair_name) = K_raw;
            pair_ufo_counts_eff.(pair_name) = K_eff;
            pair_scores.(pair_name) = NaN;
            pair_N_eff.(pair_name) = 0;
            save_partial_pair_results(save_path, pair_pvals, pair_scores, pair_null_means, pair_null_stds, ...
                pair_ufo_counts, pair_ufo_counts_raw, pair_ufo_counts_eff, pair_N_eff, ...
                skip_no_usable_events, skip_out_of_bounds, skip_degenerate, skip_other, N, desc_final);
            continue;
        end
        
        % Get macro channel for this pair
        macro_ch = pair_ufos(1).macro_channel;
        mac_mef = macro_ch;
        if any(strcmp(patient, {'pat061','pat072','pat079'})) && startsWith(macro_ch,'B')
            mac_mef = sprintf('B''%s', macro_ch(2:end));
        end
        
        % Get recording bounds for this specific macro channel from cache
        if bounds_cache.isKey(mac_mef)
            bounds = bounds_cache(mac_mef);
            t_min = bounds(1);
            t_max = bounds(2);
        else
            fprintf('  WARNING: No bounds cached for %s (MEF: %s). Using global fallback.\n', macro_ch, mac_mef);
            t_min = global_t_min_fallback;
            t_max = global_t_max_fallback;
        end
        pre_us = round(cfg.pre_ms * 1e3);
        post_us = round(cfg.post_ms * 1e3);
        
        % Ensure boundaries allow at least one window
        % Window spans [null_time - pre_us, null_time + post_us]
        min_null_time = t_min + pre_us;
        max_null_time = t_max - post_us;
        if max_null_time <= min_null_time
            skip_out_of_bounds = skip_out_of_bounds + 1;
            fprintf('  WARNING: Recording boundaries too narrow for pair %s. Skipping nulls.\n', pair_key);
            pair_pvals.(pair_name) = NaN;
            pair_null_means.(pair_name) = NaN;
            pair_null_stds.(pair_name) = NaN;
            pair_ufo_counts.(pair_name) = K_raw;
            pair_ufo_counts_raw.(pair_name) = K_raw;
            pair_ufo_counts_eff.(pair_name) = K_eff;
            pair_scores.(pair_name) = NaN;
            pair_N_eff.(pair_name) = 0;
            save_partial_pair_results(save_path, pair_pvals, pair_scores, pair_null_means, pair_null_stds, ...
                pair_ufo_counts, pair_ufo_counts_raw, pair_ufo_counts_eff, pair_N_eff, ...
                skip_no_usable_events, skip_out_of_bounds, skip_degenerate, skip_other, N, desc_final);
            continue;
        end
        
        % Check if this pair already has results (when loading from cache)
        if skip_null_sampling
            if isfield(pair_pvals, pair_name) && isfinite(pair_pvals.(pair_name))
                fprintf('      Using cached null results for pair %s\n', pair_key);
            else
                % This pair wasn't in cache, need to compute it
                skip_null_sampling = false;
            end
        end
        
        if skip_null_sampling && isfield(pair_pvals, pair_name) && isfinite(pair_pvals.(pair_name))
            continue;
        end
        
        % Generate null scores (use K_eff to match empirical statistic)
        null_scores = NaN(N, K_eff);
        max_tries_per_iter = 100;  % Maximum tries per null sample
        
        % Try to load partial progress for this pair (same N resume, or extend e.g. 1000 -> 10000)
        checkpoint_var_name = sprintf('null_scores_%s', pair_name);
        n_completed = 0;
        if exist(save_path, 'file')
            warning('off', 'MATLAB:load:variableNotFound');
            S_checkpoint = load(save_path, checkpoint_var_name, 'N', 'desc_final');
            warning('on', 'MATLAB:load:variableNotFound');
            if isfield(S_checkpoint, checkpoint_var_name) && ...
               isfield(S_checkpoint, 'desc_final') && strcmp(S_checkpoint.desc_final, desc_final)
                loaded_null_scores = S_checkpoint.(checkpoint_var_name);
                n_loaded = size(loaded_null_scores, 1);
                min_cols = min(size(loaded_null_scores, 2), K_eff);
                if min_cols > 0 && n_loaded > 0
                    % Copy up to min(n_loaded, N) rows; allows extending (e.g. keep 1000, do 9000 more)
                    n_copy = min(n_loaded, N);
                    null_scores(1:n_copy, 1:min_cols) = loaded_null_scores(1:n_copy, 1:min_cols);
                    n_completed = sum(any(isfinite(null_scores(:, 1:min_cols)), 2));
                    if n_completed > 0
                        if n_loaded >= N
                            fprintf('      Resuming null distribution: %d/%d iterations already completed (K_eff: %d->%d)\n', ...
                                n_completed, N, size(loaded_null_scores, 2), K_eff);
                        else
                            fprintf('      Extending null distribution: reusing %d iterations, computing %d more (target N=%d)\n', ...
                                n_completed, N - n_completed, N);
                        end
                    end
                end
            else
                if exist(save_path, 'file') && (~isfield(S_checkpoint, checkpoint_var_name) || ~isfield(S_checkpoint, 'desc_final') || ~strcmp(S_checkpoint.desc_final, desc_final))
                    fprintf('      No checkpoint found for pair %s (or parameter mismatch). Starting fresh.\n', pair_key);
                end
            end
        end
        
        if n_completed < N
            fprintf('      Generating null distribution: %d remaining (N=%d, K_eff=%d)...\n', N - n_completed, N, K_eff);
            
            % Initialize parallel pool if not already running
            p = gcp('nocreate');
            if isempty(p)
                want = max(1, feature('numcores') - 2);
                parpool('local', want);
                fprintf('      Initialized parallel pool with %d workers\n', want);
            else
                fprintf('      Using existing parallel pool (%d workers)\n', p.NumWorkers);
            end
            
            last_progress = 0;
            last_save_iter = 0;
            checkpoint_interval = max(100, floor(N / 20));  % Save every 5% or every 100 iterations, whichever is larger
            
            % Pre-compute random seeds for each iteration (for reproducibility)
            base_seed_local = randi(2^31-1);  % Random base seed for this pair
            iter_seeds = mod(base_seed_local + (1:N), 2^31-1);
            
            % Broadcast variables needed in parfor (MATLAB requirement)
        merged_starts_par = merged_starts;
        merged_ends_par = merged_ends;
        t_min_par = t_min;
        t_max_par = t_max;
        pre_us_par = pre_us;
        post_us_par = post_us;
        min_null_time_par = min_null_time;
        max_null_time_par = max_null_time;
        mac_mef_par = mac_mef;
        mef_path_par = mef_path;
        pw_par = pw;
        cfg_par = cfg;
        force_python_mef_par = force_python_mef;
        
            % Parallelize the outer loop - only run remaining iterations (keep n_completed)
        iter_range = (n_completed + 1):N;
        parfor_ok = false;
        fprintf('      Running parfor (no progress output until complete; may take several minutes)...\n');
        try
            parfor iter = iter_range
                % Set random seed for this iteration (for reproducibility)
                rng(iter_seeds(iter), 'twister');
                
                for k = 1:K_eff
                % Sample a random time, excluding ALL UFO windows from JSON
                valid_time = false;
                null_time = [];
                tried = 0;
                
                while tried < max_tries_per_iter && ~valid_time
                    % Sample uniformly from valid range
                    % null_time must be in [t_min + pre_us, t_max - post_us]
                    span = double(max_null_time_par - min_null_time_par);
                    if span <= 1
                        tried = tried + 1;
                        continue;
                    end
                    % Use rand() instead of randi() to avoid overflow for large spans (multi-hour recordings)
                    null_time = min_null_time_par + floor(rand() * span);
                    
                    % Define null window matching scored window: [null_time - pre_us, null_time + post_us]
                    null_start = null_time - pre_us_par;
                    null_end = null_time + post_us_par;
                    
                    % Note: Can't use fprintf inside parfor, so sanity printout moved outside
                    
                    % Ensure we don't exceed boundaries
                    if null_start < t_min_par || null_end > t_max_par
                        tried = tried + 1;
                        continue;
                    end
                    
                    % Check overlap with merged UFO exclusion intervals (fast binary search)
                    overlaps = false;
                    if ~isempty(merged_starts_par)
                        % Binary search: find rightmost interval where start <= null_end
                        left = 1;
                        right = numel(merged_starts_par);
                        idx = 0;
                        while left <= right
                            mid = floor((left + right) / 2);
                            if merged_starts_par(mid) <= null_end
                                idx = mid;
                                left = mid + 1;
                            else
                                right = mid - 1;
                            end
                        end
                        % Check candidate interval(s) for overlap
                        % Overlap condition: null_start <= merged_ends AND null_end >= merged_starts
                        if idx > 0
                            % Check interval idx
                            if null_start <= merged_ends_par(idx)
                                overlaps = true;
                            end
                            % Also check interval idx+1 (could overlap even if idx doesn't)
                            if ~overlaps && (idx < numel(merged_starts_par))
                                if merged_starts_par(idx+1) <= null_end && null_start <= merged_ends_par(idx+1)
                                    overlaps = true;
                                end
                            end
                        end
                        % idx==0 means all starts > null_end, so no overlap (overlaps remains false)
                    end
                    
                    if ~overlaps
                        valid_time = true;
                    else
                        tried = tried + 1;
                    end
                end
                
                if ~valid_time
                    % Could not find valid time after max tries - leave as NaN
                    continue;
                end
                
                if valid_time
                    try
                        read_before = max(cfg.pre_ms + 2, 10);
                        read_after = max(cfg.post_ms + 2, 10);
                        
                        % Clamp read window to MEF recording bounds to avoid warnings
                        read_start = max(t_min_par, null_time - round(read_before*1e3));
                        read_end = min(t_max_par, null_time + round(read_after*1e3));
                        
                        if force_python_mef_par
                            [meta, tr] = readMef3_python(mef_path_par, pw_par, mac_mef_par, 'time', read_start, read_end);
                            if isempty(tr)
                                [~, tr] = readMef3(mef_path_par, pw_par, mac_mef_par, 'time', read_start, read_end);
                            end
                            if ~isempty(tr) && ~isempty(meta) && isfield(meta, 'sampling_frequency')
                                fs = meta.sampling_frequency;
                            elseif ~isempty(tr)
                                fs = double(numel(tr)) / (double((read_end - read_start) / 1e6));
                            else
                                fs = NaN;
                            end
                        else
                            [~, tr] = readMef3(mef_path_par, pw_par, mac_mef_par, 'time', read_start, read_end);
                            fs = double(numel(tr)) / (double((read_end - read_start) / 1e6));
                        end
                        
                        if ~isempty(tr)
                            % Calculate ufo_start_idx based on actual position of null_time in read data
                            % (accounting for potential clamping of read window)
                            ufo_start_idx = round((null_time - read_start) * fs / 1e6) + 1;
                            pre_samples = round(cfg_par.pre_ms * fs / 1000);
                            post_samples = round(cfg_par.post_ms * fs / 1000);
                            
                            % Exact-length extraction: skip if out of bounds
                            idx_start = ufo_start_idx - pre_samples + 1;
                            idx_end = ufo_start_idx + post_samples;
                            
                            if idx_start < 1 || idx_end > numel(tr)
                                % Out of bounds - skip this null sample
                                continue;
                            end
                            
                            idx_range = idx_start:idx_end;
                            
                            null_scores(iter, k) = greens_envelope_score(tr(idx_range), fs, cfg_par);
                        end
                    catch
                        % Skip if read fails
                    end
                end
            end
            end  % End parfor
            parfor_ok = true;
        catch ME_par
            fprintf('      WARNING: parfor failed (%s). Retrying null sampling serially...\n', ME_par.message);
        end
        if ~parfor_ok
            % Fallback: run same loop serially (no parallel pool needed)
            n_serial = numel(iter_range);
            progress_interval = max(500, floor(n_serial / 20));
            for ii = 1:n_serial
                iter = iter_range(ii);
                if ii > 1 && mod(ii - 1, progress_interval) == 0
                    fprintf('      Serial progress: %d / %d iterations...\n', ii - 1, n_serial);
                end
                rng(iter_seeds(iter), 'twister');
                for k = 1:K_eff
                    valid_time = false;
                    null_time = [];
                    tried = 0;
                    while tried < max_tries_per_iter && ~valid_time
                        span = double(max_null_time_par - min_null_time_par);
                        if span <= 1, tried = tried + 1; continue; end
                        null_time = min_null_time_par + floor(rand() * span);
                        null_start = null_time - pre_us_par;
                        null_end = null_time + post_us_par;
                        if null_start < t_min_par || null_end > t_max_par, tried = tried + 1; continue; end
                        overlaps = false;
                        if ~isempty(merged_starts_par)
                            left = 1; right = numel(merged_starts_par); idx = 0;
                            while left <= right
                                mid = floor((left + right) / 2);
                                if merged_starts_par(mid) <= null_end, idx = mid; left = mid + 1; else right = mid - 1; end
                            end
                            if idx > 0 && null_start <= merged_ends_par(idx), overlaps = true; end
                            if ~overlaps && idx < numel(merged_starts_par) && merged_starts_par(idx+1) <= null_end && null_start <= merged_ends_par(idx+1), overlaps = true; end
                        end
                        if ~overlaps, valid_time = true; else tried = tried + 1; end
                    end
                    if ~valid_time, continue; end
                    try
                        read_before = max(cfg.pre_ms + 2, 10); read_after = max(cfg.post_ms + 2, 10);
                        read_start = max(t_min_par, null_time - round(read_before*1e3)); read_end = min(t_max_par, null_time + round(read_after*1e3));
                        if force_python_mef_par
                            [meta, tr] = readMef3_python(mef_path_par, pw_par, mac_mef_par, 'time', read_start, read_end);
                            if isempty(tr), [~, tr] = readMef3(mef_path_par, pw_par, mac_mef_par, 'time', read_start, read_end); end
                            if ~isempty(tr) && ~isempty(meta) && isfield(meta, 'sampling_frequency'), fs = meta.sampling_frequency; elseif ~isempty(tr), fs = double(numel(tr)) / (double((read_end - read_start) / 1e6)); else fs = NaN; end
                        else
                            [~, tr] = readMef3(mef_path_par, pw_par, mac_mef_par, 'time', read_start, read_end);
                            fs = double(numel(tr)) / (double((read_end - read_start) / 1e6));
                        end
                        if ~isempty(tr)
                            ufo_start_idx = round((null_time - read_start) * fs / 1e6) + 1;
                            pre_samples = round(cfg_par.pre_ms * fs / 1000); post_samples = round(cfg_par.post_ms * fs / 1000);
                            idx_start = ufo_start_idx - pre_samples + 1; idx_end = ufo_start_idx + post_samples;
                            if idx_start >= 1 && idx_end <= numel(tr)
                                null_scores(iter, k) = greens_envelope_score(tr(idx_start:idx_end), fs, cfg_par);
                            end
                        end
                    catch
                    end
                end
            end
            fprintf('      Serial null sampling completed (%d iterations).\n', n_serial);
        end
        end  % End if n_completed < N
        
        % Sanity printout for first pair (after parallel section)
        if p == 1
            fprintf('  Null window definition (first pair): pre_us=%d, post_us=%d\n', pre_us, post_us);
            if ~isempty(null_scores) && any(isfinite(null_scores(:)))
                % Find first valid null sample for example
                [first_iter, first_k] = find(isfinite(null_scores), 1);
                if ~isempty(first_iter)
                    % Reconstruct example (approximate)
                    example_null_time = min_null_time + floor(rand() * double(max_null_time - min_null_time));
                    fprintf('  Example null window (first pair): null_start=%.1f, null_time=%.1f, null_end=%.1f (seconds)\n', ...
                        (example_null_time - pre_us)/1e6, example_null_time/1e6, (example_null_time + post_us)/1e6);
                end
            end
        end
        
        % Periodic checkpointing: save progress after parallel section
        % (Checkpointing inside parfor is not thread-safe, so we do it after)
        % NOTE: If MATLAB crashes during the parfor loop, progress for this pair will be lost
        % and it will restart from the beginning. This is a limitation of parfor.
        try
            % Save partial progress for this pair
            checkpoint_data = struct(checkpoint_var_name, null_scores, 'N', N, 'desc_final', desc_final);
            if exist(save_path, 'file')
                % Append to existing file
                S_existing = load(save_path);
                S_existing.(checkpoint_var_name) = null_scores;
                S_existing.N = N;
                S_existing.desc_final = desc_final;
                save(save_path, '-struct', 'S_existing', '-v7.3');
            else
                % Create new file
                save(save_path, '-struct', 'checkpoint_data', '-v7.3');
            end
            n_finite = sum(any(isfinite(null_scores), 2));
            fprintf('        Checkpoint saved: %d/%d iterations complete (%d finite)\n', n_finite, N, n_finite);
        catch ME_checkpoint
            fprintf('        WARNING: Failed to save checkpoint: %s\n', ME_checkpoint.message);
        end
        fprintf('      Null distribution complete (N_eff will be calculated)\n');
        
        % Calculate empirical mean using valid scores (matches null sampling with K_eff)
        emp_mean = mean(valid_pair_scores);
        pair_scores.(pair_name) = emp_mean;  % Store for output
        pair_ufo_counts.(pair_name) = K_eff;  % Store effective count
        pair_ufo_counts_raw.(pair_name) = K_raw;  % Store raw count
        pair_ufo_counts_eff.(pair_name) = K_eff;  % Store effective count
        
        % Calculate null mean per iteration
        null_means = mean(null_scores, 2, 'omitnan');
        null_means = null_means(isfinite(null_means));
        N_eff = numel(null_means);
        
        if ~isempty(null_means)
            pair_null_means.(pair_name) = mean(null_means);
            pair_null_stds.(pair_name) = std(null_means);
            pair_N_eff.(pair_name) = N_eff;
            
            % Monte Carlo p-value with +1 correction: p = (1 + #(null >= emp)) / (N_null + 1)
            if cfg.one_sided
                p_val = (1 + sum(null_means >= emp_mean)) / (N_eff + 1);
            else
                p_hi = (1 + sum(null_means >= emp_mean)) / (N_eff + 1);
                p_lo = (1 + sum(null_means <= emp_mean)) / (N_eff + 1);
                p_val = min(1, 2*min(p_hi, p_lo));
            end
            pair_pvals.(pair_name) = p_val;
        else
            pair_pvals.(pair_name) = NaN;
            pair_null_means.(pair_name) = NaN;
            pair_null_stds.(pair_name) = NaN;
            pair_N_eff.(pair_name) = 0;
            skip_degenerate = skip_degenerate + 1;
        end
        % Save after each pair so interrupt resumes at next pair
        save_partial_pair_results(save_path, pair_pvals, pair_scores, pair_null_means, pair_null_stds, ...
            pair_ufo_counts, pair_ufo_counts_raw, pair_ufo_counts_eff, pair_N_eff, ...
            skip_no_usable_events, skip_out_of_bounds, skip_degenerate, skip_other, N, desc_final);
    end  % End of for p = 1:numel(unique_pairs) loop (null sampling)
    end  % End of if ~skip_null_sampling block
    
    % When we loaded from cache we don't have skip_* from this run; load from file if present
    if skip_null_sampling && exist(save_path, 'file')
        warning('off', 'MATLAB:load:variableNotFound');
        S_skip = load(save_path, 'skip_no_usable_events', 'skip_out_of_bounds', 'skip_degenerate', 'skip_other');
        warning('on', 'MATLAB:load:variableNotFound');
        if isfield(S_skip, 'skip_no_usable_events'), skip_no_usable_events = S_skip.skip_no_usable_events; end
        if isfield(S_skip, 'skip_out_of_bounds'),   skip_out_of_bounds   = S_skip.skip_out_of_bounds; end
        if isfield(S_skip, 'skip_degenerate'),      skip_degenerate      = S_skip.skip_degenerate; end
        if isfield(S_skip, 'skip_other'),           skip_other           = S_skip.skip_other; end
    end
    
    % Calculate summary statistics per pair
    n_significant = 0;
    n_total_pairs = 0;
    pair_pval_list = [];
    
    for p = 1:numel(unique_pairs)
        pair_name = matlab.lang.makeValidName(unique_pairs{p});
        if isfield(pair_pvals, pair_name) && isfinite(pair_pvals.(pair_name))
            n_total_pairs = n_total_pairs + 1;
            pair_pval_list(end+1) = pair_pvals.(pair_name);
            if pair_pvals.(pair_name) < 0.05
                n_significant = n_significant + 1;
            end
        end
    end
    
    N_pairs_tested = n_total_pairs;
    N_pairs_testable = N_pairs_attempted - skip_no_usable_events - skip_out_of_bounds;
    N_pairs_skipped = N_pairs_attempted - N_pairs_tested;
    
    % Pair accounting (exact format for audit)
    fprintf('\nPair accounting:\n');
    fprintf('  Candidate pairs (attempted): %d\n', N_pairs_attempted);
    fprintf('  Testable pairs (after filtering): %d\n', max(0, N_pairs_testable));
    fprintf('  Tested pairs (p-values computed): %d\n', N_pairs_tested);
    fprintf('  Skipped pairs: %d\n', N_pairs_skipped);
    fprintf('    - No usable empirical events: %d\n', skip_no_usable_events);
    fprintf('    - Missing / out-of-bounds macro data: %d\n', skip_out_of_bounds);
    fprintf('    - Degenerate empirical score: %d\n', skip_degenerate);
    fprintf('    - Other: %d\n', skip_other);
    
    % FDR: within patient-selection across tested pairs (default and only supported scope)
    fdr_scope = 'within_patient_selection';
    fdr_alpha = 0.05;
    fdr_corrected = [];
    if ~isempty(pair_pval_list) && numel(pair_pval_list) > 1
        fdr_corrected = fdr_bh_local(pair_pval_list, fdr_alpha);
        n_fdr_significant = sum(fdr_corrected < fdr_alpha);
    else
        n_fdr_significant = 0;
    end
    
    % Count NaN scores (degenerate events at UFO level)
    n_nan_scores = sum(isnan(empirical_scores));
    n_valid_scores = n_sel - n_nan_scores;
    
    % Invariant checks
    if N_pairs_attempted < N_pairs_tested
        fprintf('\n*** ERROR: Invariant failed: N_pairs_attempted (%d) < N_pairs_tested (%d) ***\n', N_pairs_attempted, N_pairs_tested);
    end
    if N_pairs_tested ~= numel(pair_pval_list)
        fprintf('\n*** ERROR: Invariant failed: N_pairs_tested (%d) ~= length(p_values) (%d) ***\n', N_pairs_tested, numel(pair_pval_list));
    end
    if ~isempty(fdr_corrected) && (numel(fdr_corrected) ~= numel(pair_pval_list))
        fprintf('\n*** ERROR: Invariant failed: length(q_values) (%d) ~= length(p_values) (%d) ***\n', numel(fdr_corrected), numel(pair_pval_list));
    end
    if n_fdr_significant > n_significant || n_significant > N_pairs_tested
        fprintf('\n*** ERROR: Invariant failed: N_FDR (%d) <= N_Sig (%d) <= N_pairs_tested (%d) required ***\n', n_fdr_significant, n_significant, N_pairs_tested);
    end
    if n_valid_scores + n_nan_scores ~= n_sel
        fprintf('\n*** ERROR: Invariant failed: N_UFOs_valid (%d) + N_UFOs_degenerate (%d) ~= N_UFOs_total (%d) ***\n', n_valid_scores, n_nan_scores, n_sel);
    end
    
    % Print results
    fprintf('\n=== RESULTS: %s ===\n', desc_final);
    fprintf('Configuration: band=[%.0f %.0f] Hz, f0=%.0f Hz, gamma=%.0f (1/s)\n', ...
        cfg.band_hz(1), cfg.band_hz(2), cfg.f0_hz, cfg.gamma);
    fprintf('Windows: pre=%.0f ms, post=%.0f ms\n', cfg.pre_ms, cfg.post_ms);
    fprintf('Monte Carlo: N_null = %d. Min resolvable p ≈ 1/(N_null+1) = %.6f\n', N, 1/(N+1));
    fprintf('Total UFOs: %d (valid scores: %d, degenerate: %d)\n', n_sel, n_valid_scores, n_nan_scores);
    % QC note: valid/missing/degenerate scores (visible at a glance)
    n_missing = n_sel - numel(empirical_scores);
    fprintf('QC: valid scores: %d, missing: %d, degenerate: %d\n', n_valid_scores, n_missing, n_nan_scores);
    fprintf('Note: Statistical tests are per micro-macro pair. N_UFOs is the number of contributing UFO events used to compute pair-level statistics; it is not the number of tests.\n');
    fprintf('Total micro-macro pairs (tested): %d\n', N_pairs_tested);
    fprintf('Significant pairs (p<0.05): %d / %d\n', n_significant, N_pairs_tested);
    fprintf('Multiple-comparisons correction:\n');
    fprintf('  Method: Benjamini-Hochberg (q=%.2f)\n', fdr_alpha);
    fprintf('  Scope: %s (m = %d tested pairs)\n', fdr_scope, N_pairs_tested);
    if ~isempty(fdr_corrected)
        fprintf('FDR-significant pairs (q<0.05): %d / %d\n', n_fdr_significant, N_pairs_tested);
    end
    if N_pairs_tested <= 1
        fprintf('Warning: Very few tested pairs (m = %d). This selection provides limited inferential power for this patient.\n', N_pairs_tested);
    end
    fprintf('K_raw: number of UFO events for this pair before any filtering. K_eff: number of UFO events actually used in the statistic (after exclusions).\n');
    fprintf('\nPer-pair statistics:\n');
    fprintf('%-30s %8s %8s %8s %12s %12s %10s\n', 'Pair', 'K_raw', 'K_eff', 'N_eff', 'Emp_Mean', 'Null_Mean±SD', 'p-value');
    fprintf('%s\n', repmat('-', 1, 100));
    for p = 1:numel(unique_pairs)
        pair_name = matlab.lang.makeValidName(unique_pairs{p});
        if isfield(pair_pvals, pair_name) && isfinite(pair_pvals.(pair_name))
            if isfield(pair_ufo_counts_raw, pair_name)
                k_raw = pair_ufo_counts_raw.(pair_name);
            else
                k_raw = NaN;
            end
            if isfield(pair_ufo_counts_eff, pair_name)
                k_eff = pair_ufo_counts_eff.(pair_name);
            else
                k_eff = NaN;
            end
            if isfield(pair_N_eff, pair_name)
                n_eff = pair_N_eff.(pair_name);
            else
                n_eff = 0;
            end
            emp_mean = pair_scores.(pair_name);
            null_mean = pair_null_means.(pair_name);
            null_std = pair_null_stds.(pair_name);
            p_val = pair_pvals.(pair_name);
            sig_marker = '';
            if p_val < 0.05, sig_marker = '*'; end
            warning_marker = '';
            if n_eff < 0.8 * N
                warning_marker = ' [low N_eff]';
            end
            fprintf('%-30s %8d %8d %8d %12.4f %12.4f±%.4f %10.4f%s%s\n', ...
                unique_pairs{p}, k_raw, k_eff, n_eff, emp_mean, null_mean, null_std, p_val, sig_marker, warning_marker);
        end
    end
    
    % Save results (keep null_scores_* so future runs can extend e.g. 1k -> 10k)
    % Load existing file to preserve checkpoint data structure
    if exist(save_path, 'file')
        S_existing = load(save_path);
        % Do not remove null_scores_* — they allow extending N on a later run
        % Merge with new results
        S_existing.empirical_scores = empirical_scores;
        S_existing.emp_avg = emp_avg;
        S_existing.desc_final = desc_final;
        S_existing.n_sel = n_sel;
        S_existing.n_valid_scores = n_valid_scores;
        S_existing.n_nan_scores = n_nan_scores;
        S_existing.N = N;
        S_existing.pair_scores = pair_scores;
        S_existing.pair_pvals = pair_pvals;
        S_existing.pair_null_means = pair_null_means;
        S_existing.pair_null_stds = pair_null_stds;
        S_existing.pair_ufo_counts = pair_ufo_counts;
        S_existing.pair_ufo_counts_raw = pair_ufo_counts_raw;
        S_existing.pair_ufo_counts_eff = pair_ufo_counts_eff;
        S_existing.pair_N_eff = pair_N_eff;
        S_existing.n_significant = n_significant;
        S_existing.n_total_pairs = n_total_pairs;
        S_existing.fdr_corrected = fdr_corrected;
        S_existing.n_fdr_significant = n_fdr_significant;
        S_existing.skip_no_usable_events = skip_no_usable_events;
        S_existing.skip_out_of_bounds = skip_out_of_bounds;
        S_existing.skip_degenerate = skip_degenerate;
        S_existing.skip_other = skip_other;
        S_existing.ANALYSIS_CONFIG = cfg;
        save(save_path, '-struct', 'S_existing', '-v7.3');
    else
        % Create new file
        save(save_path, 'empirical_scores', 'emp_avg', 'desc_final', 'n_sel', 'n_valid_scores', 'n_nan_scores', 'N', ...
            'pair_scores', 'pair_pvals', 'pair_null_means', 'pair_null_stds', 'pair_ufo_counts', ...
            'pair_ufo_counts_raw', 'pair_ufo_counts_eff', 'pair_N_eff', ...
            'n_significant', 'n_total_pairs', 'fdr_corrected', 'n_fdr_significant', ...
            'skip_no_usable_events', 'skip_out_of_bounds', 'skip_degenerate', 'skip_other', ...
            'ANALYSIS_CONFIG', '-v7.3');
    end
    
    % Store summary: patient, selection, n_ufos, n_pairs, n_sig, n_fdr_sig
    summary_results{end+1} = {patient, desc_final, n_sel, n_total_pairs, n_significant, n_fdr_significant};
end

end % end main

function score = greens_envelope_score(x, fs, cfg)
% GREENS_ENVELOPE_SCORE  Frequency-domain Green's envelope score (damped oscillator template)
%
% Implements Eq. 14 style damped oscillator template:
%   T(f) = 1 ./ ((w0^2 - w.^2).^2 + (2*gamma*w).^2)
%
% Inputs:
%   x - signal vector [baseline | post-UFO]
%   fs - sampling frequency (Hz)
%   cfg - config struct with: pre_ms, post_ms, band_hz, f0_hz, gamma, nperseg_ms, overlap, nfft_min, detrend
%
% Output:
%   score - log10(A_post / A_base) where A is the gain fit to the template

x = double(x(:));

if any(isnan(x)) || all(x == 0)
    score = NaN;
    return;
end

pre_samples  = round(cfg.pre_ms  * fs / 1000);
post_samples = round(cfg.post_ms * fs / 1000);

if numel(x) < pre_samples + post_samples
    score = NaN;
    return;
end

baseline = x(1:pre_samples);
segment = x(pre_samples+1:pre_samples+post_samples);

if all(baseline == 0) || all(segment == 0)
    score = NaN;
    return;
end

% Preprocessing: demean + detrend
baseline = baseline - mean(baseline);
segment = segment - mean(segment);
if cfg.detrend
    if numel(baseline) > 1, baseline = detrend(baseline); end
    if numel(segment) > 1,  segment  = detrend(segment);  end
end

% PSD estimation with defensible parameters
nperseg = max(round(cfg.nperseg_ms * fs / 1000), 64);  % At least 64 samples
noverlap = round(cfg.overlap * nperseg);
nfft = max(cfg.nfft_min, 2^nextpow2(nperseg));  % At least nfft_min, or next power of 2

[P_post, f] = pwelch(segment, nperseg, noverlap, nfft, fs);
[P_base, ~] = pwelch(baseline, nperseg, noverlap, nfft, fs);

% Restrict to analysis band
mask = (f >= cfg.band_hz(1) & f <= cfg.band_hz(2));
if ~any(mask)
    score = NaN;
    return;
end

% Build damped oscillator template (Eq. 14 style)
w = 2*pi*f(mask);           % Frequency in rad/s
w0 = 2*pi*cfg.f0_hz;        % Natural frequency in rad/s
gamma = cfg.gamma;          % Damping (1/s)

% Template: T(f) = 1 ./ ((w0^2 - w.^2).^2 + (2*gamma*w).^2)
T = 1 ./ ((w0^2 - w.^2).^2 + (2*gamma*w).^2);
T = real(T);  % Ensure real (should be, but be safe)
T(~isfinite(T)) = 0;

% Normalize template to unit energy
T_norm = sqrt(sum(T.^2));
if T_norm > eps
    T = T / T_norm;
else
    score = NaN;
    return;
end

% Fit gain via projection (nonnegative)
% A = max(0, dot(T, P) / dot(T, T))
% Since T is normalized, dot(T, T) = 1
A_post = max(0, T' * P_post(mask));
A_base = max(0, T' * P_base(mask));

% Score = log10(A_post / A_base) with safe eps handling
% Return NaN for degenerate amplitudes (excluded by omitnan in means)
if ~(isfinite(A_base) && isfinite(A_post)) || A_base <= 0 || A_post <= 0
    score = NaN;
    return;
end
score = log10((A_post + eps) / (A_base + eps));
score = real(score);

end % end greens_envelope_score

function save_partial_pair_results(sp, pp, ps, pnm, pns, puc, pucr, puce, pne, snue, soob, sd, so, N_val, desc)
% SAVE_PARTIAL_PAIR_RESULTS  Persist pair results and skip counts so interrupt resumes at next pair.
    if ~exist(sp, 'file'), return; end
    try
        S = load(sp);
        S.pair_pvals = pp;
        S.pair_scores = ps;
        S.pair_null_means = pnm;
        S.pair_null_stds = pns;
        S.pair_ufo_counts = puc;
        S.pair_ufo_counts_raw = pucr;
        S.pair_ufo_counts_eff = puce;
        S.pair_N_eff = pne;
        S.skip_no_usable_events = snue;
        S.skip_out_of_bounds = soob;
        S.skip_degenerate = sd;
        S.skip_other = so;
        S.N = N_val;
        S.desc_final = desc;
        save(sp, '-struct', 'S', '-v7.3');
    catch ME
        fprintf('        WARNING: Failed to save partial pair results: %s\n', ME.message);
    end
end

function q = fdr_bh_local(p, alpha)
% FDR_BH_LOCAL  Benjamini-Hochberg FDR correction
% Simple implementation of BH FDR procedure
if isempty(p)
    q = [];
    return;
end

p = p(:);
m = numel(p);
[~, idx] = sort(p);
q = zeros(size(p));

for i = m:-1:1
    q(idx(i)) = min(1, p(idx(i)) * m / i);
    if i < m
        q(idx(i)) = min(q(idx(i)), q(idx(i+1)));
    end
end

end % end fdr_bh_local

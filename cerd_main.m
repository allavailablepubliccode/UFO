function cerd_main
% ========================= CERD MAIN (THEORY ONLY) =========================
% Theory-driven matched-filter detector for cross-scale UFO echoes
%
% Implements cross-scale wave propagation theory from Notes.txt:
%   - UFOs occur at microscale (2-8 kHz) and propagate to macroscale
%   - Macroscale echo described by Green's function (Eq. 20):
%     G = e^(-γt)/v * J₀[v√((ω²-γ²)(t²-Δu²))] * H(t-Δu)
%   - Echo frequency: f_echo = √(ω²-γ²)/(2π) where ω = 2πf0
%   - Echo arrives with delay τ = Δu/v (cross-scale propagation time)
%
% The detector searches for Bessel-shaped, exponentially-decaying echoes
% in macroelectrode signals following microscale UFO events.
%
% Options:
%   CERD_THEORY_ENABLE   (logical)   default: true
%   CERD_THEORY_F0_HZ    [fL fH]     default: [60 90]     % macro natural frequency range (Hz)
%   CERD_THEORY_DELAY_MS [dL dH]     default: [2 10]      % arrival delay search window (ms)
%   CERD_THEORY_GAMMA    (1/s)       default: 80          % damping rate; tau≈1/gamma (s)
%   CERD_THEORY_NFREQ    (int)       default: 7           % grid points over f0
%   CERD_THEORY_NDELAY   (int)       default: 9           % grid points over delay
%
% Metric:
%   THEORY : max normalized correlation between macro snippet and
%            theoretical Green's function template over [delay × f0 × gamma] grid.
%
% Author: drop-in for EDF pipeline
% ====================================================================================

clearvars -except CERD_*; %#ok<CLVAR>

warning('off','all'); addpath('/Users/erik/matmef');

% --------- Top-level switches (existing) ----------
mode    = lower(string(getin('CERD_MODE','g')));   % 'g' | 'a' | 's'
patient = getin('CERD_PATIENT', 4);


% --------- Core knobs ----------
post_ufo_ms     = getin('CERD_POST_UFO_MS', 300);
gap_ufo_ms      = getin('CERD_GAP_UFO_MS',  [10 200]);
Nnull           = getin('CERD_NNULL',       10000);
alpha           = getin('CERD_ALPHA',       0.05);
patient         = getin('CERD_PATIENT',     4);
local_win_sec   = getin('CERD_LOCAL_WIN_SEC',10800);
max_redraw      = getin('CERD_MAX_REDRAW',  50);
min_ufo_freq_hz = getin('CERD_MIN_UFO_FREQ_HZ',[]);
min_ufo_dur_ms  = getin('CERD_MIN_UFO_DUR_MS',[]);
deterministic   = getin('CERD_DETERMINISTIC',false);
base_seed       = getin('CERD_BASE_SEED',   54833);

% --------- NEW: Theory-detector knobs ----------
theory_enable   = getin('CERD_THEORY_ENABLE',   true);
theory_f0_hz    = getin('CERD_THEORY_F0_HZ',    [60 90]);
theory_delay_ms = getin('CERD_THEORY_DELAY_MS', [2 10]);
theory_gamma    = getin('CERD_THEORY_GAMMA',    80);      % 1/s  (~12.5 ms)
theory_nfreq    = getin('CERD_THEORY_NFREQ',    7);
theory_ndelay   = getin('CERD_THEORY_NDELAY',   9);
% Note: theory_template_mode was removed - it was a legacy remnant from forward/backward template modes
% The template is always generated using the forward (causal) Green's function
metric_mode       = "greens";

spect_cfg = struct();
spect_cfg.pre_ms         = double(getin('CERD_SPECT_PRE_MS', 10));
spect_cfg.post_ms        = double(getin('CERD_SPECT_POST_MS', max(5, post_ufo_ms - spect_cfg.pre_ms)));
spect_cfg.band_hz        = getin('CERD_SPECT_BAND_HZ', []);
spect_cfg.bandwidth_hz   = double(getin('CERD_SPECT_BANDWIDTH_HZ', 600));
spect_cfg.window_ms      = double(getin('CERD_SPECT_WINDOW_MS', 5));
spect_cfg.overlap        = double(getin('CERD_SPECT_OVERLAP', 0.5));
spect_cfg.nfft           = double(getin('CERD_SPECT_NFFT', 512));
spect_cfg.eps            = double(getin('CERD_SPECT_EPS', 1e-12));
spect_cfg.detrend        = logical(getin('CERD_SPECT_DETREND', true));
spect_cfg.min_samples    = double(getin('CERD_SPECT_MIN_SAMPLES', 64));
spect_cfg.use_log_ratio  = logical(getin('CERD_SPECT_USE_LOG_RATIO', true));
spect_cfg.fallback_band  = double(getin('CERD_SPECT_FALLBACK_BAND', [500 3500]));
spect_cfg.max_freq_hz    = double(getin('CERD_SPECT_MAX_FREQ_HZ', 6000));

if spect_cfg.pre_ms < 0,  spect_cfg.pre_ms = 0; end
if spect_cfg.post_ms <= 0, spect_cfg.post_ms = max(5, post_ufo_ms - spect_cfg.pre_ms); end
spect_cfg.total_ms = spect_cfg.pre_ms + spect_cfg.post_ms;
% For greens mode, ensure spectral config is properly set
if ~strcmp(metric_mode,'greens')
    error('Only ''greens'' metric mode is supported.');
end
event_permute       = logical(getin('CERD_EVENT_PERMUTE', false));


% ——— Prep & UFOs (with optional caps) –––––
Dat    = prep_basics(patient);

mima_use = logical(Dat.mima);
if ~any(mima_use,'all')
    warning('Dat.mima is empty for %s – no micro ↔ macro assignments.', Dat.ID);
end

% --- Load manually edited UFOs (REQUIRED) ---
% Only use UFOs that have been visually inspected and curated
UFOdat = load_edited_ufos(Dat.ID);

if isempty(UFOdat)
    error('cerd_main: No edited UFOs found for %s. Please run edit_ufos_interactive(''%s'') first to review and edit UFOs.', Dat.ID, Dat.ID);
end

fprintf('[ufo-load] Using manually edited UFOs for %s\n', Dat.ID);

% No capping: use all inspected UFOs as-is
u_counts = cellfun(@(u) size(u,1), UFOdat);
fprintf('[ufo-load] kept_total=%d | per-micro min/med/max=%d/%d/%d\n', ...
    sum(u_counts), min([u_counts,0]), median([u_counts,0]), max([u_counts,0]));

if ~any(u_counts)
    error('cerd_main: No UFO events found in edited UFOs for %s. Please review UFOs using edit_ufos_interactive.', Dat.ID);
end

if event_permute
    UFOdat = permute_ufo_events(Dat, UFOdat, post_ufo_ms, deterministic, base_seed);
    fprintf('[control] Permuted UFO event times for patient %s.\n', Dat.ID);
end

Exclude = exclusion_zones(Dat, int64(gap_ufo_ms*1e3), mima_use);

Results_dir = fullfile(Dat.basedir,'Results',Dat.ID); Results_dir = char(Results_dir);
if ~exist(Results_dir,'dir'), mkdir(Results_dir); end
% Simple analysis tag: just empirical vs permute (parameters are fixed/settled)
analysis_tag = ternary(event_permute, 'permute', 'empirical');
% Template mode tag (for compatibility checking - always 'greens' now)
mode_tag = 'greens';

% Check for null directory for this analysis configuration
Null_dir = fullfile(Results_dir,'Nulls', analysis_tag);
null_dir_found = false;

% Check if new directory exists and has null files
if exist(Null_dir,'dir')
    files = dir(fullfile(Null_dir,'Null_*.mat'));
    if ~isempty(files)
        null_dir_found = true;
    end
end

% Create directory if it doesn't exist (for generating new nulls)
if ~exist(Null_dir,'dir'), mkdir(Null_dir); end

mode = lower(string(getin('CERD_MODE','g')));

% ====================================================================================
% MODE: GENERATE NULLS
% ====================================================================================
if mode == "g"

    % — pool management (reuse across subjects) —
    manage_pool = getin('CERD_MANAGE_POOL', true);
    if manage_pool
        p = gcp('nocreate');
        want = max(1, feature('numcores')-2);
        if isempty(p)
            parpool('local', want);
        elseif p.NumWorkers ~= want
            delete(p);              % only if you need a different size
            parpool('local', want);
        end

    end
    fprintf('Using %d parallel workers.\n', gcp().NumWorkers);

    null_range = parse_null_range(getin('CERD_NULL_RANGE',[]), Nnull);
    if isempty(null_range)
        rngStr = input(sprintf('Enter null index range 1..%d [default=1:%d]: ', Nnull, Nnull), 's');
        if isempty(rngStr), rngStr = sprintf('1:%d', Nnull); end
        null_range = parse_null_range(rngStr, Nnull);
    end
    fprintf('Mode: GENERATE — populating Null_%s in %s\n', range_to_str(null_range), Null_dir);

    all_idx   = null_range;
    existmask = arrayfun(@(k) exist(char(fullfile(Null_dir, sprintf('Null_%d.mat',double(k)))), 'file')==2, all_idx);
    todo      = double(all_idx(~existmask));
    if isempty(todo)
        fprintf('Nothing to do: all requested null files already exist.\n'); return;
    end

    % Broadcast immutable pieces to workers
    Dat_meta       = Dat.meta_path; Dat_t0 = Dat.t0; Dat_t1 = Dat.t1; Dat_password = Dat.password;
    Dat_fs         = Dat.fs; Dat_ma = Dat.ma; Dat_mima = mima_use;
    
    % Save prep_basics.mat for visualization scripts
    prep_file = fullfile(Results_dir, 'prep_basics.mat');
    try
        save(prep_file, 'Dat', 'UFOdat', '-v7.3');
        fprintf('[prep_basics] Saved Dat and UFOdat to %s\n', prep_file);
    catch ME_prep
        warning('Could not save prep_basics.mat: %s', ME_prep.message);
    end
    UFOdat_local   = UFOdat; Exclude_local = Exclude;

    % Save paths & seeds
    nJobs = numel(todo); save_names = cell(1,nJobs); seeds = zeros(1,nJobs);
    for q=1:nJobs
        ii_d = todo(q);
        save_names{q} = char(fullfile(Null_dir, sprintf('Null_%d.mat', double(ii_d))));
        seeds(q)      = deterministic * (double(base_seed) + double(ii_d));
    end
    assert(exist(Dat_meta,'dir')==7, 'Session directory not found: %s', Dat_meta);

    % -------- parfor null generation --------
    parfor q=1:nJobs
        idx = todo(q); save_name = save_names{q}; seed_q = seeds(q);
        if exist(save_name,'file'), continue; end
        if seed_q ~= 0, rng(seed_q,'twister'); else, rng('shuffle'); end

        [theoryN, r0_log, r1_log, nan_rejects] = ...
            null_echoes_local_run( ...
            Dat_meta, Dat_t0, Dat_t1, Dat_password, Dat_fs, Dat_ma, Dat_mima, ...
            UFOdat_local, Exclude_local, post_ufo_ms, ...
            local_win_sec, max_redraw, deterministic, base_seed, idx, ...
            theory_enable, theory_f0_hz, theory_delay_ms, theory_gamma, theory_nfreq, theory_ndelay, ...
            metric_mode, spect_cfg);

        rng_info = getGlobalStreamInfo();
        params = struct('idx',double(idx), 'post_ufo_ms',post_ufo_ms, ...
            'local_win_sec',local_win_sec,'max_redraw',max_redraw, ...
            'theory_enable',logical(theory_enable), 'theory_f0_hz',theory_f0_hz, ...
            'theory_delay_ms',theory_delay_ms, 'theory_gamma',theory_gamma, ...
            'event_permute', logical(event_permute), ...
            'theory_metric', char(metric_mode), 'spect_cfg', spect_cfg);

        S = struct(); S.theoryN=theoryN; S.r0_log=r0_log; S.r1_log=r1_log; S.nan_rejects=nan_rejects;
        S.params=params; S.rng_info=rng_info;
        parsave(save_name,S);
    end

    totalNow = numel(dir(fullfile(Null_dir,'Null_*.mat')));
    fprintf('Done. Generated %d new null files (total present now: %d).\n', numel(todo), totalNow);
    return;
end

% ====================================================================================
% MODE: SCORE-ONLY (build score file if missing)
% ====================================================================================
if mode == "s"
    scorefile = fullfile(Results_dir, sprintf('score_%s.mat', analysis_tag));
    if exist(scorefile,'file')
        fprintf('✅ score file already exists (%s) — skipping generation.\n', scorefile);
        return;
    end
    fprintf('Generating score file only (no nulls)...\n');
    theory = empirical_echoes( ...
        Dat, UFOdat, mima_use, post_ufo_ms, Dat.fs, theory_enable, theory_f0_hz, theory_delay_ms, ...
        theory_gamma, theory_nfreq, theory_ndelay, metric_mode, spect_cfg);
    template_meta = struct('template_mode', mode_tag, ...
        'event_permute', logical(event_permute), ...
        'metric_mode', metric_mode);
    save(scorefile,'theory','template_meta','-v7.3');
    fprintf('✅ score file saved to %s\n', scorefile);
    
    % Save prep_basics.mat for visualization scripts
    prep_file = fullfile(Results_dir, 'prep_basics.mat');
    try
        save(prep_file, 'Dat', 'UFOdat', '-v7.3');
        fprintf('[prep_basics] Saved Dat and UFOdat to %s\n', prep_file);
    catch ME_prep
        warning('Could not save prep_basics.mat: %s', ME_prep.message);
    end
    return;
end

% Disable parallel pool for analyze mode (purely serial)
% But only if we're managing the pool (not when called from batch)
manage_pool = getin('CERD_MANAGE_POOL', true);
if mode == "a" && manage_pool
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);               % close if already running
    end
end

% ====================================================================================
% MODE: ANALYZE
% ====================================================================================
if mode == "a"
    doMean = true;
    fprintf('Analyze: averaging across events (doMean=1)\n');
    
    % Save prep_basics.mat for visualization scripts (ensure it exists)
    prep_file = fullfile(Results_dir, 'prep_basics.mat');
    if ~exist(prep_file, 'file')
        try
            save(prep_file, 'Dat', 'UFOdat', '-v7.3');
            fprintf('[prep_basics] Saved Dat and UFOdat to %s\n', prep_file);
        catch ME_prep
            warning('Could not save prep_basics.mat: %s', ME_prep.message);
        end
    end

    files = dir(fullfile(Null_dir,'Null_*.mat')); if isempty(files), error('No null files in %s. Run GENERATE first.', Null_dir); end
    % sort by index
    nums = cellfun(@(s) sscanf(s,'Null_%d.mat'), {files.name}); [~,ord] = sort(nums); files = files(ord);
    files = files(1:min(Nnull, numel(files)));

    % Probe dimensions
    probe_idx = []; for k=1:numel(files), T = load(fullfile(Null_dir, files(k).name), 'theoryN'); if isfield(T,'theoryN') && ~isempty(T.theoryN), probe_idx = k; break; end, end
    if isempty(probe_idx), error('Found null files but no theoryN arrays.'); end
    S0 = load(fullfile(Null_dir, files(probe_idx).name), 'theoryN');
    sz0 = pad_to_ndims(size(S0.theoryN),3);
    M = sz0(1); N = sz0(2); P = sz0(3);

    % Collect nulls (THEORY only)
    Th_cells={};
    template_warned = false;
    for k=1:numel(files)
        S = load(fullfile(Null_dir, files(k).name));
        template_ok = true;
        if isfield(S,'params')
            params_local = S.params;
            % Core control checks
            if isfield(params_local,'event_permute')
                if logical(params_local.event_permute) ~= event_permute, template_ok = false; end
            end
            if isfield(params_local,'theory_metric')
                if ~strcmpi(char(params_local.theory_metric), char(metric_mode)), template_ok = false; end
            else
                template_ok = false;  % Must match metric mode
            end
            % For greens mode, enforce parameter match
            if strcmpi(char(metric_mode),'greens')
                if isfield(params_local,'theory_f0_hz')
                    if ~isequal(params_local.theory_f0_hz, theory_f0_hz), template_ok = false; end
                end
                if isfield(params_local,'theory_gamma')
                    if ~isequal(double(params_local.theory_gamma), double(theory_gamma)), template_ok = false; end
                end
                if isfield(params_local,'spect_cfg')
                    pl = params_local.spect_cfg;
                    if isfield(pl,'post_ms')      && ~isequal(pl.post_ms,        spect_cfg.post_ms),        template_ok = false; end
                    if isfield(pl,'bandwidth_hz') && ~isequal(pl.bandwidth_hz,   spect_cfg.bandwidth_hz),   template_ok = false; end
                    if isfield(pl,'window_ms')    && ~isequal(pl.window_ms,      spect_cfg.window_ms),      template_ok = false; end
                    if isfield(pl,'nfft')         && ~isequal(pl.nfft,           spect_cfg.nfft),           template_ok = false; end
                    if isfield(pl,'band_hz') && ~isempty(pl.band_hz) && ~isempty(spect_cfg.band_hz)
                        if ~isequal(pl.band_hz, spect_cfg.band_hz), template_ok = false; end
                    end
                end
            end
        end
        if ~template_ok
            if ~template_warned
                fprintf('[warn] Skipping null %s due to control mismatch (expected tag %s).\n', ...
                    files(k).name, analysis_tag);
                template_warned = true;
            end
            continue;
        end
        if isfield(S,'theoryN') && ~isempty(S.theoryN)
            if isequal(pad_to_ndims(size(S.theoryN),3), [M N P])
                Th_cells{end+1} = reshape(S.theoryN, [M N P]); %#ok<AGROW>
            end
        end
    end
    if isempty(Th_cells), error('No usable theory nulls loaded.'); end

    theorynull = cat(4,Th_cells{:}); fprintf('Loaded %d/%d null replicates for THEORY.\n', size(theorynull,4), Nnull);

    % Load/create empirical
    scorefile = fullfile(Results_dir, sprintf('score_%s.mat', analysis_tag)); recompute_score=false;
    fprintf('[score] Looking for file: %s\n', scorefile);
    if exist(scorefile,'file')
        try
            S = load(scorefile);
            if ~isfield(S,'theory'), fprintf('score file missing theory field — will recompute.\n'); recompute_score=true; end
            if isfield(S,'template_meta')
                tm = S.template_meta;
                if ~strcmpi(tm.template_mode, mode_tag)
                    fprintf('score file template mode mismatch (%s vs %s) — will recompute.\n', ...
                        tm.template_mode, mode_tag);
                    recompute_score=true;
                end
                if ~isfield(tm,'event_permute') || logical(tm.event_permute) ~= event_permute
                    old_perm = -1;
                    if isfield(tm,'event_permute'), old_perm = logical(tm.event_permute); end
                    fprintf('score file event permutation mismatch (%d vs %d) — will recompute.\n', ...
                        old_perm, event_permute);
                    recompute_score=true;
                end
                old_metric = 'bessel';
                if isfield(tm,'metric_mode')
                    old_metric = tm.metric_mode;
                end
                if ~strcmpi(old_metric, metric_mode)
                    fprintf('score file metric mismatch (%s vs %s) — will recompute.\n', ...
                        old_metric, metric_mode);
                    recompute_score=true;
                end
            else
                fprintf('score file lacks template meta — will recompute.\n'); recompute_score=true;
            end
        catch, fprintf('score file unreadable — will recompute.\n'); recompute_score=true;
        end
    else
        fprintf('score file not found — will compute from scratch (%s).\n', scorefile); recompute_score=true;
    end

    if recompute_score
    theory = empirical_echoes( ...
        Dat, UFOdat, mima_use, post_ufo_ms, Dat.fs, theory_enable, theory_f0_hz, theory_delay_ms, ...
        theory_gamma, theory_nfreq, theory_ndelay, metric_mode, spect_cfg);
        % Diagnostic: check how many scores are valid
        n_valid = nnz(isfinite(theory));
        n_total = numel(theory);
        fprintf('[score] Computed %d/%d valid scores (%.1f%%)\n', n_valid, n_total, 100*n_valid/max(1,n_total));
        template_meta = struct('template_mode', mode_tag, ...
            'event_permute', logical(event_permute), ...
            'metric_mode', metric_mode);
        save(scorefile,'theory','template_meta','-v7.3');
    else
        S = load(scorefile);
        if isfield(S,'theory')
            theory = S.theory;
            % Diagnostic: check how many scores are valid in loaded file
            n_valid = nnz(isfinite(theory));
            n_total = numel(theory);
            fprintf('[score] Loaded existing file: %d/%d valid scores (%.1f%%)\n', n_valid, n_total, 100*n_valid/max(1,n_total));
        else
            error('Score file %s does not contain theory field. Please regenerate.', scorefile);
        end
    end

    % Rebuild drop mask using current filters
    [drop_cells, kept_cnt, total_cnt] = build_ufo_drop_mask(UFOdat, min_ufo_freq_hz, min_ufo_dur_ms);
    if size(theory,3) ~= P
        error(['Empirical P=%d does not match null P=%d.\n' ...
            'Delete %s and rerun (or regenerate nulls).'], size(theory,3), P, scorefile);
    end

    fprintf('[filter] kept %d / %d UFOs after freq/dur gate (min_f=%.1f Hz, min_d=%.1f ms)\n', ...
        kept_cnt, total_cnt, double(firstNonEmpty(min_ufo_freq_hz,NaN)), double(firstNonEmpty(min_ufo_dur_ms,NaN)));

    % Apply drop mask to empirical & nulls
    n_valid_before_drop = nnz(isfinite(theory));
    drop3 = false(M,N,P);
    for ii=1:M
        di = drop_cells{ii}; if isempty(di), continue; end
        nThis = min(P, numel(di));
        drop3(ii,:,1:nThis) = repmat(reshape(di(1:nThis),[1 1 nThis]),[1 N 1]);
    end
    theory(drop3)=NaN;
    n_valid_after_drop = nnz(isfinite(theory));
    if n_valid_before_drop > 0 && n_valid_after_drop == 0
        fprintf('[warn] Drop mask removed all %d valid scores!\n', n_valid_before_drop);
    elseif n_valid_before_drop > n_valid_after_drop
        fprintf('[info] Drop mask: %d -> %d valid scores (removed %d)\n', ...
            n_valid_before_drop, n_valid_after_drop, n_valid_before_drop - n_valid_after_drop);
    end

    Kth = size(theorynull,4);
    drop4_th = repmat(drop3,[1 1 1 Kth]);
    theorynull(drop4_th)=NaN;

    % Average across events (optional)
    if doMean
        theory = nanmean(theory,3);
        theorynull = nanmean(theorynull,3);
        theory = reshape(theory,[M,N,1]);
        theorynull = reshape(theorynull,[M,N,1,size(theorynull,4)]);
        P = 1;
    end

    % Connectivity mask
    connMN = mima_use'; mask3 = repmat(~connMN, [1,1,P]);
    mask4  = repmat(~connMN, [1,1,P, Kth]);
    
    n_valid_before_conn = nnz(isfinite(theory));
    theory(mask3)=NaN;
    theorynull(mask4)=NaN;
    n_valid_after_conn = nnz(isfinite(theory));
    if n_valid_before_conn > 0 && n_valid_after_conn == 0
        fprintf('[warn] Connectivity mask removed all %d valid scores!\n', n_valid_before_conn);
        fprintf('[warn] theory size: [%d %d %d], mima_use size: [%d %d], connMN size: [%d %d]\n', ...
            size(theory,1), size(theory,2), size(theory,3), size(mima_use,1), size(mima_use,2), size(connMN,1), size(connMN,2));
        fprintf('[warn] nnz(mima_use)=%d, nnz(connMN)=%d\n', nnz(mima_use), nnz(connMN));
    elseif n_valid_before_conn > n_valid_after_conn
        fprintf('[info] Connectivity mask: %d -> %d valid scores (removed %d)\n', ...
            n_valid_before_conn, n_valid_after_conn, n_valid_before_conn - n_valid_after_conn);
    end

    % FDR groups & diagnostics (global grouping only)
    topN = 25;
    G = build_fdr_groups(mima_use, size(theory), 'global');

    % Diagnostic: check if nulls have valid values where theory has valid values
    n_valid_theory = nnz(isfinite(theory));
    n_valid_nulls = nnz(isfinite(theorynull));
    Keff_map = sum(isfinite(theorynull), 4);
    valid_emp = isfinite(theory);
    valid_test = valid_emp & (Keff_map > 0);
    n_valid_test = nnz(valid_test);
    if n_valid_theory > 0 && n_valid_test == 0
        fprintf('[warn] No testable voxels: %d valid theory scores, but 0 have valid nulls!\n', n_valid_theory);
        fprintf('[warn] Keff_map: min=%d, max=%d, nnz>0=%d\n', ...
            min(Keff_map(:)), max(Keff_map(:)), nnz(Keff_map > 0));
    end

    if theory_enable && ~isempty(theorynull)
        diagT = write_metric_diagnostics(theory, theorynull, alpha, 'theory [global]', topN, doMean, true, G);
        print_diag(diagT,'THEORY');
    else
        error('Theory analysis required but theorynull is empty.');
    end

    % Extract significant UFO details (get full output from compare_to_nulls_fdr)
    [p_map, sig_mask, crit_p, ~, ~] = compare_to_nulls_fdr(theory, theorynull, alpha, G);
    sig_ufos = extract_significant_ufos(sig_mask, UFOdat, Dat, theory, p_map, mima_use);

    % Save significant UFO details
    sig_file = fullfile(Results_dir, sprintf('significant_ufos_%s.mat', analysis_tag));
    save(sig_file, 'sig_ufos', 'p_map', 'sig_mask', 'crit_p', 'theory', 'theorynull', '-v7.3');
    fprintf('\n[SAVED] Significant UFO details to %s\n', sig_file);
    
    % Optionally visualize significant UFOs
    if getin('CERD_VISUALIZE_SIG', false) && ~isempty(sig_ufos)
        vis_dir = fullfile(Results_dir, 'visualizations');
        
        % Determine which events to visualize
        vis_event_indices = getin('CERD_VIS_EVENT_INDICES', []);  % Indices into sig_ufos array (1-based)
        vis_event_numbers = getin('CERD_VIS_EVENT_NUMBERS', []);  % Actual event_idx values to match
        
        fprintf('[VIS DEBUG] Found %d significant UFOs. vis_event_indices=%s, vis_event_numbers=%s\n', ...
            numel(sig_ufos), mat2str(vis_event_indices), mat2str(vis_event_numbers));
        
        % List available events for user reference
        if numel(sig_ufos) > 0
            fprintf('[VIS DEBUG] Available events:\n');
            for s = 1:min(10, numel(sig_ufos))
                fprintf('  Index %d: micro=%s, macro=%s, event_idx=%d, p=%.4g\n', ...
                    s, sig_ufos(s).micro_channel, sig_ufos(s).macro_name, sig_ufos(s).event_idx, sig_ufos(s).p_value);
            end
            if numel(sig_ufos) > 10
                fprintf('  ... and %d more\n', numel(sig_ufos) - 10);
            end
        end
        
        if ~isempty(vis_event_indices)
            % Use specified indices
            indices_to_plot = vis_event_indices(vis_event_indices >= 1 & vis_event_indices <= numel(sig_ufos));
            if isempty(indices_to_plot)
                warning('No valid event indices specified (range: 1-%d, requested: %s)', ...
                    numel(sig_ufos), mat2str(vis_event_indices));
                indices_to_plot = [];  % Explicitly set to empty
            else
                fprintf('\n[VISUALIZE] Generating plots for %d specified events (indices: %s)...\n', ...
                    numel(indices_to_plot), mat2str(indices_to_plot));
            end
        elseif ~isempty(vis_event_numbers)
            % Match by event_idx
            all_event_numbers = [sig_ufos.event_idx];
            indices_to_plot = [];
            for evt_num = vis_event_numbers(:)'
                match_idx = find(all_event_numbers == evt_num, 1);
                if ~isempty(match_idx)
                    indices_to_plot(end+1) = match_idx;
                else
                    warning('Event number %d not found in significant UFOs', evt_num);
                end
            end
            if isempty(indices_to_plot)
                warning('No matching event numbers found');
            else
                fprintf('\n[VISUALIZE] Generating plots for %d specified events (event numbers: %s)...\n', ...
                    numel(indices_to_plot), mat2str(vis_event_numbers));
            end
        else
            % Default: plot first N events (sorted by p-value)
            max_plots = getin('CERD_VIS_MAX_PLOTS', 50);
            indices_to_plot = 1:min(numel(sig_ufos), max_plots);
            fprintf('\n[VISUALIZE] Generating plots for %d significant UFOs (first %d by p-value)...\n', ...
                numel(indices_to_plot), numel(indices_to_plot));
        end
        
        % Generate plots
        if ~isempty(indices_to_plot)
            for s = indices_to_plot
                fprintf('  Plotting event %d (sig_ufos index %d): micro=%s, macro=%s, event_idx=%d, p=%.4g\n', ...
                    s, s, sig_ufos(s).micro_channel, sig_ufos(s).macro_name, sig_ufos(s).event_idx, sig_ufos(s).p_value);
            visualize_significant_ufo(Dat, sig_ufos(s), post_ufo_ms, vis_dir);
        end
        fprintf('[VISUALIZE] Done. Saved to %s\n', vis_dir);
        else
            fprintf('[VISUALIZE] No events to plot (indices_to_plot is empty)\n');
        end
    end

    % Summary
    if event_permute
        fprintf('[info] UFO event times permuted before analysis.\n');
    end
    fprintf('[info] THEORY metric: %s\n', upper(char(metric_mode)));

    fprintf('\n=== SUMMARY @ FDR q=%.3f (GLOBAL) ===\n', alpha);
    fprintf('THEORY: %d / %d significant (K=%d, crit_p=%.4g)\n', diagT.n_sig, diagT.m, diagT.K, diagT.crit_p);
    fprintf('Significant UFO events: %d\n', numel(sig_ufos));
    
    % Print significant UFO details
    if ~isempty(sig_ufos)
        fprintf('\n=== SIGNIFICANT UFO DETAILS ===\n');
        fprintf('%-8s %-12s %-10s %-12s %-10s %-10s %-10s\n', 'Micro#', 'Macro#', 'Event#', 'Duration(ms)', 'Freq(Hz)', 'Theory', 'p-value');
        fprintf('%s\n', repmat('-',1,80));
        for s=1:min(20, numel(sig_ufos))  % Show first 20
            u = sig_ufos(s);
            fprintf('%-8d %-12s %-10d %-12.2f %-10.1f %-10.4f %-10.4g\n', ...
                u.micro_idx, u.macro_name, u.event_idx, u.duration_ms, u.frequency_hz, u.theory_score, u.p_value);
        end
        if numel(sig_ufos) > 20
            fprintf('... and %d more (see significant_ufos.mat)\n', numel(sig_ufos)-20);
        end
    end

    % Export summary
    S.patient = Dat.ID;
    S.event_permute = logical(event_permute);
    S.metric_mode = char(metric_mode);
    S.analysis_tag = analysis_tag;
    S.rows = struct('name','THEORY','m',diagT.m,'sig',diagT.n_sig, ...
        'min_p',diagT.min_p,'crit_p',diagT.crit_p,'p_floor',diagT.p_floor, ...
        'K',diagT.K,'K_needed_floor',diagT.K_needed_floor,'more',diagT.worth_more_nulls);
    assignin('base','CERD_LAST_SUMMARY', S);
    fprintf('=====================================================\n\n');
    return;
end
end
% =========================== END MAIN DISPATCH ======================================


% =========================== SUPPORTING FUNCTIONS ===================================
function theory = empirical_echoes( ...
    Dat, UFOdat, mima_use, post_ufo_ms, fs, theory_enable, theory_f0_hz, theory_delay_ms, ...
    theory_gamma, theory_nfreq, theory_ndelay, metric_mode, spect_cfg)

disp('...extracting empirical echoes (THEORY only) ...'); drawnow;

fast_meta = getin('CERD_FAST_META', false);
fast_meta_ids = getin('CERD_FAST_META_IDS', {});
if ischar(fast_meta_ids) || isstring(fast_meta_ids)
    fast_meta_ids = cellstr(fast_meta_ids);
end
use_fast = fast_meta || any(strcmpi(fast_meta_ids, 'ALL')) || ...
    any(strcmpi(fast_meta_ids, Dat.ID));

if nargin < 3 || isempty(mima_use)
    mima_use = Dat.mima;
end
mima_use = logical(mima_use);
M = size(mima_use,2); N = size(mima_use,1);
Pmax = max(cellfun(@(u) size(u,1), UFOdat));
theory = NaN(M,N,Pmax);

post_ufo_us = int64(post_ufo_ms * 1e3);
pre_us  = int64(spect_cfg.pre_ms * 1e3);

for ii=1:M
    UFO = UFOdat{ii};
    for jj=1:N
        if ~mima_use(jj,ii), continue; end
        mac = Dat.ma{jj}{1};
        for kk=1:size(UFO,1)

            % For greens mode: start AFTER UFO ends (conservative, avoids contamination)
            % Segment includes baseline (pre_ms) and post-UFO window (post_ms)
                start_from_ufo_beginning = getin('CERD_UFO_START_FROM_BEGINNING', false);
                if start_from_ufo_beginning
                event_start = int64(UFO(kk,1));  % Start when UFO begins
                else
                event_start = int64(UFO(kk,2));  % Start after UFO ends (default)
                end
            seg_start = event_start - pre_us;
            seg_end   = event_start + int64(spect_cfg.post_ms * 1e3);
            if seg_start < int64(Dat.t0) || seg_end > int64(Dat.t1)
                    continue;
            end

            if use_fast
                % FAST META: skip actual readMef3, generate synthetic macro data
                fs_local = Dat.fs;
                macro = synthetic_macro(fs_local, post_ufo_ms);
            else
                [~, macro] = readMef3(Dat.meta_path, Dat.password, mac, 'time', seg_start, seg_end);
                if isempty(macro) || any(isnan(macro),'all'), continue; end
                fs_local = local_fs(fs, macro, seg_start, seg_end);
            end

            % Theory matched filter (greens mode: power spectrum approach)
            if theory_enable
                if ~strcmp(metric_mode,'greens')
                    error('Only ''greens'' metric mode is supported. Set CERD_THEORY_METRIC=''greens''.');
                end
                    score_val = greens_spectrum_score(macro, fs_local, spect_cfg);
                    theory(ii,jj,kk) = score_val;
            end
        end
    end
end
end

function [theoryN, r0_log, r1_log, nan_rejects] = ...
    null_echoes_local_run(Dat_meta, Dat_t0, Dat_t1, Dat_password, Dat_fs, Dat_ma, Dat_mima, ...
    UFOdat, Exclude, post_ufo_ms, ...
    local_win_sec, max_redraw, deterministic, base_seed, idx_seed, ...
    theory_enable, theory_f0_hz, theory_delay_ms, theory_gamma, theory_nfreq, theory_ndelay, ...
    metric_mode, spect_cfg)

fast_meta = getin('CERD_FAST_META', false);
fast_meta_ids = getin('CERD_FAST_META_IDS', {});
if ischar(fast_meta_ids) || isstring(fast_meta_ids)
    fast_meta_ids = cellstr(fast_meta_ids);
end
pat_id = string(getin('CERD_PATIENT',''));
use_fast = fast_meta || any(strcmpi(fast_meta_ids, 'ALL')) || ...
    any(strcmpi(fast_meta_ids, pat_id));

M    = size(Dat_mima,2); N = size(Dat_mima,1);
Pmax = max(cellfun(@(u) size(u,1), UFOdat));
theoryN=NaN(M,N,Pmax);
r0_log=NaN(M,N,Pmax); r1_log=NaN(M,N,Pmax); nan_rejects=zeros(M,N,Pmax);

post_ufo_us=int64(post_ufo_ms*1e3);
pre_us  = int64(spect_cfg.pre_ms * 1e3);
local_half = int64(local_win_sec*1e6);

start_min=int64(Dat_t0); start_max=int64(Dat_t1)-post_ufo_us;
if deterministic, rng(double(base_seed)+double(idx_seed),'twister'); else, rng('shuffle'); end

for ii=1:M
    UFO = UFOdat{ii};
    for jj=1:N
        if ~Dat_mima(jj,ii), continue; end
        mac = Dat_ma{jj}{1}; Exc_use = Exclude{jj};

        for kk=1:size(UFO,1)

            u_end = int64(UFO(kk,2));
            L = max(start_min, u_end - local_half); R = min(start_max, u_end + local_half);
            if R<=L, continue; end
            base_dom = [L R];

            bad = zeros(0,2,'like',start_min);
            if ~isempty(Exc_use)
                tmp = [Exc_use(:,1)-post_ufo_us,  Exc_use(:,2)];
                tmp(:,1)=max(tmp(:,1),start_min); tmp(:,2)=min(tmp(:,2),start_max);
                tmp = tmp(tmp(:,2)>=tmp(:,1),:); bad = merge_intervals(tmp);
            end
            allow_local = intersect_allow_with_domain(complement_intervals([start_min start_max], bad), base_dom);
            if isempty(allow_local), continue; end

            tried=0; ok=false;
            while tried < max_redraw && ~ok
                r0 = sample_from_intervals(allow_local); r1 = r0 + post_ufo_us;
                if use_fast
                    % FAST META mode: skip readMef3, synthesize placeholder waveform
                    macro = synthetic_macro(Dat_fs, double(post_ufo_ms));
                    fs_local = Dat_fs;
                else
                    [~, macro] = readMef3(Dat_meta, Dat_password, mac, 'time', r0, r1);
                    tried = tried + 1;

                    if isempty(macro) || any(isnan(macro),'all')
                        nan_rejects(ii,jj,kk) = nan_rejects(ii,jj,kk)+1;
                        continue;
                    end
                    fs_local = local_fs(Dat_fs, macro, r0, r1);
                end

                if theory_enable
                    if ~strcmp(metric_mode,'greens')
                        error('Only ''greens'' metric mode is supported. Set CERD_THEORY_METRIC=''greens''.');
                    end
                        theoryN(ii,jj,kk) = greens_spectrum_score(macro, fs_local, spect_cfg);
                end

                r0_log(ii,jj,kk)=double(r0); r1_log(ii,jj,kk)=double(r1); ok=true;
            end
        end
    end
end
end



function s = theory_match_score(x, fs, f0_range_hz, delay_range_ms, gamma_1ps, win_ms, nfreq, ndelay)
% THEORY_MATCH_SCORE  Matched-filter detector for cross-scale UFO echoes
%
% Implements theory from Notes.txt Eq. (20): Green's function template matching
%   G = e^(-γt)/v * J₀[v√((ω²-γ²)(t²-Δu²))] * H(t-Δu)
%
% The template uses echo frequency ω_echo = √(ω²-γ²) and searches over:
%   - Natural frequency f0 (ω = 2πf0) in range f0_range_hz
%   - Delay τ (corresponds to Δu from theory) in range delay_range_ms  
%   - Damping γ in grid around gamma_1ps
%
% Note: For light damping (γ << ω), echo frequency ≈ f0, so band-limiting
% around f0±BW is appropriate. For heavier damping, echo frequency is lower.
%
% Processing steps:
% - AR(1) prewhitening (optional, improves SNR)
% - Band-limit around f0±BW to improve SNR
% - Sweep delay, f0, and gamma grid
% - Early-time Tukey gate emphasizes echo onset
% - Normalized correlation for template matching
%
% Keeps the same API; pulls extra knobs from base if present.

% ---- fetch optional knobs ----
BW      = getin('CERD_THEORY_BW_HZ',15);
Grng    = getin('CERD_THEORY_GAMMA_RANGE',[max(20,gamma_1ps/2) min(2000,2*gamma_1ps)]);
Ngamma  = max(1, round(getin('CERD_THEORY_NGAMMA',3)));
gate_ms = getin('CERD_THEORY_GATE_MS',40);
doPW    = getin('CERD_THEORY_PREWHITEN',true);

x = double(x(:));
x = x - mean(x);
if std(x)>0, x = x/std(x); end

L = numel(x);
t = (0:L-1)'/fs;
dL = delay_range_ms(1)/1000; dH = delay_range_ms(2)/1000;
fL = f0_range_hz(1);         fH = f0_range_hz(2);
f_grid = linspace(fL, fH, max(1,nfreq));
d_grid = linspace(dL, dH, max(1,ndelay));
if numel(Grng)==1, gam_grid = Grng; else, gam_grid = linspace(Grng(1), Grng(2), max(1,Ngamma)); end

% Prewhiten once (broadband), then band-limit per f0
if doPW
    xw = prewhiten_ar1(x);
else
    xw = x;
end

best = 0;
for fi = 1:numel(f_grid)
    f0 = f_grid(fi);
    omega = 2*pi*f0;
    
    % Note: Echo frequency from theory is f_echo = √(ω²-γ²)/(2π)
    % For light damping (γ << ω), f_echo ≈ f0, so band-limiting around f0 is appropriate
    % For heavier damping, f_echo is lower, but since we're searching over a grid
    % and the difference is typically small, using f0 for band-limiting is fine
    xl = bandlimit(xw, fs, max(1,f0-BW), f0+BW);

    % normalize post-bandlimit
    xl = xl - mean(xl);
    sd = std(xl); if sd>0, xl = xl/sd; end

    for gi = 1:numel(gam_grid)
        gamma = gam_grid(gi);
        for di = 1:numel(d_grid)
            tau = d_grid(di);

            % build template & gate the early echo
            g = theory_template_gated(t, tau, 2*pi*f0, gamma, gate_ms);
            g = g - mean(g); sg = std(g); if sg>0, g = g/sg; end

            % Normalized correlation coefficient (both xl and g are zero-mean, unit-std)
            % For properly normalized vectors: corr = dot(xl,g) / sqrt(sum(xl^2)*sum(g^2))
            % Since std=1: sum(xl^2) = N-1, so corr = dot(xl,g) / (N-1) [exact]
            % Dividing by N introduces small bias (factor ~(N-1)/N), but acceptable for large N
            % Note: For zero-mean unit-std vectors, max(|corr|) ≤ 1
            N = numel(xl);
            c = (xl.'*g) / max(eps, N - 1);  % Use N-1 for unbiased correlation estimate
            if c > best, best = c; end
        end
    end
end
s = max(0, best);
end

function g = theory_template_gated(t, tau, omega, gamma_1ps, gate_ms)
% THEORY_TEMPLATE_GATED  Generate Green's function template for cross-scale echo
% Based on Eq. (20) from Notes.txt:
%   G = e^(-γt)/v * J₀[v√((ω²-γ²)(t²-Δu²))] * H(t-Δu)
% where Δu = |u - u_micro|/v is the scale distance (tau in our code)
%
% For practical implementation:
%   - Use echo frequency ω_echo = √(ω²-γ²) 
%   - Bessel argument: ω_echo * √(max(0, t²-τ²))
%   - This matches the theoretical form when v is absorbed into the frequency scaling
%
% Inputs:
%   t: time vector (seconds)
%   tau: delay (seconds) = Δu from theory
%   omega: natural frequency (rad/s) = 2πf0
%   gamma_1ps: damping rate (1/s)
%   gate_ms: early-time gate duration (ms)

dt = t - tau;
mask = (dt >= 0);

% Calculate echo frequency from theory: ω_echo = √(ω²-γ²)
omega_sq = omega^2;
gamma_sq = gamma_1ps^2;
omega_echo = sqrt(max(0, omega_sq - gamma_sq));  % Echo frequency (rad/s)

% Build template according to theory: J₀[ω_echo * √(t²-τ²)] for t >= τ
% For t < τ, the argument would be imaginary, so we use max(0, t²-τ²)
t_sq = t.^2;
tau_sq = tau^2;
arg_sq = max(0, t_sq - tau_sq);  % Ensure non-negative for sqrt
bessel_arg = omega_echo * sqrt(arg_sq);

% Exponential decay: e^(-γ(t-τ)) for t >= τ
decay = exp(-gamma_1ps * dt);

% Bessel function with theoretical argument
bessel_val = besselj(0, bessel_arg);

% Combine: decay * Bessel, with causality mask
a = decay .* bessel_val;
g = zeros(size(t));
g(mask) = a(mask);

% Apply normalization factor 1/v from theory (absorbed into overall scaling)
% The 1/v factor is typically handled by normalization in the correlation step

% apply a short early-time gate (Tukey) to emphasize the echo onset
if isfinite(gate_ms) && gate_ms > 0
    Tgate = gate_ms / 1000;
    gate = zeros(size(t));
    on = mask & (dt <= Tgate);
    if any(on)
        u = linspace(0,1, sum(on))';
        alpha = 0.25;
        w = tukeywin(numel(u), alpha);
        gate(on) = w;
    end
    g = g .* gate;
end
end

function score = spectral_power_score(x, fs, cfg, center_freq)
% SPECTRAL_POWER_SCORE  Compute band-limited power change (post vs baseline).
% The macro snippet is assumed to contain [baseline | post-ufo] samples.

if isempty(x)
    score = NaN;
    return;
end

x = double(x);
if size(x,1) <= size(x,2)
    x = x(:);
else
    x = x(:);
end

total_samples = numel(x);
pre_samples  = max(1, round(cfg.pre_ms  * fs / 1000));
post_samples = max(1, round(cfg.post_ms * fs / 1000));

if pre_samples + post_samples > total_samples
    post_samples = total_samples - pre_samples;
end

if pre_samples < cfg.min_samples || post_samples < cfg.min_samples || post_samples <= 0
    score = NaN;
    return;
end

baseline = x(1:pre_samples);
segment  = x(pre_samples+1:pre_samples+post_samples);

if cfg.detrend
    if numel(baseline) > 1, baseline = detrend(baseline); end
    if numel(segment)  > 1, segment  = detrend(segment);  end
end

if numel(segment) < 8 || numel(baseline) < 8
    score = NaN;
    return;
end

win = max(16, round(cfg.window_ms * fs / 1000));
win = min(win, numel(segment));
if win >= numel(segment)
    win = max(8, floor(numel(segment) / 2));
end

if win < 8
    score = NaN;
    return;
end

noverlap = min(win-1, max(0, round(cfg.overlap * win)));
nfft = max(win, round(cfg.nfft));

[P_post,f] = pwelch(segment, win, noverlap, nfft, fs);
[P_base,~] = pwelch(baseline, win, noverlap, nfft, fs);

band = cfg.band_hz;
if isempty(band) || numel(band) < 2 || any(~isfinite(band))
    bw = cfg.bandwidth_hz;
    if ~isfinite(center_freq) || center_freq <= 0
        band = cfg.fallback_band;
    else
        if ~isfinite(bw) || bw <= 0
            bw = max(100, center_freq * 0.5);
        end
        half_bw = bw / 2;
        band = [max(0, center_freq - half_bw), center_freq + half_bw];
    end
end
band(1) = max(0, band(1));
band(2) = min(cfg.max_freq_hz, max(band(1) + 1, band(2)));

mask = f >= band(1) & f <= band(2);
if ~any(mask)
    [~, idx] = min(abs(f - mean(band)));
    mask = idx;
end

pow_post = mean(P_post(mask));
pow_base = mean(P_base(mask));

if cfg.use_log_ratio
    score = log10((pow_post + cfg.eps) / (pow_base + cfg.eps));
else
    score = (pow_post - pow_base) / (pow_base + cfg.eps);
end
score = real(score);
end

function score = power_bump_score(x, fs, cfg, center_freq)
% POWER_BUMP_SCORE  Gamma-free transient power bump detector.
% Looks for power increase in early window (2-30 ms) followed by decrease in late window (30-150 ms).
% Score = early_power / late_power (higher = transient bump).
% The macro snippet is assumed to contain [baseline | post-ufo] samples.

if isempty(x)
    score = NaN;
    return;
end

x = double(x);
if size(x,1) <= size(x,2)
    x = x(:);
else
    x = x(:);
end

% Time windows (ms) - now configurable
bump_onset_ms = getin('CERD_BUMP_ONSET_MS', 2);        % Earliest to look for bump (propagation delay)
bump_peak_end_ms = getin('CERD_BUMP_PEAK_END_MS', 30);     % End of early window (where power should peak)
bump_decay_start_ms = getin('CERD_BUMP_DECAY_START_MS', 30);  % Start of late window (where power should decay)
% Default decay end: min(150, post_ms) to ensure it doesn't exceed available window
bump_decay_end_default = min(150, cfg.post_ms);
bump_decay_end_ms = getin('CERD_BUMP_DECAY_END_MS', bump_decay_end_default);  % End of late window

% Convert to samples
pre_samples = max(1, round(cfg.pre_ms * fs / 1000));
total_post_samples = max(1, round(cfg.post_ms * fs / 1000));
total_samples = numel(x);

if pre_samples + total_post_samples > total_samples
    total_post_samples = total_samples - pre_samples;
end

if pre_samples < cfg.min_samples || total_post_samples < cfg.min_samples || total_post_samples <= 0
    score = NaN;
    return;
end

% Extract segments
baseline = x(1:pre_samples);
post_segment = x(pre_samples+1:pre_samples+total_post_samples);

if cfg.detrend
    if numel(baseline) > 1, baseline = detrend(baseline); end
    if numel(post_segment) > 1, post_segment = detrend(post_segment); end
end

if numel(post_segment) < 8 || numel(baseline) < 8
    score = NaN;
    return;
end

% Convert time windows to sample indices (relative to post_segment start)
bump_onset_samples = max(1, round(bump_onset_ms * fs / 1000));
bump_peak_end_samples = min(numel(post_segment), round(bump_peak_end_ms * fs / 1000));
bump_decay_start_samples = min(numel(post_segment), round(bump_decay_start_ms * fs / 1000));
bump_decay_end_samples = min(numel(post_segment), round(bump_decay_end_ms * fs / 1000));

% Ensure windows are valid
if bump_peak_end_samples < bump_onset_samples
    bump_peak_end_samples = bump_onset_samples + 1;
end
if bump_decay_end_samples <= bump_decay_start_samples
    bump_decay_end_samples = bump_decay_start_samples + 1;
end

early_window = post_segment(bump_onset_samples:bump_peak_end_samples);
late_window = post_segment(bump_decay_start_samples:bump_decay_end_samples);

if numel(early_window) < 4 || numel(late_window) < 4
    score = NaN;
    return;
end

% Compute power in each window using pwelch
win = max(16, round(cfg.window_ms * fs / 1000));
win = min(win, min(numel(early_window), numel(late_window)));
if win >= min(numel(early_window), numel(late_window))
    win = max(8, floor(min(numel(early_window), numel(late_window)) / 2));
end

if win < 8
    score = NaN;
    return;
end

noverlap = min(win-1, max(0, round(cfg.overlap * win)));
nfft = max(win, round(cfg.nfft));

[P_early, f] = pwelch(early_window, win, noverlap, nfft, fs);
[P_late, ~] = pwelch(late_window, win, noverlap, nfft, fs);
[P_base, ~] = pwelch(baseline, win, noverlap, nfft, fs);

% Determine frequency band
band = cfg.band_hz;
if isempty(band) || numel(band) < 2 || any(~isfinite(band))
    % Use conservative macro frequency range [60 90] Hz (original theory range)
    band = [60 90];
    % If UFO frequency is available, could center around f_UFO/10 or similar
    % For now, use fixed [60 90] Hz range
end
band(1) = max(0, band(1));
band(2) = min(cfg.max_freq_hz, max(band(1) + 1, band(2)));

mask = f >= band(1) & f <= band(2);
if ~any(mask)
    [~, idx] = min(abs(f - mean(band)));
    mask = idx;
end

pow_early = mean(P_early(mask));
pow_late = mean(P_late(mask));
pow_base = mean(P_base(mask));

% Score: ratio of early to late power (higher = transient bump)
% Penalize if early power < baseline (no increase)
if pow_early < pow_base
    % No increase detected - return low score
    score = pow_early / (pow_late + cfg.eps);
else
    % Normal case: early > baseline, score = early/late
    score = pow_early / (pow_late + cfg.eps);
end

% Use log ratio for better dynamic range (optional)
if cfg.use_log_ratio
    score = log10(score + cfg.eps);
end

score = real(score);
end

function score = greens_spectrum_score(x, fs, cfg)
% GREENS_SPECTRUM_SCORE  Theory-aligned spectral test using frequency-domain Green's function
%
% Based on Notes.txt Eq. (25): Frequency-domain Green's function
%   G̃(ω') = 1/(2vq) * e^(-qΔu), where q = (1/v)√(ω² - 2iγω' - ω'²)
%
% The power spectrum |G̃(ω')|² has the form of a damped oscillator transfer function:
%   |H(ω')|² = 1 / ((ω_echo² - ω'²)² + (2γω')²)
% where ω_echo = √(ω²-γ²) is the echo frequency from the theory.
%
% This function:
% - Fits |H(ω')|² template to power spectrum in frequency band
% - Compares fitted gain post-UFO vs baseline
% - Returns log10 gain ratio (post/base)
%
% Uses:
%   CERD_THEORY_F0_HZ (macro natural frequency band, default [60 90] Hz)
%   CERD_THEORY_GAMMA (damping, default 80 1/s; we override to user-provided)

if isempty(x)
    score = NaN; return;
end
x = double(x(:));

% Segment sizes
pre_samples  = max(1, round(cfg.pre_ms  * fs / 1000));
post_samples = max(1, round(cfg.post_ms * fs / 1000));
if pre_samples + post_samples > numel(x)
    post_samples = numel(x) - pre_samples;
end
if pre_samples < cfg.min_samples || post_samples < cfg.min_samples || post_samples <= 0
    score = NaN; return;
end

baseline = x(1:pre_samples);
segment  = x(pre_samples+1:pre_samples+post_samples);
if cfg.detrend
    if numel(baseline)>1, baseline = detrend(baseline); end
    if numel(segment) >1, segment  = detrend(segment);  end
end
if numel(segment) < 8 || numel(baseline) < 8
    score = NaN; return;
end

% Welch PSD
win = max(16, round(cfg.window_ms * fs / 1000));
win = min(win, numel(segment));
if win >= numel(segment), win = max(8, floor(numel(segment)/2)); end
if win < 8, score = NaN; return; end
noverlap = min(win-1, max(0, round(cfg.overlap * win)));
nfft = max(win, round(cfg.nfft));
[P_post,f] = pwelch(segment, win, noverlap, nfft, fs);
[P_base,~] = pwelch(baseline, win, noverlap, nfft, fs);

% Band selection
band = cfg.band_hz;
if isempty(band) || numel(band) < 2 || any(~isfinite(band))
    % Default to macro theory band (60–90 Hz) if not specified
    band = getin('CERD_THEORY_F0_HZ', [60 90]);
end
band(1) = max(0, band(1)); band(2) = min(cfg.max_freq_hz, max(band(1)+1, band(2)));
mask = (f >= band(1) & f <= band(2));
if ~any(mask)
    [~,idx] = min(abs(f - mean(band)));
    mask = false(size(f)); mask(idx)=true;
end

% Theoretical shape (normalized) over band
% Based on Notes.txt Eq. (25-26): Power spectrum |G̃(ω')|² has form of damped oscillator
% Theory: Echo frequency ω_echo = √(ω²-γ²), where ω is natural frequency, γ is damping
f0_band = getin('CERD_THEORY_F0_HZ', [60 90]);
f0 = mean(double(f0_band(:)));              % macro natural frequency (Hz)
omega = 2*pi*max(1e-6, f0);                 % rad/s (natural frequency ω)
gamma = double(getin('CERD_THEORY_GAMMA', 80));  % 1/s (damping γ)
omega_echo = sqrt(max(0, omega^2 - gamma^2));  % rad/s (echo frequency √(ω²-γ²))
w = 2*pi*f(:);                              % ω' grid (frequency in rad/s)
% Transfer function: |H(ω')|² = 1/((ω_echo²-ω'²)² + (2γω')²)
% This is the standard form for a damped oscillator power spectrum
H2 = 1 ./ ((omega_echo^2 - w.^2).^2 + (2*gamma.*w).^2);
H2 = real(H2);
H2(~isfinite(H2)) = 0;
t = H2(mask);
if ~any(isfinite(t)) || all(t==0)
    score = NaN; return;
end
% Normalize template to unit energy to make gain identifiable
t = t / max(eps, sqrt(sum(t.^2)));

% Fit gain (least squares, non-negative)
s_post = P_post(mask);
s_base = P_base(mask);
A_post = max(0, (t.'*s_post) / max(eps, (t.'*t)));
A_base = max(0, (t.'*s_base) / max(eps, (t.'*t)));

% Return log10 ratio (post vs base)
score = log10( (A_post + cfg.eps) / (A_base + cfg.eps) );
score = real(score);
end

function diag = compute_psd_overlay(x, fs, cfg, center_freq)
diag = struct('valid',false,'f',[],'P_base',[],'P_post',[],'theory',[],'score_db',NaN);
if isempty(x) || ~isfinite(fs)
    return;
end
x = double(x(:));
pre_samples  = max(1, round(cfg.pre_ms  * fs / 1000));
post_samples = max(1, round(cfg.post_ms * fs / 1000));
if pre_samples + post_samples > numel(x)
    return;
end
baseline = x(1:pre_samples);
segment  = x(pre_samples+1:pre_samples+post_samples);
if cfg.detrend
    if numel(baseline)>1, baseline = detrend(baseline); end
    if numel(segment)>1,  segment  = detrend(segment);  end
end
if numel(segment)<8 || numel(baseline)<8
    return;
end
win = max(16, round(cfg.window_ms * fs / 1000));
win = min(win, numel(segment));
if win >= numel(segment)
    win = max(8, floor(numel(segment)/2));
end
if win < 8, return; end
noverlap = min(win-1, max(0, round(cfg.overlap * win)));
nfft = max(win, round(cfg.nfft));
[P_post,f] = pwelch(segment, win, noverlap, nfft, fs);
[P_base,~] = pwelch(baseline, win, noverlap, nfft, fs);
band = cfg.band_hz;
if isempty(band) || numel(band) < 2 || any(~isfinite(band))
    bw = cfg.bandwidth_hz;
    if ~isfinite(center_freq) || center_freq <= 0
        band = cfg.fallback_band;
    else
        if ~isfinite(bw) || bw <= 0, bw = max(100, center_freq*0.5); end
        half_bw = bw/2;
        band = [max(0, center_freq-half_bw), center_freq+half_bw];
    end
end
band(1) = max(0, band(1));
band(2) = min(cfg.max_freq_hz, max(band(1)+1, band(2)));
mask = (f >= band(1) & f <= band(2));
if ~any(mask)
    [~, idx] = min(abs(f - mean(band)));
    mask = false(size(f)); mask(idx) = true;
end
f0_band = getin('CERD_THEORY_F0_HZ', [60 90]);
f0 = mean(double(f0_band(:)));
omega = 2*pi*max(1e-6, f0);
gamma = double(getin('CERD_THEORY_GAMMA', 80));
omega_echo = sqrt(max(0, omega^2 - gamma^2));
w = 2*pi*f(:);
H2 = 1 ./ ((omega_echo^2 - w.^2).^2 + (2*gamma.*w).^2);
H2 = real(H2);
H2(~isfinite(H2)) = 0;
t = H2(mask);
if ~any(isfinite(t)) || all(t==0)
    return;
end
t = t / max(eps, sqrt(sum(t.^2)));
s_post = P_post(mask);
s_base = P_base(mask);
A_post = max(0, (t.'*s_post) / max(eps, (t.'*t)));
A_base = max(0, (t.'*s_base) / max(eps, (t.'*t)));
theory_curve_full = nan(size(f));
band_theory = A_post * t;
theory_curve_full(mask) = band_theory;
diag.f = f;
diag.P_base = P_base;
diag.P_post = P_post;
diag.theory = theory_curve_full;
diag.f_band = f(mask);
diag.theory_band = band_theory;
diag.score_db = 10*log10((A_post + cfg.eps) / (A_base + cfg.eps));
diag.valid = true;
end

function y = prewhiten_ar1(x)
% simple AR(1) prewhitener: y[n] = x[n] - rho x[n-1]
x = x(:);
if numel(x) < 5, y = x; return; end
x0 = x - mean(x);
% lag-1 autocorr (biased)
r0 = sum(x0.*x0); r1 = sum(x0(2:end).*x0(1:end-1));
rho = 0; if r0>0, rho = max(min(r1/r0, 0.98), -0.98); end
y = filter([1 -rho], 1, x0);
sd = std(y); if sd>0, y = y/sd; end
end

function y = bandlimit(x, fs, f1, f2)
% 4th-order zero-phase Butterworth bandpass (falls back gracefully)
x = x(:);
f1 = max(0.1, min(f1, fs/2-1));
f2 = max(f1+0.1, min(f2, fs/2-0.1));
try
    [b,a] = butter(4, [f1 f2]/(fs/2), 'bandpass');
    y = filtfilt(b,a,x);
catch
    % if Signal Toolbox missing, do a crude FFT mask
    N = numel(x); X = fft(x);
    f = (0:N-1)'*(fs/N);
    mask = (f>=f1 & f<=f2) | (f>=fs-f2 & f<=fs-f1);
    X(~mask) = 0;
    y = real(ifft(X));
end
end

function g = theory_template(t, tau, omega, gamma_1ps)
% THEORY_TEMPLATE  Generate Green's function template (non-gated version)
% Based on Eq. (20) from Notes.txt:
%   G = e^(-γt)/v * J₀[v√((ω²-γ²)(t²-Δu²))] * H(t-Δu)
%
% Uses echo frequency ω_echo = √(ω²-γ²) and Bessel argument ω_echo * √(t²-τ²)
%
% Inputs:
%   t: time vector (seconds)
%   tau: delay (seconds) = Δu from theory
%   omega: natural frequency (rad/s) = 2πf0
%   gamma_1ps: damping rate (1/s)

dt = t - tau;
mask = (dt >= 0);

% Calculate echo frequency from theory: ω_echo = √(ω²-γ²)
omega_sq = omega^2;
gamma_sq = gamma_1ps^2;
omega_echo = sqrt(max(0, omega_sq - gamma_sq));  % Echo frequency (rad/s)

% Build template according to theory: J₀[ω_echo * √(t²-τ²)] for t >= τ
t_sq = t.^2;
tau_sq = tau^2;
arg_sq = max(0, t_sq - tau_sq);  % Ensure non-negative for sqrt
bessel_arg = omega_echo * sqrt(arg_sq);

% Exponential decay: e^(-γ(t-τ)) for t >= τ
decay = exp(-gamma_1ps * dt);

% Bessel function with theoretical argument
bessel_val = besselj(0, bessel_arg);

% Combine: decay * Bessel, with causality mask
a = decay .* bessel_val;
g = zeros(size(t)); 
g(mask) = a(mask);
end

% ------------------------ helpers from your current file -------------------
function fs_local = local_fs(fs_global, macro, t0, t1)
if ~isscalar(fs_global) || ~isfinite(fs_global) || fs_global <= 0
    dur_s = double(t1 - t0)/1e6;
    fs_local = numel(macro)/max(eps,dur_s);
else
    fs_local = fs_global;
end
end

function fb = local_band(lb,ub,f0,macro_band_mode)
if macro_band_mode
    fb = [lb ub];
else
    if ~isfinite(f0), f0 = (lb+ub)/2; end
    halfwin = min(min(100,(ub-lb)/2), max(1e-6, min(f0-lb, ub-f0)));
    fb = [f0 - halfwin, f0 + halfwin];
end
end

function f_ref = local_fref(lb,ub,f0,macro_band_mode)
if macro_band_mode
    f_ref = (lb+ub)/2;
else
    if ~isfinite(f0), f0 = (lb+ub)/2; end
    f_ref = f0;
end
end

% ===================== (keep the rest of your existing helpers) =====================
% Everything below here is identical to the versions you posted already:
%   - prep_basics, extract_UFOs, exclusion_zones, parse_null_range, range_to_str,
%     write_metric_diagnostics, print_diag, compare_to_nulls_fdr, build_ufo_drop_mask,
%     empirical spectral functions (estimate_power, init_emp_arrays),
%     masking utilities, merge_intervals, sample_from_intervals, etc.
%
% For completeness, I've included the ones your new code calls directly.

function [bandpow, snr, ppeak, fpeak] = init_emp_arrays(M,N,Pmax)
bandpow = NaN(M,N,Pmax);
snr     = NaN(M,N,Pmax);
ppeak   = NaN(M,N,Pmax);
fpeak   = NaN(M,N,Pmax);
end

function print_diag(D, tag)
kn = ternary(isfinite(D.K_needed_floor), sprintf('%d', D.K_needed_floor), 'Inf');
fprintf('[%s] K=%d | m=%d | p_floor=%.5g | p_needed=%.3g | crit_p=%.3g | min_p=%.5g | min_b=%d | K_needed_floor=%s | worth_more_nulls=%d\n', ...
    tag, D.K, D.m, D.p_floor, D.p_needed, D.crit_p, D.min_p, D.min_b, kn, D.worth_more_nulls); drawnow;
end

function y = ternary(cond,a,b)
if ~islogical(cond), cond = logical(cond); end
if cond, y=a; else, y=b; end
end

% ==== include your previously posted implementations for these unchanged ====
% prep_basics, extract_UFOs, exclusion_zones, parse_null_range, range_to_str,
% compare_to_nulls_fdr, write_metric_diagnostics, build_fdr_groups,
% merge_intervals, intersect_allow_with_domain, complement_intervals,
% sample_from_intervals, parsave, getGlobalStreamInfo, estimate_power, getin
% (Use the versions from your existing file; they don't need edits for THEORY.)
% ============================================================================

function v = getin(name, default)
if evalin('base', sprintf('exist(''%s'',''var'')', name))
    v = evalin('base', name);
else
    v = default;
end
end

function D = prep_basics(pID)

% --- General Setup ---
D.basedir   = '/Users/erik/Library/CloudStorage/Dropbox/Work/UFO/';
D.password  = 'bemena';
patientIDs  = {'pat001','pat003','pat004','pat061','pat072','pat079','pat448','pat453'};

fast_meta = getin('CERD_FAST_META', false);          % set true for problematic subjects (e.g., pat448)
fast_meta_ids = getin('CERD_FAST_META_IDS', {});     % cellstr of IDs to force fast path

% --- Patient ID Handling ---
if isnumeric(pID)
    if pID < 1 || pID > numel(patientIDs)
        error('prep_basics: patient index %d out of range 1..%d', pID, numel(patientIDs));
    end
    D.ID = patientIDs{pID};
else
    pID = string(pID);
    ix = find(strcmp(patientIDs, pID), 1);
    if isempty(ix), error('Unknown patient ID: %s', pID); end
    D.ID = patientIDs{ix};
end

% --- Read Metadata ---
default_path = fullfile(D.basedir, 'Data', D.ID, 'sub.mefd');
local_root   = '/Users/erik/Desktop/tmp_mef';
local_copy   = fullfile(local_root, D.ID, 'sub.mefd');

if exist(local_copy,'dir') == 7
    D.meta_path = local_copy;
else
    D.meta_path = default_path;
end
fprintf('[prep_basics] Using MEF for %s at: %s\n', D.ID, D.meta_path);

assert(exist(D.meta_path,'dir') == 7, ...
    'prep_basics: MEF folder not found for %s at %s', D.ID, D.meta_path);

% --- Fast metadata mode handling ---
fast_meta_ids = getin('CERD_FAST_META_IDS', {});
if ischar(fast_meta_ids) || isstring(fast_meta_ids)
    fast_meta_ids = cellstr(fast_meta_ids);
end
use_fast = fast_meta || any(strcmpi(fast_meta_ids, 'ALL')) || ...
    any(strcmpi(fast_meta_ids, D.ID));

if use_fast
    fprintf('[prep_basics] FAST META mode active for %s — skipping full MEF read.\n', D.ID);
    D.meta = struct();  % placeholder – not used by the readers below
else
    D.meta = readMef3(D.meta_path, D.password);
end

% --- Patch for FAST META: create minimal fake channel metadata ---
if use_fast
    chfile = fullfile(D.basedir, 'Data', D.ID, ...
        sprintf('sub-%s_ses-01_task-ufos_run-00_channels.txt', upper(D.ID)));
    if exist(chfile,'file')
        T = readtable(chfile, 'FileType', 'text', 'Delimiter', '\t');
        if any(strcmpi(T.Properties.VariableNames, 'name'))
            chnames = string(T{:, strcmpi(T.Properties.VariableNames, 'name')});
        else
            chnames = [];
        end
        D.meta.time_series_channels = struct('name', num2cell(chnames));
        if isempty(D.meta.time_series_channels)
            warning('prep_basics: FAST META created empty channel struct for %s. Check channel file path.', D.ID);
        end
        if any(strcmpi(T.Properties.VariableNames, 'sampling_frequency'))
            fsvals = T{:, strcmpi(T.Properties.VariableNames, 'sampling_frequency')};
            for k = 1:numel(D.meta.time_series_channels)
                D.meta.time_series_channels(k).metadata.sampling_frequency = fsvals(k);
            end
            D.fs = mode(fsvals(fsvals>0));
        end
    else
        warning('prep_basics: FAST META but no channel file found for %s.', D.ID);
        D.meta.time_series_channels = struct([]);
    end
end


% --- Load UFOs ---
ufo_json = fullfile(D.basedir, 'Data', D.ID, 'UFO.json');
if exist(ufo_json, 'file')
    try
        D.UFOs = jsondecode(fileread(ufo_json));
    catch
        warning('prep_basics: UFO.json unreadable for %s; proceeding without.', D.ID);
        D.UFOs = struct([]);
    end
else
    D.UFOs = struct([]);
end

% --- Sampling Frequency ---
if isfield(D.meta, 'time_series_metadata') && ...
        isfield(D.meta.time_series_metadata, 'section_2') && ...
        isfield(D.meta.time_series_metadata.section_2, 'sampling_frequency')
    % Standard Brno format
    fs_meta = D.meta.time_series_metadata.section_2.sampling_frequency;
elseif isfield(D.meta, 'time_series_channels') && ...
        isstruct(D.meta.time_series_channels) && ...
        ~isempty(D.meta.time_series_channels) && ...
        isfield(D.meta.time_series_channels(1), 'metadata') && ...
        isfield(D.meta.time_series_channels(1).metadata, 'sampling_frequency')
    % Mayo format (per channel)
    fs_values = arrayfun(@(c) c.metadata.sampling_frequency, D.meta.time_series_channels);
    fs_meta = mode(fs_values(fs_values>0));  % mode of valid freqs
else
    warning('Sampling frequency field not found; setting fs_meta = NaN');
    fs_meta = NaN;
end

if fs_meta > 0
    D.fs = fs_meta;
else
    chfile = fullfile(D.basedir, 'Data', D.ID, ...
        sprintf('sub-%s_ses-01_task-ufos_run-00_channels.txt', upper(D.ID)));
    if exist(chfile, 'file')
        T = readtable(chfile, 'FileType', 'text', 'Delimiter', '\t');
        if any(strcmpi(T.Properties.VariableNames, 'sampling_frequency'))
            fs_vals = T{:, strcmpi(T.Properties.VariableNames, 'sampling_frequency')};
            uniq_fs = unique(fs_vals);
            if numel(uniq_fs) == 1
                D.fs = uniq_fs;
            else
                warning('prep_basics: multiple fs values in %s; fs set to -1.', D.ID);
                D.fs = -1;
            end
        else
            warning('prep_basics: no sampling_frequency column in channel file for %s.', D.ID);
            D.fs = -1;
        end
    else
        warning('prep_basics: channel file not found for %s; fs set to -1.', D.ID);
        D.fs = -1;
    end
end

% --- Time Bounds ---
if ~isempty(D.meta) && isfield(D.meta,'earliest_start_time') && isfield(D.meta,'latest_end_time') ...
        && isfinite(double(D.meta.earliest_start_time)) && isfinite(double(D.meta.latest_end_time))
    D.t0 = D.meta.earliest_start_time;
    D.t1 = D.meta.latest_end_time;
else
    [t0_est, t1_est] = time_bounds_from_ufo_struct(D.UFOs);
    if isfinite(double(t0_est)) && isfinite(double(t1_est)) && double(t1_est) > double(t0_est)
        D.t0 = t0_est;
        D.t1 = t1_est;
    else
        warning('prep_basics: could not infer t0/t1; using [0, inf).');
        D.t0 = int64(0);
        D.t1 = int64(2^62);
    end
end

% --- Patient-Specific Channel Mapping ---
switch D.ID
    case 'pat001'
        [mi_groups, ma_names, allmi, mima] = auto_map_wroclaw(D);
        if isempty(allmi) || isempty(ma_names)
            error(['prep_basics: could not auto-map channels for %s. ' ...
                'Ensure Data/%s/channels.tsv exists.'], D.ID, D.ID);
        end
        D.allmi = allmi;
        for j = 1:numel(mi_groups), D.mi{j} = mi_groups{j}; end
        for j = 1:numel(ma_names),  D.ma{j} = {ma_names{j}};  end
        D.mima = mima;

    case 'pat003'
        D.mi{1} = {'CSC64','CSC65','CSC66','CSC67','CSC68','CSC69','CSC70','CSC71'};
        D.mi{2} = {'CSC72','CSC73','CSC74','CSC75','CSC76','CSC77','CSC78','CSC79'};
        D.ma{1} = {'SMAAL1'};
        D.ma{2} = {'SMAL1'};
        D.allmi = [D.mi{1}, D.mi{2}];
        D.mima = eye(2);

    case 'pat004'
        D.mi{1} = {'CSC72','CSC73','CSC74','CSC75','CSC76','CSC77','CSC78','CSC79'};
        D.mi{2} = {'CSC96','CSC97','CSC98','CSC99','CSC100','CSC101','CSC102','CSC103','CSC104','CSC105'};
        D.mi{3} = {'CSC80','CSC81','CSC82','CSC83','CSC84','CSC85','CSC86','CSC87'};
        D.mi{4} = {'CSC64','CSC65','CSC66','CSC67','CSC68','CSC69','CSC70','CSC71'};
        D.ma{1} = {'HH1'};
        D.ma{2} = {'TP1'};
        D.ma{3} = {'HB1'};
        D.ma{4} = {'AMG1'};
        D.allmi = [D.mi{1}, D.mi{2}, D.mi{3}, D.mi{4}];
        D.mima = eye(4);

    case 'pat061'
        D.mi{1} = {'Bm1a','Bm2a','Bm3a'};
        D.mi{2} = {'Bm1a','Bm2a','Bm3a','Bm4a','Bm5a'};
        D.mi{3} = {'Bm4a','Bm5a'};
        D.ma{1} = {'B1'};
        D.ma{2} = {'B2'};
        D.ma{3} = {'B3'};
        D.allmi = {'Bm1a','Bm2a','Bm3a','Bm4a','Bm5a'};
        D.mima  = [1,1,1,0,0 ; 1,1,1,1,1 ; 0,0,0,1,1];

    case {'pat072','pat079'}
        D.mi{1} = {'Bm1a','Bm2a','Bm3a'};
        D.mi{2} = {'Bm1a','Bm2a','Bm3a','Bm4a','Bm5a','Bm6a'};
        D.mi{3} = {'Bm4a','Bm5a','Bm6a'};
        D.ma{1} = {'B''1'};
        D.ma{2} = {'B''2'};
        D.ma{3} = {'B''3'};
        D.allmi = {'Bm1a','Bm2a','Bm3a','Bm4a','Bm5a','Bm6a'};
        D.mima  = [1 1 1 0 0 0; 1 1 1 1 1 1; 0 0 0 1 1 1];

    case 'pat448'
        D.mi{1} = {'ADShaft_01','ADShaft_02','ADShaft_03','ADShaft_04'};
        D.mi{2} = {'ADShaft_01','ADShaft_02','ADShaft_03','ADShaft_04', ...
            'ADShaft_05','ADShaft_06','ADShaft_07','ADShaft_08'};
        D.mi{3} = {'ADShaft_05','ADShaft_06','ADShaft_07','ADShaft_08', ...
            'ADShaft_09','ADShaft_10','ADShaft_11','ADShaft_12'};
        D.mi{4} = {'ADShaft_09','ADShaft_10','ADShaft_11','ADShaft_12', ...
            'ADShaft_13','ADShaft_14','ADShaft_15','ADShaft_16'};
        D.ma{1} = {'PDMacro_01'};
        D.ma{2} = {'PDMacro_02'};
        D.ma{3} = {'PDMacro_03'};
        D.ma{4} = {'PDMacro_04'};
        D.allmi = {'ADShaft_01','ADShaft_02','ADShaft_03','ADShaft_04', ...
            'ADShaft_05','ADShaft_06','ADShaft_07','ADShaft_08', ...
            'ADShaft_09','ADShaft_10','ADShaft_11','ADShaft_12', ...
            'ADShaft_13','ADShaft_14','ADShaft_15','ADShaft_16'};
        D.mima = zeros(4,16);
        D.mima(1,1:4)   = 1;
        D.mima(2,1:8)   = 1;
        D.mima(3,5:12)  = 1;
        D.mima(4,9:16)  = 1;

    case 'pat453'
        D.mi{1} = {'RBundle_01','RBundle_02','RBundle_03','RBundle_04', ...
            'RBundle_05','RBundle_06','RBundle_07','RBundle_08'};
        D.ma{1} = {'RMacro_01'};
        D.allmi = D.mi{1};
        D.mima  = ones(1,8);
end

% Set allma = ma for compatibility with visualization scripts
D.allma = D.ma;

end


function [mi_groups, ma_names, allmi, mima] = auto_map_wroclaw(D)
chtsv = fullfile(D.basedir,'Data',D.ID,'channels.tsv');
have_tsv = exist(chtsv,'file')==2;
names = {};
fscol = [];
if have_tsv
    T = readtable(chtsv, 'FileType','text', 'Delimiter','\t');
    nameVar = local_find_var(T, {'name','channel','label','ch_name','chan_name'});
    if isempty(nameVar), error('auto_map_wroclaw: no channel-name column in channels.tsv for %s.', D.ID); end
    names = T.(nameVar);
    if isstring(names), names = cellstr(names); end
    if ischar(names),   names = {names}; end
    if ~iscellstr(names), error('auto_map_wroclaw: channel names column is not text for %s.', D.ID); end
    fsVar = local_find_var(T, {'sampling_frequency','fs','freq','samplingRate','rate'});
    if ~isempty(fsVar), fscol = T.(fsVar); end
end

uCSC = {};
if ~isempty(names)
    mask = strncmpi(names,'CSC',3);
    if ~isempty(fscol) && isnumeric(fscol)
        mask = mask | (fscol >= 30000);
    end
    uCSC = unique(names(mask));
else
    U = D.UFOs;
    chans = {};
    if ~isempty(U) && isstruct(U)
        fld = '';
        if isfield(U,'channels'), fld = 'channels'; end
        if isempty(fld) && isfield(U,'channel'), fld = 'channel'; end
        if ~isempty(fld)
            for k = 1:numel(U)
                c = U(k).(fld);
                if ischar(c)
                    chans{end+1} = c; %#ok<AGROW>
                elseif isstring(c)
                    chans = [chans, cellstr(c)']; %#ok<AGROW>
                elseif iscell(c)
                    tmp = cellfun(@(x)string(x), c, 'UniformOutput', false);
                    tmp = cellfun(@char, tmp, 'UniformOutput', false);
                    chans = [chans, tmp]; %#ok<AGROW>
                end
            end
        end
    end
    u = unique(chans);
    uCSC = u(strncmp(u,'CSC',3));
end

if isempty(uCSC)
    mi_groups = {}; ma_names = {}; allmi = {}; mima = zeros(0,0);
    return
end

pref = cell(size(uCSC));
keep = false(size(uCSC));
for k = 1:numel(uCSC)
    s = uCSC{k};
    tokA = regexp(s, '^CSC_([A-Za-z]+)_?0*([0-9]+)$', 'tokens', 'once');
    tokB = regexp(s, '^CSC0*([0-9]+)_([A-Za-z]+)$', 'tokens', 'once');
    if ~isempty(tokA)
        pref{k} = upper(tokA{1});
        keep(k) = true;
    elseif ~isempty(tokB)
        pref{k} = upper(tokB{2});
        keep(k) = true;
    end
end
uCSC = uCSC(keep);
pref = pref(keep);

if isempty(uCSC)
    mi_groups = {}; ma_names = {}; allmi = {}; mima = zeros(0,0);
    return
end

families = unique(pref, 'stable');

ma_names = cell(1, numel(families));
for j = 1:numel(families)
    fam = families{j};
    candidates = {[fam '1'], [fam '01'], upper([fam '1']), lower([fam '1'])};
    chosen = '';
    if ~isempty(names)
        for c = 1:numel(candidates)
            hit = find(strcmpi(names, candidates{c}), 1);
            if ~isempty(hit), chosen = names{hit}; break; end
        end
    end
    if isempty(chosen), chosen = [fam '1']; end
    ma_names{j} = chosen;
end

mi_groups = cell(1, numel(families));
for j = 1:numel(families)
    mi_groups{j} = uCSC(strcmp(pref, families{j}));
end
allmi = unique(vertcat(mi_groups{:}), 'stable');
K = numel(ma_names);
M = numel(allmi);
mima = zeros(K, M);
for j = 1:K
    hits = ismember(allmi, mi_groups{j});
    mima(j, hits) = 1;
end
end


function varname = local_find_var(T, candidates)
varname = '';
V = string(T.Properties.VariableNames);
for c = 1:numel(candidates)
    ix = find(strcmpi(V, candidates{c}), 1);
    if ~isempty(ix)
        varname = T.Properties.VariableNames{ix};
        return
    end
end
end


function UFOd = extract_UFOs(Dat)
M = numel(Dat.allmi);
UFOd = cell(1,M);

U = Dat.UFOs;
if isempty(U) || ~isstruct(U)
    for ii = 1:M, UFOd{ii} = zeros(0,3); end
    return;
end

chanField = '';
if isfield(U, 'channels'), chanField = 'channels';
elseif isfield(U, 'channel'), chanField = 'channel';
else
    warning('extract_UFOs: no channel(s) field found in UFO.json; returning empties.');
    for ii = 1:M, UFOd{ii} = zeros(0,3); end
    return;
end

leftField  = local_pick_field(U, {'uutc_left','left','t0','start','start_uutc'});
rightField = local_pick_field(U, {'uutc_right','right','t1','end','end_uutc'});
freqField  = local_pick_field(U, {'frequency','freq','f0','peak_frequency_hz','peak_hz'});

if isempty(leftField) || isempty(rightField) || isempty(freqField)
    warning('extract_UFOs: missing one of left/right/freq fields; returning empties.');
    for ii = 1:M, UFOd{ii} = zeros(0,3); end
    return;
end

for ii = 1:M
    mic = Dat.allmi{ii};
    keep = false(1, numel(U));
    for k = 1:numel(U)
        ch = U(k).(chanField);
        keep(k) = local_has_label(ch, mic);
    end
    if ~any(keep)
        UFOd{ii} = zeros(0,3);
        continue;
    end

    Uk = U(keep);
    n  = numel(Uk);
    LRF = nan(n,3);
    for j = 1:n
        a = Uk(j).(leftField);
        b = Uk(j).(rightField);
        f = Uk(j).(freqField);
        if isempty(a) || isempty(b) || ~isfinite(double(a)) || ~isfinite(double(b)) || double(b) <= double(a)
            continue;
        end
        LRF(j,1) = double(a);
        LRF(j,2) = double(b);
        LRF(j,3) = double(f);
    end

    LRF = LRF(all(isfinite(LRF),2), :);
    if isempty(LRF), UFOd{ii} = zeros(0,3); continue; end
    if isfinite(double(Dat.t0)) && isfinite(double(Dat.t1))
        inb = LRF(:,1) >= double(Dat.t0) & LRF(:,2) <= double(Dat.t1);
        LRF = LRF(inb, :);
    end
    if isempty(LRF), UFOd{ii} = zeros(0,3); continue; end

    [~,ord] = sort(LRF(:,2), 'ascend');
    LRF = LRF(ord, :);
    UFOd{ii} = LRF;
end
end


function fld = local_pick_field(S, candidates)
fld = '';
for c = 1:numel(candidates)
    if isfield(S, candidates{c}), fld = candidates{c}; return; end
end
end



function tf = local_has_label(ch, mic)
if ischar(ch) || (isstring(ch) && isscalar(ch))
    tf = strcmp(char(ch), mic);
elseif isstring(ch) || iscellstr(ch) || iscell(ch)
    if isstring(ch), ch = cellstr(ch); end
    ch = ch(:);
    for i = 1:numel(ch), if ~ischar(ch{i}), ch{i} = char(string(ch{i})); end, end
    tf = any(strcmp(ch, mic));
else
    tf = false;
end
end



function Exclude = exclusion_zones(Dat, gap_ufo_us, mima_override)
if isempty(gap_ufo_us)
    pre = int64(0); post = int64(0);
elseif numel(gap_ufo_us) == 1
    pre  = int64(gap_ufo_us);
    post = int64(gap_ufo_us);
else
    pre  = int64(gap_ufo_us(1));
    post = int64(gap_ufo_us(2));
end

UFOd = extract_UFOs(Dat);
M    = numel(UFOd);
t0   = int64(Dat.t0);
t1   = int64(Dat.t1);

exc = cell(1, M);
for ii = 1:M
    U = UFOd{ii};
    if isempty(U)
        exc{ii} = zeros(0,2,'int64');
        continue;
    end
    L = int64(U(:,1)) - pre;
    R = int64(U(:,2)) + post;
    L = max(L, t0);
    R = min(R, t1);
    ok = R >= L;
    if ~any(ok)
        exc{ii} = zeros(0,2,'int64');
    else
        I = [L(ok) R(ok)];
        exc{ii} = merge_intervals(I);
    end
end

if nargin < 3 || isempty(mima_override)
    mima_override = Dat.mima;
end

Nmac = size(mima_override, 1);
Exclude = cell(1, Nmac);
for jj = 1:Nmac
    micIdx = find(mima_override(jj, :) ~= 0);
    if isempty(micIdx)
        Exclude{jj} = zeros(0,2,'int64');
        continue;
    end
    U = vertcat(exc{micIdx});
    if isempty(U)
        Exclude{jj} = zeros(0,2,'int64');
    else
        Exclude{jj} = merge_intervals(U);
    end
end
end


function M = merge_intervals(I, gap_tol)
if nargin < 2 || isempty(gap_tol), gap_tol = 1; end
if isempty(I)
    M = reshape(I, 0, 2);
    return;
end
if size(I,2) ~= 2
    error('merge_intervals: I must be n x 2.');
end

cls = class(I);
tol = cast(gap_tol, cls);

ok = all(isfinite(double(I)), 2);
I  = I(ok, :);
if isempty(I)
    M = zeros(0,2, cls);
    return;
end

swap = I(:,1) > I(:,2);
if any(swap)
    I(swap,:) = [I(swap,2) I(swap,1)];
end

I = sortrows(I, [1 2]);

M = I(1,:);
for k = 2:size(I,1)
    if I(k,1) <= M(end,2) + tol
        if I(k,2) > M(end,2)
            M(end,2) = I(k,2);
        end
    else
        M(end+1,:) = I(k,:); %#ok<AGROW>
    end
end
end


function [bp, f_peak, prom_lin, score_val, snr_val] = estimate_power( ...
    mac, fs, freq_band, min_prom_dB, use_multitaper, f0, sigma_hz, compute_snr)

x = double(mac(:));
N = numel(x);
bp = NaN; f_peak = NaN; prom_lin = NaN; score_val = NaN; snr_val = NaN;
if N < 8, return; end

persistent fs0 b a
if isempty(fs0) || fs0 ~= fs
    fs0 = fs;
    [b,a] = butter(2, 1/(fs/2), 'high');
end
x = filtfilt(b,a,x);

if use_multitaper && exist('pmtm','file')==2
    [Pxx, f] = pmtm(x, 3.5, max(256, 2^nextpow2(N)), fs);
else
    if exist('hann','builtin')==5 || exist('hann','file')==2
        win = hann(N,'periodic');
    else
        win = hanning(N);
    end
    winE = sum(win.^2);
    nfft = 2^nextpow2(N);
    X = fft(x.*win, nfft);
    P = (abs(X).^2) / (fs * winE);
    keep = 1:floor(nfft/2)+1;
    Pxx = P(keep);
    if numel(keep)>2, Pxx(2:end-1) = 2*Pxx(2:end-1); end
    f = (keep-1)' * (fs/nfft);
end

lb = freq_band(1); ub = freq_band(2);
ib  = (f >= lb) & (f <= ub);
f_b = f(ib);  P_b = Pxx(ib);
if isempty(f_b)
    f_b = f; P_b = Pxx; lb = f(1); ub = f(end);
end

bp = trapz(f_b, P_b);

dB = 10*log10(P_b + eps);
[pk_heights_dB, locs, ~, prom_dB] = findpeaks(dB, f_b, 'MinPeakProminence', min_prom_dB);
if isempty(pk_heights_dB)
    [~, imax] = max(dB);
    f_peak = f_b(imax);
    prom_dB_eff = max(0, dB(imax) - median(dB));
else
    [~, i_max] = max(prom_dB);
    f_peak = locs(i_max);
    prom_dB_eff = prom_dB(i_max);
end

[~, i_loc] = min(abs(f_b - f_peak));
if i_loc > 1 && i_loc < numel(f_b)
    f0p = f_b(i_loc-1:i_loc+1);
    d0p = dB(i_loc-1:i_loc+1);
    pp  = polyfit(f0p, d0p, 2);
    if isfinite(pp(1)) && pp(1) ~= 0
        f_peak = -pp(2)/(2*pp(1));
        f_peak = min(max(f_peak, lb), ub);
    end
end

prom_lin = 10^(prom_dB_eff / 10);

if ~exist('f0','var') || ~isfinite(f0), f0 = (lb + ub)/2; end
if ~exist('sigma_hz','var') || ~isfinite(sigma_hz) || sigma_hz <= 0
    sigma_hz = (ub - lb)/6;
end
if bp <= 0
    score_val = 0;
else
    score_val = bp * exp(-0.5 * ((f_peak - f0)/sigma_hz)^2);
end

if compute_snr
    gap = 5; flw = 20;
    L1 = max(min(f), lb - (gap+flw)); L2 = max(min(f), lb - gap);
    R1 = min(max(f), ub + gap);       R2 = min(max(f), ub + (gap+flw));
    iL = (f >= L1) & (f <= L2);
    iR = (f >= R1) & (f <= R2);
    flank_power = 0;
    if any(iL), flank_power = flank_power + trapz(f(iL), Pxx(iL)); end
    if any(iR), flank_power = flank_power + trapz(f(iR), Pxx(iR)); end
    nseg = double(any(iL)) + double(any(iR));
    if nseg > 0, flank_power = flank_power / nseg; end
    snr_val = bp / max(eps, flank_power);
end
end


function rngv = parse_null_range(val, Nnull)
if isempty(val) || (isstring(val) && strlength(val)==0)
    rngv = int64(1:Nnull);
    return;
end

if isnumeric(val)
    v = val(:);
    if any(abs(v - round(v)) > 1e-9)
        error('CERD_NULL_RANGE numeric vector must contain integers.');
    end
    v = int64(round(v));
    v = v(v>=1 & v<=Nnull);
    rngv = unique(v,'sorted');
    return;
end

if ischar(val) || isstring(val)
    s = char(val);
    s = strtrim(s);

    if ~isempty(s) && any(s(1)==['[','(']) && any(s(end)==[']',')'])
        s = strtrim(s(2:end-1));
    end

    tok3 = regexp(s, '^\s*(\d+)\s*:\s*(\d+)\s*:\s*(\d+)\s*$', 'tokens', 'once');
    if ~isempty(tok3)
        a = str2double(tok3{1}); step = str2double(tok3{2}); b = str2double(tok3{3});
        if step == 0, error('CERD_NULL_RANGE step cannot be zero.'); end
        if (b >= a && step < 0) || (b <= a && step > 0)
            v = [];
        else
            v = a:step:b;
        end
        v = int64(v(:));
        v = v(v>=1 & v<=Nnull);
        rngv = unique(v,'sorted');
        return;
    end

    tok2 = regexp(s, '^\s*(\d+)\s*:\s*(\d+)\s*$', 'tokens', 'once');
    if ~isempty(tok2)
        a = str2double(tok2{1}); b = str2double(tok2{2});
        if b >= a
            v = a:1:b;
        else
            v = a:-1:b;
        end
        v = int64(v(:));
        v = v(v>=1 & v<=Nnull);
        rngv = unique(v,'sorted');
        return;
    end

    s = regexprep(s, '[,;]+', ' ');
    parts = regexp(s, '\s+', 'split');
    parts = parts(~cellfun('isempty',parts));
    if isempty(parts)
        rngv = int64(1:Nnull);
        return;
    end
    nums = str2double(parts);
    if any(~isfinite(nums))
        error('Could not parse CERD_NULL_RANGE="%s".', char(val));
    end
    if any(abs(nums - round(nums)) > 1e-9)
        error('CERD_NULL_RANGE list must contain integers.');
    end
    v = int64(round(nums(:)));
    v = v(v>=1 & v<=Nnull);
    rngv = unique(v,'sorted');
    return;
end

error('CERD_NULL_RANGE must be string or numeric vector.');
end


function s = range_to_str(v)
v = v(:)';
if isempty(v), s = '(empty)'; return; end
if numel(v)>1 && all(diff(v)==1)
    s = sprintf('%d:%d', v(1), v(end));
else
    s = sprintf('%d,', v); s(end) = [];
end
end


function allow = complement_intervals(domain, bad)
a = domain(1); b = domain(2);
if isempty(bad)
    if b >= a, allow = [a b]; else, allow = zeros(0,2,'like',a); end
    return
end
bad = bad(bad(:,2) >= bad(:,1), :);
if isempty(bad), allow = [a b]; return; end
bad = sortrows(bad,1);
allow = zeros(0,2,'like',a);
cur = a;
for t = 1:size(bad,1)
    if bad(t,1) > cur, allow(end+1,:) = [cur, bad(t,1)-1]; end
    cur = max(cur, bad(t,2)+1);
    if cur > b, break; end
end
if cur <= b, allow(end+1,:) = [cur, b]; end
end


function allow2 = intersect_allow_with_domain(allow, dom)
if isempty(allow)
    allow2 = reshape(allow, 0, 2);
    return;
end
if size(allow,2) ~= 2
    error('intersect_allow_with_domain: allow must be n x 2.');
end
if numel(dom) ~= 2
    error('intersect_allow_with_domain: dom must be 1 x 2.');
end

cls = class(allow);
dom  = cast(dom(:).', cls);
L = dom(1); R = dom(2);
if R < L
    tmp = L; L = R; R = tmp;
end

A = allow;

if ~isinteger(A)
    ok = all(isfinite(A), 2);
    A  = A(ok, :);
    if isempty(A)
        allow2 = zeros(0,2, cls);
        return;
    end
end

swap = A(:,1) > A(:,2);
if any(swap)
    A(swap,:) = [A(swap,2) A(swap,1)];
end

a = max(A(:,1), L);
b = min(A(:,2), R);

keep = (b >= a);
if any(keep)
    allow2 = [a(keep) b(keep)];
    allow2 = sortrows(allow2, [1 2]);
else
    allow2 = zeros(0,2, cls);
end
end


function r0 = sample_from_intervals(allow)
if isempty(allow)
    error('sample_from_intervals: empty allowed set');
end

if size(allow,2) ~= 2
    error('sample_from_intervals: allow must be n x 2');
end

allow = int64(allow);
swap = allow(:,1) > allow(:,2);
if any(swap)
    allow(swap,:) = [allow(swap,2) allow(swap,1)];
end

lens = double(allow(:,2) - allow(:,1)) + 1;
keep = isfinite(lens) & (lens > 0);
allow = allow(keep,:);
lens  = lens(keep);

if isempty(allow)
    error('sample_from_intervals: no valid intervals after cleaning');
end

[allow, ord] = sortrows(allow,1);
lens = lens(ord);

total = sum(lens);
if ~(isfinite(total) && total > 0)
    error('sample_from_intervals: total lattice size is nonpositive');
end
u = 1 + floor(rand * total);

cum  = cumsum(lens);
idx  = find(u <= cum, 1, 'first');
prev = 0;
if idx > 1, prev = cum(idx-1); end

offset = u - prev;
r0     = allow(idx,1) + int64(offset - 1);
end


function info = getGlobalStreamInfo()
s = RandStream.getGlobalStream();
info = struct();
info.Type   = s.Type;
info.Seed   = s.Seed;
info.Substream = s.Substream;
end


function parsave(fname, S)
% make sure filename is a text scalar (char)
if isstring(fname)
    if numel(fname) ~= 1
        error('parsave:BadName','Filename must be a single string scalar.');
    end
    fname = char(fname);
elseif ~ischar(fname)
    error('parsave:BadName','Filename must be char or string scalar.');
end
% ensure directory exists
fdir = fileparts(fname);
if ~isempty(fdir) && ~exist(fdir,'dir')
    mkdir(fdir);
end

% save struct
if isstruct(S)
    save(fname,'-struct','S','-v7.3');
else
    save(fname,'S','-v7.3');
end
end


function [drop_cells, kept_cnt, total_cnt] = build_ufo_drop_mask(UFOdat, min_f_hz, min_d_ms)
if ~isempty(min_f_hz)
    if ~isscalar(min_f_hz) || ~isfinite(min_f_hz) || min_f_hz < 0
        error('build_ufo_drop_mask: min_f_hz must be a non-negative finite scalar or [].');
    end
end
if ~isempty(min_d_ms)
    if ~isscalar(min_d_ms) || ~isfinite(min_d_ms) || min_d_ms < 0
        error('build_ufo_drop_mask: min_d_ms must be a non-negative finite scalar or [].');
    end
end

M = numel(UFOdat);
drop_cells = cell(1, M);
kept_cnt   = 0;
total_cnt  = 0;

for ii = 1:M
    U = UFOdat{ii};
    if isempty(U)
        drop_cells{ii} = false(0,1);
        continue;
    end

    if size(U,2) < 3
        error('build_ufo_drop_mask: UFOdat{%d} must be n x 3.', ii);
    end

    n = size(U,1);
    total_cnt = total_cnt + n;

    bad = U(:,2) < U(:,1);
    if any(bad)
        U(bad, [1 2]) = U(bad, [2 1]);
    end

    f = double(U(:,3));
    dur_ms = double(U(:,2) - U(:,1)) / 1e3;

    keep = true(n,1);
    if ~isempty(min_f_hz)
        keep = keep & (f >= min_f_hz);
    end
    if ~isempty(min_d_ms)
        keep = keep & (dur_ms >= min_d_ms);
    end

    keep = keep & isfinite(f) & isfinite(dur_ms);

    drop_cells{ii} = ~keep;
    kept_cnt = kept_cnt + nnz(keep);
end
end


function G = build_fdr_groups(mima, szEmp, mode)
% Build FDR grouping matrix for global correction across all micro-macro pairs.
% Only 'global' mode is supported (all valid connections grouped together).
% 
% Inputs:
%   mima: logical matrix (N x M) indicating valid micro-macro connections
%   szEmp: size vector of empirical scores array [M x N x P ...]
%   mode: ignored (kept for compatibility, always uses 'global')
%
% Output:
%   G: grouping matrix same size as empirical scores, with 1 for valid connections, 0 otherwise

M = size(mima, 2);
N = size(mima, 1);
if ~(numel(szEmp) >= 2 && szEmp(1)==M && szEmp(2)==N)
    error('build_fdr_groups: size mismatch. Expected emp size (%d x %d x ...), got %s.', ...
        M, N, mat2str(szEmp));
end
if numel(szEmp) < 3, P = 1; else, P = szEmp(3); end

connMN = logical(mima.');
if P == 1
    conn3 = connMN;
else
    conn3 = repmat(connMN, [1 1 P]);
end

G = zeros(szEmp, 'like', zeros(1));
G(conn3) = 1;  % Global grouping: all valid connections in one group
end


function OUT = write_metric_diagnostics(emp, nulls, alpha, name, topN, did_mean, is_primary, groups)
if nargin < 8, groups = []; end
if nargin < 5 || isempty(topN), topN = 25; end

[p_map, sig_mask, crit_p, p_needed, per_group] = compare_to_nulls_fdr(emp, nulls, alpha, groups);

K = size(nulls,4);
Keff_map   = sum(isfinite(nulls), 4);
valid_emp  = isfinite(emp);
valid_test = valid_emp & (Keff_map > 0);

p_floor_nominal = 1 / (K + 1);
if any(valid_test(:))
    eff_floor = median( 1 ./ (Keff_map(valid_test) + 1), 'omitnan' );
else
    eff_floor = NaN;
end

if any(isfinite(p_map(:)))
    plist = p_map; plist(~isfinite(plist)) = Inf;
    [min_p, lin_min] = min(plist(:));
    Keff_min = Keff_map(lin_min);
    if isfinite(Keff_min) && Keff_min > 0
        min_b = max(0, ceil(min_p * (Keff_min + 1) - 1 - 1e-12));
    else
        min_b = NaN;
    end
else
    min_p = NaN;
    min_b = NaN;
end

m_test = nnz(valid_test);
n_sig  = nnz(sig_mask);

if isfinite(p_needed) && p_needed > 0
    K_needed_floor = max(0, ceil(1/p_needed) - 1);
else
    K_needed_floor = Inf;
end
worth_more_nulls = isfinite(p_floor_nominal) && isfinite(p_needed) && (p_floor_nominal > p_needed + 1e-15);

fprintf('\n[%s] Diagnostics:\n', upper(name));
fprintf(['  K = %d nulls | MC floor nominal = %.4g (median eff over tested = %.4g) | ', ...
    'BH crit_p = %.4g | min_p = %.4g (min_b = %s)\n'], ...
    K, p_floor_nominal, eff_floor, crit_p, min_p, ...
    ternary(isfinite(min_b), sprintf('%d', min_b), 'NaN'));
fprintf('  Significant voxels = %d / %d (%.2f%%) at FDR q = %.3f\n', ...
    n_sig, m_test, 100*n_sig/max(1,m_test), alpha);
fprintf('  Worth generating more nulls? %s\n', string(worth_more_nulls));

if ~isempty(per_group)
    fprintf('  Grouped FDR: %d group(s)\n', numel(per_group));
    fprintf('      %5s %8s %10s %8s\n','group','m','crit_p','sig');
    for t = 1:numel(per_group)
        Gt = per_group(t);
        fprintf('      %5g %8d %10.3g %8d\n', Gt.gid, Gt.m, Gt.crit_p, Gt.n_sig);
    end
end

[~, order] = sort(p_map(:), 'ascend', 'MissingPlacement','last');
order = order(isfinite(p_map(order)));
order = order(1:min(topN, numel(order)));
[I,J,P] = ind2sub(size(emp), order);

have_groups = ~isempty(groups) && isequal(size(groups), size(emp));

fprintf('\n  Top %d entries:\n', numel(order));
if have_groups
    fprintf('%5s %5s %5s | %10s %10s %10s %10s %10s %10s %6s %6s\n', ...
        'i','j','p','p_value','emp_value','null_mean','null_std','null_median','null_min','null_max','Keff','gid');
else
    fprintf('%5s %5s %5s | %10s %10s %10s %10s %10s %10s %6s\n', ...
        'i','j','p','p_value','emp_value','null_mean','null_std','null_median','null_min','null_max','Keff');
end

for t = 1:numel(order)
    ii = I(t); jj = J(t); pp = P(t);
    ni = squeeze(nulls(ii,jj,pp,:));
    ok = isfinite(ni);
    Keff_here = sum(ok);
    if Keff_here > 0
        mu  = mean(ni(ok));
        sd  = std(ni(ok));
        med = median(ni(ok));
        mn  = min(ni(ok));
        mx  = max(ni(ok));
    else
        mu = NaN; sd = NaN; med = NaN; mn = NaN; mx = NaN;
    end

    if have_groups
        gid = groups(ii,jj,pp);
        fprintf('%5d %5d %5d | %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %6d %6g\n', ...
            ii, jj, pp, ...
            p_map(order(t)), emp(order(t)), ...
            mu, sd, med, mn, mx, Keff_here, gid);
    else
        fprintf('%5d %5d %5d | %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %6d\n', ...
            ii, jj, pp, ...
            p_map(order(t)), emp(order(t)), ...
            mu, sd, med, mn, mx, Keff_here);
    end
end

OUT = struct();
OUT.metric = name;
OUT.primary = logical(is_primary);
OUT.did_mean = logical(did_mean);
OUT.K = K;
OUT.m = m_test;
OUT.alpha = alpha;
OUT.crit_p = crit_p;
OUT.p_needed = p_needed;
OUT.p_floor = p_floor_nominal;
OUT.min_p = min_p;
OUT.min_b = min_b;
OUT.K_needed_floor = K_needed_floor;
OUT.worth_more_nulls = logical(worth_more_nulls);
OUT.per_group = per_group;
OUT.n_sig = n_sig;
fprintf('\n'); drawnow;
end


function [p_map, sig_mask, crit_p, p_needed, per_group] = compare_to_nulls_fdr(orig, nulls, alpha, groups)
if nargin < 3 || isempty(alpha), alpha = 0.05; end
if nargin < 4, groups = []; end

dimK = ndims(nulls);
K    = size(nulls, dimK);

sizO = size(orig);
sizN = size(nulls);
if numel(sizO) < dimK-1
    sizO = [sizO, ones(1, (dimK-1) - numel(sizO))];
end
if any(sizO(1:dimK-1) ~= sizN(1:dimK-1))
    error('Size mismatch: orig %s vs nulls %s (leading dims).', mat2str(size(orig)), mat2str(size(nulls)));
end

valid   = isfinite(orig);
rep_sz  = ones(1, ndims(nulls));
rep_sz(dimK) = K;
validK  = repmat(valid, rep_sz);

finite_nulls = isfinite(nulls) & validK;
ge_mask      = (nulls >= orig) & finite_nulls;
b            = sum(ge_mask, dimK);
Keff         = sum(finite_nulls, dimK);

p_vals = (b + 1) ./ (Keff + 1);

p_map    = NaN(sizO);
sig_mask = false(sizO);

valid_p = valid & (Keff > 0);
p_map(valid_p) = p_vals(valid_p);

if ~isempty(groups)
    szg = size(groups);
    if numel(szg) < numel(sizO), szg = [szg, ones(1, numel(sizO)-numel(szg))]; end
    if ~isequal(szg(1:numel(sizO)), sizO)
        if numel(groups) == prod(sizO)
            groups = reshape(groups, sizO);
        else
            error('groups must have the same size as the metric (orig).');
        end
    end
end

if isempty(groups)
    p = p_map(valid_p);
    [p_sorted, ~] = sort(p(:));
    m = numel(p_sorted);
    if m == 0
        crit_p = 0;  p_needed = NaN;  sig_mask(:) = false;
        per_group = struct('gid',{},'m',{},'crit_p',{},'p_needed',{},'n_sig',{});
        return;
    end
    thresholds = (1:m)'/m * alpha;
    idx = find(p_sorted <= thresholds);
    p_needed  = alpha / m;
    if isempty(idx)
        crit_p = p_needed;  % Show what threshold would be needed (not 0, which is confusing)
        sig_mask(:) = false;
    else
        crit_p = p_sorted(max(idx));
        sig_mask(valid_p) = p_map(valid_p) <= crit_p;
    end
    per_group = struct('gid',{},'m',{},'crit_p',{},'p_needed',{},'n_sig',{});
else
    g = groups;
    g(~valid_p) = 0;
    gids = unique(g(:));
    gids = gids(isfinite(gids) & gids > 0);

    sig_mask(:) = false;
    per_group = struct('gid', {}, 'm', {}, 'crit_p', {}, 'p_needed', {}, 'n_sig', {});
    crits = []; pneed = [];

    for kk = 1:numel(gids)
        gid = gids(kk);
        gm  = (g == gid);
        p   = p_map(gm);  p = p(isfinite(p));
        m   = numel(p);
        if m == 0
            per_group(end+1) = struct('gid',gid,'m',0,'crit_p',0,'p_needed',NaN,'n_sig',0); %#ok<AGROW>
            continue;
        end
        [p_sorted, ~] = sort(p(:));
        thresholds = (1:m)'/m * alpha;
        idx = find(p_sorted <= thresholds);
        p_needed_g = alpha / m;
        if isempty(idx)
            crit_p_g = p_needed_g;  % Show what threshold would be needed (not 0, which is confusing)
            n_sig_g  = 0;
        else
            crit_p_g = p_sorted(max(idx));
            sig_mask(gm) = isfinite(p_map(gm)) & (p_map(gm) <= crit_p_g);
            n_sig_g      = nnz(sig_mask(gm));
        end
        per_group(end+1) = struct('gid',gid,'m',m,'crit_p',crit_p_g,'p_needed',p_needed_g,'n_sig',n_sig_g); %#ok<AGROW>
        crits(end+1) = crit_p_g; %#ok<AGROW>
        pneed(end+1) = p_needed_g; %#ok<AGROW>
    end
    if isempty(crits)
        crit_p = NaN;  % No groups had valid tests
    else
        crit_p = max(crits);  % Max across groups (not including 0)
    end
    p_needed = ternary(isempty(pneed), NaN, min(pneed));
end
end

function v = firstNonEmpty(x, fallback)
if isempty(x), v = fallback; else, v = x(1); end
end

function s = theory_spec_score(x, fs, f0_range_hz, BW)
% Max z-scored band power in f0±BW across the f0 grid (simple, fast).
x = double(x(:));
x = x - mean(x);
sd = std(x); if sd>0, x = x/sd; end

fL = f0_range_hz(1); fH = f0_range_hz(2);
f_grid = linspace(fL, fH, max(1, getin('CERD_THEORY_NFREQ',7)));

% Welch PSD
nfft = 2^nextpow2(max(256, numel(x)));
win  = round(0.250*fs); if win<64, win=64; end
nover= round(0.5*win);
[pxx,f] = pwelch(x, win, nover, nfft, fs);

bestz = -Inf;
for f0 = f_grid
    band = (f >= max(0,f0-BW)) & (f <= min(fs/2, f0+BW));
    ref  = (f >= max(0,f0-4*BW)) & (f <= min(fs/2, f0+4*BW)) & ~band; % local reference
    pb   = mean(pxx(band)); pr = mean(pxx(ref)); ps = std(pxx(ref));
    z = (pb - pr) / max(ps, eps);
    if z > bestz, bestz = z; end
end
s = max(0,bestz);
end



function [UFOdat_out, info] = limit_ufo_events(UFOdat_in, per_micro, total_cap, mode, deterministic, base_seed, Dat)
% Limit UFOs per micro and/or total per patient, optionally deterministically.
% mode: 'random' | 'first'. If 'first', keeps earliest-by-right-time.
% Deterministic sampling uses RandStream with seeds derived from base_seed and micro index.
UFOdat_out = UFOdat_in;
M = numel(UFOdat_in);

% normalize limits
if isempty(per_micro) || ~isfinite(double(firstNonEmpty(per_micro,NaN))) || per_micro<=0, per_micro = []; end
if isempty(total_cap) || ~isfinite(double(firstNonEmpty(total_cap,NaN))) || total_cap<=0, total_cap = []; end
if ~(strcmpi(mode,'random') || strcmpi(mode,'first')), mode = 'random'; end

n_before = cellfun(@(u) size(u,1), UFOdat_in);
dropped_pm = 0;

% — cap per micro —
if ~isempty(per_micro)
    for ii = 1:M
        U = UFOdat_out{ii};
        n = size(U,1);
        if n > per_micro
            if strcmpi(mode,'first')
                idx_keep = 1:per_micro;
            else
                seed = mod(double(base_seed) + 10007*ii + 3*sum(double(char(Dat.ID))), 2^31-1);
                if ~deterministic, seed = mod(seed + round(1e6*now), 2^31-1); end
                rs = RandStream('mt19937ar','Seed', max(1, seed));
                idx = randperm(rs, n);
                idx_keep = sort(idx(1:per_micro)); % keep chronological order (input already sorted)
            end
            UFOdat_out{ii} = U(idx_keep, :);
            dropped_pm = dropped_pm + (n - numel(idx_keep));
        end
    end
end

n_after_pm = cellfun(@(u) size(u,1), UFOdat_out);
dropped_total_extra = 0;

% — cap total across the patient —
if ~isempty(total_cap)
    T = sum(n_after_pm);
    if T > total_cap
        % collect (micro_idx, within_micro_idx, right_time)
        idx_trip = zeros(T,3);
        pos = 1;
        for ii = 1:M
            U = UFOdat_out{ii};
            n = size(U,1);
            if n==0, continue; end
            idx_trip(pos:pos+n-1, :) = [repmat(ii,n,1), (1:n)', U(:,2)];
            pos = pos + n;
        end
        idx_trip = idx_trip(1:pos-1, :);

        if strcmpi(mode,'first')
            [~, ord] = sort(idx_trip(:,3), 'ascend');
            keep_lin = ord(1:total_cap);
        else
            seed = mod(double(base_seed) + 777 + 7*sum(double(char(Dat.ID))), 2^31-1);
            if ~deterministic, seed = mod(seed + round(1e6*now), 2^31-1); end
            rs = RandStream('mt19937ar','Seed', max(1, seed));
            rp = randperm(rs, size(idx_trip,1));
            keep_lin = rp(1:total_cap);
        end

        keep_mask = false(size(idx_trip,1),1);
        keep_mask(keep_lin) = true;

        UFOdat_new = cell(1,M);
        for ii = 1:M
            U = UFOdat_out{ii}; n = size(U,1);
            if n==0, UFOdat_new{ii}=U; continue; end
            rows = find(idx_trip(:,1)==ii & keep_mask);
            idx_keep_micro = idx_trip(rows, 2);
            % preserve chronological order by right_time
            if ~isempty(idx_keep_micro)
                [~, o2] = sort(U(idx_keep_micro,2), 'ascend');
                idx_keep_micro = idx_keep_micro(o2);
                UFOdat_new{ii} = U(idx_keep_micro, :);
            else
                UFOdat_new{ii} = zeros(0,3, 'like', U);
            end
        end
        dropped_total_extra = size(idx_trip,1) - total_cap;
        UFOdat_out = UFOdat_new;
    end

end

n_after = cellfun(@(u) size(u,1), UFOdat_out);
info = struct();
info.before_per_micro    = n_before;
info.after_per_micro     = n_after;
info.dropped_per_micro   = dropped_pm;
info.dropped_total_extra = dropped_total_extra;
info.dropped_total       = sum(n_before) - sum(n_after);
end

function UFOdat_perm = permute_ufo_events(Dat, UFOdat, post_ufo_ms, deterministic, base_seed)
% Randomly reposition UFO events within the recording (control analysis)
if isempty(UFOdat)
    UFOdat_perm = UFOdat;
    return;
end

rng_state = rng;
if deterministic
    rng(double(base_seed) + 12345, 'twister');
else
    rng('shuffle');
end

post_ufo_us = int64(post_ufo_ms * 1e3);
t0 = int64(Dat.t0);
t1 = int64(Dat.t1);
if ~isfinite(double(t0)) || ~isfinite(double(t1)) || t1 <= t0
    warning('permute_ufo_events: invalid recording bounds; returning original UFOdat.');
    UFOdat_perm = UFOdat;
    rng(rng_state);
    return;
end

UFOdat_perm = cell(size(UFOdat));
for ii = 1:numel(UFOdat)
    U = UFOdat{ii};
    if isempty(U)
        UFOdat_perm{ii} = U;
        continue;
    end
    newU = zeros(size(U));
    for kk = 1:size(U,1)
        dur = int64(U(kk,2) - U(kk,1));
        freq = U(kk,3);
        max_start = double(t1 - post_ufo_us - dur);
        min_start = double(t0);
        if max_start <= min_start
            new_start = t0;
        else
            offset = int64(floor(rand() * (max_start - min_start + 1)));
            new_start = t0 + offset;
        end
        new_end = new_start + dur;
        if new_end + post_ufo_us > t1
            new_end = t1 - post_ufo_us;
            new_start = max(t0, new_end - dur);
        end
        newU(kk,1) = double(new_start);
        newU(kk,2) = double(new_end);
        newU(kk,3) = freq;
    end
    UFOdat_perm{ii} = newU;
end

rng(rng_state);
end


function [t0, t1] = time_bounds_from_ufo_struct(U)
% Robustly compute [min left, max right] across UFO.json
t0 = int64(0); t1 = int64(0);
if isempty(U) || ~isstruct(U), return; end
leftField  = local_pick_field(U, {'uutc_left','left','t0','start','start_uutc'});
rightField = local_pick_field(U, {'uutc_right','right','t1','end','end_uutc'});
if isempty(leftField) || isempty(rightField), return; end
L = []; R = [];
for k = 1:numel(U)
    a = U(k).(leftField);
    b = U(k).(rightField);
    if isfinite(double(a)) && isfinite(double(b)) && double(b) > double(a)
        L(end+1) = double(a); %#ok
        R(end+1) = double(b); %#ok
    end
end
if ~isempty(L)
    t0 = int64(min(L));
    t1 = int64(max(R));
end
end

function macro = synthetic_macro(fs, duration_ms)
% Generate white Gaussian placeholder macro (so later code still runs)
n_samples = max(1, round((duration_ms/1000) * fs));
macro = randn(n_samples,1) * 1e-6; % arbitrary µV-scale noise
end

function print_concise_summary(items, alpha, primary_label)
% Human-readable one-liners per metric, placed at the very end of each patient printout.
% items: {{'LABEL', diagStruct}, ...}; diagStruct has fields K, m, n_sig, min_p, crit_p, p_floor, worth_more_nulls, K_needed_floor
labels   = cellfun(@(c)c{1}, items, 'uni', false);
isPrimIx = find(strcmp(labels, primary_label), 1);
if isempty(isPrimIx), isPrimIx = 1; end

fprintf('\n--- END-OF-PATIENT SUMMARY (concise) ---\n');
fprintf('At FDR q=%.3f (GLOBAL). Primary metric: %s.\n', alpha, labels{isPrimIx});

anySig = false;
for t=1:numel(items), anySig = anySig || (items{t}{2}.n_sig > 0); end

for t=1:numel(items)
    L = items{t}{1};
    D = items{t}{2};
    if D.m == 0
        fprintf('• %-9s — none tested (m=0).\n', L);
        continue;
    end
    if D.n_sig > 0
        fprintf('• %-9s — %d/%d significant; BH crit_p=%s; min_p=%s; K=%d.\n', ...
            L, D.n_sig, D.m, fmt_p(D.crit_p), fmt_p(D.min_p), D.K);
    else
        needs = (isfield(D,'worth_more_nulls') && D.worth_more_nulls);
        needK = (isfield(D,'K_needed_floor') && isfinite(D.K_needed_floor));
        addK  = '';
        if needs && needK
            addK = sprintf(' Need K≥%d to beat floor.', D.K_needed_floor);
        end
        fprintf('• %-9s — 0/%d significant; min_p=%s, crit_p=%s; MC floor=%s (K=%d).%s\n', ...
            L, D.m, fmt_p(D.min_p), fmt_p(D.crit_p), fmt_p(D.p_floor), D.K, addK);
    end
end

if ~anySig
    % Highlight why nothing popped on the primary metric
    Dp = items{isPrimIx}{2};
    if Dp.m==0
        why = 'no valid tests (m=0).';
    elseif ~isfinite(Dp.min_p) || ~isfinite(Dp.crit_p)
        why = 'insufficient data to compute p or threshold.';
    elseif Dp.min_p > max(Dp.crit_p, realmin)
        why = 'min_p above BH threshold.';
    else
        why = 'conservative FDR threshold or coarse MC floor.';
    end
    fprintf('Overall: no discoveries. Primary (%s): %s\n', labels{isPrimIx}, why);
else
    sigLabs = labels(cellfun(@(c) c{2}.n_sig>0, items));
    fprintf('Overall: discoveries in %s.\n', strjoin(sigLabs, ', '));
end
end

function s = fmt_p(x)
% Pretty 3-sig-fig p formatting with safe NaN/Inf handling
if ~isfinite(x), s = 'NaN'; return; end
if x >= 0.001, s = sprintf('%.3f', x); else, s = sprintf('%.3g', x); end
end

function sz = pad_to_ndims(sz, n)
if numel(sz) < n
    sz = [sz, ones(1, n-numel(sz))];
end
end

function sig_ufos = extract_significant_ufos(sig_mask, UFOdat, Dat, theory, p_map, mima_use)
% Extract details about significant UFO events
% Returns struct array with: micro_idx, macro_idx, macro_name, event_idx, 
%   duration_ms, frequency_hz, theory_score, p_value, start_time, end_time

sig_ufos = struct('micro_idx',{}, 'macro_idx',{}, 'macro_name',{}, 'event_idx',{}, ...
    'duration_ms',{}, 'frequency_hz',{}, 'theory_score',{}, 'p_value',{}, ...
    'start_time',{}, 'end_time',{}, 'micro_channel',{});

[M, N, P] = size(sig_mask);

for ii = 1:M
    UFO = UFOdat{ii};
    if isempty(UFO), continue; end
    micro_ch = Dat.allmi{ii};
    
    for jj = 1:N
        if nargin >= 6 && ~isempty(mima_use)
            if ~mima_use(jj,ii), continue; end
        else
            if ~Dat.mima(jj,ii), continue; end
        end
        macro_name = Dat.ma{jj}{1};
        if iscell(macro_name), macro_name = macro_name{1}; end
        
        for kk = 1:min(P, size(UFO,1))
            if sig_mask(ii, jj, kk)
                ufo_event = UFO(kk, :);
                start_time = ufo_event(1);
                end_time = ufo_event(2);
                frequency = ufo_event(3);
                duration_ms = double(end_time - start_time) / 1e3;
                
                sig_ufos(end+1).micro_idx = ii; %#ok<AGROW>
                sig_ufos(end).macro_idx = jj;
                sig_ufos(end).macro_name = char(macro_name);
                sig_ufos(end).event_idx = kk;
                sig_ufos(end).duration_ms = duration_ms;
                sig_ufos(end).frequency_hz = frequency;
                sig_ufos(end).theory_score = theory(ii, jj, kk);
                sig_ufos(end).p_value = p_map(ii, jj, kk);
                sig_ufos(end).start_time = start_time;
                sig_ufos(end).end_time = end_time;
                sig_ufos(end).micro_channel = char(micro_ch);
            end
        end
    end
end

% Sort by p-value (most significant first)
if ~isempty(sig_ufos)
    p_vals = [sig_ufos.p_value];
    [~, ord] = sort(p_vals, 'ascend');
    sig_ufos = sig_ufos(ord);
end
end

function fig = visualize_significant_ufo(Dat, sig_ufo, post_ufo_ms, save_dir)
% Visualize a single significant UFO with micro, raw macro, filtered macro + theory overlay
fig = [];
if nargin < 4, save_dir = []; end

macro_name = sig_ufo.macro_name;
if iscell(macro_name), macro_name = macro_name{1}; end
micro_name = sig_ufo.micro_channel;

metric_mode = lower(string(getin('CERD_THEORY_METRIC','greens')));
% Visualization for greens mode (power spectrum approach)

viz_pre_ms      = getin('CERD_VIZ_PRE_MS', 5);
spect_pre_ms    = double(getin('CERD_SPECT_PRE_MS', max(5, viz_pre_ms)));
spect_post_ms   = double(getin('CERD_SPECT_POST_MS', post_ufo_ms));
pre_ms          = max(viz_pre_ms, spect_pre_ms);
macro_window_ms = min(post_ufo_ms, max(getin('CERD_VIZ_MACRO_MS', 120), spect_post_ms));
micro_window_ms = max(sig_ufo.duration_ms + 10, getin('CERD_VIZ_MICRO_MS', 30));
% plot_xlim       = [-30 300];
plot_xlim       = [-sig_ufo.duration_ms 2*sig_ufo.duration_ms];
view_xlim       = plot_xlim;
scale_mode      = lower(string(getin('CERD_VIZ_SCALE_MODE', 'match'))); % 'match'|'peak'|'rms'|'plot'

win_start = max(int64(sig_ufo.start_time) - int64(pre_ms*1e3), int64(Dat.t0));
win_end   = min(int64(sig_ufo.start_time) + int64(max(macro_window_ms, micro_window_ms)*1e3), int64(Dat.t1));

try
    % ---- micro signal ----
    [~, micro] = readMef3(Dat.meta_path, Dat.password, micro_name, 'time', win_start, win_end);
    if isempty(micro) || any(isnan(micro),'all'), warning('No micro signal for %s', micro_name); return; end
    micro = double(micro(:));
    fs_micro = local_fs(Dat.fs, micro, win_start, win_end);
    t_micro_ms = (0:numel(micro)-1)'/fs_micro*1000 - double(sig_ufo.start_time - win_start)/1e3;
    highlight_mask = t_micro_ms >= 0 & t_micro_ms <= sig_ufo.duration_ms;

    % ---- macro signal (raw) ----
    [~, macro] = readMef3(Dat.meta_path, Dat.password, macro_name, 'time', win_start, win_end);
    if isempty(macro) || any(isnan(macro),'all'), warning('No macro signal for %s', macro_name); return; end
    macro = double(macro(:));
    fs_macro = local_fs(Dat.fs, macro, win_start, win_end);
    t_macro_ms = (0:numel(macro)-1)'/fs_macro*1000 - double(sig_ufo.start_time - win_start)/1e3;

    % ---- filtered macro used by matcher ----
    macro_proc = macro - mean(macro);
    if std(macro_proc)>0, macro_proc = macro_proc/std(macro_proc); end
    if getin('CERD_THEORY_PREWHITEN',true)
        macro_proc = prewhiten_ar1(macro_proc);
    end
    BW = getin('CERD_THEORY_BW_HZ',15);

    [template_ms_raw, template_sig, meta] = best_theory_template(macro, fs_macro);
    if isempty(template_sig), warning('Theory template unavailable for %s', macro_name); return; end

    % Diagnostic: Verify template parameters and check oscillatory pattern
    omega = 2*pi*meta.f0_hz;
    omega_echo = sqrt(max(0, omega^2 - meta.gamma_1ps^2));
    f_echo_calc = omega_echo / (2*pi);
    period_ms = 1000 / f_echo_calc;  % Period in milliseconds
    fprintf('[TEMPLATE DIAG] %s→%s: Template params: f0=%.2f Hz, γ=%.1f s⁻¹, delay=%.2f ms, f_echo=%.2f Hz, period=%.2f ms\n', ...
        micro_name, macro_name, meta.f0_hz, meta.gamma_1ps, meta.delay_ms, f_echo_calc, period_ms);
    
    % Check template oscillatory content: count zero-crossings in the gated window (first 40 ms after delay)
    gate_ms = getin('CERD_THEORY_GATE_MS', 40);
    delay_samples = round(meta.delay_ms / 1000 * fs_macro);
    gate_samples = round(gate_ms / 1000 * fs_macro);
    if delay_samples < numel(template_sig)
        % Only look at the gated portion (where template is non-zero)
        end_sample = min(delay_samples + gate_samples, numel(template_sig));
        template_gated = template_sig(delay_samples+1:end_sample);
        if numel(template_gated) > 2
            % Count zero-crossings relative to mean
            template_centered = template_gated - mean(template_gated);
            zero_crossings = sum(diff(sign(template_centered)) ~= 0);
            time_gated_ms = (numel(template_gated) / fs_macro) * 1000;
            if time_gated_ms > 0 && zero_crossings > 0
                apparent_freq_hz = (zero_crossings / 2) / (time_gated_ms / 1000);
                fprintf('  Template oscillatory pattern: %d zero-crossings in %.1f ms (gated window) → apparent freq=%.2f Hz (expected %.2f Hz)\n', ...
                    zero_crossings, time_gated_ms, apparent_freq_hz, f_echo_calc);
                if abs(apparent_freq_hz - f_echo_calc) > 10
                    warning('Template frequency mismatch: apparent (%.2f Hz) vs expected (%.2f Hz) - check damping/gate', ...
                        apparent_freq_hz, f_echo_calc);
                end
            end
        end
    end
    
    % Template time vector should match macro time vector exactly (same signal, same fs)
    % best_theory_template returns time relative to window start, but we need it relative to UFO onset
    % Use t_macro_ms directly since both are from the same macro signal
    template_ms = t_macro_ms;

    % For visualization, use a narrower band to show oscillations clearly
    % Use separate viz bandwidth if specified, otherwise use min(BW, f0*0.2) to ensure reasonable band
    BW_viz = getin('CERD_VIZ_BW_HZ', []);
    if isempty(BW_viz)
        % Default: use 15 Hz or 20% of f0, whichever is smaller, but at least 10 Hz
        BW_viz = min(BW, max(10, meta.f0_hz * 0.2));
    end

    % Use visualization bandwidth for filtering (narrower band shows oscillations better)
    f_low = max(10, meta.f0_hz - BW_viz);  % Minimum 10 Hz to avoid DC and very low frequencies
    f_high = meta.f0_hz + BW_viz;
    macro_bw_raw = bandlimit(macro_proc, fs_macro, f_low, f_high);
    macro_bw_raw = macro_bw_raw - mean(macro_bw_raw);
    std_bw = max(std(macro_bw_raw), eps);
    macro_bw = macro_bw_raw / std_bw;

    % Diagnostic: Check frequency content of filtered signal
    if numel(macro_bw) > 64  % Need sufficient samples for FFT
        N_fft = min(2^nextpow2(numel(macro_bw)), 4096);
        X_fft = fft(macro_bw, N_fft);
        f_fft = (0:N_fft-1)' * (fs_macro / N_fft);
        P_fft = abs(X_fft).^2;
        % Find dominant frequency (excluding DC)
        [~, max_idx] = max(P_fft(2:min(floor(N_fft/2), numel(P_fft))));
        dom_freq_hz = f_fft(max_idx + 1);
        % Calculate power in expected band (using actual filter band)
        band_mask = (f_fft >= f_low) & (f_fft <= f_high);
        power_in_band = sum(P_fft(band_mask));
        total_power = sum(P_fft(2:floor(N_fft/2))); % Exclude DC
        power_ratio = power_in_band / max(total_power, eps);
        % Check if signal looks oscillatory (high frequency content relative to low)
        low_freq_mask = (f_fft >= 1) & (f_fft < f_low);
        high_freq_mask = (f_fft >= f_low) & (f_fft <= f_high);
        low_power = sum(P_fft(low_freq_mask));
        high_power = sum(P_fft(high_freq_mask));
        osc_ratio = high_power / max(low_power, eps);
        fprintf('[FILTER DIAG] %s→%s: Filter band=[%.1f, %.1f] Hz (f0=%.1f, BW_viz=%.1f, BW_analysis=%.1f)\n', ...
            micro_name, macro_name, f_low, f_high, meta.f0_hz, BW_viz, BW);
        fprintf('  Dominant freq: %.2f Hz | Power in band: %.1f%% | Osc ratio (high/low): %.2f\n', ...
            dom_freq_hz, power_ratio*100, osc_ratio);
        if dom_freq_hz < f_low || dom_freq_hz > f_high
            warning('Filtered signal dominant frequency (%.2f Hz) outside filter band [%.1f, %.1f] Hz', ...
                dom_freq_hz, f_low, f_high);
        end
        if osc_ratio < 0.5
            warning('Filtered signal has low oscillatory content (osc_ratio=%.2f) - may appear smooth', osc_ratio);
        end
    end

    % For visualization, generate an ungated template to show full theoretical decay
    % The gated template (template_sig) is used for detection, but we want to see the full decay
    % Use the same time vector as the original template (convert from ms to seconds)
    t_template = template_ms_raw / 1000;  % Convert milliseconds to seconds
    template_ungated = theory_template_gated(t_template, meta.delay_ms/1000, 2*pi*meta.f0_hz, meta.gamma_1ps, Inf);
    
    % Diagnostic: Check template decay at different time points
    if numel(template_ungated) > 0
        t_check_ms = [50, 100, 150, 200];
        for t_check = t_check_ms
            idx = find(template_ms_raw >= t_check, 1);
            if ~isempty(idx) && idx <= numel(template_ungated)
                val = template_ungated(idx);
                fprintf('[TEMPLATE DECAY] At t=%.0f ms: template value = %.6f\n', t_check, val);
            end
        end
    end
    
    % Use ungated template for visualization (full decay), but keep same normalization as gated version
    % Normalize only using non-zero values (after delay) to avoid DC offset issues
    delay_idx = find(template_ms_raw >= meta.delay_ms, 1);
    if isempty(delay_idx), delay_idx = 1; end
    template_active = template_ungated(delay_idx:end);
    if numel(template_active) > 0
        template_mean = mean(template_active);
        template_std = std(template_active);
    else
        template_mean = mean(template_ungated);
        template_std = std(template_ungated);
    end
    template_norm = template_ungated - template_mean;
    template_norm = template_norm / max(template_std, eps);
    scale_coeff = (macro_bw' * template_norm) / (template_norm' * template_norm + eps);
    overlay_bw = scale_coeff * template_norm;

    % Display scaling
    overlay_disp_bw = overlay_bw;
    overlay_disp_raw = overlay_bw * std_bw + mean(macro);
    annotation_txt = 'Scale: match (correlation)';
    switch scale_mode
        case "peak"
            peak_macro = max(abs(macro_bw));
            peak_tpl   = max(abs(overlay_bw));
            if peak_macro>0 && peak_tpl>0
                overlay_disp_bw = overlay_bw * (peak_macro / peak_tpl);
                overlay_disp_raw = overlay_disp_bw * std_bw + mean(macro);
                annotation_txt = 'Scale: matched peak amplitude';
            end
        case "rms"
            rms_macro = sqrt(mean(macro_bw.^2));
            rms_tpl   = sqrt(mean(overlay_bw.^2));
            if rms_macro>0 && rms_tpl>0
                overlay_disp_bw = overlay_bw * (rms_macro / rms_tpl);
                overlay_disp_raw = overlay_disp_bw * std_bw + mean(macro);
                annotation_txt = 'Scale: matched RMS';
            end
        case "plot"
            % Scale to match macro mean absolute amplitude for presentation
            meanabs_macro = mean(abs(macro_bw));
            meanabs_tpl   = mean(abs(overlay_bw));
            if meanabs_macro>0 && meanabs_tpl>0
                overlay_disp_bw = overlay_bw * (meanabs_macro / meanabs_tpl);
                overlay_disp_raw = overlay_disp_bw * std_bw + mean(macro);
                annotation_txt = 'Scale: matched mean |amp|';
            end
        otherwise
            % keep correlation scale
    end

    % Mask template values for times before UFO onset (0 ms) - the echo can only occur after the UFO starts
    % The template represents the predicted echo, which should only be visible at or after t=0
    % Note: template_ms is already aligned with t_macro_ms (relative to UFO onset)
    template_before_ufo = template_ms < 0;
    if any(template_before_ufo)
        % Debug: check template values before masking
        n_before = sum(template_before_ufo);
        max_before = max(abs(overlay_disp_raw(template_before_ufo)));
        fprintf('[TEMPLATE MASK] Masking %d points before t=0 (max abs value before masking: %.4f)\n', n_before, max_before);
        overlay_disp_raw(template_before_ufo) = NaN;
        overlay_disp_bw(template_before_ufo) = NaN;
        % Debug: check template values after t=0
        template_after_ufo = template_ms >= 0;
        n_after = sum(template_after_ufo);
        max_after = max(abs(overlay_disp_raw(template_after_ufo)));
        fprintf('[TEMPLATE MASK] After t=0: %d points remain (max abs value: %.4f)\n', n_after, max_after);
    else
        fprintf('[TEMPLATE MASK] No points before t=0 to mask\n');
    end

    vis_mask = template_ms >= plot_xlim(1) & template_ms <= plot_xlim(2);
    vis_times = template_ms(vis_mask);
    if isempty(vis_times)
        vis_times = template_ms;
        vis_mask = true(size(template_ms));
    end

    % Align overlay amplitude with macro traces in visualization window
    % Only use non-NaN template values for scaling
    macro_at_template = interp1(t_macro_ms, macro, vis_times, 'linear', 'extrap');
    tpl_vis_nonnan = vis_mask & ~isnan(overlay_disp_raw);
    if any(tpl_vis_nonnan)
        tpl_amp = max(abs(overlay_disp_raw(tpl_vis_nonnan)));
    else
        tpl_amp = 0;
    end
    mac_amp = max(abs(macro_at_template));
    if mac_amp>0 && tpl_amp>0
        overlay_disp_raw = overlay_disp_raw * (mac_amp / tpl_amp);
    end
    % Match DC level in window for raw macro
    macro_mean = mean(macro_at_template);
    tpl_mean = mean(overlay_disp_raw(tpl_vis_nonnan), 'omitnan');
    if ~isnan(tpl_mean) && any(tpl_vis_nonnan)
        overlay_disp_raw = overlay_disp_raw - tpl_mean + macro_mean;
    end

    macro_bw_at_template = interp1(t_macro_ms, macro_bw, vis_times, 'linear', 'extrap');
    tpl_bw_vis_nonnan = vis_mask & ~isnan(overlay_disp_bw);
    if any(tpl_bw_vis_nonnan)
        tpl_bw_amp = max(abs(overlay_disp_bw(tpl_bw_vis_nonnan)));
    else
        tpl_bw_amp = 0;
    end
    mac_bw_amp = max(abs(macro_bw_at_template));
    if mac_bw_amp>0 && tpl_bw_amp>0
        overlay_disp_bw = overlay_disp_bw * (mac_bw_amp / tpl_bw_amp);
    end
    % Force zero mean alignment for filtered traces
    tpl_bw_mean = mean(overlay_disp_bw(tpl_bw_vis_nonnan), 'omitnan');
    if ~isnan(tpl_bw_mean) && any(tpl_bw_vis_nonnan)
        overlay_disp_bw = overlay_disp_bw - tpl_bw_mean + mean(macro_bw_at_template);
    end

    % PSD diagnostics (reuse spectral config)
    psd_cfg = struct();
    psd_cfg.pre_ms        = spect_pre_ms;
    psd_cfg.post_ms       = spect_post_ms;
    psd_cfg.band_hz       = getin('CERD_SPECT_BAND_HZ', []);
    psd_cfg.bandwidth_hz  = double(getin('CERD_SPECT_BANDWIDTH_HZ', 600));
    psd_cfg.window_ms     = double(getin('CERD_SPECT_WINDOW_MS', 10));
    psd_cfg.overlap       = double(getin('CERD_SPECT_OVERLAP', 0.6));
    psd_cfg.nfft          = double(getin('CERD_SPECT_NFFT', 4096));
    psd_cfg.eps           = double(getin('CERD_SPECT_EPS', 1e-12));
    psd_cfg.detrend       = logical(getin('CERD_SPECT_DETREND', true));
    psd_cfg.min_samples   = double(getin('CERD_SPECT_MIN_SAMPLES', 64));
    psd_cfg.fallback_band = double(getin('CERD_SPECT_FALLBACK_BAND', [50 1000]));
    psd_cfg.max_freq_hz   = double(getin('CERD_SPECT_MAX_FREQ_HZ', 6000));

    psd_diag = compute_psd_overlay(macro, fs_macro, psd_cfg, double(sig_ufo.frequency_hz));
    % Debug: Check PSD diagnostic validity
    if ~psd_diag.valid
        fprintf('[PSD DIAG] PSD computation failed for %s→%s: valid=%d\n', micro_name, macro_name, psd_diag.valid);
    elseif isempty(psd_diag.f) || isempty(psd_diag.P_base) || isempty(psd_diag.P_post)
        fprintf('[PSD DIAG] PSD data empty for %s→%s: f=%d, P_base=%d, P_post=%d\n', ...
            micro_name, macro_name, isempty(psd_diag.f), isempty(psd_diag.P_base), isempty(psd_diag.P_post));
    end
    fig_visible = ternary(isempty(save_dir), 'on', 'off');
    fig = figure('Visible', fig_visible, 'Position', [100, 100, 1400, 800]);
    tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact','Padding','compact');

    % Micro panel (subplot 1)
    ax1 = nexttile(tl);
    plot(t_micro_ms, micro, 'Color', [0.6 0.6 0.6], 'LineWidth', 1); hold on;
    if any(highlight_mask)
        plot(t_micro_ms(highlight_mask), micro(highlight_mask), 'Color', [0.85 0 0], 'LineWidth', 1.6);
    end
    xline(0, '--k');
    xlabel('Time from UFO onset (ms)');
    ylabel(sprintf('Micro %s (\muV)', micro_name));
    title(sprintf('Micro %s | Event %d | Duration %.2f ms | Freq %.1f Hz', micro_name, sig_ufo.event_idx, sig_ufo.duration_ms, sig_ufo.frequency_hz));
    grid on; box on;
    xlim(plot_xlim);

    % Macro panel (subplot 2) - macro only, no template
    zoom_mask = t_macro_ms >= plot_xlim(1) & t_macro_ms <= plot_xlim(2);
    ax2 = nexttile(tl);
    plot(t_macro_ms(zoom_mask), macro(zoom_mask), 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);
    xline(0, '--k');
    xlabel('Time from UFO onset (ms)');
    ylabel(sprintf('Macro %s (\muV)', macro_name));
    title(sprintf('Macro %s | Event %d', macro_name, sig_ufo.event_idx));
    grid on; box on;
    xlim(plot_xlim);
    
    % Set y-axis limits based on macro signal only
    macro_window = macro(zoom_mask);
    if ~isempty(macro_window) && any(isfinite(macro_window))
        pad_raw = max(1e-6, range(macro_window))*0.1;
        ylim(ax2, [min(macro_window)-pad_raw, max(macro_window)+pad_raw]);
    end
    
    % Raw correlation not computed for macro-only plot (template shown in subplot 3)
    raw_corr = NaN;
    
    % Diagnostic: Check if empirical signal frequency matches theoretical
    if numel(macro_window) > 64 && any(isfinite(macro_window))
        % Compute dominant frequency of empirical signal in the visible window
        N_fft_raw = min(2^nextpow2(numel(macro_window)), 4096);
        X_fft_raw = fft(macro_window - mean(macro_window), N_fft_raw);
        f_fft_raw = (0:N_fft_raw-1)' * (fs_macro / N_fft_raw);
        P_fft_raw = abs(X_fft_raw).^2;
        % Find dominant frequency (excluding DC and very low frequencies)
        valid_range = (f_fft_raw >= 10) & (f_fft_raw <= 200);
        if any(valid_range)
            [~, max_idx_raw] = max(P_fft_raw(valid_range));
            dom_freq_empirical = f_fft_raw(valid_range);
            dom_freq_empirical = dom_freq_empirical(max_idx_raw);
            % Expected frequency from theory (echo frequency)
            omega = 2*pi*meta.f0_hz;
            omega_echo = sqrt(max(0, omega^2 - meta.gamma_1ps^2));
            f_echo_theory = omega_echo / (2*pi);
            fprintf('[ALIGN DIAG] %s→%s: Empirical dom freq=%.2f Hz | Theory f_echo=%.2f Hz (f0=%.1f, γ=%.1f) | Delay=%.2f ms\n', ...
                micro_name, macro_name, dom_freq_empirical, f_echo_theory, meta.f0_hz, meta.gamma_1ps, meta.delay_ms);
            if abs(dom_freq_empirical - f_echo_theory) > 10
                warning('Frequency mismatch: empirical (%.2f Hz) vs theory (%.2f Hz) - check alignment or parameters', ...
                    dom_freq_empirical, f_echo_theory);
            end
        end
    end

    % Filtered macro panel
    ax3 = nexttile(tl);
    bw_mask = t_macro_ms >= plot_xlim(1) & t_macro_ms <= plot_xlim(2);
    tpl_win_mask_bw = template_ms >= plot_xlim(1) & template_ms <= plot_xlim(2);
    if ~any(tpl_win_mask_bw), tpl_win_mask_bw = true(size(template_ms)); end
    plot(t_macro_ms(bw_mask), macro_bw(bw_mask), 'k-', 'LineWidth', 1.2); hold on;
    plot(template_ms, overlay_disp_bw, 'r-', 'LineWidth', 1.5);
    xline(0, '--k');
    xlabel('Time from UFO onset (ms)');
    ylabel('Normalized amplitude');
    title(sprintf('Macro %s (filtered) | Theory %.4f | p=%.4g', macro_name, sig_ufo.theory_score, sig_ufo.p_value));
    legend({'Filtered macro','Bessel template'}, 'Location','best');
    grid on; box on;
    text(0.02, 0.92, annotation_txt, 'Units','normalized','HorizontalAlignment','left', 'FontSize',10,'Color',[0.3 0.3 0.3]);
    text(0.99, 0.85, sprintf('Best f0=%.1f Hz | delay=%.2f ms | \gamma=%.1f s^{-1}', meta.f0_hz, meta.delay_ms, meta.gamma_1ps), ...
        'Units','normalized', 'HorizontalAlignment','right', 'FontSize', 10, 'Color', [0.35 0.35 0.35]);
    xlim(plot_xlim);
    % Don't force template min/max to match macro - preserve template's theoretical shape
    % Amplitude scaling is already done above (peak/RMS/correlation matching)
    macro_bw_window = macro_bw(bw_mask);
    tpl_window_bw = overlay_disp_bw(tpl_win_mask_bw);

    tpl_interp_bw = interp1(template_ms, overlay_disp_bw, t_macro_ms(bw_mask), 'linear', NaN);
    bw_valid = isfinite(macro_bw_window) & isfinite(tpl_interp_bw);
    bw_corr = NaN;
    if any(bw_valid)
        Cbw = corrcoef(macro_bw_window(bw_valid), tpl_interp_bw(bw_valid));
        if numel(Cbw) >= 4
            bw_corr = Cbw(1,2);
        end
    end

    filt_vals = [macro_bw_window; tpl_window_bw];
    if ~isempty(filt_vals) && any(isfinite(filt_vals))
        pad_bw = max(1e-6, range(filt_vals))*0.1;
        ylim(ax3, [min(filt_vals)-pad_bw, max(filt_vals)+pad_bw]);
    end

    fprintf(['[VIS] %s | %s→%s (event %d): raw corr=%.3f, filtered corr=%.3f, ' ...
        'theory score=%.3f, p=%.4g\n'], Dat.ID, micro_name, macro_name, sig_ufo.event_idx, ...
        raw_corr, bw_corr, sig_ufo.theory_score, sig_ufo.p_value);

    ax4 = nexttile(tl);
    % Check if PSD data is valid and has content
    has_psd_data = psd_diag.valid && ...
                   ~isempty(psd_diag.f) && numel(psd_diag.f) > 0 && ...
                   ~isempty(psd_diag.P_base) && numel(psd_diag.P_base) > 0 && ...
                   ~isempty(psd_diag.P_post) && numel(psd_diag.P_post) > 0 && ...
                   numel(psd_diag.f) == numel(psd_diag.P_base) && ...
                   numel(psd_diag.f) == numel(psd_diag.P_post);
    
    % Always print diagnostic info for debugging
    fprintf('[PSD DIAG] %s→%s: valid=%d, has_data=%d, f_size=%d, P_base_size=%d, P_post_size=%d\n', ...
        micro_name, macro_name, psd_diag.valid, has_psd_data, ...
        numel(psd_diag.f), numel(psd_diag.P_base), numel(psd_diag.P_post));
    
    if has_psd_data
        % Check if data has valid (finite) values
        base_valid = isfinite(psd_diag.P_base) & isfinite(psd_diag.f);
        post_valid = isfinite(psd_diag.P_post) & isfinite(psd_diag.f);
        
        % Set x-axis range first (before plotting) to determine what data to show
        % Use configurable PSD xlim if specified, otherwise default to 0-200 Hz
        psd_xlim = getin('CERD_VIZ_PSD_XLIM', []);
        if isempty(psd_xlim) || numel(psd_xlim) ~= 2
            % Default: 0-200 Hz to see analysis band clearly
            freq_xlim = min([200, psd_cfg.max_freq_hz, fs_macro/2]);
            psd_xlim = [0, freq_xlim];
            fprintf('[PSD DIAG] %s→%s: Using default xlim [0 %.1f] Hz\n', micro_name, macro_name, freq_xlim);
        else
            freq_xlim = psd_xlim(2);  % Update for use in theory line scaling
            fprintf('[PSD DIAG] %s→%s: Using user-specified xlim [%.1f %.1f] Hz\n', micro_name, macro_name, psd_xlim(1), psd_xlim(2));
        end
        
        if any(base_valid) && any(post_valid)
            % Filter data to only plot points within the xlim range (plus a small margin for visibility)
            f_margin = (psd_xlim(2) - psd_xlim(1)) * 0.1; % 10% margin
            f_plot_range = [max(0, psd_xlim(1) - f_margin), psd_xlim(2) + f_margin];
            f_in_range = (psd_diag.f >= f_plot_range(1)) & (psd_diag.f <= f_plot_range(2));
            
            base_plot = base_valid & f_in_range;
            post_plot = post_valid & f_in_range;
            
            if any(base_plot) && any(post_plot)
                % Plot the data within the visible range
                % Use markers if there are few points to ensure visibility
                n_points = sum(base_plot);
                if n_points < 20
                    % Sparse data - use markers and lines
                    h1 = plot(ax4, psd_diag.f(base_plot), 10*log10(psd_diag.P_base(base_plot) + psd_cfg.eps), ...
                        'Color',[0.5 0.5 0.5], 'LineWidth',1.1, 'Marker','o', 'MarkerSize',4, 'MarkerFaceColor',[0.5 0.5 0.5]); 
                    hold(ax4,'on');
                    h2 = plot(ax4, psd_diag.f(post_plot), 10*log10(psd_diag.P_post(post_plot) + psd_cfg.eps), ...
                        'Color',[0.1 0.4 0.8], 'LineWidth',1.3, 'Marker','s', 'MarkerSize',4, 'MarkerFaceColor',[0.1 0.4 0.8]);
                else
                    % Dense data - use lines only
                    h1 = plot(ax4, psd_diag.f(base_plot), 10*log10(psd_diag.P_base(base_plot) + psd_cfg.eps), ...
                        'Color',[0.5 0.5 0.5], 'LineWidth',1.1); 
                    hold(ax4,'on');
                    h2 = plot(ax4, psd_diag.f(post_plot), 10*log10(psd_diag.P_post(post_plot) + psd_cfg.eps), ...
                        'Color',[0.1 0.4 0.8], 'LineWidth',1.3);
                end
                
                % Verify plots were created and check their properties
                if isempty(h1) || isempty(h2) || ~isvalid(h1) || ~isvalid(h2)
                    fprintf('[PSD DIAG] %s→%s: Plot handles invalid - h1=%d, h2=%d\n', ...
                        micro_name, macro_name, ~isempty(h1) && isvalid(h1), ~isempty(h2) && isvalid(h2));
                else
                    % Check plot data
                    xdata1 = get(h1, 'XData');
                    ydata1 = get(h1, 'YData');
                    xdata2 = get(h2, 'XData');
                    ydata2 = get(h2, 'YData');
                    fprintf('[PSD DIAG] %s→%s: Plot handles valid - h1 has %d points (x=[%.1f %.1f], y=[%.2f %.2f]), h2 has %d points (x=[%.1f %.1f], y=[%.2f %.2f])\n', ...
                        micro_name, macro_name, numel(xdata1), min(xdata1), max(xdata1), min(ydata1), max(ydata1), ...
                        numel(xdata2), min(xdata2), max(xdata2), min(ydata2), max(ydata2));
                    % Force redraw
                    drawnow;
                end
                
                % Check data ranges
                f_range = [min(psd_diag.f(base_plot)), max(psd_diag.f(base_plot))];
                base_range = [min(10*log10(psd_diag.P_base(base_plot) + psd_cfg.eps)), max(10*log10(psd_diag.P_base(base_plot) + psd_cfg.eps))];
                post_range = [min(10*log10(psd_diag.P_post(post_plot) + psd_cfg.eps)), max(10*log10(psd_diag.P_post(post_plot) + psd_cfg.eps))];
                fprintf('[PSD DIAG] %s→%s: Plotted PSD with %d base points and %d post points (in range [%.1f %.1f] Hz)\n', ...
                    micro_name, macro_name, sum(base_plot), sum(post_plot), f_plot_range(1), f_plot_range(2));
                fprintf('[PSD DIAG] %s→%s: f_range=[%.1f %.1f] Hz, base_range=[%.2f %.2f] dB, post_range=[%.2f %.2f] dB\n', ...
                    micro_name, macro_name, f_range(1), f_range(2), base_range(1), base_range(2), post_range(1), post_range(2));
                
                % Set y-axis limits based on the actual plotted data
                y_data = [10*log10(psd_diag.P_base(base_plot) + psd_cfg.eps); 10*log10(psd_diag.P_post(post_plot) + psd_cfg.eps)];
                y_finite = y_data(isfinite(y_data));
                if ~isempty(y_finite)
                    y_pad = max(1, range(y_finite)) * 0.15; % 15% padding
                    y_auto = [min(y_finite) - y_pad, max(y_finite) + y_pad];
                    ylim(ax4, y_auto);
                    fprintf('[PSD DIAG] %s→%s: Set ylim to [%.2f %.2f] dB (based on plotted data)\n', micro_name, macro_name, y_auto(1), y_auto(2));
                else
                    fprintf('[PSD DIAG] %s→%s: No finite y-data in plotted range to set ylim\n', micro_name, macro_name);
                end
            else
                fprintf('[PSD DIAG] %s→%s: No data points in xlim range [%.1f %.1f] Hz - base_plot=%d, post_plot=%d\n', ...
                    micro_name, macro_name, f_plot_range(1), f_plot_range(2), sum(base_plot), sum(post_plot));
                has_psd_data = false; % Force else block
            end
        else
            fprintf('[PSD DIAG] %s→%s: PSD data has no finite values - base_valid=%d, post_valid=%d\n', ...
                micro_name, macro_name, sum(base_valid), sum(post_valid));
            has_psd_data = false; % Force else block
        end
        
        % Set x-axis limits (do this AFTER plotting to ensure data is visible)
        xlim(ax4, psd_xlim);
        
        % Force axes to update and ensure data is visible
        if has_psd_data && any(base_plot) && any(post_plot)
            % Verify axes children exist
            ax4_children = get(ax4, 'Children');
            fprintf('[PSD DIAG] %s→%s: ax4 has %d children after plotting\n', micro_name, macro_name, numel(ax4_children));
            drawnow;
        end
        
        % Compute theory line (but plot it AFTER shaded region so it's on top)
        theory_plotted = false;
        h_theory = [];
        theory_f_plot = [];
        theory_y_plot = [];
        theory_template_normalized = []; % Store normalized template for weighting
        
        % Get theory parameters
        f0_band = getin('CERD_THEORY_F0_HZ', [60 90]);
        f0 = mean(double(f0_band(:)));
        omega = 2*pi*max(1e-6, f0);
        gamma = double(getin('CERD_THEORY_GAMMA', 10));
        omega_echo = sqrt(max(0, omega^2 - gamma^2));
        
        % Determine frequency band for theory curve - MUST use macro analysis band, not micro UFO frequency!
        % Use the macro theory band (f0_band, typically 55-85 Hz), not psd_diag.f_band which may use micro frequency
        band = psd_cfg.band_hz;
        if isempty(band) || numel(band) < 2 || any(~isfinite(band))
            band = f0_band; % Use macro analysis band (e.g., 55-85 Hz)
        end
        
        % Ensure band is within visible range
        band(1) = max(0, band(1));
        band(2) = min(freq_xlim, max(band(1)+1, band(2)));
        
        % Create frequency grid from psd_diag.f that falls within the macro analysis band
        if isfield(psd_diag,'f') && ~isempty(psd_diag.f)
            f_mask = (psd_diag.f >= band(1)) & (psd_diag.f <= band(2));
            if any(f_mask)
                f_theory = psd_diag.f(f_mask);
            else
                % Create a frequency grid in the macro band
                f_theory = linspace(band(1), band(2), 100)';
            end
        else
            % Create a frequency grid in the macro band
            f_theory = linspace(band(1), band(2), 100)';
        end
        
        % Ensure f_theory is always valid and within visible range
        if isempty(f_theory) || numel(f_theory) == 0 || max(f_theory) > freq_xlim
            f_theory = linspace(max(0, f0-20), min(freq_xlim, f0+20), 100)';
        end
        
        if ~isempty(f_theory) && numel(f_theory) > 0
            % Compute theoretical transfer function |H(ω')|²
            w_theory = 2*pi*f_theory(:);
            H2_theory = 1 ./ ((omega_echo^2 - w_theory.^2).^2 + (2*gamma.*w_theory).^2);
            H2_theory = real(H2_theory);
            H2_theory(~isfinite(H2_theory) | H2_theory < 0) = 0;
            
            % Scale theory template to show what portion of post-UFO PSD it detects
            % Use the actual power level in the analysis band (where detection happens)
            post_mask = (psd_diag.f >= min(f_theory)) & (psd_diag.f <= max(f_theory));
            if any(post_mask)
                post_band_vals = psd_diag.P_post(post_mask);
                post_vals_pos = post_band_vals(post_band_vals > 0);
                if ~isempty(post_vals_pos)
                    % Use median power in band (more representative than max)
                    post_scale = median(post_vals_pos);
                else
                    post_scale = [];
                end
            else
                post_scale = [];
            end
            
            % Fallback: use median of post-UFO in visible range
            if isempty(post_scale) || post_scale <= 0
                post_visible = psd_diag.P_post(psd_diag.f <= freq_xlim & psd_diag.P_post > 0);
                if ~isempty(post_visible)
                    post_scale = median(post_visible);
                else
                    post_scale = median(psd_diag.P_post(psd_diag.P_post > 0));
                end
            end
            
            H2_max = max(H2_theory(H2_theory > 0));
            if H2_max > 0
                % Store normalized template for weighting (unit energy)
                theory_template_normalized = H2_theory / H2_max;
                
                if ~isempty(post_scale) && post_scale > 0
                    % Scale template to match median power level in analysis band
                    H2_scaled = H2_theory * (post_scale / H2_max);
                else
                    % If can't scale, just normalize to unit max and use a reasonable dB level
                    H2_scaled = H2_theory / H2_max;
                    % Scale to be visible (around -5 to 5 dB)
                    H2_scaled = H2_scaled * 10^(5/10); % scale to ~5 dB peak
                end
                theory_y = 10*log10(H2_scaled + psd_cfg.eps);
                
                % Always try to plot - be very lenient with validation
                valid_f = (f_theory >= 0) & (f_theory <= freq_xlim);
                valid_y = isfinite(theory_y);
                both_valid = valid_f & valid_y;
                
                if any(both_valid)
                    % Store for plotting after shaded region
                    theory_f_plot = f_theory(both_valid);
                    theory_y_plot = theory_y(both_valid);
                    theory_template_normalized = theory_template_normalized(both_valid);
                    theory_plotted = true;
                elseif numel(f_theory) > 0
                    % Last resort: plot whatever we have, even if partially invalid
                    finite_mask = isfinite(theory_y) & isfinite(f_theory);
                    if any(finite_mask)
                        theory_f_plot = f_theory(finite_mask);
                        theory_y_plot = theory_y(finite_mask);
                        theory_template_normalized = theory_template_normalized(finite_mask);
                        theory_plotted = true;
                    end
                end
            end
        end
        
        % DEBUG: If theory wasn't plotted, try a simple test plot to verify mechanism
        if ~theory_plotted
            % Try plotting a simple line to verify plotting works
            try
                hold(ax4, 'on');
                % Plot a simple test line at a known location - make it VERY visible
                test_f = linspace(60, 90, 50);
                test_y = 5 * ones(size(test_f)); % 5 dB - should be clearly visible
                h_test = plot(ax4, test_f, test_y, 'Color',[1 0 0], 'LineWidth',4.0, 'LineStyle','--', ...
                    'Marker','s', 'MarkerSize',6, 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor',[1 0 0], ...
                    'DisplayName','Theory (test)');
                uistack(h_test, 'top');
                theory_plotted = true; % Mark as plotted so legend includes it
                drawnow;
                % Print warning that we're using test plot
                fprintf('[DEBUG] Theory computation failed - using test plot at 5 dB\n');
                warning('Theory computation failed - using test plot in PSD subplot');
            catch ME
                % If even test plot fails, something is wrong with axes
                fprintf('[DEBUG] Failed to plot test line: %s\n', ME.message);
                warning('Failed to plot theory line in PSD subplot: %s', ME.message);
            end
        end
        
        % Add shaded region for analysis band (after setting xlim so ylim is established)
        % BUT send it to back so it doesn't cover the data
        has_band = isfield(psd_diag,'f_band') && ~isempty(psd_diag.f_band) && numel(psd_diag.f_band) > 0;
        if has_band
            band_min = min(psd_diag.f_band);
            band_max = max(psd_diag.f_band);
            y_curr = ylim(ax4);
            if all(isfinite(y_curr))
                h_fill = fill(ax4, [band_min, band_max, band_max, band_min], [y_curr(1), y_curr(1), y_curr(2), y_curr(2)], ...
                    [0.9 0.95 1], 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off');
                % Send shaded region to back so data lines are on top
                uistack(h_fill, 'bottom');
                % Vertical lines marking band edges
                h_vline1 = xline(ax4, band_min, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'HandleVisibility','off');
                h_vline2 = xline(ax4, band_max, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'HandleVisibility','off');
                uistack([h_vline1, h_vline2], 'bottom');
            end
        end
        
        % NOW plot theory line on top of everything (after shaded region)
        % Also compute and show weighted PSD (what's actually being detected)
        if theory_plotted && ~isempty(theory_f_plot) && ~isempty(theory_y_plot)
            hold(ax4, 'on');
            % Make it VERY visible: bright red, thick, dashed, with markers
            h_theory = plot(ax4, theory_f_plot, theory_y_plot, 'Color',[1 0 0], 'LineWidth',3.5, 'LineStyle','--', ...
                'Marker','o', 'MarkerSize',4, 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor',[1 0 0], ...
                'DisplayName','Detection template');
            uistack(h_theory, 'top');
            
            % Compute weighted PSD data (will plot in separate subplot)
            weighted_post = [];
            weighted_base = [];
            template_on_psd = [];
            if isfield(psd_diag,'f') && ~isempty(psd_diag.f) && numel(psd_diag.f) > 0 && ...
               ~isempty(theory_template_normalized) && numel(theory_template_normalized) == numel(theory_f_plot)
                % Interpolate normalized template to PSD frequency grid
                template_on_psd = interp1(theory_f_plot, theory_template_normalized, psd_diag.f, 'linear', 0);
                
                % Weighted PSD = PSD * normalized template (what matched-filter detects)
                weighted_post = psd_diag.P_post .* template_on_psd;
                weighted_base = psd_diag.P_base .* template_on_psd;
            end
            
            drawnow;
            % Debug: print range of values
            fprintf('[DEBUG] Theory line plotted: f_range=[%.1f %.1f] Hz, y_range=[%.2f %.2f] dB, n_points=%d\n', ...
                min(theory_f_plot), max(theory_f_plot), min(theory_y_plot), max(theory_y_plot), numel(theory_f_plot));
        else
            fprintf('[DEBUG] Theory line NOT plotted: theory_plotted=%d, f_empty=%d, y_empty=%d\n', ...
                theory_plotted, isempty(theory_f_plot), isempty(theory_y_plot));
        end
        
        xlabel(ax4, 'Frequency (Hz)');
        ylabel(ax4, 'PSD (dB)');
        title(ax4, 'Raw Power Spectral Density');
        
        % Verify axes limits are set correctly
        xlim_curr = xlim(ax4);
        ylim_curr = ylim(ax4);
        fprintf('[PSD DIAG] %s→%s: Final axes limits - xlim=[%.1f %.1f] Hz, ylim=[%.2f %.2f] dB\n', ...
            micro_name, macro_name, xlim_curr(1), xlim_curr(2), ylim_curr(1), ylim_curr(2));
        
        % Set y-axis limits if specified (do this after all plots to override auto-scaling)
        % BUT only if the user-specified range actually makes sense for the data
        psd_ylim = getin('CERD_VIZ_PSD_YLIM', []);
        if ~isempty(psd_ylim) && numel(psd_ylim) == 2
            % Check if user-specified ylim would hide the data
            if has_psd_data && any(base_plot) && any(post_plot)
                y_data = [10*log10(psd_diag.P_base(base_plot) + psd_cfg.eps); 10*log10(psd_diag.P_post(post_plot) + psd_cfg.eps)];
                y_finite = y_data(isfinite(y_data));
                if ~isempty(y_finite)
                    y_data_min = min(y_finite);
                    y_data_max = max(y_finite);
                    % Only use user-specified ylim if it actually contains the data
                    if psd_ylim(1) <= y_data_min && psd_ylim(2) >= y_data_max
                        ylim(ax4, psd_ylim);
                        fprintf('[PSD DIAG] %s→%s: Using user-specified ylim [%.2f %.2f] dB (contains data range [%.2f %.2f])\n', ...
                            micro_name, macro_name, psd_ylim(1), psd_ylim(2), y_data_min, y_data_max);
                    else
                        fprintf('[PSD DIAG] %s→%s: User-specified ylim [%.2f %.2f] dB would hide data [%.2f %.2f] - keeping auto ylim\n', ...
                            micro_name, macro_name, psd_ylim(1), psd_ylim(2), y_data_min, y_data_max);
                    end
                end
            else
                % No data to check, just use user-specified ylim
                ylim(ax4, psd_ylim);
            end
        end
        
        if theory_plotted
            legend(ax4, {'Baseline (pre-UFO)','Post-UFO (empirical)','Detection template'}, ...
                'Location','best', 'FontSize',8);
        else
            legend(ax4, {'Baseline (pre-UFO)','Post-UFO (empirical)'}, 'Location','best');
        end
        
        % Subplot 5: Weighted PSD (what's actually detected)
        if ~isempty(weighted_post) && ~isempty(weighted_base) && ~isempty(template_on_psd)
            ax5 = nexttile(tl);
            template_mask = template_on_psd > 0.01 * max(template_on_psd);
            if any(template_mask)
                % Get transform type for accentuating differences
                % Options: 'log' (10*log10, default), 'sqrt' (square root), 'linear', or 'power' (with exponent)
                transform_type = getin('CERD_VIZ_WEIGHTED_PSD_TRANSFORM', 'power');
                transform_power = getin('CERD_VIZ_WEIGHTED_PSD_POWER', 2); % For 'power' transform
                
                % Apply transform to weighted PSDs
                w_base_vals = weighted_base(template_mask) + psd_cfg.eps;
                w_post_vals = weighted_post(template_mask) + psd_cfg.eps;
                
                switch lower(transform_type)
                    case 'sqrt'
                        % Square root transform: accentuates differences, compresses large values
                        y_base = sqrt(w_base_vals);
                        y_post = sqrt(w_post_vals);
                        ylabel_str = 'Weighted PSD (√)';
                    case 'linear'
                        % Linear scale: no transform
                        y_base = w_base_vals;
                        y_post = w_post_vals;
                        ylabel_str = 'Weighted PSD';
                    case 'power'
                        % Power transform with configurable exponent
                        y_base = w_base_vals .^ transform_power;
                        y_post = w_post_vals .^ transform_power;
                        ylabel_str = sprintf('Weighted PSD (^%.2f)', transform_power);
                    case 'log'
                        % Log transform (original)
                        y_base = 10*log10(w_base_vals);
                        y_post = 10*log10(w_post_vals);
                        ylabel_str = 'Weighted PSD (dB)';
                    otherwise
                        % Default to sqrt
                        y_base = sqrt(w_base_vals);
                        y_post = sqrt(w_post_vals);
                        ylabel_str = 'Weighted PSD (√)';
                end
                
                % Plot weighted PSDs with transform
                plot(ax5, psd_diag.f(template_mask), y_base, ...
                    'Color',[0.5 0.5 0.5], 'LineWidth',2.5, 'DisplayName','Weighted baseline');
                hold(ax5, 'on');
                plot(ax5, psd_diag.f(template_mask), y_post, ...
                    'Color',[0.1 0.4 0.8], 'LineWidth',2.5, 'DisplayName','Weighted post-UFO');
                
                % Add detection template as reference (scaled to show frequency band emphasis)
                if theory_plotted && ~isempty(theory_f_plot) && ~isempty(theory_template_normalized)
                    % Interpolate template to match the frequency grid used in weighted plot
                    template_ref = interp1(theory_f_plot, theory_template_normalized, psd_diag.f(template_mask), 'linear', 0);
                    % Scale template to be visible as a reference (normalize to max of weighted PSDs for visibility)
                    max_weighted = max([max(weighted_post(template_mask)), max(weighted_base(template_mask))]);
                    if max_weighted > 0
                        template_scaled = template_ref * max_weighted;
                        % Apply same transform to template
                        template_vals = template_scaled + psd_cfg.eps;
                        switch lower(transform_type)
                            case 'sqrt'
                                y_template = sqrt(template_vals);
                            case 'linear'
                                y_template = template_vals;
                            case 'power'
                                y_template = template_vals .^ transform_power;
                            case 'log'
                                y_template = 10*log10(template_vals);
                            otherwise
                                y_template = sqrt(template_vals);
                        end
                        plot(ax5, psd_diag.f(template_mask), y_template, ...
                            'Color',[1 0 0], 'LineWidth',2.0, 'LineStyle','--', ...
                            'Marker','o', 'MarkerSize',3, 'DisplayName','Detection template (ref)');
                    end
                end
                
                % Set x-axis to focus on analysis band
                f_visible = psd_diag.f(template_mask);
                xlim(ax5, [max(0, min(f_visible)-5), min(200, max(f_visible)+5)]);
                
                % Set y-axis to show the difference clearly
                y_vals = [y_base; y_post];
                y_vals = y_vals(isfinite(y_vals));
                if ~isempty(y_vals)
                    y_pad = max(0.5, range(y_vals) * 0.1);
                    ylim(ax5, [min(y_vals) - y_pad, max(y_vals) + y_pad]);
                end
                
                xlabel(ax5, 'Frequency (Hz)');
                ylabel(ax5, ylabel_str);
                title(ax5, sprintf('Weighted PSD (PSD × template) | Detection signal'));
                if theory_plotted && ~isempty(theory_f_plot)
                    legend(ax5, {'Weighted baseline','Weighted post-UFO','Detection template (ref)'}, 'Location','best', 'FontSize',9);
                else
                    legend(ax5, {'Weighted baseline','Weighted post-UFO'}, 'Location','best', 'FontSize',9);
                end
                grid(ax5, 'on');
                box(ax5, 'on');
                
                % Add shaded region for analysis band
                has_band = isfield(psd_diag,'f_band') && ~isempty(psd_diag.f_band) && numel(psd_diag.f_band) > 0;
                if has_band
                    band_min = min(psd_diag.f_band);
                    band_max = max(psd_diag.f_band);
                    y_curr = ylim(ax5);
                    if all(isfinite(y_curr))
                        fill(ax5, [band_min, band_max, band_max, band_min], [y_curr(1), y_curr(1), y_curr(2), y_curr(2)], ...
                            [0.9 0.95 1], 'EdgeColor','none', 'FaceAlpha',0.2, 'HandleVisibility','off');
                        xline(ax5, band_min, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'HandleVisibility','off');
                        xline(ax5, band_max, '--', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'HandleVisibility','off');
                    end
                end
                
                % Add annotation explaining what this shows
                text(ax5, 0.02, 0.95, sprintf('Detection = weighted post > weighted baseline\nBlue line should be above grey line'), ...
                    'Units','normalized', 'HorizontalAlignment','left', 'FontSize',9, ...
                    'Color',[0.2 0.2 0.2], 'BackgroundColor',[1 1 1 0.9], ...
                    'EdgeColor',[0.5 0.5 0.5], 'LineWidth',1, 'VerticalAlignment','top');
            end
        end
        grid(ax4,'on'); box(ax4,'on');
    else
        % PSD data invalid or empty - show message with diagnostic info
        fprintf('[PSD DIAG] PSD plot failed for %s→%s: valid=%d, f_empty=%d, P_base_empty=%d, P_post_empty=%d\n', ...
            micro_name, macro_name, psd_diag.valid, isempty(psd_diag.f), isempty(psd_diag.P_base), isempty(psd_diag.P_post));
        if ~isempty(psd_diag.f) && ~isempty(psd_diag.P_base) && ~isempty(psd_diag.P_post)
            fprintf('[PSD DIAG] Array sizes: f=%d, P_base=%d, P_post=%d\n', ...
                numel(psd_diag.f), numel(psd_diag.P_base), numel(psd_diag.P_post));
        end
        text(ax4, 0.5, 0.5, 'PSD unavailable (insufficient samples or invalid data)', 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle','FontSize',11,'Color',[0.4 0.4 0.4]);
        xlabel(ax4, 'Frequency (Hz)');
        ylabel(ax4, 'PSD (dB)');
        title(ax4, 'PSD vs Theory');
        grid(ax4,'on'); box(ax4,'on');
        xlim(ax4, [0, min(psd_cfg.max_freq_hz, fs_macro/2)]);
        ylim(ax4, [-20, 20]); % Set a default y-range so the text is visible
    end

    linkaxes([ax1, ax2, ax3],'x');

    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir'), mkdir(save_dir); end
        fname = sprintf('ufo_%s_macro_%s_evt_%02d.pdf', sig_ufo.micro_channel, macro_name, sig_ufo.event_idx);
        fname = regexprep(fname, '[^A-Za-z0-9_\.]', '_');
        fname = fullfile(save_dir, fname);
        print(fig, fname, '-dpdf', '-vector');
        fprintf('Saved: %s\n', fname);
    end

catch ME
    warning('Error visualizing UFO (%s->%s): %s', micro_name, macro_name, ME.message);
    if ~isempty(fig) && ishghandle(fig), close(fig); end
    fig = [];
end
end

function fig = visualize_significant_ufo_spectrum(Dat, sig_ufo, save_dir)
% Alternative visualization for spectral metric: micro trace + macro PSD bump.
fig = [];
if nargin < 3, save_dir = []; end

cfg.pre_ms      = double(getin('CERD_SPECT_PRE_MS', 10));
cfg.post_ms     = double(getin('CERD_SPECT_POST_MS', 20));
cfg.band_hz     = getin('CERD_SPECT_BAND_HZ', []);
cfg.bandwidth_hz= double(getin('CERD_SPECT_BANDWIDTH_HZ', 600));
cfg.window_ms   = double(getin('CERD_SPECT_WINDOW_MS', 5));
cfg.overlap     = double(getin('CERD_SPECT_OVERLAP', 0.5));
cfg.nfft        = double(getin('CERD_SPECT_NFFT', 1024));
cfg.eps         = double(getin('CERD_SPECT_EPS', 1e-12));
cfg.detrend     = logical(getin('CERD_SPECT_DETREND', true));
cfg.min_samples = double(getin('CERD_SPECT_MIN_SAMPLES', 64));
cfg.use_log_ratio = logical(getin('CERD_SPECT_USE_LOG_RATIO', true));
cfg.fallback_band = double(getin('CERD_SPECT_FALLBACK_BAND', [500 3500]));
cfg.max_freq_hz   = double(getin('CERD_SPECT_MAX_FREQ_HZ', 6000));
cfg.total_ms      = cfg.pre_ms + cfg.post_ms;

macro_name = char(sig_ufo.macro_name);
micro_name = char(sig_ufo.micro_channel);

spec_start = int64(sig_ufo.start_time) - int64(cfg.pre_ms * 1e3);
spec_end   = spec_start + int64(cfg.total_ms * 1e3);
if spec_start < int64(Dat.t0) || spec_end > int64(Dat.t1)
    warning('visualize_spectrum: window exceeds recording bounds for %s', Dat.ID);
    return;
end

try
    [~, micro] = readMef3(Dat.meta_path, Dat.password, micro_name, 'time', spec_start, spec_end);
    [~, macro] = readMef3(Dat.meta_path, Dat.password, macro_name, 'time', spec_start, spec_end);
catch ME
    warning('visualize_spectrum: failed to read MEF3 data (%s)', ME.message);
    return;
end

if isempty(micro) || isempty(macro)
    warning('visualize_spectrum: empty micro/macro data for %s', Dat.ID);
    return;
end
if any(isnan(micro),'all') || any(isnan(macro),'all')
    warning('visualize_spectrum: NaNs detected for %s', Dat.ID);
    return;
end

micro = double(micro(:));
macro = double(macro(:));

fs_micro = local_fs(Dat.fs, micro, spec_start, spec_end);
fs_macro = local_fs(Dat.fs, macro, spec_start, spec_end);

pre_samples  = max(1, round(cfg.pre_ms  * fs_macro / 1000));
post_samples = max(1, round(cfg.post_ms * fs_macro / 1000));
if pre_samples + post_samples > numel(macro)
    post_samples = numel(macro) - pre_samples;
end
if pre_samples < cfg.min_samples || post_samples < cfg.min_samples
    warning('visualize_spectrum: insufficient samples for PSD (pre=%d, post=%d).', pre_samples, post_samples);
    return;
end

baseline = macro(1:pre_samples);
segment  = macro(pre_samples+1:pre_samples+post_samples);

score = spectral_power_score(macro, fs_macro, cfg, double(sig_ufo.frequency_hz));
score_db = 10*score;

win = max(16, round(cfg.window_ms * fs_macro / 1000));
win = min(win, numel(segment));
if win >= numel(segment)
    win = max(8, floor(numel(segment)/2));
end
noverlap = min(win-1, max(0, round(cfg.overlap * win)));
nfft = max(win, round(cfg.nfft));
[P_post,f] = pwelch(segment, win, noverlap, nfft, fs_macro);
[P_base,~] = pwelch(baseline, win, noverlap, nfft, fs_macro);

band = cfg.band_hz;
center_freq = double(sig_ufo.frequency_hz);
if isempty(band) || numel(band)<2 || any(~isfinite(band))
    bw = cfg.bandwidth_hz;
    if ~isfinite(center_freq) || center_freq <= 0
        band = cfg.fallback_band;
    else
        if ~isfinite(bw) || bw <= 0, bw = max(100, center_freq*0.5); end
        half_bw = bw/2;
        band = [max(0, center_freq-half_bw), center_freq+half_bw];
    end
end
band(1) = max(0, band(1));
band(2) = min(cfg.max_freq_hz, max(band(1)+1, band(2)));
band_mask = f >= band(1) & f <= band(2);
if ~any(band_mask)
    [~, idx] = min(abs(f - mean(band)));
    band_mask = idx;
end
pow_base_band = mean(P_base(band_mask));
pow_post_band = mean(P_post(band_mask));
ratio_db = 10*log10((pow_post_band + cfg.eps)/(pow_base_band + cfg.eps));

t_rel_ms = (0:numel(macro)-1)'/fs_macro*1000 - cfg.pre_ms;
t_micro_ms = (0:numel(micro)-1)'/fs_micro*1000 - cfg.pre_ms;

fig = figure('Visible','off','Color','w','Position',[100 100 800 700]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

% Micro trace
nexttile(tl);
plot(t_micro_ms, micro, 'Color',[0.2 0.2 0.2], 'LineWidth',0.8); hold on;
mask = t_micro_ms >= 0 & t_micro_ms <= sig_ufo.duration_ms;
plot(t_micro_ms(mask), micro(mask), 'Color',[0.85 0.1 0.1], 'LineWidth',1.2);
xlabel('Time (ms)'); ylabel('Micro (\muV)');
title(sprintf('%s | Micro %s (Event %d)', Dat.ID, micro_name, sig_ufo.event_idx));
xline(0,'k:'); xline(sig_ufo.duration_ms,'k:');
grid on;

% Macro trace with baseline/post shading
nexttile(tl);
plot(t_rel_ms, macro, 'Color',[0.1 0.3 0.6], 'LineWidth',0.9); hold on;
yl = ylim;
patch([-cfg.pre_ms, 0, 0, -cfg.pre_ms], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.9], ...
    'FaceAlpha',0.3, 'EdgeColor','none');
patch([0, cfg.post_ms, cfg.post_ms, 0], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.95 1], ...
    'FaceAlpha',0.3, 'EdgeColor','none');
plot(t_rel_ms, macro, 'Color',[0.1 0.3 0.6], 'LineWidth',0.9);
ylim(yl);
xlabel('Time (ms)'); ylabel('Macro (\muV)');
title(sprintf('Macro %s | Score = %.2f dB', macro_name, score_db));
xline(0,'k:','LineWidth',1.0);
grid on;

% Power spectral density
xax = nexttile(tl);
plot(f, 10*log10(P_base + cfg.eps), 'Color',[0.5 0.5 0.5], 'LineWidth',1.1); hold on;
plot(f, 10*log10(P_post + cfg.eps), 'Color',[0.1 0.4 0.8], 'LineWidth',1.3);
yl_psd = ylim;
patch([band(1) band(2) band(2) band(1)], [yl_psd(1) yl_psd(1) yl_psd(2) yl_psd(2)], [0.8 0.9 1], ...
    'FaceAlpha',0.2, 'EdgeColor','none');
uistack(findobj(xax,'Type','patch'),'bottom');
xlabel('Frequency (Hz)'); ylabel('PSD (dB/\mathrm{Hz})');
legend({'Baseline','Post-UFO'},'Location','northeast');
title(sprintf('Band [%0.0f, %0.0f] Hz | Ratio = %.2f dB', band(1), band(2), ratio_db));
xlim([0 min(cfg.max_freq_hz, fs_macro/2)]);
grid on;

sgtitle(sprintf('%s | Micro %s → Macro %s | Event %d', Dat.ID, micro_name, macro_name, sig_ufo.event_idx), 'Interpreter','none');

if ~isempty(save_dir)
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    fname = sprintf('%s_%s_to_%s_ev%03d_spectrum.pdf', Dat.ID, micro_name, macro_name, sig_ufo.event_idx);
    fname = regexprep(fname,'[^A-Za-z0-9_\.]+','_');
    saveas(fig, fullfile(save_dir, fname));
end
end

function [t_ms, template, meta] = best_theory_template(x, fs)
% Compute best-fit theory template mirroring theory_match_score search
x = double(x(:));
if isempty(x)
    t_ms = []; template = []; meta = struct();
    return;
end

f0_range_hz    = getin('CERD_THEORY_F0_HZ',    [60 90]);
delay_ms_range = getin('CERD_THEORY_DELAY_MS', [2 10]);
gamma_1ps      = getin('CERD_THEORY_GAMMA',    80);
nfreq          = getin('CERD_THEORY_NFREQ',    7);
ndelay         = getin('CERD_THEORY_NDELAY',   9);
BW             = getin('CERD_THEORY_BW_HZ',    15);
Grng           = getin('CERD_THEORY_GAMMA_RANGE', [max(20,gamma_1ps/2) min(2000,2*gamma_1ps)]);
Ngamma         = max(1, round(getin('CERD_THEORY_NGAMMA',3)));
gate_ms        = getin('CERD_THEORY_GATE_MS',  40);
doPW           = getin('CERD_THEORY_PREWHITEN',true);

L = numel(x);
t = (0:L-1)'/fs;

x = x - mean(x);
if std(x)>0, x = x/std(x); end
if doPW, x = prewhiten_ar1(x); end

f_grid = linspace(f0_range_hz(1), f0_range_hz(2), max(1,nfreq));
d_grid = linspace(delay_ms_range(1)/1000, delay_ms_range(2)/1000, max(1,ndelay));
if numel(Grng)==1
    gamma_grid = Grng;
else
    gamma_grid = linspace(Grng(1), Grng(2), max(1,Ngamma));
end

best_score = -Inf; best_template = []; best_meta = struct('f0_hz',NaN,'delay_ms',NaN,'gamma_1ps',NaN);
for fi = 1:numel(f_grid)
    f0 = f_grid(fi);
    xl = bandlimit(x, fs, max(1,f0-BW), f0+BW);
    xl = xl - mean(xl);
    if std(xl)>0, xl = xl/std(xl); end
    for gi = 1:numel(gamma_grid)
        gamma = gamma_grid(gi);
        for di = 1:numel(d_grid)
            tau = d_grid(di);
            g = theory_template_gated(t, tau, 2*pi*f0, gamma, gate_ms);
            g = g - mean(g);
            sg = std(g); if sg<=0, continue; end
            g = g/sg;
            score = (xl.'*g)/max(eps, numel(xl));
            if score > best_score
                best_score = score;
                best_template = g;
                best_meta = struct('f0_hz', f0, 'delay_ms', tau*1000, 'gamma_1ps', gamma, 'score', score);
            end
        end
    end
end

t_ms = t*1000;
template = best_template;
meta = best_meta;
end
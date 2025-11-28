function cerd_batch
% ========================= CERD BATCH (all 8 patients) =========================
% Runs SCORE → GENERATE → ANALYZE for each patient in one go, non-interactive.
% - Opens the parallel pool once here; cerd_main won't touch it.
% - Avoids the interactive null-range prompt via CERD_NULL_RANGE.
% - Optional UFO caps (total / per-micro) with deterministic sampling.
% ==============================================================================

% ---- Patients ----
default_patients = {'pat001','pat003','pat004','pat061','pat072','pat079'};
if evalin('base','exist(''CERD_BATCH_PATIENTS'',''var'')')
    patientIDs = evalin('base','CERD_BATCH_PATIENTS');
else
    patientIDs = default_patients;
end

% ---- Steps ----
STEPS.score    = false;     % build/update score file
STEPS.generate = false;     % generate nulls (K = NNULL over NULL_RANGE_STR)
STEPS.analyze  = true;      % run FDR analysis and print diagnostics
if evalin('base','exist(''CERD_BATCH_STEPS'',''var'')')
    STEPS_override = evalin('base','CERD_BATCH_STEPS');
    if isfield(STEPS_override,'score'),    STEPS.score    = logical(STEPS_override.score);    end
    if isfield(STEPS_override,'generate'), STEPS.generate = logical(STEPS_override.generate); end
    if isfield(STEPS_override,'analyze'),  STEPS.analyze  = logical(STEPS_override.analyze);  end
end

% ---- Stats / FDR ----
if evalin('base','exist(''CERD_NNULL'',''var'')')
    NNULL = evalin('base','CERD_NNULL');
else
    NNULL = 10000;                     % total nulls to target
    assignin('base','CERD_NNULL', NNULL);
end
if evalin('base','exist(''CERD_NULL_RANGE'',''var'')')
    NULL_RANGE_STR = evalin('base','CERD_NULL_RANGE');
else
    NULL_RANGE_STR = sprintf('1:%d', NNULL);   % prevents prompt in 'g'
    assignin('base','CERD_NULL_RANGE', NULL_RANGE_STR);
end
if evalin('base','exist(''CERD_ALPHA'',''var'')')
    ALPHA = evalin('base','CERD_ALPHA');
else
    ALPHA = 0.05;
    assignin('base','CERD_ALPHA', ALPHA);
end
% FDR correction uses global grouping (hardcoded)
DETERMINISTIC       = true;
BASE_SEED           = 54833;

% ---- Spectral/theory knobs (safe defaults) ----
POST_UFO_MS     = 300;
if evalin('base','exist(''CERD_POST_UFO_MS'',''var'')')
    POST_UFO_MS = evalin('base','CERD_POST_UFO_MS');
end
if evalin('base','exist(''CERD_GAP_UFO_MS'',''var'')')
    GAP_UFO_MS = evalin('base','CERD_GAP_UFO_MS');
else
    GAP_UFO_MS = [30 300];
end
FREQ_BAND_HZ    = [60 300];
MIN_PROM_DB     = 3;
USE_MULTITAPER  = exist('pmtm','file')==2;

THEORY_ENABLE   = true;        % supported by your cerd_main
% Frequency range for macro echo detection
% Theory (Notes.txt Eq. 20): Echo frequency f_echo = √(ω²-γ²)/(2π) where ω = 2πf0
% For typical macroscale natural frequencies (60-90 Hz) and light damping (γ=10 1/s):
%   f_echo ≈ f0 (since γ << ω), so [60 90] Hz is the standard range
% Expanded range [100 1000] Hz may be used for:
%   - Higher natural frequencies at different scales
%   - Conservative search when exact frequency is unknown
%   - Exploring intermediate frequencies between micro UFOs (2-8 kHz) and low frequencies
if evalin('base','exist(''CERD_THEORY_F0_HZ'',''var'')')
    THEORY_F0_HZ = evalin('base','CERD_THEORY_F0_HZ');
else
    THEORY_F0_HZ    = [100 1000];  % Hz - expanded range (override with [60 90] for standard macroscale)
end
if evalin('base','exist(''CERD_THEORY_DELAY_MS'',''var'')')
    THEORY_DELAY_MS = evalin('base','CERD_THEORY_DELAY_MS');
else
    THEORY_DELAY_MS = [2 10];
end
if evalin('base','exist(''CERD_THEORY_GAMMA'',''var'')')
    THEORY_GAMMA = evalin('base','CERD_THEORY_GAMMA');
else
    THEORY_GAMMA    = 10;          % 1/s - very light damping (essentially undamped: echo freq ≈ omega)
end
THEORY_NFREQ    = 15;          % increased from 7 to better sample wider range
THEORY_NDELAY   = 9;
THEORY_BW_HZ    = 100;         % increased from 15 to match wider frequency range
THEORY_GAMMA_RANGE = [5 20];   % narrow range around very light damping
THEORY_NGAMMA   = 3;           % fewer gamma points since range is narrow
THEORY_GATE_MS  = 40;
THEORY_PREWHITEN= true;
if evalin('base','exist(''CERD_THEORY_METRIC'',''var'')')
    METRIC_MODE = evalin('base','CERD_THEORY_METRIC');
else
    METRIC_MODE = 'greens';  % Default: power spectrum approach (implements Notes.txt Eq. 26)
    assignin('base','CERD_THEORY_METRIC', METRIC_MODE);
end
fprintf('[batch] Theory metric: %s\n', upper(string(METRIC_MODE)));

if ~evalin('base','exist(''CERD_EVENT_PERMUTE'',''var'')')
    assignin('base','CERD_EVENT_PERMUTE', false);
end

% ---- Optional empirical UFO filters (empty = use all) ----
MIN_UFO_FREQ_HZ = [];          % e.g., 2000hz
MIN_UFO_DUR_MS  = [];          % e.g., 5ms

% ---- Paths / MEF copy (must match your layout) ----
drop_root  = '/Users/erik/Library/CloudStorage/Dropbox/Work/UFO';
local_root = '/Users/erik/Desktop/tmp_mef';

% ==============================================================================
% Parallel pool: only if we will GENERATE nulls
% ==============================================================================
% Check if caller wants to manage the pool (e.g., parameter sweep script)
% If not set, default to true (we manage it)
if evalin('base','exist(''CERD_MANAGE_POOL'',''var'')')
    caller_manages_pool = ~evalin('base','CERD_MANAGE_POOL');
else
    caller_manages_pool = false;  % Default: we manage it
end
assignin('base','CERD_MANAGE_POOL', false);   % cerd_main won't create/kill pools
if STEPS.generate
    p = gcp('nocreate');
    want = max(1, feature('numcores')-2);
    if isempty(p), parpool('local', want); end
end

% Common knobs (not patient-specific)
assignin('base','CERD_ALPHA',        ALPHA);
assignin('base','CERD_NNULL',        NNULL);
assignin('base','CERD_NULL_RANGE',   NULL_RANGE_STR);    % avoids prompt.  [oai_citation:2‡cerd_main.txt](sediment://file_00000000b5047246aa24aaf2fb39252c)

assignin('base','CERD_POST_UFO_MS',  POST_UFO_MS);
assignin('base','CERD_GAP_UFO_MS',   GAP_UFO_MS);
assignin('base','CERD_FREQ_BAND_HZ', FREQ_BAND_HZ);
assignin('base','CERD_MIN_PROM_DB',  MIN_PROM_DB);
assignin('base','CERD_MULTITAPER',   USE_MULTITAPER);

assignin('base','CERD_DETERMINISTIC',    DETERMINISTIC);
assignin('base','CERD_BASE_SEED',        BASE_SEED);

% Theory knobs used by cerd_main/theory_match_score
assignin('base','CERD_THEORY_ENABLE',     THEORY_ENABLE);
assignin('base','CERD_THEORY_F0_HZ',      THEORY_F0_HZ);
assignin('base','CERD_THEORY_DELAY_MS',   THEORY_DELAY_MS);
assignin('base','CERD_THEORY_GAMMA',      THEORY_GAMMA);
assignin('base','CERD_THEORY_NFREQ',      THEORY_NFREQ);
assignin('base','CERD_THEORY_NDELAY',     THEORY_NDELAY);
assignin('base','CERD_THEORY_BW_HZ',      THEORY_BW_HZ);
assignin('base','CERD_THEORY_GAMMA_RANGE',THEORY_GAMMA_RANGE);
assignin('base','CERD_THEORY_NGAMMA',     THEORY_NGAMMA);
assignin('base','CERD_THEORY_GATE_MS',    THEORY_GATE_MS);
assignin('base','CERD_THEORY_PREWHITEN',  THEORY_PREWHITEN);

assignin('base','CERD_FAST_META',false);
assignin('base','CERD_FAST_META_IDS',{});

% Optional empirical filters
if ~isempty(MIN_UFO_FREQ_HZ), assignin('base','CERD_MIN_UFO_FREQ_HZ', MIN_UFO_FREQ_HZ); end
if ~isempty(MIN_UFO_DUR_MS),  assignin('base','CERD_MIN_UFO_DUR_MS',  MIN_UFO_DUR_MS);  end

% ==============================================================================
% Loop patients
% ==============================================================================
rng(BASE_SEED,'twister'); warning('off','all'); addpath('/Users/erik/matmef');
allS = struct([]);
for pIx = 1:numel(patientIDs)
    pat_id = patientIDs{pIx};
    fprintf('\n==================== %s (%d/%d) ====================\n', pat_id, pIx, numel(patientIDs));

    % Ensure local MEF mirror; prep_basics prefers it automatically.  [oai_citation:3‡cerd_main.txt](sediment://file_00000000b5047246aa24aaf2fb39252c)
    ensure_local_mef(drop_root, local_root, pat_id);

    % Tell cerd_main which patient to operate on
    assignin('base','CERD_PATIENT', pat_id);

    % Steps
    if STEPS.score
        assignin('base','CERD_MODE','s');  % build score file only (skips if present).  [oai_citation:4‡cerd_main.txt](sediment://file_00000000b5047246aa24aaf2fb39252c)
        cerd_main;
    end
    if STEPS.generate
        assignin('base','CERD_MODE','g');  % generate nulls using CERD_NULL_RANGE.  [oai_citation:5‡cerd_main.txt](sediment://file_00000000b5047246aa24aaf2fb39252c)
        cerd_main;
    end
    if STEPS.analyze
        assignin('base','CERD_MODE','a');  % analyze at FDR q=ALPHA.  [oai_citation:6‡cerd_main.txt](sediment://file_00000000b5047246aa24aaf2fb39252c)
        cerd_main;
    end
    if evalin('base','exist(''CERD_LAST_SUMMARY'',''var'')')
        S = evalin('base','CERD_LAST_SUMMARY');
        if isempty(allS)
            allS = S;                         % adopt full schema on first insert
        else
            [allS, S] = local_align_struct_fields(allS, S);
            allS(end+1) = S;                  % safe append thereafter
        end
    end

    fprintf('-------------------- %s done --------------------\n', pat_id);
    fprintf('Finished %s — summaries retained above.\n', pat_id);
end

% –– Across-patient summary (concise) ––
if ~isempty(allS)
    fprintf('\n=== ACROSS-PATIENT SUMMARY (primary = THEORY) ===\n');
    fprintf('[Metric] %s\n', upper(string(METRIC_MODE)));
    % collect all metric labels seen
    labs = {};
    for s=1:numel(allS), for r=1:numel(allS(s).rows), labs{end+1}=allS(s).rows(r).name; end, end %#ok
    ulabs = unique(labs,'stable');
    % print per metric
    for L = ulabs
        lab = L{1};
        fprintf('\n[%s]\n', lab);
        for s=1:numel(allS)
            % find this metric in S
            ix = find(strcmp({allS(s).rows.name}, lab), 1);
            if isempty(ix)
                fprintf('  %-8s — (not run)\n', allS(s).patient);
            else
                R = allS(s).rows(ix);
                if R.m==0
                    fprintf('  %-8s — none tested (m=0)\n', allS(s).patient);
                else
                    fprintf('  %-8s — %d/%d sig; min_p=%s; crit_p=%s; K=%d\n', ...
                        allS(s).patient, R.sig, R.m, ...
                        ternary(isfinite(R.min_p), sprintf('%.3g',R.min_p), 'NaN'), ...
                        ternary(isfinite(R.crit_p), sprintf('%.3g',R.crit_p), 'NaN'), ...
                        R.K);
                end
            end
        end
    end
    fprintf('\n');
end
fprintf('\nALL PATIENTS COMPLETE.\n');

% Optionally close once at the end (only if we're managing the pool)
% If caller_manages_pool is true, the caller is managing the pool, so don't delete it
if ~caller_manages_pool
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
end
end

% ============================ helpers ============================
function [A, B] = local_align_struct_fields(A, B)
% Ensure struct arrays A and B share identical field sets (recursively).
if isempty(A)
    % nothing to align
    return;
end

fieldsA = fieldnames(A);
fieldsB = fieldnames(B);
allFields = unique([fieldsA; fieldsB]);

for k = 1:numel(allFields)
    fname = allFields{k};
    if ~isfield(A, fname)
        [A(1:numel(A)).(fname)] = deal([]);
    end
    if ~isfield(B, fname)
        B.(fname) = [];
    end
end

A = orderfields(A, allFields);
B = orderfields(B, allFields);

% Recurse into known nested struct arrays (e.g., 'rows') if present.
if isfield(A, 'rows') && isfield(B, 'rows') && ~isempty(A(1).rows) && ~isempty(B(1).rows)
    [A_rows, B_rows] = local_align_struct_fields([A.rows], B.rows);
    % write back reshaped rows
    startIdx = 1;
    for idx = 1:numel(A)
        n = numel(A(idx).rows);
        if n>0
            A(idx).rows = A_rows(startIdx:startIdx+n-1);
            startIdx = startIdx + n;
        else
            A(idx).rows = struct([]);
        end
    end
    B.rows = B_rows;
end
end

function ensure_local_mef(drop_root, local_root, pat_id)
src_mef_dir  = fullfile(drop_root, 'Data', pat_id, 'sub.mefd');
dest_parent  = fullfile(local_root, pat_id);
dest_mef_dir = fullfile(dest_parent, 'sub.mefd');

if ~exist(src_mef_dir,'dir')
    error('Source MEF folder not found: %s', src_mef_dir);
end
if ~exist(dest_parent,'dir'), mkdir(dest_parent); end

% fix mistaken nested copy
if exist(fullfile(dest_mef_dir,'sub.mefd'),'dir')
    warning('Found nested sub.mefd at %s — fixing.', fullfile(dest_mef_dir,'sub.mefd'));
    rmdir(dest_mef_dir,'s');
end

need_copy = ~exist(dest_mef_dir,'dir') || isempty(dir(fullfile(dest_mef_dir,'')));
if need_copy
    copyfile(src_mef_dir, dest_mef_dir);
else
    d = dir(fullfile(dest_mef_dir,''));
    if numel(d) <= 2
        rmdir(dest_mef_dir,'s');
        copyfile(src_mef_dir, dest_mef_dir);
    end
end
fprintf('[batch] Using MEF at: %s\n', dest_mef_dir);
end

function y = ternary(cond,a,b)
if ischar(cond) || isstring(cond)
    cond = lower(string(cond));
    cond = (cond=="y" | cond=="yes" | cond=="true" | cond=="1");
end
if cond, y=a; else, y=b; end
end
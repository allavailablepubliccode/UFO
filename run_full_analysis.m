% RUN_FULL_ANALYSIS  End-to-end UFO cross-scale analysis + visualization.
%
% Quick start (see README.md for details):
%   1. Make sure MEF mirrors live under ~/Desktop/tmp_mef/<patient>/sub.mefd.
%   2. Set CERD_NNULL in the base workspace (e.g., 10 or 10000).
%   3. Optionally set CERD_BATCH_PATIENTS to limit the cohort.
%   4. Run this script:  >> run_full_analysis
%
% The script:
%   • Applies the theory/PSD settings used in the paper.
%   • Runs cerd_batch twice (empirical timing + event-permute control).
%   • Leaves all outputs in Results/<patient>/.
%
clear;clc;close all;
warning('off','all');

% --- Configure analysis mode (baseline NONE, GREENS) ---
metric_mode   = 'greens';

spect_pre_ms   = 30;   % Baseline window (pre-UFO) - short to stay local and avoid contamination
spect_post_ms  = 300;  % Post-UFO window - long to capture full echo decay (~200-300ms empirically)
spect_bandwidth_hz = 600;  % Spectral bandwidth

% --- Analysis combinations to run sequentially ---
% 1) Empirical UFO timing (baseline)
% 2) Event-permuted timing (control)
analysis_combos = {
    struct('event_permute', false, 'desc', 'BASELINE: empirical timing');
    struct('event_permute', true,  'desc', 'CONTROL: event permute');
};

% Core analysis knobs
assignin('base','CERD_THEORY_METRIC',        metric_mode);
assignin('base','CERD_DETERMINISTIC',        true);
assignin('base','CERD_MANAGE_POOL',          false); % We'll manage the pool ourselves to keep it alive

% --- GREENS parameters (BEST COMBO from parameter sweep) ---
assignin('base','CERD_THEORY_F0_HZ',         [55 85]);  % BEST COMBO: Wider band, lower start
assignin('base','CERD_THEORY_NFREQ',         15);
assignin('base','CERD_THEORY_BW_HZ',         100);
assignin('base','CERD_THEORY_DELAY_MS',      [2 10]);
assignin('base','CERD_THEORY_GAMMA',         10);  % BEST COMBO: gamma=10
assignin('base','CERD_THEORY_GAMMA_RANGE',   [5 20]);
assignin('base','CERD_THEORY_NGAMMA',        3);
assignin('base','CERD_THEORY_GATE_MS',       40);
assignin('base','CERD_THEORY_PREWHITEN',     true);

% Spectral config used by GREENS PSD step
assignin('base','CERD_SPECT_PRE_MS',         spect_pre_ms);
assignin('base','CERD_SPECT_POST_MS',        spect_post_ms);
assignin('base','CERD_SPECT_BANDWIDTH_HZ',   spect_bandwidth_hz);
assignin('base','CERD_SPECT_BAND_HZ',        []);           % auto-center using UFO freq
assignin('base','CERD_SPECT_WINDOW_MS',      10);   % Larger window for better frequency resolution
assignin('base','CERD_SPECT_OVERLAP',        0.6);  % More overlap for smoother estimates
assignin('base','CERD_SPECT_NFFT',           4096);  % Higher resolution
assignin('base','CERD_SPECT_MIN_SAMPLES',    15);    % Reduced threshold to allow more valid scores
assignin('base','CERD_SPECT_FALLBACK_BAND',  [50 1000]);  % Focus on lower frequencies where macro signals are stronger

% Post-UFO window for spectral analysis
assignin('base','CERD_POST_UFO_MS', spect_pre_ms + spect_post_ms);  % Total window = pre + post

% --- Visualization control ---
DO_VISUALIZATION = false;  % Set to false to skip all visualization

% --- Visualization knobs (tweak as desired) ---
assignin('base','CERD_VISUALIZE_SIG',        false);   % auto-save plots during analysis
assignin('base','CERD_VIS_MAX_PLOTS',        50);     % per-patient cap during cerd_main
assignin('base','CERD_VIZ_PRE_MS',           5);      % baseline before UFO onset (ms)
assignin('base','CERD_VIZ_MICRO_MS',         30);     % min micro window (ms)
assignin('base','CERD_VIZ_MACRO_MS',         120);    % macro window cap (ms)
assignin('base','CERD_VIZ_XLIM_MS',          [-5 25]);% x-axis for filtered macro panel
assignin('base','CERD_POST_UFO_MS',          300);    % macro window after UFO (ms)

% Run full pipeline (score → generate → analyze). Force regeneration for new params.
assignin('base','CERD_BATCH_STEPS', struct('score', true, 'generate', true, 'analyze', true));

% Check if we need parallel pool (only for score/generate, not analyze-only)
batch_steps = struct('score', false, 'generate', false, 'analyze', true);
try
    batch_steps = evalin('base', 'CERD_BATCH_STEPS');
catch
    % Use default if variable doesn't exist
end
need_pool = (batch_steps.score || batch_steps.generate);

% --- Initialize parallel pool only if needed (for score/generate steps) ---
if need_pool
    fprintf('\n========================================\n');
    fprintf('INITIALIZING PARALLEL POOL (needed for score/generate)\n');
    fprintf('========================================\n');
    p = gcp('nocreate');
    if isempty(p)
        want = max(1, feature('numcores')-2);
        parpool('local', want);
        fprintf('Parallel pool initialized with %d workers.\n', want);
    else
        fprintf('Using existing parallel pool (%d workers).\n', p.NumWorkers);
    end
    fprintf('========================================\n\n');
else
    fprintf('\n[Note] Analyze-only mode: skipping parallel pool initialization\n\n');
end

% --- Run all analysis combinations sequentially ---
fprintf('\n========================================\n');
fprintf('RUNNING %d ANALYSIS COMBINATIONS SEQUENTIALLY\n', length(analysis_combos));
fprintf('========================================\n\n');

for combo_idx = 1:length(analysis_combos)
    combo = analysis_combos{combo_idx};
    
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('COMBINATION %d/%d: %s\n', combo_idx, length(analysis_combos), combo.desc);
    fprintf('========================================\n');
    fprintf('  event_permute: %d\n', combo.event_permute);
    
    % Set parameters for this combination
    assignin('base','CERD_EVENT_PERMUTE',        combo.event_permute);
    
    fprintf('\n>>> Running cerd_batch for combination %d/%d...\n', combo_idx, length(analysis_combos));
    
    try
        cerd_batch;
        fprintf('\n✓ Combination %d/%d completed successfully.\n', combo_idx, length(analysis_combos));
    catch ME
        fprintf('\n✗ ERROR in combination %d/%d: %s\n', combo_idx, length(analysis_combos), ME.message);
        fprintf('  Continuing with next combination...\n');
        % Continue to next combination even if this one fails
    end
    
    fprintf('\n');
end

fprintf('\n========================================\n');
fprintf('ALL COMBINATIONS COMPLETE\n');
fprintf('========================================\n');

% --- Close parallel pool at the end (only if we created it) ---
if need_pool
    p = gcp('nocreate');
    if ~isempty(p)
        fprintf('\nClosing parallel pool...\n');
        delete(p);
        fprintf('Parallel pool closed.\n');
    end
end

% --- Visualization pass ---
% patientIDs = {'pat001','pat003','pat004','pat061','pat072','pat079'};
patientIDs = {'pat001'};
maxPlots   = 100;                                   % hard cap per patient for export
% Use ~/Desktop/plots explicitly
home_dir = getenv('HOME');
if isempty(home_dir)
    home_dir = '~';
end
outRoot = fullfile(home_dir, 'Desktop', 'plots');
if ~exist(outRoot, 'dir')
    mkdir(outRoot);
    fprintf('Created output directory: %s\n', outRoot);
end
post_ms    = 300;

% --- Visualization pass (optional) ---
if DO_VISUALIZATION
    fprintf('>>> Generating plots for significant UFOs...\n');
    if exist('plot_all_significant_ufos','file') == 2
        plot_all_significant_ufos(patientIDs, maxPlots, outRoot, post_ms);
    else
        warning('plot_all_significant_ufos.m not found. Skipping visualization export.');
    end

    % Also visualize top candidates (even if not significant) for inspection
    fprintf('>>> Generating plots for top candidates (for inspection)...\n');
    if exist('visualize_top_candidates','file') == 2
        for p = 1:numel(patientIDs)
            try
                fprintf('  Processing %s...\n', patientIDs{p});
                visualize_top_candidates(patientIDs{p}, 10, [], [], []);
                fprintf('  ✓ %s done\n', patientIDs{p});
            catch ME
                fprintf('  ✗ %s failed: %s\n', patientIDs{p}, ME.message);
                fprintf('    Stack trace:\n');
                for k = 1:numel(ME.stack)
                    fprintf('      %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
                end
            end
        end
    else
        warning('visualize_top_candidates.m not found. Skipping top candidate visualization.');
    end

    fprintf('All done. Plots available in %s\n', outRoot);
    fprintf('Top candidate plots saved in Results/<patient>/top_candidates/\n');
else
    fprintf('>>> Visualization skipped (DO_VISUALIZATION = false)\n');
    fprintf('All done.\n');
end

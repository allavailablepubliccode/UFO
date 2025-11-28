% Visualize significant UFOs using cerd_main analyze-only mode
% This script sets the visualization flags and re-runs cerd_main so that
% saved significant UFOs are plotted without recomputing scores/nulls.

clear;

% Patients with significant echoes (adjust as needed)
patients = {'pat003', 'pat079'};

% Maximum plots per patient (raise if more sig UFOs appear)
% Set to [] to use CERD_VIS_EVENT_INDICES or CERD_VIS_EVENT_NUMBERS instead
vis_max_plots = 5;

% Optional: Select specific events to visualize
% Option 1: By index in sig_ufos array (1-based, sorted by p-value)
%   assignin('base','CERD_VIS_EVENT_INDICES', [1, 2, 3]);  % Plot first 3 most significant
% Option 2: By actual event_idx number
%   assignin('base','CERD_VIS_EVENT_NUMBERS', [1, 5, 10]);  % Plot events with these event_idx values
% If both are empty, uses vis_max_plots to plot first N events

fprintf('=== visualize_sig_ufos_direct ===\n');
fprintf('Re-running cerd_main analyze-only with visualization enabled.\n\n');

% --- Match the analysis parameters used for the saved results ---
assignin('base','CERD_THEORY_METRIC',   'greens');
assignin('base','CERD_DETERMINISTIC',   true);
assignin('base','CERD_MANAGE_POOL',     false);

assignin('base','CERD_THEORY_F0_HZ',    [55 85]);
assignin('base','CERD_THEORY_NFREQ',    15);
assignin('base','CERD_THEORY_BW_HZ',    100);
assignin('base','CERD_THEORY_DELAY_MS', [2 10]);
assignin('base','CERD_THEORY_GAMMA',    10);
assignin('base','CERD_THEORY_GAMMA_RANGE', [5 20]);
assignin('base','CERD_THEORY_NGAMMA',   3);
assignin('base','CERD_THEORY_GATE_MS',  40);
assignin('base','CERD_THEORY_PREWHITEN', true);

spect_pre_ms  = 30;
spect_post_ms = 300;
spect_bw_hz   = 600;

assignin('base','CERD_SPECT_PRE_MS',        spect_pre_ms);
assignin('base','CERD_SPECT_POST_MS',       spect_post_ms);
assignin('base','CERD_SPECT_BANDWIDTH_HZ',  spect_bw_hz);
assignin('base','CERD_SPECT_BAND_HZ',       []);
assignin('base','CERD_SPECT_WINDOW_MS',     10);
assignin('base','CERD_SPECT_OVERLAP',       0.6);
assignin('base','CERD_SPECT_NFFT',          4096);
assignin('base','CERD_SPECT_MIN_SAMPLES',   15);
assignin('base','CERD_SPECT_FALLBACK_BAND', [50 1000]);

assignin('base','CERD_POST_UFO_MS', spect_post_ms);  % match analyze-only setting (300 ms)

% Visualization knobs (match run_full_analysis defaults)
assignin('base','CERD_VIZ_PRE_MS',    5);
assignin('base','CERD_VIZ_MICRO_MS',  30);
assignin('base','CERD_VIZ_MACRO_MS',  120);
assignin('base','CERD_VIZ_XLIM_MS',   [-5 25]);
assignin('base','CERD_EVENT_PERMUTE',        false);

% PSD plot axis limits
% x-axis: set to [55 85] to focus on analysis band, or [] for default 0-200 Hz
assignin('base','CERD_VIZ_PSD_XLIM',  [55-30 85+30]);
% y-axis: set to [min max] in dB, or [] for auto-scaling
assignin('base','CERD_VIZ_PSD_YLIM',  []);  % e.g., [-20 10] for specific range

for p = 1:numel(patients)
    patient = patients{p};
    fprintf('--- %s ---\n', patient);
    
    % Configure flags in the base workspace for this patient
    assignin('base','CERD_PATIENT', patient);
    assignin('base','CERD_MODE', 'a');               % analyze-only
    assignin('base','CERD_VISUALIZE_SIG', true);
    assignin('base','CERD_VIS_MAX_PLOTS', vis_max_plots);
    
    % Optional: Select specific events to visualize for this patient
    % Uncomment and modify one of these options:
    % assignin('base','CERD_VIS_EVENT_INDICES', [2]);  % Plot by index (1=most significant)
    % assignin('base','CERD_VIS_EVENT_NUMBERS', [2]); % Plot by event_idx number
    assignin('base','CERD_NNULL', 1000);
    assignin('base','CERD_NULL_RANGE', '1:1000');
    assignin('base','CERD_BATCH_STEPS', struct('score', false, 'generate', false, 'analyze', true));
    
    try
        cerd_main;
        fprintf('  ✓ Visualizations generated for %s\n\n', patient);
    catch ME
        warning('  ✗ Failed to visualize %s: %s', patient, ME.message);
    end
end

fprintf('All done. Check Results/<patient>/visualizations for outputs.\n');


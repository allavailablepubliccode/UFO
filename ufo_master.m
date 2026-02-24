function ufo_master()
% UFO_MASTER  Pipeline controller for UFO inspection and merging.
%
% Usage:
%   ufo_master()
%
% Lets you choose patients and which steps to run:
%   Step 1: step1_inspect (Interactive labelling)
%   Step 2: step2_merge   (Auto-merge overlaps and refresh figures)
%
% Directory structure (configure in ufo_config.m):
%   BASE_DIR/
%     Data/<patient>/           - MEF3 data (sub.mefd or *.mefd), UFO.json
%     Data/<patient>/UFO_edited/ - Edited UFOs (UFOdat_edited.mat)
%     Results/<patient>/         - ufo_inspect_log.mat, prep_basics.mat, UFO_figures/

clc; close all;

patIDs = {'pat001','pat003','pat004','pat061','pat072','pat079','pat038','pat397','pat448','pat452','pat453'};
% Patient ID mapping (internal → paper):
%   pat001 = WROC00001    pat038 = MAYO00038    pat397 = MAYO00397
%   pat003 = WROC00003    pat061 = FNUS00061    pat448 = MAYO00448
%   pat004 = WROC00004    pat072 = FNUS00072    pat452 = MAYO00452
%                         pat079 = FNUS00079    pat453 = MAYO00453

fprintf('=== UFO PIPELINE MASTER ===\n\n');
fprintf('Available patients:\n');
for i = 1:numel(patIDs)
    fprintf('  %d) %s\n', i, patIDs{i});
end

p_choice = input('\nSelect patient index (or indices like [1 4]): ');
if isempty(p_choice), return; end

fprintf('\nAvailable Steps:\n');
fprintf('  1) Inspect & Classify (step1_inspect)\n');
fprintf('  2) Merge Overlaps & Refresh Figs (step2_merge)\n');
fprintf('  3) Statistical Analysis (step3_statistics_fixed — fixed frequency-domain Green''s envelope)\n');

s_choice = input('\nSelect steps to run (e.g. 1 or [1 2]): ');
if isempty(s_choice), return; end

% If step 3 is selected, prompt for selection and N once (applies to all patients)
stats_choices = [];
stats_N = [];
stats_equalize_options = [];
if ismember(3, s_choice)
    fprintf('\n=== Statistics Configuration (applies to all patients) ===\n');
    fprintf('Using fixed analysis parameters:\n');
    fprintf('  Frequency band: [60, 100] Hz\n');
    fprintf('  Macro natural frequency: f0 = 75 Hz\n');
    fprintf('  Damping: gamma = 40 (1/s)\n');
    fprintf('  Baseline window: pre_ms = 200 ms\n');
    fprintf('  Post window: post_ms = 200 ms\n');
    fprintf('  Statistical test: One p-value per micro-macro pair (one-sided, greater)\n\n');

    fprintf('Selection Options:\n');
    fprintf('  1) Type 1 only, all frequencies\n');
    fprintf('  2) Type 2 only, all frequencies\n');
    fprintf('  3) Both Type 1 & 2, all frequencies\n');
    fprintf('  4) Type 1 only, only >2kHz\n');
    fprintf('  5) Type 2 only, only >2kHz\n');
    fprintf('  6) Both Type 1 & 2, only >2kHz\n');
    fprintf('  7) ALL UFOs (original detections, ignore post-hoc labels)\n');
    stats_choices = input('\nSelect option [1-7] (can be a vector, e.g. [4 5 7]): ');
    if isempty(stats_choices), return; end
    stats_N = input('Enter number of null/bootstrap iterations (e.g. 1000): ');
    if isempty(stats_N) || stats_N < 1, return; end

    % Allow both y and n for equalize
    eq_resp = lower(strtrim(input('Equalize sample sizes? [y/n, or [y n] for both, default n]: ','s')));
    if isempty(eq_resp), eq_resp = 'n'; end
    % Parse response - could be 'y', 'n', 'y n', '[y n]', '[y,n]', etc.
    eq_resp_clean = strrep(eq_resp, '[', '');
    eq_resp_clean = strrep(eq_resp_clean, ']', '');
    eq_resp_clean = strrep(eq_resp_clean, ',', ' ');
    eq_parts = strsplit(strtrim(eq_resp_clean));

    if numel(eq_parts) > 1 || (contains(eq_resp, 'y') && contains(eq_resp, 'n'))
        % Multiple options
        stats_equalize_options = [true, false];
        fprintf('  Will run both: with and without equalization.\n');
    elseif any(strcmp(eq_parts, 'y')) || any(strcmp(eq_parts, 'yes'))
        stats_equalize_options = true;
        fprintf('  Will equalize sample sizes (use minimum N across all selections).\n');
    else
        stats_equalize_options = false;
    end

    % Minimum duration filter
    stats_min_duration = input('Minimum UFO duration (ms) to include [default: 0 = no filter]: ', 's');
    if isempty(stats_min_duration)
        stats_min_duration_ms = 0;
    else
        stats_min_duration_ms = str2double(stats_min_duration);
        if isnan(stats_min_duration_ms) || stats_min_duration_ms < 0
            fprintf('  Invalid input. Using default: no filter.\n');
            stats_min_duration_ms = 0;
        elseif stats_min_duration_ms > 0
            fprintf('  Will only analyze UFOs with duration >= %.1f ms.\n', stats_min_duration_ms);
        end
    end

end

% Collect results for summary table
% Structure: {patient, selection, n_ufos, emp_avg, p_val, version, equalize}
summary_results = {};

fprintf('\n>>> Processing %d patient(s)...\n', numel(p_choice));
p_count = 0;
for p_idx = p_choice
    p_count = p_count + 1;
    patient = patIDs{p_idx};
    fprintf('\n========================================\n');
    fprintf('Patient %d/%d: %s\n', p_count, numel(p_choice), patient);
    fprintf('========================================\n');

    for s_idx = s_choice
        if s_idx == 1
            fprintf('\n--- Step 1: Inspecting %s ---\n', patient);
            step1_inspect(patient);
        elseif s_idx == 2
            fprintf('\n--- Step 2: Merging %s ---\n', patient);
            step2_merge(patient);
        elseif s_idx == 3
            % Loop over equalize combinations
            for e_idx = 1:numel(stats_equalize_options)
                stats_equalize = stats_equalize_options(e_idx);
                eq_str = 'equalized';
                if ~stats_equalize, eq_str = 'not-equalized'; end

                patient_results = step3_statistics_fixed(patient, stats_choices, stats_N, stats_equalize, stats_min_duration_ms);
                if ~isempty(patient_results)
                    % Add equalize metadata to each result
                    for r_idx = 1:numel(patient_results)
                        r = patient_results{r_idx};
                        % Result structure: {patient, desc_final, n_sel, n_total_pairs, n_significant, n_fdr_significant}
                        % Append equalize info as 7th element
                        r{7} = stats_equalize;
                        summary_results{end+1} = r;
                    end
                end
            end
        else
            fprintf('\nUnknown step index: %d\n', s_idx);
        end
    end
end

% Display summary table if we have results
if ismember(3, s_choice)
    if isempty(summary_results)
        fprintf('\n\n=== NO RESULTS TO DISPLAY ===\n');
        fprintf('No valid results were generated. This may occur if:\n');
        fprintf('  - All selections had NaN empirical averages\n');
        fprintf('  - All selections were skipped (no UFOs found)\n');
        fprintf('  - An error occurred during processing\n');
        fprintf('\nCheck the individual patient processing for details.\n\n');
    else
        % Group results by equalize
        equalize_flags = unique(cellfun(@(x) x{7}, summary_results));

        fprintf('\n\n=== SUMMARY TABLES (Grouped by Analysis Type) ===\n\n');

        % Display tables grouped by equalize
        for eq = equalize_flags(:)'
            eq_str = 'Equalized';
            if ~eq, eq_str = 'Not-Equalized'; end
            fprintf('--- Fixed Parameters, %s ---\n', eq_str);
            fprintf('%-10s %-20s %-8s %-8s %-8s %-8s\n', ...
                'Patient', 'Selection', 'N_UFOs', 'N_Pairs', 'N_Sig', 'N_FDR');
            fprintf('%s\n', repmat('-', 1, 70));

            % Filter results for this combination
            filtered = {};
            for i = 1:numel(summary_results)
                r = summary_results{i};
                % Equalize is at index 7
                if numel(r) >= 7 && r{7} == eq
                    filtered{end+1} = r;
                end
            end

            % Display results (already sorted)
            if ~isempty(filtered)
                % Define selection type order (for consistent sorting)
                selection_order = {'AllUFOs_OrigJSON', 'Type12_AllFreq', 'Type1_AllFreq', 'Type2_AllFreq', ...
                    'Type12_HighFreq', 'Type1_HighFreq', 'Type2_HighFreq'};

                % Extract selection base names and patient numbers for sorting
                sel_indices = zeros(numel(filtered), 1);
                pat_nums = zeros(numel(filtered), 1);
                for i = 1:numel(filtered)
                    sel_name = filtered{i}{2};
                    % Strip _Equalized suffix if present
                    if length(sel_name) >= 11 && strcmp(sel_name(end-10:end), '_Equalized')
                        sel_name = sel_name(1:end-11);
                    end
                    % Strip _minDur*.ms suffix if present
                    idx = strfind(sel_name, '_minDur');
                    if ~isempty(idx)
                        sel_name = sel_name(1:idx(end)-1);
                    end

                    % Find selection index
                    sel_idx = find(strcmp(selection_order, sel_name), 1);
                    if isempty(sel_idx), sel_idx = 999; end
                    sel_indices(i) = sel_idx;

                    % Extract patient number
                    pat_str = filtered{i}{1};
                    if length(pat_str) >= 4 && strcmp(pat_str(1:3), 'pat')
                        pat_num_str = pat_str(4:end);
                        pat_nums(i) = str2double(pat_num_str);
                        if isnan(pat_nums(i)), pat_nums(i) = 9999; end
                    else
                        pat_nums(i) = 9999;
                    end
                end

                % Sort by selection index, then patient number
                [~, sort_idx] = sortrows([sel_indices, pat_nums]);
                filtered = filtered(sort_idx);

                % Display results
                for i = 1:numel(filtered)
                    r = filtered{i};
                    % Structure: {patient, desc_final, n_sel, n_total_pairs, n_significant, n_fdr_significant, equalize}
                    fprintf('%-10s %-20s %8d %8d %8d %8d\n', ...
                        r{1}, r{2}, r{3}, r{4}, r{5}, r{6});
                end
            end
        end
        fprintf('\n');

        % Summary of most significant results
        fprintf('=== ANALYSIS TYPE SUMMARY ===\n\n');

        % Count significant results by equalize type
        sig_counts = struct();
        total_counts = struct();
        sig_patients = struct();
        fdr_sig_counts = struct();
        fdr_sig_patients = struct();

        n_combinations_uncorrected = struct();
        n_combinations_fdr = struct();
        aggregate_pairs = struct();
        aggregate_sig = struct();
        aggregate_fdr = struct();
        for eq = equalize_flags(:)'
            eq_str_field = 'Equalized';
            if ~eq, eq_str_field = 'NotEqualized'; end
            field_name = eq_str_field;
            sig_counts.(field_name) = 0;
            total_counts.(field_name) = 0;
            sig_patients.(field_name) = {};
            fdr_sig_counts.(field_name) = 0;
            fdr_sig_patients.(field_name) = {};
            n_combinations_uncorrected.(field_name) = 0;
            n_combinations_fdr.(field_name) = 0;
            aggregate_pairs.(field_name) = 0;
            aggregate_sig.(field_name) = 0;
            aggregate_fdr.(field_name) = 0;

            for i = 1:numel(summary_results)
                r = summary_results{i};
                % Equalize is at index 7
                if numel(r) >= 7 && r{7} == eq
                    total_counts.(field_name) = total_counts.(field_name) + 1;
                    if r{5} >= 1
                        n_combinations_uncorrected.(field_name) = n_combinations_uncorrected.(field_name) + 1;
                    end
                    if r{6} >= 1
                        n_combinations_fdr.(field_name) = n_combinations_fdr.(field_name) + 1;
                    end
                    if numel(r) >= 6 && ~isnan(r{4})
                        aggregate_pairs.(field_name) = aggregate_pairs.(field_name) + r{4};
                        aggregate_sig.(field_name) = aggregate_sig.(field_name) + r{5};
                        aggregate_fdr.(field_name) = aggregate_fdr.(field_name) + r{6};
                    end
                    if r{5} > 0
                        sig_counts.(field_name) = sig_counts.(field_name) + r{5};
                        pat_id = r{1};
                        if ~any(strcmp(sig_patients.(field_name), pat_id))
                            sig_patients.(field_name){end+1} = pat_id;
                        end
                    end
                    if r{6} > 0
                        fdr_sig_counts.(field_name) = fdr_sig_counts.(field_name) + r{6};
                        pat_id = r{1};
                        if ~any(strcmp(fdr_sig_patients.(field_name), pat_id))
                            fdr_sig_patients.(field_name){end+1} = pat_id;
                        end
                    end
                end
            end
        end

        % FDR scope (explicit one-line statement)
        fprintf('FDR correction is applied within each patient/selection across that patient''s micro-macro pairs (not pooled across patients).\n\n');

        % Display summary (by patient-selection combinations)
        fprintf('=== ANALYSIS TYPE SUMMARY (by patient-selection combinations) ===\n\n');
        fprintf('%-30s %-28s %-28s %-45s\n', 'Analysis Type', 'Patients with >=1 uncorrected hit', 'Patients with >=1 FDR hit', 'Total patients evaluated (patient-selection combinations)');
        fprintf('%s\n', repmat('-', 1, 135));
        for eq = equalize_flags(:)'
            eq_str_field = 'Equalized';
            if ~eq, eq_str_field = 'NotEqualized'; end
            field_name = eq_str_field;
            eq_str_display = 'Equalized';
            if ~eq, eq_str_display = 'Not-Equalized'; end
            display_name = sprintf('Fixed Parameters, %s', eq_str_display);
            n_comb = total_counts.(field_name);
            n_pat_unc = numel(sig_patients.(field_name));
            n_pat_fdr = numel(fdr_sig_patients.(field_name));
            fprintf('%-30s %-28d %-28d %-45d\n', display_name, n_pat_unc, n_pat_fdr, n_comb);
        end
        fprintf('\n');

        % Pair-level summary (total_pairs = sum of N_Pairs, hit-rates)
        fprintf('Pair-level summary (per selection, sum over patients):\n');
        for eq = equalize_flags(:)'
            eq_str_field = 'Equalized';
            if ~eq, eq_str_field = 'NotEqualized'; end
            field_name = eq_str_field;
            eq_str_display = 'Equalized';
            if ~eq, eq_str_display = 'Not-Equalized'; end
            display_name = sprintf('Fixed Parameters, %s', eq_str_display);
            total_pairs = aggregate_pairs.(field_name);
            total_sig_pairs = aggregate_sig.(field_name);
            total_fdr_pairs = aggregate_fdr.(field_name);
            fprintf('  %s: total_pairs = %d, total_sig_pairs = %d, total_fdr_pairs = %d\n', ...
                display_name, total_pairs, total_sig_pairs, total_fdr_pairs);
            if total_pairs > 0
                fprintf('    Pair-level uncorrected hit-rate: %d / %d (%.2f%%)\n', total_sig_pairs, total_pairs, 100 * total_sig_pairs / total_pairs);
                fprintf('    Pair-level FDR hit-rate: %d / %d (%.2f%%)\n', total_fdr_pairs, total_pairs, 100 * total_fdr_pairs / total_pairs);
            end
        end
        fprintf('\n');

        % Aggregate pairwise hypothesis counts (same totals, labeled for clarity)
        fprintf('Aggregate pairwise hypothesis counts:\n');
        for eq = equalize_flags(:)'
            eq_str_field = 'Equalized';
            if ~eq, eq_str_field = 'NotEqualized'; end
            field_name = eq_str_field;
            eq_str_display = 'Equalized';
            if ~eq, eq_str_display = 'Not-Equalized'; end
            display_name = sprintf('Fixed Parameters, %s', eq_str_display);
            M = aggregate_pairs.(field_name);
            S = aggregate_sig.(field_name);
            F = aggregate_fdr.(field_name);
            fprintf('  %s: Total pairwise tests (sum over selections): %d; Uncorrected significant pairs: %d; FDR significant pairs: %d\n', ...
                display_name, M, S, F);
        end
        fprintf('\nImportant: N_FDR counts are based on the above scope; they are not comparable across different scopes.\n\n');

        % Aggregate statistics across all patients
        fprintf('=== AGGREGATE STATISTICS ACROSS ALL PATIENTS ===\n\n');

        % Calculate overall averages
        all_n_ufos = [];
        all_n_pairs = [];
        all_n_sig = [];
        all_n_fdr_sig = [];

        for i = 1:numel(summary_results)
            r = summary_results{i};
            if numel(r) >= 6
                all_n_ufos(end+1) = r{3};      % N_UFOs
                all_n_pairs(end+1) = r{4};     % N_Pairs
                all_n_sig(end+1) = r{5};       % N_Significant
                all_n_fdr_sig(end+1) = r{6};   % N_FDR_Significant
            end
        end

        fprintf('Overall Statistics:\n');
        fprintf('  Total patient-selection combinations: %d\n', numel(summary_results));
        fprintf('  Average UFOs per combination: %.1f (range: %d-%d)\n', ...
            mean(all_n_ufos), min(all_n_ufos), max(all_n_ufos));
        fprintf('  Average pairs per combination: %.1f (range: %d-%d)\n', ...
            mean(all_n_pairs), min(all_n_pairs), max(all_n_pairs));
        fprintf('  Total significant pairs (p<0.05): %d\n', sum(all_n_sig));
        fprintf('  Total FDR-significant pairs (q<0.05): %d\n', sum(all_n_fdr_sig));

        % Statistics by selection type
        fprintf('\nStatistics by Selection Type:\n');
        selection_types = unique(cellfun(@(x) x{2}, summary_results, 'UniformOutput', false));
        for sel = selection_types(:)'
            sel_name = sel{1};
            % Strip suffixes for grouping
            sel_base = sel_name;
            if contains(sel_base, '_Equalized')
                sel_base = sel_base(1:strfind(sel_base, '_Equalized')-1);
            end
            if contains(sel_base, '_minDur')
                idx = strfind(sel_base, '_minDur');
                sel_base = sel_base(1:idx-1);
            end

            sel_results = summary_results(cellfun(@(x) strcmp(x{2}, sel_name), summary_results));
            if ~isempty(sel_results)
                sel_n_ufos = cellfun(@(x) x{3}, sel_results);
                sel_n_pairs = cellfun(@(x) x{4}, sel_results);
                sel_n_sig = cellfun(@(x) x{5}, sel_results);
                sel_n_fdr_sig = cellfun(@(x) x{6}, sel_results);

                fprintf('  %s:\n', sel_base);
                fprintf('    Combinations: %d | Avg UFOs: %.1f | Avg Pairs: %.1f | Sig Pairs: %d | FDR-Sig Pairs: %d\n', ...
                    numel(sel_results), mean(sel_n_ufos), mean(sel_n_pairs), sum(sel_n_sig), sum(sel_n_fdr_sig));
            end
        end
        fprintf('\n');
    end  % Close the if isempty(summary_results) block (includes else)
end  % Close the if ismember(3, s_choice) block

fprintf('\nPipeline master finished.\n');



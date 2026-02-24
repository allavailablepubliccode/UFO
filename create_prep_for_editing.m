function create_prep_for_editing(patient)
% CREATE_PREP_FOR_EDITING  Create prep_basics.mat for UFO editing without full analysis
%
% Usage:
%   create_prep_for_editing('pat448')
%
% This creates prep_basics.mat by running the prep steps from cerd_main
% without running the full analysis or null generation.
%
% After running this, you can use: edit_ufos_interactive('pat448')

if nargin < 1 || isempty(patient)
    error('create_prep_for_editing: patient ID required');
end

fprintf('=== Creating prep_basics.mat for %s ===\n\n', patient);

[base_dir, ~] = ufo_config();
results_dir = fullfile(base_dir, 'Results', patient);
prep_file = fullfile(results_dir, 'prep_basics.mat');

if exist(prep_file, 'file')
    resp = input(sprintf('prep_basics.mat already exists. Overwrite? [y/n]: '), 's');
    if ~strcmpi(strtrim(resp), 'y')
        fprintf('Aborted.\n');
        return;
    end
end

% Ensure results directory exists
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('Created Results directory: %s\n', results_dir);
end

fprintf('Running prep steps (this will load patient metadata and UFO.json)...\n\n');

% Since cerd_main requires edited UFOs to exist, we'll create a temporary
% empty edited UFO file, let cerd_main create prep_basics.mat, then clean up.
% This is a workaround to bypass the edited UFO check.

% Check if cerd_main exists
if ~exist('cerd_main', 'file')
    error('cerd_main.m not found in path. Please ensure it is accessible.');
end

% Create temporary edited UFO directory and file
ufo_edited_dir = fullfile(base_dir, 'Data', patient, 'UFO_edited');
if ~exist(ufo_edited_dir, 'dir')
    mkdir(ufo_edited_dir);
end

temp_ufo_file = fullfile(ufo_edited_dir, 'UFOdat_edited.mat');
temp_ufo_exists = exist(temp_ufo_file, 'file');

try
    % Create a minimal empty UFOdat_edited to satisfy cerd_main's check
    % We'll use the original UFOs from UFO.json via extract_UFOs
    fprintf('Creating temporary edited UFO file (will be replaced during editing)...\n');
    
    % Save current settings
    old_mode = [];
    old_patient = [];
    old_event_permute = [];
    if evalin('base', 'exist(''CERD_MODE'', ''var'')')
        old_mode = evalin('base', 'CERD_MODE');
    end
    if evalin('base', 'exist(''CERD_PATIENT'', ''var'')')
        old_patient = evalin('base', 'CERD_PATIENT');
    end
    if evalin('base', 'exist(''CERD_EVENT_PERMUTE'', ''var'')')
        old_event_permute = evalin('base', 'CERD_EVENT_PERMUTE');
    end
    
    % Set parameters for cerd_main
    evalin('base', 'CERD_MODE = ''s'';');  % Use score mode - creates prep_basics.mat without generating nulls
    evalin('base', sprintf('CERD_PATIENT = ''%s'';', patient));
    evalin('base', 'CERD_EVENT_PERMUTE = false;');  % Must be false to create prep_basics.mat
    
    % We need to create a valid UFOdat_edited structure
    % The simplest: create it by running prep_basics first, then extract_UFOs
    % But we can't call those directly. Instead, we'll create a minimal version
    % that cerd_main can use, then it will create the full prep_basics.mat
    
    % Actually, better: Let cerd_main do the work, but we need to provide
    % edited UFOs. We can create them from the original UFO.json by
    % calling the extraction logic. But since we can't call extract_UFOs...
    
    % We need to create a valid NON-EMPTY UFOdat_edited from the original UFO.json
    % This requires calling prep_basics first to get Dat, then extract_UFOs.
    % Since we can't call those directly, we'll manually extract UFOs from UFO.json.
    
    fprintf('  Extracting UFOs from UFO.json to create edited UFO file...\n');
    
    % Load original UFO.json
    ufo_json = fullfile(base_dir, 'Data', patient, 'UFO.json');
    
    if ~exist(ufo_json, 'file')
        error('UFO.json not found at %s', ufo_json);
    end
    
    try
        U = jsondecode(fileread(ufo_json));
    catch ME_json
        error('Failed to read UFO.json: %s', ME_json.message);
    end
    
    if isempty(U) || ~isstruct(U)
        error('UFO.json is empty or invalid for %s', patient);
    end
    
    % We need to know the micro channels to properly extract UFOs
    % For now, create a simple structure that matches what extract_UFOs would create
    % We'll extract all UFOs and group them by channel
    
    % Find channel field
    if isfield(U, 'channels')
        chanField = 'channels';
    elseif isfield(U, 'channel')
        chanField = 'channel';
    else
        error('UFO.json missing channel(s) field');
    end
    
    % Find time and frequency fields
    leftField = '';
    for f = {'uutc_left','left','t0','start','start_uutc'}
        if isfield(U, f{1}), leftField = f{1}; break; end
    end
    rightField = '';
    for f = {'uutc_right','right','t1','end','end_uutc'}
        if isfield(U, f{1}), rightField = f{1}; break; end
    end
    freqField = '';
    for f = {'frequency','freq','f0','peak_frequency_hz','peak_hz'}
        if isfield(U, f{1}), freqField = f{1}; break; end
    end
    
    if isempty(leftField) || isempty(rightField) || isempty(freqField)
        error('UFO.json missing required time/frequency fields');
    end
    
    % Group UFOs by channel (we'll create one cell per unique channel)
    % For simplicity, create a structure that cerd_main can use
    % We'll let cerd_main's prep_basics figure out the channels
    
    % Actually, simpler: Create UFOdat_edited by copying all UFOs into a single cell
    % cerd_main will then properly extract them per channel via extract_UFOs
    % But load_edited_ufos expects the structure to already be per-channel.
    
    % Best approach: Create a minimal valid structure that will pass the check
    % We'll create one UFO per channel found in the JSON
    channels = {};
    for k = 1:numel(U)
        ch = U(k).(chanField);
        if ischar(ch) || (isstring(ch) && isscalar(ch))
            channels{end+1} = char(ch);
        elseif iscell(ch) || (isstring(ch) && ~isscalar(ch))
            if isstring(ch), ch = cellstr(ch); end
            channels(end+1:end+numel(ch)) = cellstr(ch);
        end
    end
    channels = unique(channels);
    
    % Create UFOdat_edited: one cell per channel
    UFOdat_edited = cell(1, numel(channels));
    for ii = 1:numel(channels)
        % Find all UFOs for this channel
        ufos_for_ch = [];
        for k = 1:numel(U)
            ch = U(k).(chanField);
            matches = false;
            if ischar(ch) || (isstring(ch) && isscalar(ch))
                matches = strcmp(char(ch), channels{ii});
            elseif iscell(ch) || (isstring(ch) && ~isscalar(ch))
                if isstring(ch), ch = cellstr(ch); end
                matches = any(strcmp(ch, channels{ii}));
            end
            if matches
                a = U(k).(leftField);
                b = U(k).(rightField);
                f = U(k).(freqField);
                if ~isempty(a) && ~isempty(b) && isfinite(double(a)) && isfinite(double(b)) && double(b) > double(a)
                    a_us = double(a);
                    b_us = double(b);
                    % Convert to microseconds if needed
                    % If field name contains 'uutc' or 'us', assume already in microseconds
                    % Check magnitude to determine units:
                    % - Unix timestamps: ~1e9 (seconds since epoch) -> convert *1e6
                    % - Small values < 1e6: likely seconds -> convert *1e6  
                    % - Very large > 1e12: likely nanoseconds -> convert /1e3
                    % - Values 1e6 to 1e9: could be microseconds or milliseconds
                    %   (microseconds if reasonable duration, milliseconds if very long)
                    if ~contains(lower(leftField), 'uutc') && ~contains(lower(leftField), 'us')
                        % Check if values look like Unix timestamps (seconds since epoch)
                        % Unix timestamps are typically 1e9 to 1e10 (years 2001-2286)
                        if a_us >= 1e9 && a_us < 1e10
                            % Likely Unix timestamp in seconds, convert to microseconds
                            a_us = a_us * 1e6;
                            b_us = b_us * 1e6;
                        elseif a_us < 1e6 && b_us < 1e6
                            % Likely in seconds (reasonable range: 0 to 1e6 seconds = ~11 days)
                            a_us = a_us * 1e6;
                            b_us = b_us * 1e6;
                        elseif a_us > 1e12 || b_us > 1e12
                            % Likely in nanoseconds, convert to microseconds
                            a_us = a_us / 1e3;
                            b_us = b_us / 1e3;
                        elseif a_us >= 1e6 && a_us < 1e9
                            % Could be microseconds or milliseconds
                            % Check duration: if duration > 1000 ms, likely milliseconds
                            dur_ms = (b_us - a_us) / 1e3;  % assuming microseconds
                            if dur_ms > 1000
                                % Likely milliseconds, convert to microseconds
                                a_us = a_us * 1e3;
                                b_us = b_us * 1e3;
                            end
                            % Otherwise assume already in microseconds
                        end
                    end
                    ufos_for_ch(end+1, :) = [a_us, b_us, double(f)];
                end
            end
        end
        if isempty(ufos_for_ch)
            UFOdat_edited{ii} = zeros(0, 3, 'int64');
        else
            UFOdat_edited{ii} = int64(ufos_for_ch);
        end
    end
    
    % Save the edited UFO file
    Dat_temp = struct('ID', patient);
    save(temp_ufo_file, 'UFOdat_edited', 'Dat_temp', '-v7.3');
    fprintf('  Created edited UFO file with %d channels, %d total UFOs\n', ...
        numel(channels), sum(cellfun(@(u) size(u,1), UFOdat_edited)));
    
    fprintf('Calling cerd_main to create prep_basics.mat...\n');
    fprintf('(Using score mode - will create prep_basics.mat without generating nulls)\n\n');
    
    % Run cerd_main - it will create prep_basics.mat and generate nulls
    try
        cerd_main;
        fprintf('\n✓ cerd_main completed successfully\n');
    catch ME_cerd
        % Check if prep_basics.mat was created despite the error
        if exist(prep_file, 'file')
            fprintf('\n✓ prep_basics.mat was created successfully!\n');
            fprintf('  (Note: cerd_main encountered an error, but prep_basics.mat was created)\n\n');
        else
            % If prep_basics.mat wasn't created, show the error
            fprintf('\nError: prep_basics.mat was not created. Error was:\n');
            fprintf('  %s\n', ME_cerd.message);
            rethrow(ME_cerd);
        end
    end
    
    % Restore original settings
    if ~isempty(old_mode)
        evalin('base', sprintf('CERD_MODE = ''%s'';', old_mode));
    else
        evalin('base', 'clear CERD_MODE');
    end
    if ~isempty(old_patient)
        evalin('base', sprintf('CERD_PATIENT = ''%s'';', old_patient));
    else
        evalin('base', 'clear CERD_PATIENT');
    end
    if ~isempty(old_event_permute)
        evalin('base', sprintf('CERD_EVENT_PERMUTE = %d;', old_event_permute));
    else
        evalin('base', 'clear CERD_EVENT_PERMUTE');
    end
    
    % Check if prep_basics.mat was created
    if exist(prep_file, 'file')
        fprintf('✓ Successfully created prep_basics.mat at:\n  %s\n\n', prep_file);
        fprintf('You can now run: edit_ufos_interactive(''%s'')\n', patient);
        
        % Clean up temporary file if it didn't exist before
        if ~temp_ufo_exists
            delete(temp_ufo_file);
            fprintf('(Removed temporary edited UFO file)\n');
        end
    else
        error('prep_basics.mat was not created. Check errors above.');
    end
    
catch ME
    % Restore original settings
    if ~isempty(old_mode)
        evalin('base', sprintf('CERD_MODE = ''%s'';', old_mode));
    else
        evalin('base', 'clear CERD_MODE');
    end
    if ~isempty(old_patient)
        evalin('base', sprintf('CERD_PATIENT = ''%s'';', old_patient));
    else
        evalin('base', 'clear CERD_PATIENT');
    end
    if ~isempty(old_event_permute)
        evalin('base', sprintf('CERD_EVENT_PERMUTE = %d;', old_event_permute));
    else
        evalin('base', 'clear CERD_EVENT_PERMUTE');
    end
    
    % Check if prep_basics.mat was created despite the error
    if exist(prep_file, 'file')
        fprintf('\n✓ prep_basics.mat was created despite the error.\n');
        fprintf('You can now run: edit_ufos_interactive(''%s'')\n', patient);
    else
        rethrow(ME);
    end
end

% Restore original settings
if ~isempty(old_mode)
    evalin('base', sprintf('CERD_MODE = ''%s'';', old_mode));
else
    evalin('base', 'clear CERD_MODE');
end
if ~isempty(old_patient)
    evalin('base', sprintf('CERD_PATIENT = ''%s'';', old_patient));
else
    evalin('base', 'clear CERD_PATIENT');
end
if ~isempty(old_event_permute)
    evalin('base', sprintf('CERD_EVENT_PERMUTE = %d;', old_event_permute));
else
    evalin('base', 'clear CERD_EVENT_PERMUTE');
end

end


function UFOdat = load_edited_ufos(patient, data_dir)
% LOAD_EDITED_UFOS  Load edited UFOs for use in the analysis pipeline
%
% Usage:
%   UFOdat = load_edited_ufos('pat001')
%   UFOdat = load_edited_ufos('pat001', data_dir)
%
% Returns UFOdat in the same format as extract_UFOs() would return,
% ready to use in cerd_main.m pipeline.
%
% If edited UFOs don't exist, returns empty and warns.

if nargin < 1 || isempty(patient)
    error('load_edited_ufos: patient ID required');
end

if nargin < 2 || isempty(data_dir)
    [base_dir, ~] = ufo_config();
    data_dir = fullfile(base_dir, 'Data', patient, 'UFO_edited');
end

edited_file = fullfile(data_dir, 'UFOdat_edited.mat');

if ~exist(edited_file, 'file')
    warning('load_edited_ufos: Edited UFO file not found at %s', edited_file);
    warning('load_edited_ufos: Returning empty. Will use original UFOs from Dat.UFOs');
    UFOdat = [];
    return;
end

try
    Edited = load(edited_file);
    if ~isfield(Edited, 'UFOdat_edited')
        error('UFOdat_edited field not found in %s', edited_file);
    end
    UFOdat = Edited.UFOdat_edited;
    
    % Verify structure
    if ~iscell(UFOdat)
        error('UFOdat_edited must be a cell array');
    end
    
    % Verify each cell contains Nx3 or Nx4 matrix
    % Nx4 format: [start_us, end_us, freq_hz, type] where type: 0=unclassified, 1=Type1<2kHz, 2=Type1>2kHz, 3=Type2<2kHz, 4=Type2>2kHz
    for ii = 1:numel(UFOdat)
        if ~isempty(UFOdat{ii})
            n_cols = size(UFOdat{ii}, 2);
            if n_cols ~= 3 && n_cols ~= 4
                error('UFOdat_edited{%d} must have 3 or 4 columns [start_us, end_us, freq_hz, type]', ii);
            end
        end
    end
    
    % If 3 columns, add 4th column (type, default 0 = unclassified) for backward compatibility
    for ii = 1:numel(UFOdat)
        if ~isempty(UFOdat{ii}) && size(UFOdat{ii}, 2) == 3
            if isinteger(UFOdat{ii})
                type_col = zeros(size(UFOdat{ii}, 1), 1, 'like', UFOdat{ii});
            else
                type_col = zeros(size(UFOdat{ii}, 1), 1, 'int8');
            end
            UFOdat{ii} = [UFOdat{ii}, type_col];
        end
    end
    
    fprintf('load_edited_ufos: Loaded edited UFOs from %s\n', edited_file);
    fprintf('load_edited_ufos: %d micro channels, %d total UFOs\n', ...
        numel(UFOdat), sum(cellfun(@(u) size(u,1), UFOdat)));
    
catch ME
    error('load_edited_ufos: Failed to load edited UFOs: %s', ME.message);
end

end


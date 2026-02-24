function [base_dir, matmef_path] = ufo_config()
% UFO_CONFIG  Central configuration for paths. Edit this file to match your setup.
%
% Expected directory structure (relative to BASE_DIR):
%   BASE_DIR/
%     Data/<patient>/           - MEF3 data (sub.mefd or *.mefd), UFO.json
%     Data/<patient>/UFO_edited/ - Edited UFOs (UFOdat_edited.mat)
%     Results/<patient>/         - ufo_inspect_log.mat, prep_basics.mat, UFO_figures/
%
% Returns:
%   base_dir    - Root directory for Data/ and Results/
%   matmef_path - Path to matmef library (for readMef3)

% === USER CONFIGURATION ===
% Default: parent of Code directory. Change to your data root.
base_dir = fullfile(fileparts(mfilename('fullpath')), '..');
matmef_path = fullfile(getenv('HOME'), 'matmef');  % Or e.g. '/path/to/matmef'

end

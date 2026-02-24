function [meta, data] = readMef3_python(mef_path, password, channel, varargin)
% READMEF3_PYTHON  Read MEF files using Python/pymef (faster for pat448)
%
% This is a drop-in replacement for readMef3 that uses Python for pat448.
% For other patients, it falls back to the original readMef3.
%
% Usage:
%   [meta, data] = readMef3_python(mef_path, password, channel, 'time', start_us, end_us)
%
% This function is specifically designed for pat448 where readMef3 is very slow.
% It calls a Python script that uses pymef library for fast MEF reading.

% Parse arguments
time_mode = false;
start_us = [];
end_us = [];

if nargin < 3
    error('readMef3_python: requires at least mef_path, password, and channel');
end

% Parse varargin for 'time' mode
if numel(varargin) >= 3 && strcmp(varargin{1}, 'time')
    time_mode = true;
    start_us = varargin{2};
    end_us = varargin{3};
end

% Get script directory
script_dir = fileparts(mfilename('fullpath'));
python_script = fullfile(script_dir, 'read_mef_python.py');
venv_python = fullfile(script_dir, 'venv_mef', 'bin', 'python3');

% Check if venv exists, otherwise use system python3
if exist(venv_python, 'file')
    python_cmd = venv_python;
else
    python_cmd = 'python3';
end

% Check if Python script exists
if ~exist(python_script, 'file')
    error('readMef3_python: Python script not found at %s', python_script);
end

if ~time_mode
    error('readMef3_python: Only ''time'' mode is currently supported');
end

% Call Python script
% Convert to int64 to ensure proper integer formatting (no scientific notation)
start_us_int = int64(start_us);
end_us_int = int64(end_us);
cmd = sprintf('%s "%s" "%s" "%s" %ld %ld "%s"', ...
    python_cmd, python_script, mef_path, channel, start_us_int, end_us_int, password);

[status, output] = system(cmd);

if status ~= 0
    warning('readMef3_python: Python script failed. Output: %s', output);
    meta = [];
    data = [];
    return;
end

% Parse output
lines = strsplit(strtrim(output), '\n');
if isempty(lines) || isempty(strtrim(lines{1}))
    warning('readMef3_python: Unexpected output format from Python script');
    meta = [];
    data = [];
    return;
end

% Bounds-only request (start_us==0, end_us==0): output is two lines = session_start, session_end
if start_us == 0 && end_us == 0
    try
        t_min = str2double(strtrim(lines{1}));
        t_max = str2double(strtrim(lines{2}));
        if isnan(t_min) || isnan(t_max)
            meta = [];
            data = [];
            return;
        end
        meta = struct();
        meta.session_start_time = t_min;
        meta.session_end_time = t_max;
        data = [];
        return;
    catch
        meta = [];
        data = [];
        return;
    end
end

% Data request: first line is sampling frequency, second line is data
if numel(lines) < 2
    warning('readMef3_python: Unexpected output format from Python script');
    meta = [];
    data = [];
    return;
end

try
    fs = str2double(lines{1});
    data_str = lines{2};
    data = str2num(data_str);  %#ok<ST2NM>
    
    if isempty(data)
        meta = [];
        data = [];
        return;
    end
    
    % Convert to column vector
    data = data(:);
    
    % Create minimal meta structure (compatible with readMef3 output)
    meta = struct();
    meta.sampling_frequency = fs;
    meta.start_time = start_us;
    meta.end_time = end_us;
    meta.channel = channel;
    meta.num_samples = numel(data);
    
catch ME
    warning('readMef3_python: Error parsing Python output: %s', ME.message);
    meta = [];
    data = [];
end

end


#!/usr/bin/env python3
"""
Python function to read MEF files - callable from MATLAB.

This script reads MEF data using pymef (much faster than MATLAB's readMef3 for pat448).

Usage from MATLAB:
    [status, data] = system(sprintf('python3 read_mef_python.py "%s" "%s" %d %d', ...
        mef_path, channel, start_us, end_us));
    if status == 0
        data = str2num(data);  % Parse the output
    end

Or use as standalone:
    python read_mef_python.py <mef_path> <channel> <start_us> <end_us> [password]
"""

import sys
import os
import json
import numpy as np

# Add the venv to path if it exists
script_dir = os.path.dirname(os.path.abspath(__file__))
venv_python = os.path.join(script_dir, 'venv_mef', 'bin', 'python3')
if os.path.exists(venv_python):
    # If running from system Python, we need to use venv's Python
    # But if pymef is already available, use it
    pass

try:
    from pymef.mef_session import MefSession
except ImportError:
    print("ERROR: pymef not found. Make sure venv_mef is activated or pymef is installed.", file=sys.stderr)
    sys.exit(1)

def read_mef_channel(mef_path, channel_name, start_us, end_us, password='bemena'):
    """
    Read MEF channel data for a specific time range.
    
    Args:
        mef_path: Path to .mefd directory
        channel_name: Channel name (e.g., 'ADMacro_01')
        start_us: Start time in microseconds (uutc)
        end_us: End time in microseconds (uutc)
        password: MEF password (default: 'bemena')
    
    Returns:
        data: numpy array of samples, or None on error
        fs: sampling frequency in Hz, or None on error
    """
    try:
        # Initialize session
        ms = MefSession(mef_path, password)
        
        # Read data
        data = ms.read_ts_channels_uutc(channel_name, [[int(start_us), int(end_us)]])
        
        if data is None or len(data) == 0:
            return None, None
        
        # Extract samples (data is a list of arrays, one per channel)
        samples = data[0] if isinstance(data, (list, tuple)) and len(data) > 0 else data
        
        if not isinstance(samples, np.ndarray) or len(samples) == 0:
            return None, None
        
        # Calculate sampling frequency
        duration_sec = (end_us - start_us) / 1e6
        if duration_sec > 0:
            fs = len(samples) / duration_sec
        else:
            fs = None
        
        return samples, fs
        
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return None, None

def get_session_bounds(mef_path, password='bemena'):
    """
    Return (session_start_us, session_end_us) from session metadata.
    """
    ms = MefSession(mef_path, password)
    md = ms.session_md['session_specific_metadata']
    t_min = int(md['earliest_start_time'][0])
    t_max = int(md['latest_end_time'][0])
    return t_min, t_max


def main():
    """Main entry point for command-line usage."""
    if len(sys.argv) < 5:
        print("Usage: python read_mef_python.py <mef_path> <channel> <start_us> <end_us> [password]", file=sys.stderr)
        sys.exit(1)
    
    mef_path = sys.argv[1]
    channel_name = sys.argv[2]
    # Handle scientific notation by converting to float first, then int
    start_us = int(float(sys.argv[3]))
    end_us = int(float(sys.argv[4]))
    password = sys.argv[5] if len(sys.argv) > 5 else 'bemena'
    
    # Bounds-only request: return session start/end (two lines) and exit
    if start_us == 0 and end_us == 0:
        try:
            t_min, t_max = get_session_bounds(mef_path, password)
            print(t_min)
            print(t_max)
            sys.exit(0)
        except Exception as e:
            print(f"ERROR: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc(file=sys.stderr)
            sys.exit(1)
    
    data, fs = read_mef_channel(mef_path, channel_name, start_us, end_us, password)
    
    if data is None:
        # Output empty result (MATLAB will parse as empty)
        print("[]")
        sys.exit(1)
    
    # Output data as space-separated values (MATLAB can parse with str2num)
    # Format: first line is sampling frequency, second line is data
    print(f"{fs:.6f}")
    # Print data in a format MATLAB can easily parse
    # Use scientific notation for precision
    data_str = ' '.join([f"{x:.6e}" for x in data])
    print(data_str)
    
    sys.exit(0)

if __name__ == "__main__":
    main()


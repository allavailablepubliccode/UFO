# UFO Analysis Pipeline

This code accompanies the paper **"Cross-Scale Echoes of Ultra-Fast Oscillations in the Human Brain."**

The pipeline has one entry point: **`ufo_master.m`**. It presents a menu to select patients and steps:

## Pipeline Steps

- **Step 1** (`step1_inspect.m`): Interactive UFO inspection and classification. Loads MEF3 data, displays each detected UFO, and allows manual classification into Type I (rhythmic) or Type II (artifact). Helper functions `read_ufo_window_simple` and `edit_single_ufo_reedit` are defined locally within this file.

- **Step 2** (`step2_merge.m`): Merges overlapping UFO detections within each microchannel and refreshes figures.

- **Step 3** (`step3_statistics_fixed.m`): Statistical analysis using fixed frequency-domain Green's envelope template matching. Contains the core detection function `greens_envelope_score` and FDR correction function `fdr_bh_local` as local functions.

## Supporting Files

- `ufo_config.m` — central path configuration (edit `base_dir` and `matmef_path`)
- `create_prep_for_editing.m` — prepares data structures for UFO editing
- `load_edited_ufos.m` — loads manually edited UFO classifications
- `edit_ufos_interactive.m` — interactive UFO review/classification tool
- `ufo_edit_master.m` — higher-level editing workflow with logging
- `export_ufo_figs.m` — exports UFO visualisations
- `export_fdr_pair_ufo_figs.m` — exports figures for FDR-significant pairs
- `readMef3_python.m` + `read_mef_python.py` — MEF3 data reading (requires MEF Python library)

## Analysis Parameters (paper defaults)

Set in `step3_statistics_fixed.m`:

| Parameter | Value |
|-----------|-------|
| Frequency band | [60, 100] Hz |
| Template centre frequency | f₀ = 75 Hz |
| Damping | γ = 40 s⁻¹ |
| Baseline window | 200 ms pre-UFO |
| Post-UFO window | 200 ms |
| Welch PSD | 100 ms segments, 50% overlap, minimum 512-point FFT |
| Null iterations | 10,000 |
| Multiple comparisons | Benjamini-Hochberg FDR (q = 0.05) |
| Test | One-sided (post > baseline) |
| UFO selection | Type I only, peak frequency > 2 kHz |

## Requirements

- MATLAB (tested on R2023b+)
- Signal Processing Toolbox
- Parallel Computing Toolbox (optional, for faster null generation)
- Python 3 with `numpy` and `pymef` for MEF3 reading

## Configuration

Edit `ufo_config.m` to set `base_dir` (data root) and `matmef_path` (MATLAB MEF library). The default `base_dir` is the parent of the Code directory. In `step3_statistics_fixed.m`, set `ANALYSIS_CONFIG.force_python_mef = true` if MATLAB readMef3 fails on your platform.

## .gitignore

The repository excludes:

- `venv_mef/`
- `.DS_Store`
- `__MACOSX/`
- `*.mat`
- `Results/`

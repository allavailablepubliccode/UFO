function export_fdr_pair_ufo_figs()
% EXPORT_FDR_PAIR_UFO_FIGS  Export manuscript-style figures for FDR-surviving pairs.
%
% Exports "micro on top / macro on bottom + PSD inset (baseline vs post-UFO) + template inset"
% figures to PNG for ALL Type1_HighFreq UFOs that survived FDR in step3_statistics_fixed.
%
% FDR-surviving pairs (hardcoded):
%   pat003: micro CSC73 with macro SMAL1
%   pat004: micro CSC65 with macro AMG1
%
% Output:
%   BASE_DIR/Results/UFO_figures/Type1_HighFreq_FDRpairs/
%
% Filename format: pat003_CSC73_SMAL1_ufoIdx0001.png
%
% Usage:
%   export_fdr_pair_ufo_figs
%
% Requires: readMef3 and ufo_inspect_log.mat. Uses MEF directly; does not require prep_basics.mat.

assignin('base', 'CERD_MODE', 'export_fdr');
cerd_main;

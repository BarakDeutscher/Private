function [bpTable, diagnostics] = generatePSDBreakpointTable(f, psd, varargin)
%GENERATEPSDBREAKPOINTTABLE Reduce a linear-spaced PSD to a breakpoint table.
%
%   [bpTable, diagnostics] = GENERATEPSDBREAKPOINTTABLE(f, psd) reduces a
%   linearly frequency-spaced PSD vibration spectrum (f in Hz, psd in
%   units^2/Hz) to a compact breakpoint (frequency, PSD) table of the
%   kind used to define random-vibration test/tolerance curves per
%   MIL-STD-810 - a small number of "breakpoints" joined by straight
%   line segments on log-log axes. The new curve is built to ENVELOPE
%   the original spectrum (at or above it everywhere, plus an optional
%   margin) using as few points as possible, subject to a user-set
%   maximum point count.
%
%   [...] = GENERATEPSDBREAKPOINTTABLE(f, psd, 'Name', Value, ...)
%
%   Name-Value options:
%     'MaxPoints'      Maximum number of breakpoints allowed (default 12).
%     'MarginDB'       Coverage margin applied above the original PSD, in
%                       dB, before the envelope is fit (default 3 dB; use
%                       0 dB for a tight envelope with no headroom).
%     'TargetRmsRatio' Desired upper bound on grms_new/grms_orig (default
%                       1.4). This is not enforced by construction -
%                       diagnostics.meetsRmsTarget reports whether it was
%                       actually achieved, with a warning if not (reduce
%                       MarginDB and/or raise MaxPoints to fix it).
%     'Interactive'    If true, opens a GUI table editor after the
%                       automatic table is built, letting you drag
%                       breakpoints, add/delete rows, and see the Grms
%                       ratio and coverage margin update live
%                       (default false).
%     'Plot'           If true (default), plots original vs breakpoint
%                       curve on log-log axes.
%
%   Outputs:
%     bpTable      MATLAB table, variables Frequency_Hz and PSD, sorted
%                  by increasing frequency.
%     diagnostics  Struct: grmsOriginal, grmsBreakpoint, rmsRatio,
%                  targetRmsRatio, meetsRmsTarget, marginDB, minMarginDB,
%                  numPoints, maxPoints.
%
%   To edit the table by hand after the fact (or re-open the editor on a
%   table you already exported), call:
%       bpTable = editBreakpointTableGUI(f, psd, bpTable, targetRmsRatio, marginDB);
%       diagnostics = evaluateBreakpointTable(f, psd, bpTable, targetRmsRatio, marginDB);
%
%   Background: MIL-STD-810 vibration test tolerance/breakpoint curves
%   are always specified as straight-line segments on log(PSD) vs
%   log(frequency) axes. This function builds a curve in exactly that
%   form, and computes Grms with the standard closed-form log-log segment
%   integral (see segmentMeanSquare.m) rather than linear trapezoidal
%   integration, which would misstate the area under log-log segments.

p = inputParser;
addRequired(p, 'f', @(v) isvector(v) && isnumeric(v));
addRequired(p, 'psd', @(v) isvector(v) && isnumeric(v));
addParameter(p, 'MaxPoints', 12, @(v) isscalar(v) && v >= 2);
addParameter(p, 'MarginDB', 3, @(v) isscalar(v) && isnumeric(v));
addParameter(p, 'TargetRmsRatio', 1.4, @(v) isscalar(v) && v > 1);
addParameter(p, 'Interactive', false, @(v) isscalar(v));
addParameter(p, 'Plot', true, @(v) isscalar(v));
parse(p, f, psd, varargin{:});

f = p.Results.f(:);
psd = p.Results.psd(:);
maxPoints = round(p.Results.MaxPoints);
marginDB = p.Results.MarginDB;
targetRmsRatio = p.Results.TargetRmsRatio;
doInteractive = logical(p.Results.Interactive);
doPlot = logical(p.Results.Plot);

if numel(f) ~= numel(psd)
    error('generatePSDBreakpointTable:sizeMismatch', ...
        'f and psd must have the same number of elements.');
end
if numel(f) < 3
    error('generatePSDBreakpointTable:tooFewPoints', ...
        'Need at least 3 points in the input spectrum.');
end

[f, order] = sort(f);
psd = psd(order);
if any(diff(f) <= 0)
    error('generatePSDBreakpointTable:duplicateFrequencies', ...
        'Frequency vector must have strictly increasing, unique values.');
end
if any(psd <= 0)
    error('generatePSDBreakpointTable:nonPositivePSD', ...
        'PSD values must be > 0 (breakpoint curves live on a log axis).');
end

dF = diff(f);
if (max(dF) - min(dF)) / mean(dF) > 0.05
    warning('generatePSDBreakpointTable:notLinearSpacing', ...
        ['Input frequency vector does not look linearly spaced (spacing ' ...
         'varies by more than 5%%). The reduction still works, but it ' ...
         'assumes the given samples are dense/representative everywhere.']);
end

% ---- original composite Grms ------------------------------------------
grmsOriginal = computeOriginalGrms(f, psd);

% ---- margined target and minimal enveloping breakpoint set ------------
marginFactor = 10^(marginDB/10);
psdTarget = psd * marginFactor;

% A uniform dB margin alone raises grms by exactly sqrt(marginFactor),
% since scaling a PSD by a constant scales its integral by that constant.
% That is therefore a hard floor on the achievable ratio for ANY curve
% that fully envelopes psdTarget (any point reduction only pushes the
% ratio up further, never below this floor). Warn early if the request
% is self-contradictory.
minAchievableRatio = sqrt(marginFactor);
if targetRmsRatio < minAchievableRatio
    maxMarginForTarget = 20*log10(targetRmsRatio);
    warning('generatePSDBreakpointTable:targetUnreachable', ...
        ['MarginDB = %.2f dB alone forces grms ratio >= %.3f (a uniform dB ' ...
         'margin scales grms by sqrt(10^(MarginDB/10))), which already ' ...
         'exceeds TargetRmsRatio = %.3f - no number of breakpoints can ' ...
         'reach the target while keeping full coverage. Lower MarginDB ' ...
         'below %.2f dB, or raise TargetRmsRatio, or accept a table that ' ...
         'does not fully envelope the original.'], ...
        marginDB, minAchievableRatio, targetRmsRatio, maxMarginForTarget);
end

hullIdx = upperConvexHullLogLog(f, psdTarget);

if numel(hullIdx) > maxPoints
    warning('generatePSDBreakpointTable:pointBudgetExceeded', ...
        ['A full envelope needs %d breakpoints but MaxPoints = %d was ' ...
         'requested; reducing the point count will leave small regions ' ...
         'of the original spectrum uncovered by the requested margin.'], ...
        numel(hullIdx), maxPoints);
    x = log10(f);
    y = log10(psdTarget);
    hullIdx = reduceBreakpointsVW(x, y, hullIdx, maxPoints);
end

bpFreq = f(hullIdx);
bpPsd  = psdTarget(hullIdx);
bpTable = table(bpFreq, bpPsd, 'VariableNames', {'Frequency_Hz', 'PSD'});

% ---- diagnostics (grms ratio + worst-case coverage margin) ------------
diagnostics = evaluateBreakpointTable(f, psd, bpTable, targetRmsRatio, marginDB);
diagnostics.maxPoints = maxPoints;

if ~diagnostics.meetsRmsTarget
    warning('generatePSDBreakpointTable:rmsRatioExceeded', ...
        ['Grms ratio (new/original) = %.3f exceeds TargetRmsRatio = %.3f. ' ...
         'Try reducing MarginDB, increasing MaxPoints, or edit the table ' ...
         'manually.'], diagnostics.rmsRatio, targetRmsRatio);
end
if diagnostics.minMarginDB < -1e-6
    warning('generatePSDBreakpointTable:coverageViolated', ...
        ['Breakpoint curve dips %.2f dB BELOW the original spectrum at ' ...
         'at least one frequency (point budget too small for full ' ...
         'coverage). Increase MaxPoints for a guaranteed envelope.'], ...
        -diagnostics.minMarginDB);
end

fprintf('--- PSD Breakpoint Reduction Summary ---\n');
fprintf('  Original Grms         : %.4f\n', diagnostics.grmsOriginal);
fprintf('  Breakpoint Grms       : %.4f\n', diagnostics.grmsBreakpoint);
fprintf('  Grms ratio (new/orig) : %.4f  (target <= %.3f)  %s\n', ...
    diagnostics.rmsRatio, targetRmsRatio, tern(diagnostics.meetsRmsTarget, 'OK', 'EXCEEDED'));
fprintf('  Min coverage margin   : %.2f dB (requested %.2f dB)\n', ...
    diagnostics.minMarginDB, marginDB);
fprintf('  Breakpoints used      : %d / %d max\n', diagnostics.numPoints, maxPoints);

if doPlot
    plotPSDBreakpoints(f, psd, bpTable, diagnostics);
end

if doInteractive
    bpTable = editBreakpointTableGUI(f, psd, bpTable, targetRmsRatio, marginDB);
    diagnostics = evaluateBreakpointTable(f, psd, bpTable, targetRmsRatio, marginDB);
    diagnostics.maxPoints = maxPoints;
end
end

function s = tern(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end

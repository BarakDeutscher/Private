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
%     'MarginMode'     'uniform' (default): every point gets the same
%                       MarginDB. 'adaptive': the margin is NOT a flat
%                       factor everywhere - it varies linearly (in dB,
%                       vs. normalized log-amplitude) between MarginDB at
%                       the spectrum's quietest point and 0 dB at its
%                       loudest peak. Since the loudest peak dominates
%                       the Grms integral the most, leaving it un-inflated
%                       cuts the ratio a lot more than uniformly shaving
%                       margin everywhere would, while the quiet floor -
%                       which barely affects Grms regardless of its
%                       margin - keeps the full requested headroom.
%     'PeakBracketing' If true (default false), significant narrowband
%                       peaks are represented by exactly two flat, equal-
%                       amplitude breakpoints straddling each one, at
%                       PeakFreqMarginFraction on either side of the peak
%                       frequency, held at the peak's own RAW psd value
%                       (no amplitude margin there at all) - the standard
%                       way a resonance is broadened in derived vibration
%                       test breakpoint tables to cover uncertainty in
%                       exactly where it sits, rather than by amplitude
%                       inflation. Applied before the normal margin
%                       scheme's hull/trim/refine steps, so it can reduce
%                       point count a lot around real resonances, at the
%                       cost of some excess area in the peak's flanks
%                       (a flat cap across a wide band overshoots more
%                       than tightly tracing the peak's true shape would)
%                       - if TargetRmsRatio isn't met, narrow
%                       PeakFreqMarginFraction or raise MaxPoints.
%     'PeakFreqMarginFraction'  Half-width of each peak bracket as a
%                       fraction of its frequency (default 0.10, i.e.
%                       +/-10%%). Only used when PeakBracketing is true.
%     'PeakMinProminenceDB'  Minimum topographic prominence, in dB, for a
%                       local maximum to count as a "peak" worth
%                       bracketing (default 6 dB). Only used when
%                       PeakBracketing is true.
%     'TargetRmsRatio' Desired upper bound on grms_new/grms_orig (default
%                       1.4). If the minimal fully-enveloping breakpoint
%                       set already fits within MaxPoints but overshoots
%                       this target, spare points are spent greedily
%                       (see refineBreakpointsGreedy.m) narrowing the gap
%                       until the target is met, MaxPoints runs out, or
%                       no further insertion can reduce the ratio (a true
%                       local optimum - not every MarginDB/MaxPoints
%                       combination can reach every target).
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
%     'OutputFolder'   If non-empty, saves the FINAL results there (the
%                       post-edit table/plot if 'Interactive' was used):
%                       breakpoint_table.csv and
%                       psd_breakpoint_comparison.png (original PSD vs
%                       breakpoint curve, log-log). The folder is created
%                       if it does not exist. Default '' (nothing saved).
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
addParameter(p, 'MarginMode', 'uniform', @(v) any(strcmpi(v, {'uniform','adaptive'})));
addParameter(p, 'PeakBracketing', false, @(v) isscalar(v));
addParameter(p, 'PeakFreqMarginFraction', 0.10, @(v) isscalar(v) && v > 0 && v < 1);
addParameter(p, 'PeakMinProminenceDB', 6, @(v) isscalar(v) && v > 0);
addParameter(p, 'TargetRmsRatio', 1.4, @(v) isscalar(v) && v > 1);
addParameter(p, 'Interactive', false, @(v) isscalar(v));
addParameter(p, 'Plot', true, @(v) isscalar(v));
addParameter(p, 'OutputFolder', '', @(v) ischar(v) || isstring(v));
parse(p, f, psd, varargin{:});

f = p.Results.f(:);
psd = p.Results.psd(:);
maxPoints = round(p.Results.MaxPoints);
marginDB = p.Results.MarginDB;
marginMode = lower(char(p.Results.MarginMode));
doPeakBracketing = logical(p.Results.PeakBracketing);
peakFreqMarginFraction = p.Results.PeakFreqMarginFraction;
peakMinProminenceDB = p.Results.PeakMinProminenceDB;
targetRmsRatio = p.Results.TargetRmsRatio;
doInteractive = logical(p.Results.Interactive);
doPlot = logical(p.Results.Plot);
outputFolder = char(p.Results.OutputFolder);

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
if strcmp(marginMode, 'uniform')
    marginFactor = 10^(marginDB/10);
    psdTarget = psd * marginFactor;

    % A uniform dB margin alone raises grms by exactly sqrt(marginFactor),
    % since scaling a PSD by a constant scales its integral by that
    % constant. That is therefore a hard floor on the achievable ratio
    % for ANY curve that fully envelopes psdTarget (any point reduction
    % only pushes the ratio up further, never below this floor). Warn
    % early if the request is self-contradictory.
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
else
    % Adaptive margin: linear in dB vs. normalized log-amplitude, from
    % MarginDB at the quietest sample (normalizedLevel=0) down to 0 dB at
    % the loudest peak (normalizedLevel=1). No simple closed-form floor
    % exists here (unlike the uniform case above), since the scaling is
    % no longer a single constant - it depends on how much of the
    % spectrum's Grms actually comes from near-peak vs. near-floor
    % content.
    psdDB = 10*log10(psd);
    minDB = min(psdDB);
    maxDB = max(psdDB);
    if maxDB > minDB
        normalizedLevel = (psdDB - minDB) / (maxDB - minDB);
    else
        normalizedLevel = ones(size(psd)); % flat spectrum: nothing to grade
    end
    marginDBLocal = marginDB * (1 - normalizedLevel);
    psdTarget = psd .* 10.^(marginDBLocal/10);
end

% ---- optional peak frequency bracketing --------------------------------
% Broadens significant narrowband peaks in FREQUENCY (a flat, unmargined
% plateau spanning +/-PeakFreqMarginFraction around each) rather than in
% amplitude, before the general hull/trim/refine machinery runs on
% whatever is left. Adds new synthetic anchor frequencies at the bracket
% edges, so fHull/psdTargetHull can have more points than the original
% (f, psd) - diagnostics are always computed against the true original
% spectrum regardless.
if doPeakBracketing
    [fHull, psdTargetHull] = applyPeakFrequencyBracketing(f, psd, psdTarget, ...
        peakFreqMarginFraction, peakMinProminenceDB);
else
    fHull = f;
    psdTargetHull = psdTarget;
end

hullIdx = upperConvexHullLogLog(fHull, psdTargetHull);

if numel(hullIdx) > maxPoints
    warning('generatePSDBreakpointTable:pointBudgetExceeded', ...
        ['A full envelope needs %d breakpoints but MaxPoints = %d was ' ...
         'requested; reducing the point count will leave small regions ' ...
         'of the original spectrum uncovered by the requested margin.'], ...
        numel(hullIdx), maxPoints);
    x = log10(fHull);
    y = log10(psdTargetHull);
    hullIdx = reduceBreakpointsVW(x, y, hullIdx, maxPoints);
elseif numel(hullIdx) < maxPoints
    % The minimal fully-enveloping hull already fits the point budget -
    % but fewest-points-for-coverage is NOT the same as lowest-area (a
    % point the hull leaves out because it's dominated by its neighbors'
    % straight line can still be worth adding: it lets the curve dip down
    % and track a real local valley between two much higher, distant
    % features, cutting Grms, as long as the two new sub-segments it
    % creates still cover everything in their own sub-ranges). Spend the
    % spare point budget on this before giving up on TargetRmsRatio.
    targetGrms = targetRmsRatio * grmsOriginal;
    hullIdx = refineBreakpointsGreedy(fHull, psdTargetHull, hullIdx, maxPoints, targetGrms);
end

bpFreq = fHull(hullIdx);
bpPsd  = psdTargetHull(hullIdx);
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
if strcmp(marginMode, 'uniform')
    fprintf('  Min coverage margin   : %.2f dB (requested %.2f dB)\n', ...
        diagnostics.minMarginDB, marginDB);
else
    fprintf(['  Min coverage margin   : %.2f dB (adaptive: %.2f dB at the ' ...
        'quiet floor tapering to 0 dB at the peak)\n'], diagnostics.minMarginDB, marginDB);
end
fprintf('  Breakpoints used      : %d / %d max\n', diagnostics.numPoints, maxPoints);

if doPlot
    plotPSDBreakpoints(f, psd, bpTable, diagnostics);
end

if doInteractive
    bpTable = editBreakpointTableGUI(f, psd, bpTable, targetRmsRatio, marginDB);
    diagnostics = evaluateBreakpointTable(f, psd, bpTable, targetRmsRatio, marginDB);
    diagnostics.maxPoints = maxPoints;
end

if ~isempty(outputFolder)
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    tablePath = fullfile(outputFolder, 'breakpoint_table.csv');
    writetable(bpTable, tablePath);

    % Always render a fresh plot of the FINAL table (post-edit, if any)
    % for saving, independent of the 'Plot' display option above.
    saveFig = plotPSDBreakpoints(f, psd, bpTable, diagnostics);
    plotPath = fullfile(outputFolder, 'psd_breakpoint_comparison.png');
    print(saveFig, plotPath, '-dpng', '-r150');
    close(saveFig);

    fprintf('  Saved breakpoint table -> %s\n', tablePath);
    fprintf('  Saved comparison plot  -> %s\n', plotPath);
end
end

function s = tern(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end

%% Sanity tests for the psdBreakpoints library (run directly, no toolbox needed)
% Checks the closed-form log-log segment integral against numerical
% integration, the convex-hull envelope property, and the
% Visvalingam-Whyatt point-budget reduction. Run with:
%   run('tests/test_psd_breakpoints.m')

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'psdBreakpoints'));
tol = 1e-6;

%% 1. computeOriginalGrms: constant PSD has an exact analytic answer
f = linspace(10, 110, 1001)';
psd = 0.05 * ones(size(f));
grms = computeOriginalGrms(f, psd);
assert(abs(grms - sqrt(0.05*(110-10))) < tol, 'computeOriginalGrms failed on constant PSD');

%% 2. segmentMeanSquare vs. numerical integration of the underlying power law
cases = [10 0.01 80 0.08; 100 0.08 100.001 0.0800001; 50 0.02 51 0.019; 1 1 1000 1];
for i = 1:size(cases,1)
    f1 = cases(i,1); w1 = cases(i,2); f2 = cases(i,3); w2 = cases(i,4);
    ms = segmentMeanSquare(f1, w1, f2, w2);
    n = log(w2/w1) / log(f2/f1);
    powerlaw = @(ff) w1 * (ff/f1).^n;
    msNumeric = integral(powerlaw, f1, f2);
    assert(abs(ms - msNumeric) < 1e-4*max(abs(msNumeric), 1e-9) + 1e-8, ...
        sprintf('segmentMeanSquare mismatch on case %d', i));
end

%% 3. computeBreakpointGrms: multi-segment table vs. numerical integration
bpF = [10 80 350 2000]';
bpW = [0.01 0.08 0.08 0.01]';
grmsBP = computeBreakpointGrms(bpF, bpW);
msTotal = 0;
for k = 1:numel(bpF)-1
    f1 = bpF(k); w1 = bpW(k); f2 = bpF(k+1); w2 = bpW(k+1);
    n = log(w2/w1)/log(f2/f1);
    powerlaw = @(ff) w1 * (ff/f1).^n;
    msTotal = msTotal + integral(powerlaw, f1, f2);
end
assert(abs(grmsBP - sqrt(msTotal)) < 1e-4*sqrt(msTotal), 'computeBreakpointGrms mismatch');

%% 4. upperConvexHullLogLog: envelope property and peak preservation
xh = [0 1 2]; yh = [0 5 0];               % strong interior peak must survive
hidx = upperConvexHullLogLog(10.^xh, 10.^yh);
assert(numel(hidx) == 3 && hidx(2) == 2, 'hull dropped a required peak');

xh2 = [0 1 2]; yh2 = [0 0.5 2];           % middle point below the chord: redundant
hidx2 = upperConvexHullLogLog(10.^xh2, 10.^yh2);
assert(numel(hidx2) == 2, 'hull kept a redundant point');

%% 5. Full-hull envelope covers every margined point, and reduction keeps endpoints
fFull = (10:1:2000)';
psdFull = 0.01 + 0.07*exp(-((fFull-150).^2)/(2*30^2));
marginDB = 3;
psdTarget = psdFull * 10^(marginDB/10);
hullIdx = upperConvexHullLogLog(fFull, psdTarget);
curveAtOrig = interpBreakpointCurve(fFull(hullIdx), psdTarget(hullIdx), fFull);
assert(min(curveAtOrig - psdTarget) >= -1e-9*max(psdTarget), 'full hull failed to envelope all points');

maxPoints = 6;
reduced = reduceBreakpointsVW(log10(fFull), log10(psdTarget), hullIdx, maxPoints);
assert(numel(reduced) == maxPoints, 'reduceBreakpointsVW did not respect MaxPoints');
assert(reduced(1) == hullIdx(1) && reduced(end) == hullIdx(end), ...
    'reduceBreakpointsVW must never drop the first/last breakpoint');

%% 6. refineBreakpointsGreedy: spare point budget must reduce Grms ratio
% Regression test for a real finding: the minimal fully-enveloping hull
% minimizes POINT COUNT, not AREA. A sharp, distant peak sitting far
% above a quiet floor forces the hull's bridging segments to badly
% overshoot the floor - even though every floor sample is individually
% "dominated" (below the segment endpoints' straight line, so the hull
% correctly omits it), voluntarily adding some of those dominated points
% back in lets the curve track the floor and can cut Grms substantially,
% as long as coverage is never lost.
fSharp = (1:1:500)';
psdSharp = 1e-3 + 0.2*exp(-((fSharp-250).^2)/(2*3^2)); % floor ~1e-3, sharp peak ~0.2 at f=250
marginDB = 1;
psdTargetSharp = psdSharp * 10^(marginDB/10);
grmsOrigSharp = computeOriginalGrms(fSharp, psdSharp);

hullSharp = upperConvexHullLogLog(fSharp, psdTargetSharp);
ratioHull = computeBreakpointGrms(fSharp(hullSharp), psdTargetSharp(hullSharp)) / grmsOrigSharp;
assert(ratioHull > 1.4, 'test setup assumption failed: expected the minimal hull to badly overshoot here');

refined = refineBreakpointsGreedy(fSharp, psdTargetSharp, hullSharp, 20, 1.05*grmsOrigSharp);
ratioRefined = computeBreakpointGrms(fSharp(refined), psdTargetSharp(refined)) / grmsOrigSharp;

assert(ratioRefined < ratioHull, 'refineBreakpointsGreedy did not improve on the plain hull');
assert(numel(refined) <= 20, 'refineBreakpointsGreedy exceeded its point budget');

curveRefined = interpBreakpointCurve(fSharp(refined), psdTargetSharp(refined), fSharp);
assert(min(curveRefined - psdTargetSharp) >= -1e-9*max(psdTargetSharp), ...
    'refineBreakpointsGreedy lost full coverage while adding points');

%% 7. Adaptive margin: 0 dB at the peak, full MarginDB at the quiet floor
% Mirrors the 'adaptive' MarginMode logic in generatePSDBreakpointTable.m.
fAdapt = (1:1:500)';
psdAdapt = 1e-3 + 0.2*exp(-((fAdapt-250).^2)/(2*3^2)); % same shape as test 6
marginDBMax = 3;

psdDB = 10*log10(psdAdapt);
normalizedLevel = (psdDB - min(psdDB)) / (max(psdDB) - min(psdDB));
marginDBLocal = marginDBMax * (1 - normalizedLevel);
psdTargetAdapt = psdAdapt .* 10.^(marginDBLocal/10);

[~, peakPos] = max(psdAdapt);
[~, floorPos] = min(psdAdapt);
assert(abs(marginDBLocal(peakPos)) < 1e-9, 'adaptive margin at the peak should be exactly 0 dB');
assert(abs(marginDBLocal(floorPos) - marginDBMax) < 1e-9, ...
    'adaptive margin at the quietest point should equal the requested MarginDB');

% A uniform-margin hull always inflates the peak; the adaptive one should
% not, so it must achieve a strictly lower ratio for the same MarginDB.
grmsOrigAdapt = computeOriginalGrms(fAdapt, psdAdapt);
hullUniform = upperConvexHullLogLog(fAdapt, psdAdapt * 10^(marginDBMax/10));
ratioUniform = computeBreakpointGrms(fAdapt(hullUniform), (psdAdapt(hullUniform))*10^(marginDBMax/10)) / grmsOrigAdapt;
hullAdapt = upperConvexHullLogLog(fAdapt, psdTargetAdapt);
ratioAdapt = computeBreakpointGrms(fAdapt(hullAdapt), psdTargetAdapt(hullAdapt)) / grmsOrigAdapt;
assert(ratioAdapt < ratioUniform, 'adaptive margin should beat uniform margin on a peaky spectrum');

% Full coverage of the adaptive target must still hold.
curveAdapt = interpBreakpointCurve(fAdapt(hullAdapt), psdTargetAdapt(hullAdapt), fAdapt);
assert(min(curveAdapt - psdTargetAdapt) >= -1e-9*max(psdTargetAdapt), ...
    'adaptive-margin hull failed to envelope its own target');

%% 8. applyPeakFrequencyBracketing: flat, unmargined +/-10% peak brackets
% Single isolated narrow peak: should produce exactly one bracket, at the
% peak's own raw amplitude (no amplitude margin), and the hull built on
% the augmented arrays should collapse the whole peak to just the two
% bracket edges (the peak's own sample becomes redundant/interior).
fPk = (1:1:300)';
psdPk = 1e-3*ones(size(fPk)) + 0.05*exp(-((fPk-100).^2)/(2*3^2));
psdTargetPk = psdPk; % no separate margin scheme layered on top, for clarity

[fAugPk, psdTargetAugPk, peakListPk] = applyPeakFrequencyBracketing(fPk, psdPk, psdTargetPk, 0.10, 6);
assert(numel(peakListPk) == 1, 'expected exactly one bracketed peak');
assert(abs(peakListPk(1).freqLeft - 90) < 1e-9 && abs(peakListPk(1).freqRight - 110) < 1e-9, ...
    'peak bracket edges should be at +/-10%% of the peak frequency');
assert(abs(peakListPk(1).psd - psdPk(fPk==100)) < 1e-12, ...
    'peak bracket amplitude should exactly equal the raw (unmargined) peak PSD');

hullPk = upperConvexHullLogLog(fAugPk, psdTargetAugPk);
hullFreqsPk = fAugPk(hullPk);
assert(~any(abs(hullFreqsPk - 100) < 1e-9), ...
    'the peak''s own sample should become redundant once bracketed (dominated by the flat plateau)');
assert(any(abs(hullFreqsPk - 90) < 1e-9) && any(abs(hullFreqsPk - 110) < 1e-9), ...
    'both bracket edges should survive as hull vertices');

curvePk = interpBreakpointCurve(fAugPk(hullPk), psdTargetAugPk(hullPk), fPk);
assert(min(curvePk - psdPk) >= -1e-9*max(psdPk), ...
    'peak-bracketed hull failed to envelope the true original spectrum');

%% 9. applyPeakFrequencyBracketing: two close peaks must merge into one bracket
% Two peaks close enough that their +/-10% bands overlap must not produce
% two conflicting brackets; they should merge into a single wider one
% whose flat level is the true max of the WHOLE merged band (so coverage
% can never be violated, even for the lower of the two peaks' own band).
fTwo = (1:1:300)';
psdTwo = 1e-3*ones(size(fTwo)) + 0.05*exp(-((fTwo-100).^2)/(2*4^2)) + 0.03*exp(-((fTwo-108).^2)/(2*4^2));
[~, ~, peakListTwo] = applyPeakFrequencyBracketing(fTwo, psdTwo, psdTwo, 0.10, 6);
assert(numel(peakListTwo) == 1, 'two overlapping-band peaks should merge into a single bracket');
assert(peakListTwo(1).freqLeft < 100 && peakListTwo(1).freqRight > 108, ...
    'merged bracket should span both original peaks'' individual bands');

%% 10. pruneNegligibleBreakpoints: drop points that don't pull their weight
% Build a breakpoint set with one obviously-redundant point (sitting
% exactly on the straight line between its neighbors, in log-log space -
% zero-cost to remove) alongside a genuinely necessary peak. Pruning
% should drop the redundant one and keep the peak.
fPrune = (1:1:400)';
psdPrune = 1e-3*ones(size(fPrune)) + 0.2*exp(-((fPrune-200).^2)/(2*3^2));
grmsOrigPrune = computeOriginalGrms(fPrune, psdPrune);

% breakpoints: start, an exactly-collinear midpoint (redundant by
% construction), the peak, and the end.
idxStart = 1; idxMid = 51; idxPeak = find(fPrune==200); idxEnd = numel(fPrune);
xA = log10(fPrune(idxStart)); yA = log10(psdPrune(idxStart));
xC = log10(fPrune(idxPeak));  yC = log10(psdPrune(idxPeak));
xMid = log10(fPrune(idxMid));
yMidCollinear = yA + (yC-yA)*(xMid-xA)/(xC-xA);
psdTargetPrune = psdPrune;
psdTargetPrune(idxMid) = 10^yMidCollinear; % force it exactly onto the A-peak line

bpIdxPrune = [idxStart, idxMid, idxPeak, idxEnd];
grmsBeforePrune = computeBreakpointGrms(fPrune(bpIdxPrune), psdTargetPrune(bpIdxPrune));

prunedIdx = pruneNegligibleBreakpoints(fPrune, psdTargetPrune, bpIdxPrune, grmsOrigPrune, ...
    grmsBeforePrune/grmsOrigPrune + 0.01); % generous target: only true negligibility matters here

assert(~ismember(idxMid, prunedIdx), 'pruning failed to drop an exactly-redundant (collinear) breakpoint');
assert(ismember(idxPeak, prunedIdx), 'pruning must never drop a genuinely necessary peak');
assert(ismember(idxStart, prunedIdx) && ismember(idxEnd, prunedIdx), ...
    'pruning must never drop the first/last breakpoint');

curvePrune = interpBreakpointCurve(fPrune(prunedIdx), psdTargetPrune(prunedIdx), fPrune);
assert(min(curvePrune - psdTargetPrune) >= -1e-9*max(psdTargetPrune), ...
    'pruning must never break coverage');

fprintf('All psdBreakpoints sanity tests passed.\n');

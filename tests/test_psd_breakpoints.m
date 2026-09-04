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

fprintf('All psdBreakpoints sanity tests passed.\n');

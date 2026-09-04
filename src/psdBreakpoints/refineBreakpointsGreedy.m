function bpIdx = refineBreakpointsGreedy(f, psdTarget, bpIdx, maxPoints, targetGrms)
%REFINEBREAKPOINTSGREEDY Spend spare point budget to shrink the Grms ratio.
%   bpIdx = REFINEBREAKPOINTSGREEDY(f, psdTarget, bpIdx, maxPoints, targetGrms)
%   starts from a breakpoint set that already fully envelopes psdTarget
%   (e.g. the output of upperConvexHullLogLog) and greedily inserts
%   additional real data points - one at a time, always the single
%   insertion that reduces the enclosed mean-square area the most - until
%   either the resulting Grms reaches targetGrms, the point budget
%   maxPoints is used up, or no remaining candidate can improve things
%   further (a true local optimum).
%
%   IMPORTANT: minimizing the NUMBER of breakpoints needed for full
%   coverage (what upperConvexHullLogLog does) is a DIFFERENT problem
%   from minimizing the ENCLOSED AREA (Grms) for a given point BUDGET.
%   A point that upperConvexHullLogLog correctly leaves out (because it
%   lies on or below the straight line joining its neighbors, so it is
%   not NEEDED for coverage) can still be worth adding voluntarily: doing
%   so lets the curve dip down and track a real local valley between two
%   distant, much higher features, which can substantially cut the
%   enclosed area - as long as the two new sub-segments this creates
%   still stay above every other sample in their own sub-ranges. This
%   function checks that condition (locally, with interpBreakpointCurve)
%   before accepting any insertion, so coverage of psdTarget is never
%   lost.
%
%   Each insertion only changes the ONE segment it splits, so candidate
%   points are evaluated locally (against just the points between the
%   segment's own endpoints), not against the whole curve - this keeps
%   the search cheap even for large input spectra.

f = f(:);
psdTarget = psdTarget(:);
bpIdx = sort(bpIdx(:))';

currentGrms = computeBreakpointGrms(f(bpIdx), psdTarget(bpIdx));

while currentGrms > targetGrms && numel(bpIdx) < maxPoints
    bestReduction = 0;
    bestK = -1;

    for s = 1:numel(bpIdx)-1
        A = bpIdx(s);
        C = bpIdx(s+1);
        if C - A <= 1
            continue; % no original samples strictly between A and C
        end

        oldArea = segmentMeanSquare(f(A), psdTarget(A), f(C), psdTarget(C));

        for k = A+1:C-1
            trialF = [f(A); f(k); f(C)];
            trialW = [psdTarget(A); psdTarget(k); psdTarget(C)];
            localCurve = interpBreakpointCurve(trialF, trialW, f(A:C));
            if min(localCurve - psdTarget(A:C)) < -1e-9 * psdTarget(C)
                continue; % this split would uncover some other point in (A,C)
            end

            newArea = segmentMeanSquare(f(A), psdTarget(A), f(k), psdTarget(k)) + ...
                      segmentMeanSquare(f(k), psdTarget(k), f(C), psdTarget(C));
            reduction = oldArea - newArea;
            if reduction > bestReduction
                bestReduction = reduction;
                bestK = k;
            end
        end
    end

    if bestK < 0
        break; % no remaining insertion can reduce the area further
    end

    bpIdx = sort([bpIdx, bestK]);
    currentGrms = computeBreakpointGrms(f(bpIdx), psdTarget(bpIdx));
end
end

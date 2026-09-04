function bpIdx = pruneNegligibleBreakpoints(f, psdTarget, bpIdx, grmsOriginal, targetRatio)
%PRUNENEGLIGIBLEBREAKPOINTS Remove breakpoints that don't actually matter.
%   bpIdx = PRUNENEGLIGIBLEBREAKPOINTS(f, psdTarget, bpIdx, grmsOriginal,
%   targetRatio) repeatedly finds the interior breakpoint whose removal
%   costs the LEAST (smallest increase in enclosed area) and removes it,
%   as long as doing so (a) keeps the curve fully covering psdTarget
%   everywhere, and (b) keeps grms_new/grmsOriginal at or below
%   targetRatio. Stops once no remaining point can be removed without
%   breaking one of those two conditions.
%
%   This is the natural cleanup pass after refineBreakpointsGreedy (a
%   forward, insert-only greedy search): inserting points one at a time
%   in whatever order most helps AT THAT STEP can leave some of them
%   redundant once the full set is assembled - e.g. a later insertion
%   nearby ends up doing the same job, or a point was only needed to
%   inch the ratio just under target and neighboring points already
%   cover its span just fine. Every surviving breakpoint, after this
%   pass, is one whose removal would either break coverage or blow the
%   ratio - i.e. one that is actually doing something.

f = f(:);
psdTarget = psdTarget(:);
bpIdx = sort(bpIdx(:))';
targetGrms = targetRatio * grmsOriginal;

improved = true;
while improved
    improved = false;
    m = numel(bpIdx);
    if m <= 2
        break;
    end

    bestCost = inf;
    bestPos = -1;
    for i = 2:m-1
        A = bpIdx(i-1);
        B = bpIdx(i);
        C = bpIdx(i+1);

        oldArea = segmentMeanSquare(f(A), psdTarget(A), f(B), psdTarget(B)) + ...
                  segmentMeanSquare(f(B), psdTarget(B), f(C), psdTarget(C));
        newArea = segmentMeanSquare(f(A), psdTarget(A), f(C), psdTarget(C));
        cost = newArea - oldArea; % area increase from removing B (>= 0 normally)

        localCurve = interpBreakpointCurve([f(A); f(C)], [psdTarget(A); psdTarget(C)], f(A:C));
        if min(localCurve - psdTarget(A:C)) < -1e-9 * psdTarget(C)
            continue; % removing B would uncover some other point between A and C
        end

        if cost < bestCost
            bestCost = cost;
            bestPos = i;
        end
    end

    if bestPos < 0
        break; % no removable-without-losing-coverage candidate left
    end

    trialIdx = bpIdx([1:bestPos-1, bestPos+1:end]);
    trialGrms = computeBreakpointGrms(f(trialIdx), psdTarget(trialIdx));
    if trialGrms <= targetGrms
        bpIdx = trialIdx;
        improved = true;
    end
    % else: even the cheapest safe removal would exceed the ratio target,
    % so every remaining point is load-bearing - stop.
end
end

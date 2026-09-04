function [fAug, psdTargetAug, peakList] = applyPeakFrequencyBracketing(f, psd, psdTarget, freqMarginFraction, minProminenceDB)
%APPLYPEAKFREQUENCYBRACKETING Broaden narrowband peaks in frequency, not amplitude.
%   [fAug, psdTargetAug, peakList] = APPLYPEAKFREQUENCYBRACKETING(f, psd,
%   psdTarget, freqMarginFraction, minProminenceDB) finds significant
%   local maxima in psd and, for each, replaces the usual multi-point
%   tracking of its shape with exactly two new breakpoint anchors at
%   f_peak*(1-freqMarginFraction) and f_peak*(1+freqMarginFraction), both
%   held at the peak's own RAW, unmargined PSD value.
%
%   This is the standard way a narrowband resonance is represented in a
%   derived vibration test breakpoint table: broadened in FREQUENCY to
%   cover uncertainty in exactly where the resonance sits (it can shift
%   with temperature, wear, unit-to-unit variation, ...), rather than
%   inflated in amplitude, since the flat plateau this creates already
%   envelopes the peak's own immediate neighborhood without needing any
%   extra headroom. upperConvexHullLogLog will naturally collapse a flat
%   plateau like this to just its two endpoints (interior points at the
%   same height are "on the line" between them, hence redundant) - no
%   special-case merge logic is needed downstream.
%
%   A peak's "significance" is judged by its topographic prominence, in
%   dB: from the peak, extend a line to each side until it either passes
%   a HIGHER point or reaches the end of the data, track the minimum
%   value crossed in each of the two resulting intervals, and take the
%   LARGER of those two interval minima as the reference level (this is
%   the standard prominence definition, e.g. as used by MATLAB's
%   findpeaks 'MinPeakProminence' - a peak embedded in a noisy slope
%   needs this, not just "the nearest local dip on each side", or a
%   single noise wiggle right next to a major peak would make it look
%   unprominent and it would wrongly be skipped).
%
%   Within each bracket [f_left, f_right], psdTarget is SET (not just
%   raised) to the true maximum raw psd found anywhere in that bracket -
%   which may exceed the individual peak sample if two nearby peaks'
%   brackets overlap and get merged into one. This always keeps
%   psdTarget >= psd inside the bracket (the max used is that exact
%   bracket's own worst case), while genuinely zeroing out any margin the
%   general scheme would otherwise have added there - "no amplitude
%   margin" holds for the WHOLE bracket, not just at the peak's own
%   sample.
%
%   peakList is a struct array (one entry per bracket - two nearby peaks
%   whose brackets overlap are merged into one) with fields freq, psd
%   (the representative peak location/level within that bracket),
%   freqLeft, freqRight, for diagnostics/plotting.

f = f(:);
psd = psd(:);
psdTarget = psdTarget(:);
n = numel(f);

peakList = struct('freq', {}, 'psd', {}, 'freqLeft', {}, 'freqRight', {});
if n < 3
    fAug = f;
    psdTargetAug = psdTarget;
    return;
end

psdDB = 10*log10(psd);
peakIdx = [];
for i = 2:n-1
    if ~(psd(i) > psd(i-1) && psd(i) > psd(i+1))
        continue; % not a local maximum
    end

    h = psdDB(i);

    leftMinDB = h;
    j = i - 1;
    while j >= 1 && psdDB(j) <= h
        leftMinDB = min(leftMinDB, psdDB(j));
        j = j - 1;
    end

    rightMinDB = h;
    j = i + 1;
    while j <= n && psdDB(j) <= h
        rightMinDB = min(rightMinDB, psdDB(j));
        j = j + 1;
    end

    referenceDB = max(leftMinDB, rightMinDB);
    prominenceDB = h - referenceDB;
    if prominenceDB >= minProminenceDB
        peakIdx(end+1) = i; %#ok<AGROW>
    end
end

if isempty(peakIdx)
    fAug = f;
    psdTargetAug = psdTarget;
    return;
end

% Raw candidate bands, one per detected peak, sorted by left edge.
fPeakAll = f(peakIdx);
bandLeft = fPeakAll * (1 - freqMarginFraction);
bandRight = fPeakAll * (1 + freqMarginFraction);
[bandLeft, ord] = sort(bandLeft);
bandRight = bandRight(ord);
fPeakAll = fPeakAll(ord);

% Pass 1: merge overlapping/adjacent bands into one, so two close peaks
% share a single flat bracket instead of producing conflicting ones.
% mergedList(k,:) = [bandLeftEdge, bandRightEdge, representativePeakFreq]
mergedList = zeros(0, 3);
curLeft = bandLeft(1);
curRight = bandRight(1);
curPeakF = fPeakAll(1);
for k = 2:numel(bandLeft)
    if bandLeft(k) <= curRight
        curRight = max(curRight, bandRight(k));
        if psd(f == fPeakAll(k)) > psd(f == curPeakF)
            curPeakF = fPeakAll(k);
        end
    else
        mergedList(end+1, :) = [curLeft, curRight, curPeakF]; %#ok<AGROW>
        curLeft = bandLeft(k);
        curRight = bandRight(k);
        curPeakF = fPeakAll(k);
    end
end
mergedList(end+1, :) = [curLeft, curRight, curPeakF];

% Pass 2: for each merged bracket, flatten psdTarget to that bracket's
% own true maximum raw psd (always >= every sample inside it), and add
% the two new bracket-edge anchor points at that level.
newF = zeros(0,1);
newPsd = zeros(0,1);
for k = 1:size(mergedList, 1)
    bLeft = mergedList(k,1);
    bRight = mergedList(k,2);
    bPeakF = mergedList(k,3);

    inBand = f >= bLeft & f <= bRight;
    aPeak = max(psd(inBand));
    psdTarget(inBand) = aPeak;

    newF = [newF; bLeft; bRight]; %#ok<AGROW>
    newPsd = [newPsd; aPeak; aPeak]; %#ok<AGROW>

    peakList(end+1) = struct('freq', bPeakF, 'psd', aPeak, ... %#ok<AGROW>
        'freqLeft', bLeft, 'freqRight', bRight);
end

fAll = [f; newF];
psdTargetAll = [psdTarget; newPsd];

[fSorted, order] = sort(fAll);
psdTargetSorted = psdTargetAll(order);

% De-duplicate repeated frequencies (a bracket edge landing exactly on an
% existing sample, or two peaks' brackets overlapping), keeping the
% larger target at each - coverage requirements only ever get raised
% here, never lowered.
[fAug, ~, groupId] = unique(fSorted);
psdTargetAug = accumarray(groupId, psdTargetSorted, [], @max);
end

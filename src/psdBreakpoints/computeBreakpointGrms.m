function grms = computeBreakpointGrms(bpFreq, bpPsd)
%COMPUTEBREAKPOINTGRMS Composite Grms of a breakpoint (log-log linear) table.
%   grms = COMPUTEBREAKPOINTGRMS(bpFreq, bpPsd) sums the closed-form
%   mean-square of every segment between consecutive breakpoints (see
%   segmentMeanSquare.m) and returns the square root of the total.

bpFreq = bpFreq(:);
bpPsd  = bpPsd(:);

if numel(bpFreq) ~= numel(bpPsd)
    error('computeBreakpointGrms:sizeMismatch', ...
        'Frequency and PSD vectors must be the same length.');
end
if numel(bpFreq) < 2
    error('computeBreakpointGrms:tooFewPoints', ...
        'At least two breakpoints are required.');
end
if any(diff(bpFreq) <= 0)
    error('computeBreakpointGrms:notSorted', ...
        'Breakpoint frequencies must be strictly increasing (sort/dedupe first).');
end

meanSquare = 0;
for k = 1:numel(bpFreq)-1
    meanSquare = meanSquare + segmentMeanSquare(bpFreq(k), bpPsd(k), bpFreq(k+1), bpPsd(k+1));
end
grms = sqrt(meanSquare);
end

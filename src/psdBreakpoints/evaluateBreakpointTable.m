function diagnostics = evaluateBreakpointTable(f, psd, bpTable, targetRmsRatio, marginDB)
%EVALUATEBREAKPOINTTABLE Recompute grms/coverage diagnostics for a table.
%   diagnostics = EVALUATEBREAKPOINTTABLE(f, psd, bpTable, targetRmsRatio,
%   marginDB) checks any breakpoint table (auto-generated or hand-edited)
%   against the original spectrum: the composite Grms ratio and the
%   worst-case coverage margin, in dB, over the original frequency grid.
%
%   bpTable must be a table with variables Frequency_Hz and PSD.

f = f(:);
psd = psd(:);

bpFreq = bpTable.Frequency_Hz(:);
bpPsd  = bpTable.PSD(:);
[bpFreq, order] = sort(bpFreq);
bpPsd = bpPsd(order);

if numel(unique(bpFreq)) ~= numel(bpFreq)
    error('evaluateBreakpointTable:duplicateFrequencies', ...
        'Breakpoint table has duplicate frequencies.');
end

grmsOriginal = computeOriginalGrms(f, psd);
grmsBreakpoint = computeBreakpointGrms(bpFreq, bpPsd);
rmsRatio = grmsBreakpoint / grmsOriginal;

curveAtOrig = interpBreakpointCurve(bpFreq, bpPsd, f);
marginDBActual = 10*log10(curveAtOrig ./ psd);
minMarginDB = min(marginDBActual);

diagnostics = struct( ...
    'grmsOriginal', grmsOriginal, ...
    'grmsBreakpoint', grmsBreakpoint, ...
    'rmsRatio', rmsRatio, ...
    'targetRmsRatio', targetRmsRatio, ...
    'meetsRmsTarget', rmsRatio <= targetRmsRatio, ...
    'marginDB', marginDB, ...
    'minMarginDB', minMarginDB, ...
    'numPoints', numel(bpFreq), ...
    'maxPoints', NaN);
end

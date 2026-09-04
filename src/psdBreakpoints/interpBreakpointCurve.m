function psdOut = interpBreakpointCurve(bpFreq, bpPsd, fQuery)
%INTERPBREAKPOINTCURVE Evaluate a log-log piecewise-linear breakpoint curve.
%   psdOut = INTERPBREAKPOINTCURVE(bpFreq, bpPsd, fQuery) interpolates the
%   breakpoint table the same way it is meant to be read/plotted: straight
%   lines on log10(freq) vs log10(psd) axes. Query points outside the
%   table range are clamped to the first/last breakpoint.

bpFreq = bpFreq(:);
bpPsd  = bpPsd(:);

[bpFreq, order] = sort(bpFreq);
bpPsd = bpPsd(order);

logF = log10(bpFreq);
logW = log10(bpPsd);

fQuery = fQuery(:);
fQueryClamped = min(max(fQuery, bpFreq(1)), bpFreq(end));

logOut = interp1(logF, logW, log10(fQueryClamped), 'linear');
psdOut = reshape(10.^logOut, size(fQuery));
end

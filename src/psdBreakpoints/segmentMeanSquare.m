function ms = segmentMeanSquare(f1, w1, f2, w2)
%SEGMENTMEANSQUARE Mean-square (integral) of one breakpoint segment.
%   ms = SEGMENTMEANSQUARE(f1, w1, f2, w2) returns the integral of the
%   PSD between two consecutive breakpoints (f1,w1) and (f2,w2), assuming
%   the segment is a straight line on log-log axes (a power law
%   W(f) = w1*(f/f1)^n). This is the standard closed-form expression used
%   to compute a composite Grms from a breakpoint table, and matches how
%   MIL-STD-810 style random-vibration test/tolerance curves -
%   piecewise-linear on log(PSD) vs log(frequency) axes - are integrated.

if f2 <= f1
    error('segmentMeanSquare:badFrequencies', 'f2 must be greater than f1.');
end
if w1 <= 0 || w2 <= 0
    error('segmentMeanSquare:nonPositivePSD', 'PSD values must be > 0.');
end

n = log(w2 / w1) / log(f2 / f1);

if abs(n + 1) < 1e-9
    % Special case n = -1: constant PSD*f, integral is logarithmic.
    ms = w1 * f1 * log(f2 / f1);
else
    ms = (w2 * f2 - w1 * f1) / (n + 1);
end
end

function grms = computeOriginalGrms(f, psd)
%COMPUTEORIGINALGRMS Overall RMS (Grms) of a linearly spaced PSD.
%   grms = COMPUTEORIGINALGRMS(f, psd) integrates the PSD over frequency
%   with trapezoidal integration on the supplied grid and returns the
%   square root of the area, i.e. the composite RMS value.
%
%   f, psd : equal length vectors, f strictly increasing.

f   = f(:);
psd = psd(:);

if numel(f) ~= numel(psd)
    error('computeOriginalGrms:sizeMismatch', ...
        'f and psd must have the same number of elements.');
end
if any(psd < 0)
    error('computeOriginalGrms:negativePSD', 'PSD values must be non-negative.');
end

meanSquare = trapz(f, psd);
grms = sqrt(meanSquare);
end

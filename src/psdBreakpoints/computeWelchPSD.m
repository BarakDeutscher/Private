function [f, psd, info] = computeWelchPSD(t, x, freqResolutionHz, overlapFraction)
%COMPUTEWELCHPSD One-sided Welch PSD of a uniformly sampled time history.
%   [f, psd, info] = COMPUTEWELCHPSD(t, x, freqResolutionHz, overlapFraction)
%   estimates the sample rate from the time vector (robust to
%   text/rounding jitter in exported timestamps: it uses the total time
%   span divided by the sample count, not local sample-to-sample diffs,
%   since limited-precision CSV timestamps can make truly uniform
%   sampling look jittery sample-to-sample), then computes a
%   Hann-windowed Welch PSD estimate at the requested frequency
%   resolution and segment overlap. Returns a ONE-SIDED PSD on a linearly
%   spaced frequency grid from 0 to fs/2 - exactly the input format
%   expected by generatePSDBreakpointTable.
%
%   t                   Time vector, seconds (only used to estimate fs).
%   x                   Time history samples (e.g. acceleration, g).
%   freqResolutionHz    Desired PSD frequency bin spacing, in Hz (e.g. 1).
%   overlapFraction     Segment overlap fraction, 0 <= overlapFraction < 1
%                        (e.g. 0.5 for 50%).
%
%   info fields: fs, nfft, noverlap, numSegments, grmsTimeDomain (std of
%   x, i.e. AC/vibration RMS with any DC bias removed), grmsFromPSD
%   (Grms recovered by integrating the PSD, excluding the DC bin - should
%   closely match grmsTimeDomain; a Parseval sanity check).

t = t(:);
x = x(:);
if numel(t) ~= numel(x)
    error('computeWelchPSD:sizeMismatch', 't and x must be the same length.');
end
if numel(x) < 2
    error('computeWelchPSD:tooFewSamples', 'Need at least 2 samples.');
end
if overlapFraction < 0 || overlapFraction >= 1
    error('computeWelchPSD:badOverlap', 'overlapFraction must be in [0, 1).');
end

% Sample rate from total span / sample count: robust against per-sample
% timestamp rounding in exported text (a handful of significant figures
% can make truly uniform sampling look jittery from one row to the next).
fsRaw = (numel(t) - 1) / (t(end) - t(1));
fs = round(fsRaw);
if abs(fsRaw - fs) > 1e-3 * fs
    warning('computeWelchPSD:nonIntegerSampleRate', ...
        ['Estimated sample rate %.4f Hz is not close to an integer; using ' ...
         'the estimated value as-is. If you know the true DAQ sample ' ...
         'rate, resample the data or override fs explicitly instead.'], fsRaw);
    fs = fsRaw;
end

nfft = round(fs / freqResolutionHz);
if abs(fs/nfft - freqResolutionHz) > 1e-9
    warning('computeWelchPSD:resolutionNotExact', ...
        ['Requested resolution %.4f Hz does not divide fs = %.4f Hz evenly; ' ...
         'actual resolution will be %.6f Hz.'], freqResolutionHz, fs, fs/nfft);
end
if nfft > numel(x)
    error('computeWelchPSD:segmentTooLong', ...
        ['Requested frequency resolution %.4f Hz needs %d samples per segment, ' ...
         'but only %d samples are available. Use a coarser resolution.'], ...
        freqResolutionHz, nfft, numel(x));
end

noverlap = round(nfft * overlapFraction);
window = hann(nfft);

[psd, f] = pwelch(x, window, noverlap, nfft, fs, 'onesided');

f = f(:);
psd = psd(:);

numSegments = floor((numel(x) - noverlap) / (nfft - noverlap));

grmsTimeDomain = std(x); % AC-only RMS (DC/mean removed)
validIdx = f > 0;        % PSD undefined at f=0 on the log axis used downstream
grmsFromPSD = computeOriginalGrms(f(validIdx), psd(validIdx));

info = struct('fs', fs, 'nfft', nfft, 'noverlap', noverlap, ...
    'numSegments', numSegments, 'grmsTimeDomain', grmsTimeDomain, ...
    'grmsFromPSD', grmsFromPSD);
end

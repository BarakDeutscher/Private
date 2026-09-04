%% Sanity test for computeWelchPSD (requires the Signal Processing Toolbox)
% Run with: run('tests/test_welch_psd.m')
%
% NOTE: this test relies on MATLAB's pwelch signature, where the 3rd
% argument (noverlap) is a number of SAMPLES. Octave's signal-package
% pwelch instead takes segment overlap as a FRACTION (0-1) - a genuine
% API difference between the two, not a bug - so this test is not
% Octave-compatible. The underlying Welch/Parseval math was validated
% separately against Octave's pwelch (using its fraction-based call)
% during development; this test checks the actual shipped MATLAB code.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'psdBreakpoints'));

fs = 2000;
dt = 1/fs;
T = 60;
N = round(T*fs);
t = (0:N-1)'*dt;

rng(42);
x = 0.05*randn(N,1) + 0.3*sin(2*pi*150*t);

[f, psd, info] = computeWelchPSD(t, x, 1, 0.5);

assert(info.fs == fs, 'sample rate estimate mismatch');
assert(abs(f(2)-f(1) - 1) < 1e-9, 'frequency resolution mismatch');
assert(abs(info.grmsFromPSD - info.grmsTimeDomain) < 0.01*info.grmsTimeDomain, ...
    'Parseval check failed: PSD-integrated grms does not match time-domain std');

% The 150 Hz tone should show up as the dominant PSD peak
[~, peakIdx] = max(psd);
assert(abs(f(peakIdx) - 150) <= 1, 'expected PSD peak near the 150 Hz tone');

% A non-uniform-looking (rounded/jittery) time vector with the SAME true
% span and sample count must still recover the same fs (robustness check
% for real exported CSV timestamps, which are often rounded to a handful
% of significant figures).
tJittered = t + 1e-6*sin((1:N)'); % tiny per-sample jitter, same start/end
[~, ~, infoJ] = computeWelchPSD(tJittered, x, 1, 0.5);
assert(infoJ.fs == fs, 'fs estimate not robust to small timestamp jitter');

fprintf('All computeWelchPSD sanity tests passed.\n');

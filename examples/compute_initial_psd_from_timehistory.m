%% Compute the initial PSD from a raw acceleration time history
% Reads a two-column CSV (time [s], acceleration [g], no header),
% estimates the sample rate directly from the data (robust to per-sample
% timestamp rounding in the exported text), and computes a one-sided
% Welch PSD at a user-chosen frequency resolution and segment overlap.
% The result (f, psd) is a linearly frequency-spaced PSD - ready to feed
% straight into generatePSDBreakpointTable as the "initial PSD".

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'psdBreakpoints'));

% ---- point this at your time-history CSV -------------------------------
csvPath = 'path/to/your/time_history.csv';     % <-- EDIT ME
freqResolutionHz = 1;                          % 1 Hz PSD resolution
overlapFraction  = 0.5;                        % 50% segment overlap
% -------------------------------------------------------------------------

data = readmatrix(csvPath);
t = data(:,1);
x = data(:,2);

[f, psd, info] = computeWelchPSD(t, x, freqResolutionHz, overlapFraction);

fprintf('--- Initial PSD Summary ---\n');
fprintf('  Estimated sample rate  : %.4f Hz\n', info.fs);
fprintf('  Segment length (nfft)  : %d samples (%.3f s)\n', info.nfft, info.nfft/info.fs);
fprintf('  Segment overlap        : %d samples (%.0f%%)\n', info.noverlap, 100*info.noverlap/info.nfft);
fprintf('  Welch segments averaged: %d\n', info.numSegments);
fprintf('  Grms, time domain (std): %.5f\n', info.grmsTimeDomain);
fprintf('  Grms, from PSD (f>0)   : %.5f  (Parseval check, should match closely)\n', info.grmsFromPSD);

% Exclude the DC bin (f=0): undefined on the log-log axis breakpoint
% curves are built on, and not meaningful vibration content anyway.
keep = f > 0;
f = f(keep);
psd = psd(keep);

figure('Name', 'Initial PSD from time history');
loglog(f, psd, '-', 'Color', [0.15 0.35 0.65], 'LineWidth', 1);
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD (g^2/Hz)');
title(sprintf('Welch PSD, %.0f%% overlap, %.1f Hz resolution (Grms = %.4f)', ...
    100*overlapFraction, freqResolutionHz, info.grmsFromPSD));

% Save (f, psd) so it can be reused directly as input to
% generatePSDBreakpointTable without recomputing the PSD each time.
writematrix([f, psd], fullfile(fileparts(mfilename('fullpath')), 'initial_psd.csv'));

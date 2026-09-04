%% Example: reduce a synthetic linearly-spaced PSD to a breakpoint table
% Builds a representative random-vibration PSD (ramp - plateau - ramp,
% the classic shape used in MIL-STD-810 category vibration test curves,
% plus a couple of resonance peaks), then reduces it to a small
% breakpoint table that envelopes it.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'psdBreakpoints'));

f = (10:1:2000)';                 % linearly spaced input, 1 Hz resolution
psd = zeros(size(f));

% 10-80 Hz: rising ramp from 0.01 to 0.08 g^2/Hz
seg1 = f >= 10 & f <= 80;
psd(seg1) = interp1([10 80], [0.01 0.08], f(seg1));

% 80-350 Hz: flat plateau at 0.08 g^2/Hz
seg2 = f > 80 & f <= 350;
psd(seg2) = 0.08;

% 350-2000 Hz: decaying ramp down to 0.01 g^2/Hz
seg3 = f > 350 & f <= 2000;
psd(seg3) = interp1([350 2000], [0.08 0.01], f(seg3));

% Narrowband resonance bumps, to show the algorithm captures peaks
psd = psd + 0.03  * exp(-((f-150).^2) / (2*8^2));
psd = psd + 0.015 * exp(-((f-900).^2) / (2*20^2));

% NOTE on MarginDB vs TargetRmsRatio: a uniform dB margin alone sets a
% floor on the achievable grms ratio of sqrt(10^(MarginDB/10)) (raising
% every PSD value by a constant factor raises grms by that factor's
% square root, before any point-reduction overshoot is even considered).
% For TargetRmsRatio = 1.4 that means MarginDB must stay below
% 20*log10(1.4) = 2.92 dB, so this example uses 2 dB.
[bpTable, diagnostics] = generatePSDBreakpointTable(f, psd, ...
    'MaxPoints', 10, ...
    'MarginDB', 2, ...
    'TargetRmsRatio', 1.4, ...
    'Interactive', false, ...
    'Plot', true);

disp(bpTable);
disp(diagnostics);

%% Manual edit pass (uncomment to try interactively)
% bpTable = editBreakpointTableGUI(f, psd, bpTable, 1.4, 2);
% diagnostics = evaluateBreakpointTable(f, psd, bpTable, 1.4, 2);
% disp(bpTable);
% disp(diagnostics);

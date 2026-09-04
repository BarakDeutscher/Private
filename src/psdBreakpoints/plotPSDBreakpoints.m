function plotPSDBreakpoints(f, psd, bpTable, diagnostics)
%PLOTPSDBREAKPOINTS Log-log plot of original spectrum vs breakpoint curve.

figure('Name', 'PSD Breakpoint Reduction');
loglog(f, psd, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'DisplayName', 'Original PSD');
hold on;
loglog(bpTable.Frequency_Hz, bpTable.PSD, '-o', 'Color', [0.85 0.1 0.1], ...
    'LineWidth', 1.5, 'MarkerFaceColor', [0.85 0.1 0.1], 'DisplayName', 'Breakpoint curve');
grid on;
xlabel('Frequency (Hz)');
ylabel('PSD');
title(sprintf('Grms ratio = %.3f (target \\leq %.2f), min margin = %.2f dB', ...
    diagnostics.rmsRatio, diagnostics.targetRmsRatio, diagnostics.minMarginDB));
legend('Location', 'best');
hold off;
end

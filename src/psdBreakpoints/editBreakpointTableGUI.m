function bpTableOut = editBreakpointTableGUI(f, psd, bpTableIn, targetRmsRatio, marginDB)
%EDITBREAKPOINTTABLEGUI Interactive editor for a PSD breakpoint table.
%   bpTableOut = EDITBREAKPOINTTABLEGUI(f, psd, bpTableIn, targetRmsRatio,
%   marginDB) opens a GUI with an editable breakpoint table (Frequency_Hz,
%   PSD) and a live log-log plot against the original spectrum. Edit
%   cells directly, or use Add Row / Delete Selected Row. The Grms ratio
%   and worst-case coverage margin update after every change and are
%   shown in green (target met) or red (target missed). Click Done to
%   close the GUI and return the (possibly hand-edited) table.

f = f(:);
psd = psd(:);

fig = uifigure('Name', 'Edit PSD Breakpoint Table', 'Position', [100 100 940 540]);
fig.CloseRequestFcn = @(~, ~) uiresume(fig);

ax = uiaxes(fig, 'Position', [30 70 540 440]);

tbl = uitable(fig, 'Data', bpTableIn, 'ColumnEditable', [true true], ...
    'Position', [590 160 320 350]);

statusLabel = uilabel(fig, 'Position', [590 100 320 55], ...
    'Text', '', 'FontSize', 11, 'WordWrap', 'on');

uibutton(fig, 'Text', 'Add Row', 'Position', [590 60 95 30], ...
    'ButtonPushedFcn', @(~, ~) addRow());
uibutton(fig, 'Text', 'Delete Selected', 'Position', [695 60 115 30], ...
    'ButtonPushedFcn', @(~, ~) deleteSelectedRow());
uibutton(fig, 'Text', 'Done', 'Position', [820 60 90 30], ...
    'ButtonPushedFcn', @(~, ~) uiresume(fig));

selectedRow = [];
tbl.CellEditCallback = @(~, ~) refresh();
tbl.CellSelectionCallback = @(~, evt) trackSelection(evt);

refresh();
uiwait(fig);

bpTableOut = tbl.Data;
delete(fig);

    function trackSelection(evt)
        if ~isempty(evt.Indices)
            selectedRow = evt.Indices(1, 1);
        end
    end

    function addRow()
        d = tbl.Data;
        if isempty(d)
            newF = f(1);
            newP = psd(1);
        else
            newF = d.Frequency_Hz(end) * 0.9; % placeholder, edit it in the table
            newP = d.PSD(end);
        end
        d = [d; {newF, newP}];
        tbl.Data = d;
        refresh();
    end

    function deleteSelectedRow()
        d = tbl.Data;
        if isempty(selectedRow) || height(d) <= 2
            uialert(fig, ...
                'Cannot delete: keep at least 2 breakpoints, and select a row first.', ...
                'Delete Row');
            return;
        end
        d(selectedRow, :) = [];
        tbl.Data = d;
        selectedRow = [];
        refresh();
    end

    function refresh()
        d = tbl.Data;
        try
            d = sortrows(d, 'Frequency_Hz');
            tbl.Data = d;
            diag = evaluateBreakpointTable(f, psd, d, targetRmsRatio, marginDB);

            cla(ax);
            loglog(ax, f, psd, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            hold(ax, 'on');
            loglog(ax, d.Frequency_Hz, d.PSD, '-o', 'Color', [0.85 0.1 0.1], ...
                'LineWidth', 1.5, 'MarkerFaceColor', [0.85 0.1 0.1]);
            hold(ax, 'off');
            grid(ax, 'on');
            xlabel(ax, 'Frequency (Hz)');
            ylabel(ax, 'PSD');
            legend(ax, {'Original', 'Breakpoint (edited)'}, 'Location', 'best');

            statusLabel.Text = sprintf(['Grms orig:  %.4f\nGrms new:   %.4f\n' ...
                'Ratio:      %.3f  (target <= %.2f)\n' ...
                'Min margin: %.2f dB (requested %.2f dB)\nPoints:     %d'], ...
                diag.grmsOriginal, diag.grmsBreakpoint, diag.rmsRatio, ...
                diag.targetRmsRatio, diag.minMarginDB, diag.marginDB, diag.numPoints);

            if diag.meetsRmsTarget && diag.minMarginDB >= -1e-6
                statusLabel.FontColor = [0.13 0.55 0.13]; % green: both targets met
            else
                statusLabel.FontColor = [0.70 0.13 0.13]; % red: ratio and/or margin missed
            end
        catch ME
            statusLabel.Text = sprintf('Invalid table: %s', ME.message);
            statusLabel.FontColor = [0.70 0 0];
        end
    end
end

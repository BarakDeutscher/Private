# PSD Breakpoint Table Generator

MATLAB tool that takes a linearly frequency-spaced PSD (Power Spectral
Density) vibration spectrum and reduces it to a compact **breakpoint
table** — the (frequency, PSD) pair list used to define random-vibration
test/tolerance curves, in the style of **MIL-STD-810**.

MIL-STD-810 vibration test curves are always specified as a handful of
breakpoints connected by **straight-line segments on log(PSD) vs.
log(frequency) axes** (i.e. power-law segments, often described in
dB/octave). This tool builds a curve in exactly that form: it finds the
fewest breakpoints needed so the piecewise-linear curve **envelopes**
(stays at or above) the original spectrum, subject to limits you control.

## What it does

0. (Optional) If you only have a raw acceleration **time history**, not
   a PSD yet, `computeWelchPSD` turns it into one first: a one-sided
   Welch PSD estimate (Hann-windowed, your choice of frequency
   resolution and segment overlap) on a linearly spaced frequency grid -
   exactly the input format the rest of the pipeline expects.
1. Computes the composite RMS (**Grms**) of your original PSD.
2. Applies your requested coverage margin (in dB) above the original.
3. Finds the minimal breakpoint set that envelopes the margined spectrum,
   using the **upper convex hull** in log-log space — this is the
   provably fewest-points curve that stays above every sample (it always
   keeps peaks/resonances, and drops points that are already covered by
   the line between their neighbors).
4. If that needs more points than your `MaxPoints` budget, it trims the
   least significant breakpoints first (Visvalingam–Whyatt area-based
   simplification), reporting any resulting shortfall in coverage. If
   instead there's spare budget left over, it spends it greedily
   (`refineBreakpointsGreedy.m`) narrowing the gap to `TargetRmsRatio` -
   fewest-points-for-coverage and lowest-area-for-a-point-budget are
   different problems, and a sharp, distant peak can force the minimal
   hull to badly overshoot a quiet floor even though every floor sample
   is individually "dominated" by its neighbors' straight line.
5. Computes the Grms of the new breakpoint curve using the standard
   closed-form log-log segment integral (not linear trapezoidal
   integration, which would misstate the area under log-log segments),
   and reports the ratio vs. the original Grms against your target.
6. Optionally opens a GUI so you can hand-edit the table (add/delete
   rows, retype any cell) and see the Grms ratio and coverage margin
   update live.

## Files

```
src/psdBreakpoints/
  computeWelchPSD.m          Time history -> one-sided Welch PSD estimate
  computeOriginalGrms.m      Grms of the linearly spaced input PSD (trapz)
  segmentMeanSquare.m        Closed-form integral of one log-log segment
  computeBreakpointGrms.m    Grms of a full breakpoint table
  interpBreakpointCurve.m    Evaluate the breakpoint curve at any frequency
  upperConvexHullLogLog.m    Minimal enveloping breakpoint set
  reduceBreakpointsVW.m      Point-budget trim (Visvalingam-Whyatt)
  refineBreakpointsGreedy.m  Spend spare points to shrink the Grms ratio
  pruneNegligibleBreakpoints.m  Drop any point that isn't pulling its weight
  applyPeakFrequencyBracketing.m  Flat +/-freq%% peak brackets, no amplitude margin
  evaluateBreakpointTable.m  Grms ratio + coverage margin for any table
  plotPSDBreakpoints.m       Log-log plot helper
  generatePSDBreakpointTable.m   Main entry point
  editBreakpointTableGUI.m   Interactive manual-edit GUI

examples/
  compute_initial_psd_from_timehistory.m   Time history CSV -> initial PSD
  example_generate_breakpoints.m           End-to-end breakpoint example

tests/
  test_psd_breakpoints.m     Self-contained sanity checks (no toolboxes)
  test_welch_psd.m           computeWelchPSD checks (needs Signal Proc. Toolbox)
```

## Usage

### Starting from a time history

```matlab
addpath('src/psdBreakpoints');

data = readmatrix('time_history.csv');   % 2 columns: time (s), acceleration (g)
t = data(:,1);
x = data(:,2);

[f, psd, info] = computeWelchPSD(t, x, 1, 0.5); % 1 Hz resolution, 50% overlap
f = f(2:end); psd = psd(2:end);                 % drop the f=0 (DC) bin
```

`computeWelchPSD` estimates the sample rate from the *total time span
divided by sample count*, not local sample-to-sample differences - real
exported CSV timestamps are often rounded to only a handful of
significant figures, which can make truly uniform sampling look jittery
row-to-row; averaging over the whole file cancels that out. `info` also
reports a Parseval cross-check: `info.grmsFromPSD` (Grms recovered by
integrating the PSD) should closely match `info.grmsTimeDomain` (the
time-domain `std` of the signal) - a sanity check that the PSD was
computed correctly, in the same spirit as the Grms-ratio check below.

### Starting from a PSD

```matlab
% f, psd: your linearly frequency-spaced input spectrum (Hz, units^2/Hz)
[bpTable, diagnostics] = generatePSDBreakpointTable(f, psd, ...
    'MaxPoints', 10, ...       % max breakpoints allowed
    'MarginDB', 2, ...         % coverage margin above the original, in dB
    'TargetRmsRatio', 1.4, ... % target ceiling on grms_new / grms_original
    'Interactive', false, ...  % true = open the manual-edit GUI after
    'Plot', true, ...
    'OutputFolder', 'output'); % saves table + comparison PNG here

disp(bpTable);        % table with Frequency_Hz, PSD
disp(diagnostics);     % grmsOriginal, grmsBreakpoint, rmsRatio, minMarginDB, ...
```

When `OutputFolder` is set, two files are written there (folder created if
needed): `breakpoint_table.csv` (the final table) and
`psd_breakpoint_comparison.png` (original PSD vs. the final breakpoint
curve, log-log). If `Interactive` is also `true`, both reflect the table
*after* your manual edits, not the auto-generated one.

To edit a table by hand at any later point (not just right after
generation):

```matlab
bpTable = editBreakpointTableGUI(f, psd, bpTable, 1.4, 2);
diagnostics = evaluateBreakpointTable(f, psd, bpTable, 1.4, 2);
```

Run `examples/example_generate_breakpoints.m` for a full worked example
with a synthetic ramp/plateau/ramp spectrum plus two resonance peaks —
the shape typical of a MIL-STD-810 category vibration test curve — and
run `tests/test_psd_breakpoints.m` to sanity-check the math on your
MATLAB install.

## Controls

- **`MaxPoints`** — hard cap on the number of breakpoints in the output
  table. If the exact envelope needs more, the least significant
  breakpoints are dropped first and a warning reports the worst-case
  coverage shortfall (`diagnostics.minMarginDB` — negative means the
  curve dips below the original by that many dB somewhere). If the exact
  envelope needs fewer, the spare points are spent narrowing the Grms
  ratio (see below) rather than left unused.
- **`MarginDB`** — how far above the original PSD the curve should sit
  before the envelope is fit. `0` = tight envelope (curve touches the
  original at its peaks), positive values add headroom (e.g. for test
  tolerance bands). Its meaning depends on `MarginMode`.
- **`MarginMode`** — `'uniform'` (default) applies `MarginDB` everywhere,
  flat. `'adaptive'` applies `MarginDB` only at the spectrum's quietest
  point and tapers it linearly (in dB, vs. normalized log-amplitude) down
  to **0 dB exactly at the loudest peak** - so major peaks are not
  inflated at all, and low-amplitude regions get the most headroom. Since
  the peak dominates the Grms integral, leaving it untouched typically
  cuts the ratio far more than shaving a uniform margin ever could, often
  needing fewer points too (on one real flight-vibration spectrum: 21
  points at ratio 1.40 with `'uniform'`, vs. 14 points at ratio 1.25 with
  `'adaptive'`, same `MarginDB`).
- **`PeakBracketing`** (default `false`) — represent each significant
  narrowband peak with exactly **two flat, equal-amplitude breakpoints**
  straddling it, at `±PeakFreqMarginFraction` (default 10%) of its
  frequency, held at the peak's own raw PSD value (no amplitude margin
  at all there). This is the standard way a resonance is broadened in a
  derived test spec: widen it in *frequency* to cover uncertainty in
  exactly where it sits (temperature, wear, unit-to-unit variation), not
  by inflating its amplitude. `upperConvexHullLogLog` naturally collapses
  the flat plateau this creates down to just its two edges - the peak's
  own sample becomes redundant, so tracing its exact shape (which can
  otherwise take several points) is no longer needed. Trade-off: a flat
  cap across a wide band overshoots more in the peak's flanks than
  tightly tracing its true shape would, so the achieved Grms ratio is
  usually a bit higher than without bracketing, for the same point
  count - narrow `PeakFreqMarginFraction` or raise `MaxPoints` if that
  pushes you over `TargetRmsRatio`. `PeakMinProminenceDB` (default 6 dB)
  sets how much a local maximum must stand out - via true topographic
  prominence, not just "the nearest dip on each side" - to count as a
  peak worth bracketing.
- **`TargetRmsRatio`** — the ceiling you want on `grms_new / grms_original`.
  This is a check, not a hard constraint: `diagnostics.meetsRmsTarget`
  reports whether it was actually achieved, with a warning if not. Every
  breakpoint is checked afterward (`pruneNegligibleBreakpoints.m`) and
  dropped if removing it neither breaks coverage nor pushes the ratio
  back over this target - so a tight target that needed many points
  inserted greedily doesn't leave some of them redundant once the full
  set is assembled (a later insertion nearby can end up doing an
  earlier one's job, or the last few points added to just barely clear
  target can leave neighbors over-provisioned).

### Why the ratio can overshoot even with points to spare

Minimizing the *number* of breakpoints needed for full coverage
(`upperConvexHullLogLog`) and minimizing the *enclosed area* (Grms) for a
given point budget are different problems. A sharp resonance peak sitting
far above a quiet noise floor forces the hull's bridging segment to badly
overshoot that floor in between - even though every floor sample is
individually "dominated" (below the straight line joining the peak to its
distant neighbor, so the hull correctly leaves it out of the *minimal*
set). Voluntarily adding some of those dominated points back in lets the
curve dip down and track the real floor, cutting Grms substantially - as
long as the two new segments each insertion creates still stay above
everything in their own sub-range, which `refineBreakpointsGreedy.m`
checks before accepting any insertion. This is exactly the case with real
flight vibration data carrying a strong tonal peak: the minimal hull can
need very few points (e.g. 9) but overshoot the target ratio by a lot
(e.g. 2.5x), while spending a realistic point budget (e.g. 21) on
targeted insertions gets it back under target.

### Important interaction: margin sets a floor on the ratio (uniform mode only)

Scaling every PSD value by a constant factor scales Grms by that
factor's square root. So with `MarginMode='uniform'`, `MarginDB` alone
imposes a **hard lower bound** on the achievable ratio:

```
min achievable ratio = sqrt(10^(MarginDB/10))
```

For `TargetRmsRatio = 1.4`, that means `MarginDB` must stay below
`20*log10(1.4) ≈ 2.92 dB`, or the target is mathematically unreachable
regardless of `MaxPoints`. The function checks this up front and warns
you if your `MarginDB`/`TargetRmsRatio` combination conflicts. This floor
does not apply to `MarginMode='adaptive'` - there the effective scaling
varies by frequency, so there's no simple closed form; the actual ratio
just gets computed and checked against `TargetRmsRatio` as usual.

## Manual editing

`editBreakpointTableGUI` opens a MATLAB `uifigure` with:
- An editable table of breakpoints (retype any Frequency_Hz/PSD cell).
- **Add Row** / **Delete Selected** buttons.
- A live log-log plot of the original spectrum vs. your edited curve.
- A status panel showing Grms ratio and worst-case coverage margin,
  colored green when both targets are met and red when they are not.

Click **Done** to close the GUI and get the (possibly hand-edited) table
back as a normal MATLAB `table`.

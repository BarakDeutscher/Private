function hullIdx = upperConvexHullLogLog(f, psd)
%UPPERCONVEXHULLLOGLOG Minimal breakpoint set that envelopes psd(f).
%   hullIdx = UPPERCONVEXHULLLOGLOG(f, psd) returns indices (into f/psd,
%   which must already be sorted by increasing frequency) of the points
%   that form the upper convex hull of the data in log10(f)-log10(psd)
%   space. Connecting these points with straight lines on a log-log plot
%   - exactly how MIL-STD-810 style breakpoint / test-tolerance curves
%   are drawn - produces the piecewise-linear curve with the FEWEST
%   possible breakpoints that still lies at or above every input point.
%   Any point not selected lies on or below the segment joining its
%   neighbors, so it is automatically enveloped by the hull curve.

f = f(:);
psd = psd(:);
if any(psd <= 0)
    error('upperConvexHullLogLog:nonPositivePSD', ...
        'PSD values must be > 0 to take log10 (check for zeros).');
end
if any(diff(f) <= 0)
    error('upperConvexHullLogLog:notSorted', ...
        'f must be strictly increasing.');
end

x = log10(f);
y = log10(psd);
n = numel(x);

stack = zeros(n,1);
k = 0;
for i = 1:n
    while k >= 2 && crossProd(x(stack(k-1)), y(stack(k-1)), ...
                               x(stack(k)),   y(stack(k)), ...
                               x(i),          y(i)) >= 0
        k = k - 1;
    end
    k = k + 1;
    stack(k) = i;
end
hullIdx = stack(1:k);
end

function c = crossProd(ax, ay, bx, by, cx, cy)
% Cross product of (B-A) x (C-A). c < 0 => right turn at B (keep, e.g. a
% peak); c >= 0 => B sits on/under the line A-C and is redundant.
c = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
end

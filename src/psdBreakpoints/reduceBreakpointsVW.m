function keepIdx = reduceBreakpointsVW(x, y, idx, maxPoints)
%REDUCEBREAKPOINTSVW Reduce a polyline to maxPoints (Visvalingam-Whyatt).
%   keepIdx = REDUCEBREAKPOINTSVW(x, y, idx, maxPoints) repeatedly drops
%   the interior vertex (from the point set idx, in order) that
%   contributes the smallest triangle area to the polyline shape, until
%   only maxPoints remain. The first and last points are never removed.
%   Each removal only touches the areas of its two former neighbors (a
%   linked list tracks the current chain), so the whole reduction is
%   cheap even when idx starts with hundreds of points.
%
%   Used when the exact enveloping hull (see upperConvexHullLogLog) needs
%   more breakpoints than the user allows: dropping the least significant
%   vertices first keeps the curve as close as possible to the true
%   envelope. Any resulting coverage/margin shortfall is reported
%   separately by the caller.

idx = idx(:)';
if maxPoints < 2
    error('reduceBreakpointsVW:tooFewPoints', 'maxPoints must be >= 2.');
end

m = numel(idx);
if m <= maxPoints
    keepIdx = idx;
    return;
end

prevPos = [0, 1:m-1];  % linked-list predecessor position (0 = none)
nextPos = [2:m, 0];    % linked-list successor position (0 = none)
alive = true(1, m);

areas = inf(1, m);     % endpoints (1 and m) stay inf: never dropped
for i = 2:m-1
    areas(i) = triArea(x, y, idx(i-1), idx(i), idx(i+1));
end

numAlive = m;
while numAlive > maxPoints
    [~, drop] = min(areas);
    alive(drop) = false;
    areas(drop) = inf;
    numAlive = numAlive - 1;

    pv = prevPos(drop);
    nx = nextPos(drop);
    if pv ~= 0
        nextPos(pv) = nx;
    end
    if nx ~= 0
        prevPos(nx) = pv;
    end

    if pv ~= 0 && prevPos(pv) ~= 0
        areas(pv) = triArea(x, y, idx(prevPos(pv)), idx(pv), idx(nextPos(pv)));
    end
    if nx ~= 0 && nextPos(nx) ~= 0
        areas(nx) = triArea(x, y, idx(prevPos(nx)), idx(nx), idx(nextPos(nx)));
    end
end

keepIdx = idx(alive);
end

function a = triArea(x, y, ia, ib, ic)
a = 0.5 * abs((x(ib)-x(ia))*(y(ic)-y(ia)) - (x(ic)-x(ia))*(y(ib)-y(ia)));
end

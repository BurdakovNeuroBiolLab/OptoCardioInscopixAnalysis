% pupilCircleFit.m
% Finds a DeepLabCut CSV with "zoompup" in its name in the current folder,
% removes low-confidence points (set to NaN), fits a circle to the remaining
% pupil points in every frame, removes blinks / bad fits / outliers, and
% plots radius and diameter over time (before vs after cleaning).
%
% File layout: 3 header lines, then per row:
%   col 1 = frame index / time, then repeating [x, y, likelihood] per point.
%
% Cleaning steps:
%   1. Fit only frames with >= minPoints good points
%   2. Remove blinks (frames with many low-confidence points) +/- blinkPad frames
%   3. Remove poor fits (point scatter around the circle > maxRelResid * radius)
%   4. (optional) Remove remaining single-frame spikes with a moving-median test
%   5. (optional) Fill short NaN gaps by linear interpolation

% When called from a parent script that sets keepWorkspace = true, the
% workspace is not cleared so the parent's variables survive.
if ~exist('keepWorkspace', 'var') || ~keepWorkspace
    clear; clc;
end

%% Settings
confThresh   = 0.99;  % points with likelihood < confThresh are set to NaN
minPoints    = 6;     % step 1: minimum good points needed to fit a circle
blinkMinBad  = 3;     % step 2: frame is a blink if >= this many points are low-confidence
blinkPad     = 2;     % step 2: also remove this many frames before/after each blink
maxRelResid  = 0.1;   % step 3: reject fits where RMS point scatter > this fraction of radius
doOutliers   = true;  % step 4: remove spikes that deviate from the moving median
outlierWin   = 31;    % step 4: moving-median window (frames)
outlierK     = 5;     % step 4: flag if > outlierK scaled MADs from the moving median
doFillGaps   = false; % step 5: interpolate across short NaN gaps
maxGapFill   = 10;    % step 5: longest gap (frames) to fill
fps          = [];    % camera frame rate (Hz). Leave empty to plot vs frame number.

%% Find and load the file
files = dir(fullfile(pwd, '*zoompup*.csv'));
% Ignore macOS "._" ghost files and this script's own output files
names = {files.name};
keep  = ~startsWith(names, '._') & ~endsWith(names, '_circlefit.csv');
files = files(keep);
if isempty(files)
    error('No CSV file containing "zoompup" found in %s', pwd);
end
if numel(files) > 1
    warning('%d "zoompup" files found; using the first: %s', numel(files), files(1).name);
end
fname = fullfile(files(1).folder, files(1).name);
fprintf('Loading %s\n', files(1).name);

data = readmatrix(fname, 'NumHeaderLines', 3);

frame   = data(:, 1);
nFrames = size(data, 1);
nPoints = floor((size(data, 2) - 1) / 3);

xCols = 2:3:(1 + 3*nPoints);
X    = data(:, xCols);          % frames x points
Y    = data(:, xCols + 1);
Conf = data(:, xCols + 2);

%% Replace low-confidence points with NaN
bad = Conf < confThresh;
X(bad) = NaN;
Y(bad) = NaN;
fprintf('Points set to NaN (likelihood < %g): %.1f%%\n', confThresh, 100*mean(bad(:)));

%% Fit a circle to each frame (any frame with >= 3 good points)
xc       = nan(nFrames, 1);
yc       = nan(nFrames, 1);
radius   = nan(nFrames, 1);
relResid = nan(nFrames, 1);
nValid   = sum(~bad, 2);

for i = 1:nFrames
    ok = ~isnan(X(i, :)) & ~isnan(Y(i, :));
    if nnz(ok) >= 3
        xi = X(i, ok)';
        yi = Y(i, ok)';
        [xc(i), yc(i), radius(i)] = fitCircle(xi, yi);
        relResid(i) = sqrt(mean((hypot(xi - xc(i), yi - yc(i)) - radius(i)).^2)) / radius(i);
    end
end
radiusRaw = radius;             % uncleaned trace, plotted in black

%% Step 1: require enough good points
rejFewPts = nValid < minPoints & ~isnan(radius);
radius(rejFewPts) = NaN;
fprintf('Step 1: frames removed with < %d good points: %d\n', minPoints, nnz(rejFewPts));

%% Step 2: blinks (+ padding)
blink = sum(bad, 2) >= blinkMinBad;
blink = movmax(double(blink), 2*blinkPad + 1) > 0;
rejBlink = blink & ~isnan(radius);
radius(blink) = NaN;
fprintf('Step 2: frames removed as blinks (+/- %d frames): %d\n', blinkPad, nnz(rejBlink));

%% Step 3: poor circle fits
rejResid = relResid > maxRelResid & ~isnan(radius);
radius(rejResid) = NaN;
fprintf('Step 3: frames removed for poor fit (rel. residual > %g): %d\n', ...
        maxRelResid, nnz(rejResid));

%% Step 4: remaining spikes (moving-median outlier test)
rejOutlier = false(nFrames, 1);
if doOutliers
    rejOutlier = isoutlier(radius, 'movmedian', outlierWin, ...
                           'ThresholdFactor', outlierK) &~isnan(radius);
    radius(rejOutlier) = NaN;
end
fprintf('Step 4: frames removed as moving-median outliers: %d\n', nnz(rejOutlier));

%% Step 5: fill short gaps
if doFillGaps
    nBefore = nnz(isnan(radius));
    radius  = fillmissing(radius, 'linear', 'MaxGap', maxGapFill);
    fprintf('Step 5: NaN frames filled (gaps <= %d frames): %d\n', ...
            maxGapFill, nBefore - nnz(isnan(radius)));
end

% Keep the centre consistent with the cleaned radius
xc(isnan(radius) & ~isnan(radiusRaw)) = NaN;
yc(isnan(radius) & ~isnan(radiusRaw)) = NaN;

diameterRaw = 2 * radiusRaw;
diameter    = 2 * radius;
fprintf('Total frames without a valid radius: %d of %d (%.2f%%)\n', ...
        nnz(isnan(radius)), nFrames, 100*mean(isnan(radius)));

%% Time axis
if isempty(fps)
    t = frame;       tLabel = 'Frame';
else
    t = frame / fps; tLabel = 'Time (s)';
end

%% Figure 1: radius and diameter, before (black) vs after cleaning (colour)
figure('Name', 'Pupil size over time');
ax1 = subplot(2, 1, 1);
plot(t, radiusRaw, 'k-'); hold on;
plot(t, radius, 'b-');    hold off;
ylabel('Radius (pixels)');
title(sprintf('Pupil circle fit  (likelihood < %g set to NaN)', confThresh));
legend({'before cleaning', 'after cleaning'}, 'Location', 'best');
grid on;

ax2 = subplot(2, 1, 2);
plot(t, diameterRaw, 'k-'); hold on;
plot(t, diameter, 'r-');    hold off;
ylabel('Diameter (pixels)');
xlabel(tLabel);
legend({'before cleaning', 'after cleaning'}, 'Location', 'best');
grid on;
linkaxes([ax1 ax2], 'x');

%% Figure 2: check the fit on the first valid frame
row = find(~isnan(radius), 1);
th = linspace(0, 2*pi, 200);
figure('Name', 'Circle fit check');
hold on;
plot(data(row, xCols), data(row, xCols + 1), 'o', 'Color', [0.7 0.7 0.7]); % all points
plot(X(row, :), Y(row, :), 'ko', 'MarkerFaceColor', 'k');                 % points used
plot(xc(row) + radius(row)*cos(th), yc(row) + radius(row)*sin(th), 'r-', 'LineWidth', 1.5);
plot(xc(row), yc(row), 'r+');
hold off;
axis equal; set(gca, 'YDir', 'reverse'); grid on; box on;
xlabel('x (pixels)'); ylabel('y (pixels)');
legend({'all points', 'used in fit', 'fitted circle', 'centre'}, 'Location', 'best');
title(sprintf('Frame %g: r = %.2f px, d = %.2f px', frame(row), radius(row), diameter(row)));

%% Results table (saved alongside the CSV)
results = table(frame, xc, yc, radiusRaw, diameterRaw, radius, diameter, nValid, ...
    relResid, rejFewPts, rejBlink, rejResid, rejOutlier, ...
    'VariableNames', {'frame', 'centerX', 'centerY', 'radiusRaw', 'diameterRaw', ...
    'radius', 'diameter', 'nValidPoints', 'relResidual', ...
    'removedFewPoints', 'removedBlink','removedPoorFit', 'removedOutlier'});
[~, base] = fileparts(fname);
outName = fullfile(files(1).folder, [base '_circlefit.csv']);
writetable(results, outName);
fprintf('Saved results to %s\n', outName);

%% ---------------------------------------------------------------------
function [xc, yc, r] = fitCircle(x, y)
% Algebraic least-squares circle fit (Kasa method).
% Solves x^2 + y^2 + D*x + E*y + F = 0 for D, E, F.
    A = [x, y, ones(size(x))];
    b = -(x.^2 + y.^2);
    p = A \ b;
    xc = -p(1) / 2;
    yc = -p(2) / 2;
    r  = sqrt(xc^2 + yc^2 - p(3));
end

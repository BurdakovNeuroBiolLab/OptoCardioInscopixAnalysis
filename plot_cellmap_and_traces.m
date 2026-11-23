% plot_cellmap_and_traces.m
% Run from the MATLAB current folder containing the Inscopix exports. The
% three files are found automatically (exactly one of each must be present):
%   *celltrace*.csv with "props" in the name    -> cell props
%   *celltrace*.csv without "props" in the name -> cell traces
%   (celltrace / celltraces, -props / _props all work)
%   *.tif / *.tiff         -> cell map
% Left panel:   the cell map with each ROI outline recoloured and labelled.
% Right panel:  every cell's trace, stacked (C00 at top), in the same colour
%               as its ROI on the map.
%
% The cell map TIFF is the image exported from IDPS with the ROI outlines
% drawn in white. The script works out how much the image was scaled relative
% to the cell set (e.g. 5x), then gives each white outline the colour of the
% nearest cell centroid listed in the props file.

% When called from a parent script that sets keepWorkspace = true, the
% workspace is not cleared so the parent's variables survive.
if ~exist('keepWorkspace', 'var') || ~keepWorkspace
    clear; clc;
end

%% ---------------- Settings ----------------
normalizeTraces = true;   % scale each trace to 0-1 so cells are comparable
spacing         = 1.1;    % vertical gap between traces (normalized units)
skipRejected    = true;   % leave out cells whose status is 'rejected'
imageScale      = [];     % [] = auto-detect; or set manually, e.g. 5
contourMaxFrac  = 0.6;    % outline pixels further than Size*scale*this from
                          % any centroid are drawn grey (unassigned)
lineThickness   = 1;      % 0 = 1-pixel outlines, 1 = 3 px, 2 = 5 px ...
labelFontSize   = 7;

% File-name patterns (regular expressions, case-insensitive) used to find
% each file in the current folder. Adjust if your exports are named differently.
% A file counts if it matches the pattern and does NOT match the exclude
% pattern ('' = exclude nothing).
propsPattern = 'celltrace.*props.*\.csv$'; propsExclude = '';
tracePattern = 'celltrace.*\.csv$';        traceExclude = 'props';
mapPattern   = '\.tiff?$';                 mapExclude   = '';

%% ---------------- Find files in the current folder ----------------
dataDir = pwd;

listing = dir(dataDir);
listing = listing(~[listing.isdir]);
names   = {listing.name};
names   = names(~startsWith(names, '.'));   % skip hidden / macOS '._' files

propsFile = findOneFile(dataDir, names, propsPattern, propsExclude, 'cell PROPS csv');
traceFile = findOneFile(dataDir, names, tracePattern, traceExclude, 'cell TRACES csv');
mapFile   = findOneFile(dataDir, names, mapPattern,   mapExclude,   'cell MAP tif');

[~, pn, pe] = fileparts(propsFile); pf = [pn pe];
[~, tn, te] = fileparts(traceFile); tf = [tn te];
[~, mn, me] = fileparts(mapFile);   mf = [mn me];
fprintf('Folder: %s\n  props : %s\n  traces: %s\n  map   : %s\n', dataDir, pf, tf, mf);

%% ---------------- Read props ----------------
props = readtable(propsFile);
if ~all(ismember({'Name', 'CentroidX', 'CentroidY'}, props.Properties.VariableNames))
    error('"%s" does not look like an IDPS props file (needs Name, CentroidX, CentroidY).', pf);
end
cellNameAll = strtrim(string(props.Name));
if ismember('Status', props.Properties.VariableNames)
    statusAll = lower(strtrim(string(props.Status)));
else
    statusAll = repmat("undecided", height(props), 1);
end
if ismember('Size', props.Properties.VariableNames)
    sizeAll = props.Size(:);
else
    sizeAll = 10 * ones(height(props), 1);
end
cxAll = props.CentroidX(:);
cyAll = props.CentroidY(:);
nAll  = numel(cellNameAll);

%% ---------------- Read traces ----------------
fid = fopen(traceFile, 'r');
line1 = fgetl(fid);
line2 = fgetl(fid);
fclose(fid);

hdr  = strtrim(erase(string(strsplit(line1, ',')), '"'));
row2 = strtrim(erase(string(strsplit(line2, ',')), '"'));
if strcmpi(hdr(1), "Name")
    error('"%s" looks like a props file, not a traces file. Check the file names in the folder.', tf);
end
% IDPS traces have a second header row with cell status (text, not numbers)
nHeader = 1 + any(isnan(str2double(row2(2:end))) & row2(2:end) ~= "");

data = readmatrix(traceFile, 'NumHeaderLines', nHeader);
data = data(:, 1:numel(hdr));
t = data(:, 1);                       % first column = time (s)
traceNames = hdr(2:end);

%% ---------------- Choose cells to plot ----------------
[inTraces, traceCol] = ismember(cellNameAll, traceNames);
keep = inTraces;
if skipRejected
    keep = keep & statusAll ~= "rejected";
end
plotIdx = find(keep);                 % props order: C00, C01, ...
nP = numel(plotIdx);
if nP == 0
    error('No cells in the props file match the trace column names.');
end

% Distinct colours: golden-ratio hue steps so neighbours differ
hues = mod((0:nAll-1)' * 0.61803398875, 1);
colsAll = hsv2rgb([hues, 0.85 * ones(nAll, 1), ones(nAll, 1)]);

%% ---------------- Read cell map ----------------
img = imread(mapFile);                % first page if multi-page
if isinteger(img)
    maxVal = double(intmax(class(img)));
else
    maxVal = 1;
end
img = double(img);
if size(img, 3) >= 3
    img = img(:, :, 1:3);
    isWhite = all(img == maxVal, 3);
    gray = mean(img, 3);
else
    isWhite = img == maxVal;
    gray = img;
end
[H, W] = size(gray);

% Background: hide the white outlines, then stretch contrast
bg = gray;
bg(isWhite) = localPrctile(gray(~isWhite), 50);
lo = localPrctile(bg(:), 1);
hi = localPrctile(bg(:), 99.8);
bg = min(max((bg - lo) / max(hi - lo, eps), 0), 1);

[wy, wx] = find(isWhite);

%% ---------------- Image scale (cell-set pixels -> image pixels) -------
if isempty(imageScale)
    if isempty(wx)
        imageScale = 1;
        warning(['No white ROI outlines found in the TIFF; assuming scale 1. ' ...
                 'Set imageScale manually if the labels are misplaced.']);
    else
        sub = 1:max(1, floor(numel(wx) / 20000)):numel(wx);   % limit memory
        sx = wx(sub)'; sy = wy(sub)';
        score = @(s) median(sqrt(min(((cxAll + 0.5) * s + 0.5 - sx).^2 + ...
                                     ((cyAll + 0.5) * s + 0.5 - sy).^2, [], 2)));
        coarse = 0.5:0.1:12;
        sc = arrayfun(score, coarse);
        [~, b] = min(sc);
        fine = max(0.5, coarse(b) - 0.1):0.01:(coarse(b) + 0.1);
        sf = arrayfun(score, fine);
        [~, b] = min(sf);
        imageScale = fine(b);
    end
    fprintf('Image scale: %.2fx (cell-set pixels -> TIFF pixels)\n', imageScale);
end

cxI = (cxAll + 0.5) * imageScale + 0.5;   % centroids in image coordinates
cyI = (cyAll + 0.5) * imageScale + 0.5;

%% ---------------- Recolour outlines ----------------
Rch = bg; Gch = bg; Bch = bg;
if ~isempty(wx)
    D2 = (cxI - wx').^2 + (cyI - wy').^2;      % cells x outline pixels
    [dmin2, nearest] = min(D2, [], 1);
    maxR = sizeAll * imageScale * contourMaxFrac;
    assigned = sqrt(dmin2(:)) < maxR(nearest(:));
    pixCol = 0.6 * ones(numel(wx), 3);          % unassigned: grey
    for k = 1:nAll
        sel = assigned & nearest(:) == k;
        if keep(k)
            pixCol(sel, :) = repmat(colsAll(k, :), nnz(sel), 1);
        end                                     % rejected/skipped stay grey
    end
    for dy = -lineThickness:lineThickness
        for dx = -lineThickness:lineThickness
            yy = min(max(wy + dy, 1), H);
            xx = min(max(wx + dx, 1), W);
            lin = sub2ind([H, W], yy, xx);
            Rch(lin) = pixCol(:, 1);
            Gch(lin) = pixCol(:, 2);
            Bch(lin) = pixCol(:, 3);
        end
    end
end
rgb = cat(3, Rch, Gch, Bch);

%% ---------------- Plot ----------------
fig = figure('Color', 'k', 'Name', mf, 'Position', [50 100 1600 800]);

% --- Left: map ---
axMap = axes(fig, 'Position', [0.02 0.06 0.46 0.88]);
image(axMap, rgb);
axis(axMap, 'image', 'off');
hold(axMap, 'on');
for j = 1:nP
    k = plotIdx(j);
    text(axMap, cxI(k) + 0.6 * sizeAll(k) * imageScale * 0.5, ...
         cyI(k) - 0.6 * sizeAll(k) * imageScale * 0.5, cellNameAll(k), ...
         'Color', colsAll(k, :), 'FontSize', labelFontSize, 'FontWeight', 'bold');
end
hold(axMap, 'off');
title(axMap, strrep(mf, '_', '\_'), 'Color', 'w');

% --- Right: stacked traces ---
tr = data(:, 1 + traceCol(plotIdx));
if normalizeTraces
    mn = min(tr, [], 1, 'omitnan');
    mx = max(tr, [], 1, 'omitnan');
    rg = mx - mn; rg(rg == 0) = 1;
    tr = (tr - mn) ./ rg;
    gap = spacing;
    labelY = 0.5;
else
    gap = spacing * max(max(tr, [], 'omitnan') - min(tr, [], 'omitnan'));
    labelY = 0;
end
offsets = (nP - 1:-1:0) * gap;                 % first cell at the top

axTr = axes(fig, 'Position', [0.55 0.08 0.43 0.86]);
hold(axTr, 'on');
for j = 1:nP
    k = plotIdx(j);
    plot(axTr, t, tr(:, j) + offsets(j), 'Color', colsAll(k, :), 'LineWidth', 0.7);
    text(axTr, t(1) - 0.005 * (t(end) - t(1)), offsets(j) + labelY, cellNameAll(k), ...
         'Color', colsAll(k, :), 'FontSize', labelFontSize, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
hold(axTr, 'off');
xlim(axTr, [t(1) t(end)]);
ylim(axTr, [min(offsets) - 0.2 * gap, max(offsets) + gap]);
set(axTr, 'Color', 'k', 'XColor', 'w', 'YColor', 'none', ...
    'TickDir', 'out', 'FontSize', 8, 'Box', 'off');
xlabel(axTr, 'Time (s)', 'Color', 'w');

set(fig, 'InvertHardcopy', 'off');             % keep black background when saving

%% ---------------- Helper ----------------
function v = localPrctile(x, p)
% Percentile without the Statistics Toolbox
x = sort(x(~isnan(x(:))));
if isempty(x), v = 0; return; end
v = x(max(1, min(numel(x), round(p / 100 * numel(x)))));
end

function f = findOneFile(folder, names, pattern, exclude, label)
% Return the full path of the single file in NAMES that matches PATTERN and
% does not match EXCLUDE (both case-insensitive regular expressions).
% Errors if no file or more than one file matches.
isHit = ~cellfun(@isempty, regexpi(names, pattern, 'once'));
if ~isempty(exclude)
    isHit = isHit & cellfun(@isempty, regexpi(names, exclude, 'once'));
end
hit = names(isHit);
if isempty(hit)
    csvs = names(~cellfun(@isempty, regexpi(names, '\.(csv|tiff?)$', 'once')));
    if isempty(csvs), csvs = {'(none)'}; end
    error(['No %s found in:\n  %s\n(expected a file name matching "%s"%s).\n' ...
           'CSV/TIFF files in this folder:\n  %s'], label, folder, pattern, ...
          excludeNote(exclude), strjoin(csvs, sprintf('\n  ')));
elseif numel(hit) > 1
    error('Found %d candidate %s files in:\n  %s\nKeep only one of:\n  %s', ...
          numel(hit), label, folder, strjoin(hit, sprintf('\n  ')));
end
f = fullfile(folder, hit{1});
end

function s = excludeNote(exclude)
if isempty(exclude), s = ''; else, s = sprintf(' and not "%s"', exclude); end
end

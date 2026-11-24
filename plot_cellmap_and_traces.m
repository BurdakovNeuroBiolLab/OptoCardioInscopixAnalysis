% plot_cellmap_and_traces.m
% Run from the MATLAB current folder containing the Inscopix exports. The
% three files are found automatically (exactly one of each must be present):
%   *celltrace.csv         -> cell traces
%   *celltrace-props.csv   -> cell props
%   *cellmap.tif / .tiff   -> cell map
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
% Cells to plot as traces, by position in the props file (1 = C00, 2 = C01,
% ...), e.g. 1:10 or 5:20. [] = all cells. Rejected cells are skipped.
% The remaining cells go into a second figure, "Unlikely to be cells".
% A parent script can set this before calling.
if ~exist('nCellsToShow', 'var')
    nCellsToShow = [];
end

% File-name patterns (regular expressions, case-insensitive) used to find
% each file in the current folder. Adjust if your exports are named differently.
propsPattern = 'celltrace-props\.csv$';    % e.g. m945_celltrace-props.csv
tracePattern = 'celltrace\.csv$';          % e.g. m945_celltrace.csv
mapPattern   = 'cellmap\.tiff?$';          % e.g. m945_cellmap.tif (not *_cellmap_traces.tif)

%% ---------------- Find files in the current folder ----------------
dataDir = pwd;

listing = dir(dataDir);
listing = listing(~[listing.isdir]);
names   = {listing.name};
names   = names(~startsWith(names, '.'));   % skip hidden files, incl. macOS '._' files

propsFile = findOneFile(dataDir, names, propsPattern, 'cell PROPS csv');
traceFile = findOneFile(dataDir, names, tracePattern, 'cell TRACES csv');
mapFile   = findOneFile(dataDir, names, mapPattern,   'cell MAP tif');

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
okIdx = find(keep);                   % props order: C00, C01, ...
if isempty(okIdx)
    error('No cells in the props file match the trace column names.');
end
if isempty(nCellsToShow)
    inRange = true(nAll, 1);
else
    r = nCellsToShow(:);
    if ~isnumeric(r) || any(r < 1) || any(r ~= round(r))
        error('nCellsToShow must be whole numbers >= 1 (e.g. 1:10 or 5:20), or [] for all cells.');
    end
    if any(r > nAll)
        warning('nCellsToShow goes up to %d but the props file has %d cells; ignoring the rest.', ...
                max(r), nAll);
        r = r(r <= nAll);
    end
    inRange = false(nAll, 1);
    inRange(r) = true;
end
plotIdx  = okIdx(inRange(okIdx));     % shown as traces (main figure)
otherIdx = okIdx(~inRange(okIdx));    % "Unlikely to be cells" figure
nP = numel(plotIdx);
if nP == 0
    error('None of the cells in nCellsToShow are available (all rejected or out of range).');
end
fprintf('Showing %d cells; %d in "Unlikely to be cells"\n', nP, numel(otherIdx));

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

%% ---------------- Assign outline pixels to cells ----------------
% owner(i) = cell that outline pixel i belongs to (0 = none, drawn grey)
owner = zeros(numel(wx), 1);
if ~isempty(wx)
    D2 = (cxI - wx').^2 + (cyI - wy').^2;      % cells x outline pixels
    [dmin2, nearest] = min(D2, [], 1);
    maxR = sizeAll * imageScale * contourMaxFrac;
    assigned = sqrt(dmin2(:)) < maxR(nearest(:));
    owner(assigned) = nearest(assigned);
end

P = struct('cxI', cxI, 'cyI', cyI, 'sizeAll', sizeAll, 'imageScale', imageScale, ...
           'names', cellNameAll, 'cols', colsAll, 'labelFontSize', labelFontSize, ...
           't', t, 'data', data, 'traceCol', traceCol, ...
           'normalizeTraces', normalizeTraces, 'spacing', spacing);

%% ---------------- Figure 1: selected cells ----------------
% Selected cells in colour; all other outlines grey.
rgb = makeOutlineRGB(bg, wy, wx, owner, plotIdx, [], colsAll, lineThickness);
[fig, axMap, axTr] = drawMapAndTraces(rgb, plotIdx, P, mf, strrep(mf, '_', '\_'));

%% ---------------- Figure 2: unlikely to be cells ----------------
% Cells outside nCellsToShow in colour (map and traces); the selected cells
% are drawn black on the map and not plotted as traces.
if ~isempty(otherIdx)
    rgbOther = makeOutlineRGB(bg, wy, wx, owner, otherIdx, plotIdx, colsAll, lineThickness);
    [figOther, axMapOther, axTrOther] = drawMapAndTraces(rgbOther, otherIdx, P, ...
        'Unlikely to be cells', ['Unlikely to be cells: ' strrep(mf, '_', '\_')]);
end

%% ---------------- Helpers ----------------
function rgb = makeOutlineRGB(bg, wy, wx, owner, colourIdx, blackIdx, cols, lineThickness)
% Greyscale background with the ROI outlines drawn in: cells in colourIdx in
% their own colour, cells in blackIdx black, everything else grey.
[H, W] = size(bg);
Rch = bg; Gch = bg; Bch = bg;
if ~isempty(wx)
    pixCol = 0.6 * ones(numel(wx), 3);
    isCol = ismember(owner, colourIdx);
    pixCol(isCol, :) = cols(owner(isCol), :);
    pixCol(ismember(owner, blackIdx), :) = 0;
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
end

function [fig, axMap, axTr] = drawMapAndTraces(rgb, idx, P, figName, mapTitle)
% Cell map (left, labels for the cells in idx) and their stacked traces
% (right, first cell at the top) in the cells' colours.
nP  = numel(idx);
fig = figure('Color', 'k', 'Name', figName, 'Position', [50 100 1600 800]);

% --- Left: map ---
axMap = axes(fig, 'Position', [0.02 0.06 0.46 0.88]);
image(axMap, rgb);
axis(axMap, 'image', 'off');
hold(axMap, 'on');
for j = 1:nP
    k = idx(j);
    text(axMap, P.cxI(k) + 0.6 * P.sizeAll(k) * P.imageScale * 0.5, ...
         P.cyI(k) - 0.6 * P.sizeAll(k) * P.imageScale * 0.5, P.names(k), ...
         'Color', P.cols(k, :), 'FontSize', P.labelFontSize, 'FontWeight', 'bold');
end
hold(axMap, 'off');
title(axMap, mapTitle, 'Color', 'w');

% --- Right: stacked traces ---
t  = P.t;
tr = P.data(:, 1 + P.traceCol(idx));
if P.normalizeTraces
    mn = min(tr, [], 1, 'omitnan');
    mx = max(tr, [], 1, 'omitnan');
    rg = mx - mn; rg(rg == 0) = 1;
    tr = (tr - mn) ./ rg;
    gap = P.spacing;
    labelY = 0.5;
else
    gap = P.spacing * max(max(tr, [], 'omitnan') - min(tr, [], 'omitnan'));
    labelY = 0;
end
offsets = (nP - 1:-1:0) * gap;                 % first cell at the top

axTr = axes(fig, 'Position', [0.55 0.08 0.43 0.86]);
hold(axTr, 'on');
for j = 1:nP
    k = idx(j);
    plot(axTr, t, tr(:, j) + offsets(j), 'Color', P.cols(k, :), 'LineWidth', 0.7);
    text(axTr, t(1) - 0.005 * (t(end) - t(1)), offsets(j) + labelY, P.names(k), ...
         'Color', P.cols(k, :), 'FontSize', P.labelFontSize, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end
hold(axTr, 'off');
xlim(axTr, [t(1) t(end)]);
ylim(axTr, [min(offsets) - 0.2 * gap, max(offsets) + gap]);
set(axTr, 'Color', 'k', 'XColor', 'w', 'YColor', 'none', ...
    'TickDir', 'out', 'FontSize', 8, 'Box', 'off');
xlabel(axTr, 'Time (s)', 'Color', 'w');

set(fig, 'InvertHardcopy', 'off');             % keep black background when saving
end

function v = localPrctile(x, p)
% Percentile without the Statistics Toolbox
x = sort(x(~isnan(x(:))));
if isempty(x), v = 0; return; end
v = x(max(1, min(numel(x), round(p / 100 * numel(x)))));
end

function f = findOneFile(folder, names, pattern, label)
% Return the full path of the single file in NAMES matching PATTERN.
% Errors if no file or more than one file matches.
hit = names(~cellfun(@isempty, regexpi(names, pattern, 'once')));
if isempty(hit)
    error('No %s found in:\n  %s\n(expected a file name matching "%s").', ...
          label, folder, pattern);
elseif numel(hit) > 1
    error('Found %d candidate %s files in:\n  %s\nKeep only one of:\n  %s', ...
          numel(hit), label, folder, strjoin(hit, sprintf('\n  ')));
end
f = fullfile(folder, hit{1});
end

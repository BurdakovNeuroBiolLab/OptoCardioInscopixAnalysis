%% B_AlignHekaToCellTraces
% Parent script: loads the Inscopix cell traces and the HEKA recording from
% the MATLAB current folder, then downsamples extractedMatGRINcrop so it has
% one row per cell-trace frame.
%
% Steps
%   1. inscopix/plot_cellmap_and_traces.m  -> t (s), data, cell map figure
%   2. heka/A_Extract_HekaData_GRIN.m      -> extractedMatGRINcrop (2.5 ms)
%   3. Bin the HEKA samples into the cell-trace frames
%
% The current folder must contain the files both child scripts look for:
%   *celltraces.csv, *celltraces-props.csv, one cell map .tif, and the HEKA
%   export *L.mat (A_Extract_HekaData_GRIN uses the first *L.mat it finds).
%
% Output (in the workspace, and saved to <traces name>_hekaAligned.mat)
%   hekaDS        (nFrames x 7), same columns as extractedMatGRINcrop except
%                 column 1 is the cell-trace time (s):
%                 1 time (s)  2 ECG  3 run  4 Laser  5 piezo  6 GRIN  7 Cam
%   hekaDSTitles  column names
%   t, traceNames, data   cell traces as read by plot_cellmap_and_traces
%
% Alignment
%   The GRIN synch channel is high for the whole Inscopix recording, so its
%   rising edge is the first trace frame and its falling edge the end of the
%   last frame. The trace frames are spread evenly over that period
%   (stretched slightly if the HEKA and Inscopix clocks disagree), and each
%   frame takes the HEKA samples that fall inside it. extractedMatGRINcrop
%   also requires the camera sync, so if the camera started late or stopped
%   early, frames outside the crop are NaN.

clear; clc;

%% ---------------- Settings ----------------
durationTol = 0.02;   % warn if GRIN-high time and trace length differ by more
saveAligned = true;   % save hekaDS etc. to a .mat file in the current folder

% How each HEKA column is combined within one trace frame:
%   'mean' = average (analogue; binary channels give the fraction of time high)
%   'sum'  = total (run edge events -> encoder ticks per frame)
%   'max'  = 1 if high at any point in the frame
alignMethod = {'time', 'mean', 'sum', 'mean', 'mean', 'mean', 'mean'};
%              time    ECG     run    Laser   piezo   GRIN    Cam

%% ---------------- Put the child scripts on the path ----------------
alignRepoDir = fileparts(mfilename('fullpath'));
if ~isempty(alignRepoDir)
    addpath(fullfile(alignRepoDir, 'inscopix'), fullfile(alignRepoDir, 'heka'));
end
if isempty(which('plot_cellmap_and_traces')) || isempty(which('A_Extract_HekaData_GRIN'))
    error(['Cannot find plot_cellmap_and_traces.m or A_Extract_HekaData_GRIN.m. ' ...
           'Add the repository''s inscopix and heka folders to the MATLAB path.']);
end

hekaFiles = dir('*L.mat');
hekaFiles = hekaFiles(~startsWith({hekaFiles.name}, '._'));
if isempty(hekaFiles)
    error('No HEKA export (*L.mat) found in:\n  %s', pwd);
elseif numel(hekaFiles) > 1
    error('Found %d HEKA exports (*L.mat) in:\n  %s\nKeep only one of:\n  %s', ...
          numel(hekaFiles), pwd, strjoin({hekaFiles.name}, sprintf('\n  ')));
end

%% ---------------- 1. Cell traces + map ----------------
keepWorkspace = true;          % stops the child script clearing our settings
plot_cellmap_and_traces;       % defines t, data, traceNames, tf, ...

%% ---------------- 2. HEKA ----------------
A_Extract_HekaData_GRIN;       % defines extractedMatGRIN(crop), grinRise/Fall, fileName
if ~exist('extractedMatGRINcrop', 'var')
    error('A_Extract_HekaData_GRIN did not produce extractedMatGRINcrop (no *L.mat in %s?).', pwd);
end
fprintf('HEKA file: %s\n', fileName);

%% ---------------- 3. Align and downsample ----------------
% HEKA time (s) of the GRIN-high period in the uncropped data, from the
% rise/fall indices found by A_Extract_HekaData_GRIN
hekaMs    = extractedMatGRIN(:, 1);
hekaDtMs  = median(diff(hekaMs));
grinOnS   = hekaMs(grinRise) / 1000;
grinOffS  = (hekaMs(grinFall) + hekaDtMs) / 1000;   % end of last high sample
grinDur   = grinOffS - grinOnS;

nFrames   = numel(t);
dtFrame   = median(diff(t));
traceDur  = (t(end) - t(1)) + dtFrame;
stretch   = grinDur / traceDur;

fprintf('GRIN synch high: %.2f s; cell traces: %d frames, %.2f s (stretch %.4f)\n', ...
        grinDur, nFrames, traceDur, stretch);
if abs(stretch - 1) > durationTol
    warning(['GRIN synch (%.2f s) and cell traces (%.2f s) differ by %.1f%%. ' ...
             'Check that the traces and the HEKA file are from the same recording.'], ...
            grinDur, traceDur, 100 * abs(stretch - 1));
end

% Frame start times in HEKA seconds, plus the end of the last frame
binStart = grinOnS + (t(:) - t(1)) * stretch;
binEnd   = grinOffS;
hekaT    = extractedMatGRINcrop(:, 1) / 1000;          % crop samples, HEKA s

% Frame j collects the HEKA samples in [binStart(j), binStart(j+1))
binIdx = discretize(hekaT, [binStart(:); binEnd]);
nCols  = size(extractedMatGRINcrop, 2);
hekaDS = nan(nFrames, nCols);
hekaDS(:, 1) = t(:);
for c = 2:nCols
    x  = extractedMatGRINcrop(:, c);
    ok = ~isnan(binIdx) & ~isnan(x);
    n  = accumarray(binIdx(ok), 1, [nFrames 1]);
    switch alignMethod{c}
        case 'mean'
            v = accumarray(binIdx(ok), x(ok), [nFrames 1]) ./ n;
        case 'sum'
            v = accumarray(binIdx(ok), x(ok), [nFrames 1]);
        case 'max'
            v = accumarray(binIdx(ok), x(ok), [nFrames 1], @max);
        otherwise
            error('Unknown alignMethod ''%s'' for column %d.', alignMethod{c}, c);
    end
    v(n == 0) = NaN;               % frames with no HEKA samples
    hekaDS(:, c) = v;
end
hekaDSTitles = ['time (s)', outTitles(2:end)];

nEmpty = nnz(isnan(hekaDS(:, 2)));
if nEmpty > 0
    warning(['%d of %d trace frames have no HEKA samples: they fall outside ' ...
             'extractedMatGRINcrop (camera sync shorter than GRIN synch).'], nEmpty, nFrames);
end

%% ---------------- Check figure ----------------
% HEKA at full rate (grey) and downsampled (colour) on the trace time axis,
% with the mean of the plotted cell traces underneath.
hekaTplot = (hekaT - grinOnS) / stretch + t(1);
showCols  = 2:nCols;
nShow     = numel(showCols) + 1;
figure('Name', 'HEKA aligned to cell traces', 'Position', [100 100 900 900]);
axAlign = gobjects(nShow, 1);
for p = 1:numel(showCols)
    c = showCols(p);
    axAlign(p) = subplot(nShow, 1, p);
    if ~strcmp(alignMethod{c}, 'sum')     % sums are on a different scale
        plot(hekaTplot, extractedMatGRINcrop(:, c), 'Color', [0.75 0.75 0.75]);
        hold on
    end
    plot(t, hekaDS(:, c), 'Color', [0 0.45 0.74], 'LineWidth', 1);
    hold off
    ylabel(hekaDSTitles{c});
end
axAlign(nShow) = subplot(nShow, 1, nShow);
plot(t, mean(tr, 2, 'omitnan'), 'k');
ylabel('mean trace');
xlabel('Time (s)');
linkaxes(axAlign, 'x');
xlim(axAlign(1), [t(1) t(end)]);

%% ---------------- Save ----------------
if saveAligned
    [~, tracesBase] = fileparts(traceFile);
    alignedFile = fullfile(pwd, [tracesBase '_hekaAligned.mat']);
    hekaFile = fileName;
    save(alignedFile, 'hekaDS', 'hekaDSTitles', 't', 'data', 'traceNames', ...
         'binStart', 'stretch', 'alignMethod', 'hekaFile', 'traceFile');
    fprintf('Saved %s\n', alignedFile);
end

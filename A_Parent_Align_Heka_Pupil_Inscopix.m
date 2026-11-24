%% A_Parent_Align_Heka_Pupil_Inscopix
% Parent script: loads the pupil video data, the Inscopix cell traces and the
% HEKA recording from the MATLAB current folder, puts the pupil on the HEKA
% time base, then downsamples extractedMatGRINcrop so it has one row per
% cell-trace frame.
%
% Steps
%   1. pupilCircleFit.m           -> pupil diameter per video frame
%   2. plot_cellmap_and_traces.m  -> t (s), data, cell map figure
%   3. Extract_HekaData_GRIN.m    -> extractedMatGRINcrop (2.5 ms)
%   4. alignPupilToHeka.m         -> pupil added as column 8, using
%                                    the Cam Synch high period
%   5. Bin the HEKA samples into the cell-trace frames
%   6. Heart rate and HRV from the ECG (continuous_hrv_mouse.m)
%   7. Running, heart rate, RMSSD, laser and pupil added under the traces in
%      the plot_cellmap_and_traces figure
%
% The current folder must contain exactly one of each:
%   *zoompup*.csv (DeepLabCut pupil), *celltrace.csv, *celltrace-props.csv,
%   *cellmap.tif, and the HEKA export *L.mat. Files starting with ._ are ignored.
%
% Output (in the workspace, and saved to <traces name>_hekaAligned.mat)
%   hekaDS        (nFrames x 8), same columns as extractedMatGRINcrop except
%                 column 1 is the cell-trace time (s):
%                 1 time (s)  2 ECG  3 run  4 Laser  5 piezo  6 GRIN  7 Cam
%                 8 pupil diameter (px)
%   hekaDSTitles  column names
%   runDS         running, encoder ticks per second, per frame
%   hrDS          beat-to-beat heart rate (bpm) from ECG R peaks, per frame
%   hrvDS         windowed HRV per frame (struct): HR (bpm), RMSSD (ms), CVSD,
%                 CV_IBI, HRV_resid, nIBI
%   hrvMetrics    continuous_hrv_mouse output (window centres in HEKA s)
%   alpha1Whole   DFA alpha1 over all valid IBIs of the crop
%   rT            R-peak times (HEKA s), refined to sub-sample precision
%   laserDS       fraction of each frame the laser was on
%   pupilDS       pupil diameter (px), per frame
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
nCellsToShow = [];    % cells shown as traces in the final graph, by position in
                      % the props file (1 = C00), e.g. 1:10 or 5:20; [] = all.
                      % The other cells go to a second figure, "Unlikely to be
                      % cells" (their traces in colour, these cells black on the map).
durationTol = 0.02;   % warn if GRIN-high time and trace length differ by more
saveAligned = true;   % save hekaDS etc. to a .mat file in the current folder
saveFigure  = true;   % save the main figure (cell map, selected cells' traces,
                      % HEKA/pupil panels) as .tif and .svg in the current folder
figureDpi   = 300;    % resolution of the .tif

% Heart rate from the ECG (defaults for mouse)
hrRange       = [250 900]; % bpm; beats outside this range are dropped
hrThreshFrac  = 0.5;       % R-peak threshold as a fraction of the ECG 99.5th percentile
hrBaselineWin = 0.05;      % s, moving-mean window removed from the ECG
showHRCheck   = false;     % true = extra figure with the detected R peaks
hrvWindowSec  = 10;        % s, sliding window for HRV (continuous_hrv_mouse)
hrvStepSec    = 0.25;      % s, step between HRV windows
ibiRelTol     = 0.3;       % drop IBIs > 30% away from the running median (missed /
                           % extra beats), for both HR and HRV

% How each HEKA column is combined within one trace frame:
%   'mean' = average (analogue; binary channels give the fraction of time high)
%   'sum'  = total (run edge events -> encoder ticks per frame)
%   'max'  = 1 if high at any point in the frame
alignMethod = {'time', 'mean', 'sum', 'mean', 'mean', 'mean', 'mean', 'mean'};
%              time    ECG     run    Laser   piezo   GRIN    Cam     pupil

%% ---------------- Put the child scripts on the path ----------------
alignRepoDir = fileparts(mfilename('fullpath'));
if ~isempty(alignRepoDir)
    addpath(alignRepoDir);
end
childScripts = {'pupilCircleFit', 'plot_cellmap_and_traces', ...
                'Extract_HekaData_GRIN', 'alignPupilToHeka', ...
                'continuous_hrv_mouse', 'dfa_alpha1_local', 'flag_bad_ibi'};
for s = 1:numel(childScripts)
    if isempty(which(childScripts{s}))
        error(['Cannot find %s.m. Keep all the repository''s .m files in one ' ...
               'folder on the MATLAB path.'], childScripts{s});
    end
end

pupilFiles = dir('*zoompup*.csv');
pupilFiles = pupilFiles(~startsWith({pupilFiles.name}, '._') & ...
                        ~endsWith({pupilFiles.name}, '_circlefit.csv'));
if isempty(pupilFiles)
    error('No pupil CSV (*zoompup*.csv) found in:\n  %s', pwd);
elseif numel(pupilFiles) > 1
    error('Found %d pupil CSVs (*zoompup*.csv) in:\n  %s\nKeep only one of:\n  %s', ...
          numel(pupilFiles), pwd, strjoin({pupilFiles.name}, sprintf('\n  ')));
end

hekaFiles = dir('*L.mat');
hekaFiles = hekaFiles(~startsWith({hekaFiles.name}, '._'));
if isempty(hekaFiles)
    error('No HEKA export (*L.mat) found in:\n  %s', pwd);
elseif numel(hekaFiles) > 1
    error('Found %d HEKA exports (*L.mat) in:\n  %s\nKeep only one of:\n  %s', ...
          numel(hekaFiles), pwd, strjoin({hekaFiles.name}, sprintf('\n  ')));
end

keepWorkspace = true;          % stops the child scripts clearing our settings

%% ---------------- 1. Pupil ----------------
% Runs first: it uses variable names (data, t, ...) that the cell-trace
% script reuses, so the results are copied to pupil* names straight away.
pupilCircleFit;                % defines diameter (cleaned, px), fname
pupilDiameter = diameter;
pupilFile     = fname;

%% ---------------- 2. Cell traces + map ----------------
plot_cellmap_and_traces;       % defines t, data, traceNames, tf, ...

%% ---------------- 3. HEKA ----------------
Extract_HekaData_GRIN;       % defines extractedMatGRIN(crop), grin/camRise/Fall, fileName
if ~exist('extractedMatGRINcrop', 'var')
    error('Extract_HekaData_GRIN did not produce extractedMatGRINcrop (no *L.mat in %s?).', pwd);
end
fprintf('HEKA file: %s\n', fileName);

%% ---------------- 4. Pupil onto the HEKA time base ----------------
alignPupilToHeka;              % adds column 8 'pupil diameter'

%% ---------------- 5. Align and downsample ----------------
% HEKA time (s) of the GRIN-high period in the uncropped data, from the
% rise/fall indices found by Extract_HekaData_GRIN
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

%% ---------------- 6. Heart rate and HRV from the ECG ----------------
% R peaks are found in the full-rate ECG of the crop, refined to sub-sample
% precision, turned into beat-to-beat heart rate and windowed HRV, and
% sampled at the centre of each trace frame.
hekaFs   = 1000 / hekaDtMs;                                % Hz
ecg      = extractedMatGRINcrop(:, 2);
ecgHP    = ecg - movmean(ecg, max(3, round(hrBaselineWin * hekaFs)), 'omitnan');
if abs(localPrctile(ecgHP, 0.1)) > abs(localPrctile(ecgHP, 99.9))
    ecgHP = -ecgHP;                                        % R peaks point down
end
hrThresh = hrThreshFrac * localPrctile(ecgHP, 99.5);
cand = find(ecgHP(2:end-1) > hrThresh & ecgHP(2:end-1) >= ecgHP(1:end-2) & ...
            ecgHP(2:end-1) > ecgHP(3:end)) + 1;
% Refractory period: of two peaks closer than 60/maxBpm, keep the larger
minGap = round(60 / hrRange(2) * hekaFs);
rIdx = zeros(numel(cand), 1); nR = 0;
for i = 1:numel(cand)
    if nR > 0 && cand(i) - rIdx(nR) < minGap
        if ecgHP(cand(i)) > ecgHP(rIdx(nR)), rIdx(nR) = cand(i); end
    else
        nR = nR + 1; rIdx(nR) = cand(i);
    end
end
rIdx  = rIdx(1:nR);
% Parabolic interpolation through each peak and its two neighbours: the peak
% time is only known to +/- half a sample (1.25 ms at 400 Hz), which would
% otherwise add several bpm of jitter to beat-to-beat HR and a floor of
% ~1.8 ms to RMSSD.
pkShift = zeros(nR, 1);
inner   = rIdx > 1 & rIdx < numel(ecgHP);
ym1 = ecgHP(rIdx(inner) - 1); y0 = ecgHP(rIdx(inner)); yp1 = ecgHP(rIdx(inner) + 1);
den = ym1 - 2 * y0 + yp1;
sh  = 0.5 * (ym1 - yp1) ./ den;
sh(~isfinite(sh)) = 0;
pkShift(inner) = max(-0.5, min(0.5, sh));
rT    = hekaT(rIdx) + pkShift / hekaFs;                     % HEKA s
beatT = (rT(1:end-1) + rT(2:end)) / 2;                     % HEKA s
beatHR = 60 ./ diff(rT);                                   % bpm
ibiRangeSec = sort(60 ./ hrRange);
okHR  = ~flag_bad_ibi(diff(rT), ibiRangeSec, ibiRelTol, 21);
frameEnd = [binStart(2:end); binEnd];
frameMid = (binStart + frameEnd) / 2;
if nnz(okHR) >= 2
    hrDS = interp1(beatT(okHR), beatHR(okHR), frameMid, 'linear', NaN);
else
    hrDS = nan(nFrames, 1);
    warning('Fewer than 2 heart beats detected; check hrThreshFrac / hrRange.');
end
hrDS(isnan(hekaDS(:, 2))) = NaN;                           % outside the crop
fprintf('Heart rate: %d R peaks, median %.0f bpm (%d IBIs dropped: outside hrRange or >%g%% from local median)\n', ...
        nR, median(beatHR(okHR)), nnz(~okHR), 100 * ibiRelTol);

% Windowed HRV. The same IBIs as above (missed / extra beats) are excluded.
hrvFields   = {'HR', 'RMSSD', 'CVSD', 'CV_IBI', 'HRV_resid', 'nIBI'};
hrvDS       = struct();
for f = 1:numel(hrvFields), hrvDS.(hrvFields{f}) = nan(nFrames, 1); end
hrvMetrics  = struct();
if nR >= 3 && rT(end) - rT(1) > hrvWindowSec
    hrvMetrics = continuous_hrv_mouse(rT, hrvWindowSec, hrvStepSec, ibiRangeSec, ibiRelTol);
    for f = 1:numel(hrvFields)
        v = interp1(hrvMetrics.t, hrvMetrics.(hrvFields{f}), frameMid, 'linear', NaN);
        v(isnan(hekaDS(:, 2))) = NaN;                       % outside the crop
        hrvDS.(hrvFields{f}) = v;
    end
    hrvDS.RMSSD = 1000 * hrvDS.RMSSD;                       % s -> ms
else
    warning('Recording shorter than hrvWindowSec (%g s); HRV not computed.', hrvWindowSec);
end

% DFA alpha1 once over the whole crop (too few beats per window to be
% reliable). Uses the valid IBIs only, joined across rejected ones.
ibiAll = diff(rT);
ibiAll = ibiAll(~flag_bad_ibi(ibiAll, ibiRangeSec, ibiRelTol, 21));
alpha1Whole = dfa_alpha1_local(ibiAll);
fprintf('HRV: %g s windows, median RMSSD %.2f ms; DFA alpha1 (whole crop) %.2f\n', ...
        hrvWindowSec, median(hrvDS.RMSSD(~isnan(hrvDS.RMSSD))), alpha1Whole);

runDS   = hekaDS(:, 3) ./ ((frameEnd - binStart) / stretch); % encoder ticks per s
laserDS = hekaDS(:, 4);                                    % fraction of frame on
pupilDS = hekaDS(:, 8);                                    % pupil diameter (px)

%% ---------------- 7. Add running, HR, RMSSD, laser and pupil to the figure
% Same figure as plot_cellmap_and_traces; the traces panel is shortened and
% five panels on the same time axis are added beneath it.
trPos  = get(axTr, 'Position');                % [left bottom width height]
panelH = 0.05; panelGap = 0.012;
nBelow = 5;
belowTop = trPos(2) + nBelow * (panelH + panelGap);
set(axTr, 'Position', [trPos(1) belowTop trPos(3) trPos(2) + trPos(4) - belowTop]);
xlabel(axTr, '');
set(axTr, 'XTickLabel', []);

belowY    = {runDS, hrDS, hrvDS.RMSSD, laserDS, pupilDS};
belowName = {'Run (ticks/s)', 'HR (bpm)', 'RMSSD (ms)', 'Laser', 'Pupil (px)'};
belowCol  = {[0.3 0.8 1], [1 0.35 0.35], [1 0.55 0.8], [0.4 1 0.4], [1 0.8 0.3]};
laserPanel = 4;
axBelow   = gobjects(nBelow, 1);
for p = 1:nBelow
    y0 = belowTop - p * (panelH + panelGap);
    axBelow(p) = axes(fig, 'Position', [trPos(1) y0 trPos(3) panelH]);
    if p == laserPanel
        stairs(axBelow(p), t, belowY{p}, 'Color', belowCol{p}, 'LineWidth', 1);
        ylim(axBelow(p), [-0.1 1.1]);
    else
        plot(axBelow(p), t, belowY{p}, 'Color', belowCol{p}, 'LineWidth', 1);
    end
    set(axBelow(p), 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'TickDir', 'out', ...
        'FontSize', 7, 'Box', 'off');
    ylabel(axBelow(p), belowName{p}, 'Color', belowCol{p}, 'FontSize', 8);
    if p < nBelow
        set(axBelow(p), 'XTickLabel', []);
    end
end
xlabel(axBelow(end), 'Time (s)', 'Color', 'w');
linkaxes([axTr; axBelow], 'x');
xlim(axTr, [t(1) t(end)]);

%% ---------------- Optional: ECG peak check ----------------
if showHRCheck
    figure('Name', 'ECG R-peak detection', 'Position', [100 100 1200 400]);
    plot(hekaT - grinOnS, ecgHP, 'k'); hold on
    plot(rT - grinOnS, ecgHP(rIdx), 'r.', 'MarkerSize', 10);
    yline(hrThresh, 'b--'); hold off
    xlabel('Time from GRIN start (s)'); ylabel('ECG (baseline removed)');
end

%% ---------------- Save data and figure ----------------
if saveAligned
    [~, tracesBase] = fileparts(traceFile);
    alignedFile = fullfile(pwd, [tracesBase '_hekaAligned.mat']);
    hekaFile = fileName;
    save(alignedFile, 'hekaDS', 'hekaDSTitles', 'runDS', 'hrDS', 'laserDS', ...
         'pupilDS', 'hrvDS', 'hrvMetrics', 'alpha1Whole', 'hrvWindowSec', ...
         't', 'data', 'traceNames', 'binStart', 'stretch', ...
         'alignMethod', 'rT', 'pupilHekaT', 'hekaFile', 'traceFile', 'pupilFile');
    fprintf('Saved %s\n', alignedFile);
end

if saveFigure
    [~, tracesBase] = fileparts(traceFile);
    figBase = fullfile(pwd, [tracesBase '_cellmap_traces']);
    set(fig, 'InvertHardcopy', 'off', 'PaperPositionMode', 'auto');  % keep black background
    print(fig, [figBase '.tif'], '-dtiff', sprintf('-r%d', figureDpi));
    set(fig, 'Renderer', 'painters');      % vector lines and text in the .svg
    print(fig, [figBase '.svg'], '-dsvg');
    fprintf('Saved %s.tif and .svg\n', figBase);
end

%% ---------------- Helper ----------------
function v = localPrctile(x, p)
% Percentile without the Statistics Toolbox
x = sort(x(~isnan(x(:))));
if isempty(x), v = 0; return; end
v = x(max(1, min(numel(x), round(p / 100 * numel(x)))));
end

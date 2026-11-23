%% Extract_HekaData_GRIN
% Loads HEKA-exported sweeps (*L.mat), concatenates them per channel,
% binarises the digital channels, turns the run channel into edge events,
% and crops everything to the period where both the GRIN and camera sync
% signals are running.
%
% Outputs (photometry channels are dropped):
%   extractedMatGRIN     (nSamples x 7)
%   extractedMatGRINcrop (cropped to the GRIN/camera overlap)
%   Columns: 1 time (ms, continuous)   5 piezo (raw)
%            2 ECG (raw)               6 GRIN synch (binarised)
%            3 run (edge events)       7 Cam Synch (binarised)
%            4 Laser (binarised)
% Channel 9 of the export is not used.
% Last updated 23.11.2026

%% Settings
nChannels = 9;         % channels per sweep in the HEKA export

cutTheFront    = 0;    % 1 = drop the first cutTheFrontVal samples
cutTheFrontVal = 5000;

sampleInterval = 2.5;     % ms, gap added between sweeps

limitLength = 0;       % 1 = keep only the first superLength samples
superLength = 100000;

% Column order of superRigMatMaster (column 1 = timestamps)
rawTitles = {'time', 'ECG', 'run', 'Laser', 'phot1 Raw', 'phot2 Raw', 'piezo', 'GRIN synch', 'Cam Synch'};
colTime      = 1;   % HEKA timestamps (ms)
colECG       = 2;   % raw
colRun       = 3;   % rotary encoder -> edge events
colLaser     = 4;   % binarised
colPhot1     = 5;   % raw, not kept in output
colPhot2     = 6;   % raw, not kept in output
colPiezo     = 7;   % raw
colGRINSynch = 8;   % binarised
colCamSynch  = 9;   % binarised

% superRig columns kept in extractedMatGRIN, in output order
keepCh    = [colTime colECG colRun colLaser colPiezo colGRINSynch colCamSynch];
outTitles = rawTitles(keepCh);
outGRIN   = find(keepCh == colGRINSynch);   % column in extractedMatGRIN
outCam    = find(keepCh == colCamSynch);

%% Find files, skipping macOS '._' resource files
fileNames = dir('*L.mat');
fileNames = fileNames(~startsWith({fileNames.name}, '._'));

filesToProcess = 1;    % change to 1:numel(fileNames) to process every file

for fileLoop = filesToProcess

    fileName = fileNames(fileLoop).name;
    fprintf('Processing file %d: %s\n', fileLoop, fileName);

    %% Load sweeps (column 1 of each = timestamp, 2.5 ms/sample)
    rawCell = struct2cell(load(fileName));

    % Remove the last sweep (one cell per channel) if it is incomplete
    if size(rawCell{end}, 1) ~= size(rawCell{1}, 1)
        rawCell(end-nChannels+1:end) = [];
    end

    %% Concatenate sweeps into one column per channel
    % Timestamps are taken from the first channel of each sweep and made
    % continuous: each sweep starts one sample interval after the last
    % timestamp of the previous sweep.
    timeCol = cellfun(@(x) x(:, 1), rawCell(1:nChannels:end), 'UniformOutput', false);
    for sw = 2:numel(timeCol)
        timeCol{sw} = timeCol{sw} - timeCol{sw}(1) + timeCol{sw-1}(end) + sampleInterval;
    end
    channelCols = cell(1, nChannels);
    for ch = 1:nChannels
        channelCols{ch} = vertcat(rawCell{ch:nChannels:end});
        channelCols{ch} = channelCols{ch}(:, 2:end);
    end
    superRigMatMaster = [vertcat(timeCol{:}) channelCols{:}];

    if limitLength
        superRigMatMaster = superRigMatMaster(1:superLength, :);
    end

    if contains(fileName, 'absLevelCheck') || cutTheFront
        disp('Cutting first section')
        superRigMatMaster(1:cutTheFrontVal, :) = [];
    end

    timeVec = superRigMatMaster(:, colTime);

    %% Figure 1: all raw HEKA traces (time column not plotted)
    figure('Renderer', 'painters', 'Position', [0 350 400 700]);
    for ch = 1:nChannels
        subplot(nChannels, 1, ch);
        plot(timeVec, superRigMatMaster(:, ch + 1), '-k')
    end
    xlabel('Time (ms)')

    %% Process channels (analogue stay raw, digital are binarised)
    processedMat = superRigMatMaster(:, 1:numel(rawTitles));

    % NOTE: the threshold is half the range, not the midpoint between min
    % and max. It only equals the midpoint if the channel's minimum is ~0.
    digitalCh  = [colRun colLaser colGRINSynch colCamSynch];
    rawDigital = superRigMatMaster(:, digitalCh);
    supThresh  = (max(rawDigital) - min(rawDigital)) / 2;
    processedMat(:, digitalCh) = rawDigital > supThresh;

    % Run (rotary encoder): mark every state change
    processedMat(:, colRun) = [abs(diff(processedMat(:, colRun))); NaN];

    % Drop the photometry channels
    extractedMatGRIN = processedMat(:, keepCh);

    %% Crop to the period where both GRIN and camera sync are running
    % Start = whichever signal rises second (first rising edge)
    % End   = whichever signal falls first  (last falling edge)
    [grinRise, grinFall] = syncOnOff(extractedMatGRIN(:, outGRIN), 'GRIN synch');
    [camRise,  camFall]  = syncOnOff(extractedMatGRIN(:, outCam),  'Cam Synch');

    cropStart = max(grinRise, camRise);
    cropEnd   = min(grinFall, camFall);
    if cropStart >= cropEnd
        error('GRIN and camera sync signals do not overlap in %s', fileName);
    end
    fprintf('  Crop: samples %d to %d (%d samples)\n', ...
        cropStart, cropEnd, cropEnd - cropStart + 1);

    extractedMatGRINcrop = extractedMatGRIN(cropStart:cropEnd, :);
    timeCrop = extractedMatGRINcrop(:, 1);

    %% Figure 2: all output traces, full (black) with crop (blue) overlaid
    figure;
    plotCols = find(keepCh ~= colTime);
    nPlot    = numel(plotCols);
    ax       = gobjects(nPlot, 1);
    for p = 1:nPlot
        k = plotCols(p);
        ax(p) = subplot(nPlot, 1, p);
        plot(timeVec, extractedMatGRIN(:, k), '-k')
        hold on
        plot(timeCrop, extractedMatGRINcrop(:, k), '-b')
        hold off
        ylabel(outTitles{k})
    end
    xlabel('Time (ms)')
    linkaxes(ax, 'x')

end

%% Local functions
function [riseIdx, fallIdx] = syncOnOff(sig, sigName)
% First sample of the first rising edge and last high sample before the
% last falling edge of a binarised sync signal.
d       = diff(sig);
riseIdx = find(d == 1, 1, 'first') + 1;
fallIdx = find(d == -1, 1, 'last');
if isempty(riseIdx) || isempty(fallIdx)
    error('%s has no complete rise and fall', sigName);
end
end

function metrics = continuous_hrv_mouse(beatTimes, windowLengthSec, stepSec, ibiRangeSec, ibiRelTol)
% CONTINUOUS_HRV_MOUSE
%   Compute continuous HRV measures from mouse beat times in sliding windows.
%
% INPUTS
%   beatTimes       - vector of heartbeat (R-peak) times in seconds
%   windowLengthSec - length of sliding window in seconds
%   stepSec         - step between window starts in seconds
%   ibiRangeSec     - (optional) [min max] plausible inter-beat interval (s).
%                     IBIs outside it (missed / extra beats) are set to NaN
%                     and left out of every metric. Default: keep all.
%                     Mouse example: 60 ./ [900 250] = [0.067 0.24].
%   ibiRelTol       - (optional) also reject IBIs that differ from the running
%                     median of 21 IBIs by more than this fraction (e.g. 0.3).
%                     Catches missed beats (IBI doubles), which can still be
%                     inside ibiRangeSec. Default: off. See flag_bad_ibi.m.
%
% OUTPUT struct 'metrics' (all column vectors, one row per window)
%   metrics.t          - centre time of each window (s)
%   metrics.HR         - heart rate per window (bpm), 60 / mean IBI
%   metrics.RMSSD      - root mean square of successive IBI differences (s)
%   metrics.CVSD       - RMSSD / mean IBI
%   metrics.CV_IBI     - SDNN / mean IBI
%   metrics.alpha1_DFA - DFA alpha1 (unreliable for windows < ~60 beats)
%   metrics.logRMSSD   - log(RMSSD)
%   metrics.HRV_resid  - residual log(RMSSD) after regressing out HR
%   metrics.nIBI       - number of IBIs used per window
%
% Successive differences are only taken between two valid, consecutive IBIs,
% so a rejected IBI never contributes a difference.

    if nargin < 4 || isempty(ibiRangeSec)
        ibiRangeSec = [0 Inf];
    end
    if nargin < 5
        ibiRelTol = [];
    end
    beatTimes = beatTimes(:);

    % Flag bad IBIs once over the whole recording; IBI k is beat k -> k+1
    ibiAll = diff(beatTimes);
    badAll = flag_bad_ibi(ibiAll, ibiRangeSec, ibiRelTol, 21);

    % Window start times
    tStart = beatTimes(1);
    tEnd   = beatTimes(end);
    winStarts = tStart:stepSec:(tEnd - windowLengthSec);
    nWin = numel(winStarts);

    % Preallocate
    tCenter   = nan(nWin, 1);
    HR        = nan(nWin, 1);
    RMSSD     = nan(nWin, 1);
    CVSD      = nan(nWin, 1);
    CV_IBI    = nan(nWin, 1);
    alpha1    = nan(nWin, 1);
    logRMSSD  = nan(nWin, 1);
    nIBI      = zeros(nWin, 1);

    for i = 1:nWin
        t0 = winStarts(i);
        t1 = t0 + windowLengthSec;
        tCenter(i) = t0 + windowLengthSec / 2;

        % Beats within the moving window
        b = find(beatTimes > t0 & beatTimes <= t1);
        if numel(b) < 2
            continue
        end
        IBI = ibiAll(b(1):b(end) - 1);          % inter-beat intervals (s)
        IBI(badAll(b(1):b(end) - 1)) = NaN;
        good = IBI(~isnan(IBI));
        nIBI(i) = numel(good);
        if isempty(good)
            continue
        end

        meanIBI = mean(good);
        HR(i)   = 60 / meanIBI;                 % bpm

        diffIBI = diff(IBI);                    % NaN next to a rejected IBI
        diffIBI = diffIBI(~isnan(diffIBI));
        if ~isempty(diffIBI)
            RMSSD(i) = sqrt(mean(diffIBI .^ 2));
        end
        sdnn = sqrt(mean((good - meanIBI) .^ 2));   % population SD, as std(x, 1)

        CVSD(i)   = RMSSD(i) / meanIBI;
        CV_IBI(i) = sdnn / meanIBI;
        alpha1(i) = dfa_alpha1_local(good);
        logRMSSD(i) = log(RMSSD(i));
    end

    % Residual HRV: regress logRMSSD on HR across windows
    valid = isfinite(HR) & isfinite(logRMSSD);     % excludes log(0) = -Inf
    HRV_resid = nan(size(HR));
    if sum(valid) >= 3
        X = [ones(sum(valid), 1), HR(valid)];
        beta = X \ logRMSSD(valid);
        HRV_resid(valid) = logRMSSD(valid) - X * beta;
    end

    % Output
    metrics = struct();
    metrics.t          = tCenter;
    metrics.HR         = HR;
    metrics.CVSD       = CVSD;
    metrics.CV_IBI     = CV_IBI;
    metrics.RMSSD      = RMSSD;
    metrics.alpha1_DFA = alpha1;
    metrics.logRMSSD   = logRMSSD;
    metrics.HRV_resid  = HRV_resid;
    metrics.nIBI       = nIBI;
end

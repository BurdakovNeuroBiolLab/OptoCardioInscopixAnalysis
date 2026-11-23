%% alignPupilToHeka
% Puts the pupil diameter (one value per video frame, from pupilCircleFit)
% on the HEKA time base and adds it as an extra column to extractedMatGRIN
% and extractedMatGRINcrop.
%
% The camera sync (Cam Synch) runs from its first rising edge to its last
% falling edge for exactly the period covered by the pupil CSV. The video
% frames are spread evenly over that period; every HEKA sample takes the
% diameter of the video frame it falls in (NaN outside the camera period
% and in frames removed by pupilCircleFit, e.g. blinks).
%
% Needs in the workspace
%   from A_Extract_HekaData_GRIN: extractedMatGRIN, camRise, camFall,
%       cropStart, cropEnd, outCam, outTitles
%   pupilDiameter   cleaned diameter per video frame (pixels), nPupil x 1
%
% Adds
%   extractedMatGRIN(:, end+1), extractedMatGRINcrop (rebuilt), outTitles{end+1}
%   pupilHekaT      centre of each video frame in HEKA time (ms)
%   pupilPeriodMs   video frame period implied by the camera sync (ms)

if any(strcmp(outTitles, 'pupil diameter'))
    error('extractedMatGRIN already has a pupil column. Re-run A_Extract_HekaData_GRIN first.');
end

hekaMsP  = extractedMatGRIN(:, 1);
dtMsP    = median(diff(hekaMsP));
camOnMs  = hekaMsP(camRise);
camOffMs = hekaMsP(camFall) + dtMsP;          % end of the last high sample

pupilDiameter = pupilDiameter(:);
nPupil        = numel(pupilDiameter);
pupilPeriodMs = (camOffMs - camOnMs) / nPupil;
pupilHekaT    = camOnMs + ((1:nPupil)' - 0.5) * pupilPeriodMs;

pupilFrameIdx = floor((hekaMsP - camOnMs) / pupilPeriodMs) + 1;
inCam         = pupilFrameIdx >= 1 & pupilFrameIdx <= nPupil;
pupilCol      = nan(size(hekaMsP));
pupilCol(inCam) = pupilDiameter(pupilFrameIdx(inCam));

extractedMatGRIN(:, end + 1) = pupilCol;
extractedMatGRINcrop = extractedMatGRIN(cropStart:cropEnd, :);
outTitles{end + 1}   = 'pupil diameter';

% Diagnostic: if the camera sync pulses once per frame, the pulse count
% should match the number of rows in the pupil CSV.
camPulses = nnz(diff(extractedMatGRIN(camRise:camFall, outCam)) == 1) + 1;
fprintf(['Pupil: %d video frames over %.2f s of camera sync -> %.2f Hz ' ...
         '(%d camera sync pulses)\n'], ...
        nPupil, (camOffMs - camOnMs) / 1000, 1000 / pupilPeriodMs, camPulses);

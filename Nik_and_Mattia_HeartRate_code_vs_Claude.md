# Nik and Mattia Heart Rate code vs Claude

Comparison of heart-rate detection in Nik and Mattia's
`preprocessMultipleMice_2.m` (with `continuous_hrv_mouse.m`) and Claude's
version in `A_Parent_Align_Heka_Pupil_Inscopix.m`.

## Comparison (before the updates below)

| | Nik and Mattia's script | Claude's, in `A_Parent_Align_Heka_Pupil_Inscopix` |
|---|---|---|
| Sampling | 1 kHz (`recordedData`) | 400 Hz (HEKA `*L.mat`, 2.5 ms) |
| Filtering | 8–30 Hz bandpass | Subtract a 50 ms moving mean |
| Polarity | Set by hand (`invertECG`) | Detected automatically |
| Threshold | Fixed, 0.15 V | Relative: 0.5 × 99.5th percentile |
| Minimum gap between beats | 50 ms | 66.7 ms (60/900 bpm), keeps the larger peak |
| Bad-beat rejection | None | Drops beats outside 250–900 bpm |
| HR output | 60/mean IBI over a 3 s sliding window (smooth, about 3 s blur) | Beat-to-beat 60/RR, interpolated to each frame |
| HRV | RMSSD, CVSD, CV, DFA α1, residual | None |
| Toolboxes | Signal Processing | None |

The two should agree on average heart rate; Nik and Mattia's is smoother
because of the 3 s window. Claude's had the sampling-resolution problem worse:
at 400 Hz each beat-to-beat HR value jittered by about ±6 bpm (one standard
deviation), and HRV from those R peaks would have an RMSSD floor of about
1.8 ms.

## Updates to Claude's version

1. **Sub-sample R-peak timing.** Each R peak is refined by parabolic
   interpolation through the peak sample and its two neighbours.
2. **Bad-beat rejection.** `flag_bad_ibi.m` drops IBIs outside 250–900 bpm
   *and* IBIs more than 30% away from the running median of 21 IBIs
   (`ibiRelTol`). At 600 bpm a missed beat gives a 0.2 s IBI, which is still
   inside the 250–900 bpm range, so the range check alone does not catch it.
3. **HRV.** `continuous_hrv_mouse.m` (Nik and Mattia's function, with fixes)
   computes HR, RMSSD, CVSD, CV_IBI and the HR-regressed residual in 10 s
   windows, step 0.25 s, from the cleaned IBIs. The results are aligned to the
   cell-trace frames (`hrvDS`), and RMSSD is plotted under the cell traces.
4. **DFA α1** is computed once over the whole crop (`alpha1Whole`) instead of
   per window, and `dfa_alpha1_local.m` caps box sizes at 16 beats (the
   standard α1 range). Windows under ~64 IBIs give the same result as before.

Fixes to `continuous_hrv_mouse.m` (with 3 inputs it gives the same results
as the original on clean data):

- `RMSSD` is preallocated as a column, like the other outputs.
- `isfinite` instead of `~isnan` before the `HRV_resid` regression, so
  `log(0) = -Inf` cannot break it.
- Optional 4th and 5th inputs reject bad IBIs; successive differences are only
  taken between two valid, consecutive IBIs.
- `nIBI` output: number of IBIs used per window.

## Test on synthetic data (Octave)

3000 beats at 600 bpm with 1.5 ms random IBI variation (true RMSSD ≈ 2.14 ms),
R waves 3 ms wide sampled at 400 Hz plus noise:

| | Before | After |
|---|---|---|
| R-peak timing error (rms) | 0.73 ms | 0.11 ms |
| Beat-to-beat HR jitter | 6.2 bpm | 0.9 bpm |
| Median RMSSD (truth 2.14 ms) | 2.72 ms | 2.16 ms |
| RMSSD in the window with one missed beat (truth ≈ 2.1 ms) | 14.70 ms | 2.47 ms |

## Remaining notes on Nik and Mattia's script

- The fixed prominence (0.15 V) depends on amplifier gain, and
  `MaxPeakWidth = 0.15 s` is about 15× a mouse QRS, so it filters nothing.
- 8–30 Hz leaves out much of the mouse QRS energy, which is mostly above
  30 Hz; check the filtered trace against the raw one.
- DFA α1 on 3 s windows (~30 IBIs, box sizes 4–7) is mostly noise.
- Hard-coded `cd('...\M887\')` and `recDate = '25-Jun-2026'` ignore the
  mouse and date loops.
- `dir` results are loaded without checking for none or several matches.
- `zscore(pupilSize)` returns all NaN if the pupil trace contains NaN.

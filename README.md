# TDT Analysis App

A desktop app for analysing Temporal Discrimination Threshold (TDT) results collected with the TDT Quest Android app (v0.94), written for the Reilly Lab at Trinity College Dublin.

TDT is the shortest interval at which two sequential stimuli are perceived as separate rather than simultaneous, and it is used as an endophenotype marker in adult-onset dystonia. The Quest headset writes every session to a plain-text log; this app parses that log for a given participant, fits a psychometric function to their responses, and reports the threshold with bootstrap confidence intervals.

## What it does

- Finds `Results.txt` automatically on a mounted USB drive (the headset's export), or falls back to a file picker.
- Filters trials by eye (left/right) and paradigm (staircase/random), skipping practice runs and any test flagged `Include test in analysis: NO`.
- Estimates the TDT per trial as the first of three consecutive "different" responses, taking the median per eye.
- Fits a cumulative Gaussian by maximum likelihood (binomial log-likelihood, L-BFGS-B) to recover PSE (mean) and JND (standard deviation).
- Resamples 2000 parametric bootstrap replicates to get 95% confidence intervals on PSE, JND and TDT by the quantile method recommended in Wichmann & Hill (2001).
- Plots the fitted curve over the bootstrap envelope and exports the figure to a timestamped folder.

## Files

- **`main_gui.py`** — Tkinter front end and the entry point. Handles participant ID entry, the trial-selection checkboxes, USB discovery of the results file, and export of the summary figure.
- **`txt_parsing.py`** — Parses the Quest log. Splits it into per-test blocks, matches the participant, applies the eye/paradigm filters, and returns a `TestResults` object of per-trial ISI and response arrays. Also strips the zero-width characters the Quest app appends to participant IDs.
- **`tdt_fitting.py`** — The analysis itself: threshold extraction, the negative log-likelihood and cumulative-Gaussian fit, the bootstrap, and the plotting routines.
- **`analyse_head_movements.py`** — Optional extra: parses quaternion head-rotation data from the log and computes the change in yaw/pitch/roll between the start and end of each trial, to check whether the participant's head drifted during a test.
- **`TDT Analyser.spec`** — PyInstaller spec used to build the standalone Windows executable.
- **`examples/Results.txt`** — A sample Quest export, kept so the parser can be tried without a headset.

## Running it

```bash
python -m venv .venv
.venv\Scripts\activate        # Windows
pip install -r requirements.txt
python main_gui.py
```

Enter the participant ID exactly as it appears in `Results.txt`, tick at least one eye and one paradigm, then **Run Analysis**. **Export Results** writes the plot to `TDT results/<ID> <timestamp>/`.

## Building the standalone executable

To distribute the app to machines without Python:

```bash
pyinstaller --noconfirm "TDT Analyser.spec"
```

The executable is written to `dist/TDT Analyser.exe`. `Calculator.ico` must be present in the project root — it is bundled as data and used as the window and application icon.

## Author

Spencer O'Connell, Reilly Lab, Trinity College Dublin.

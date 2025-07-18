TDT Patient Analysis App

Author: Spencer O’Connell (Reilly Lab)

This application provides a graphical interface to analyze Temporal Discrimination Threshold (TDT) results exported from the TDT Quest Android app. It supports both staircase and random paradigms, left/right eye data selection, and bootstrap-based confidence intervals for PSE, JND, and TDT estimates.

Key features:

Automatic parsing of Results.txt files from USB sticks or manual selection

Flag-based trial filtering (left eye/right eye, staircase/random)

Robust fitting to Gaussian psychometric function

Bootstrap resampling for 95% confidence intervals

Interactive plotting of psychometric curves

Easy export of analysis results and figures


Requirements

Python 3.9+
Virtual environment (recommended)

Dependencies listed in requirements.txt:
numpy
pandas
scipy
matplotlib
scikits.bootstrap
tkinter (standard library)

Usage

Run the GUI
python main_gui.py
In the application window:
Enter Patient ID exactly as it appears in the Results.txt file.
Select eye(s) and paradigm(s) using the checkboxes.
Click Run Analysis to parse data, fit the psychometric function, and display results.
Click Export Results to save summary and plots into a timestamped folder under TDT results/.

Building a Standalone Executable

To distribute the app without requiring Python:
Ensure the .ico file is present in the project root.
From the project folder (with virtual environment activated), run:

pyinstaller \
  --noconfirm \
  --onefile \
  --windowed \
  --name "TDT Analyser" \
  --add-data "Calculator.ico;." \
  main_gui.py

The standalone executable will be in dist/TDT Analyser.exe.

Project Structure

tdt-analyser/
├── main_gui.py         # Main Tkinter application
├── tdt_fitting.py      # Analysis and fitting routines
├── txt_parsing.py      # Functions for extracting test data
├── Calculator.ico      # Application icon for builds
├── requirements.txt    # Python dependencies
├── dist/               # Output folder for PyInstaller builds
├── build/              # Temporary build artifacts
└── README.md           # This documentation



I apologise for the messiness of the files, I am new to GitHub

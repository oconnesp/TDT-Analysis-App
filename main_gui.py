"""""
GUI Main Script

for Reilly Lab TDT Analysis
Author: Spencer O'Connell

For use with the TDT Quest App v0.94

"""""

# Import required libraries
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import tdt_fitting as fit
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scikits.bootstrap as bootci
import os, sys
from datetime import datetime
from txt_parsing import extract_from_txt, TestResults
from scipy.stats import norm

no_bootstraps = 2000 # Number of bootstrap samples for CI estimation

def resource_path(relative_path):
    """
    Get absolute path to resource, works for dev and for PyInstaller.
    Used for loading icons/resources regardless of packaging.
    """
    if getattr(sys, 'frozen', False):
        # PyInstaller unpacks to a temp folder stored in _MEIPASS
        base_path = sys._MEIPASS
    else:
        base_path = os.path.dirname(__file__)
    return os.path.join(base_path, relative_path)

def find_file_on_usb(filename="Results.txt"):
    """
    Attempts to locate the Results.txt file on connected USB drives.
    Returns the absolute path if found, otherwise None.
    """
    from pathlib import Path
    potential_mounts = []

    if os.name == 'nt':  # Windows
        # Check common drive letters for USB devices
        for drive in "DEFGHIJKLMNOPQRSTUVWXYZ":
            potential_mounts.append(f"{drive}:/")
    else:  # macOS/Linux
        # Check typical mount points
        potential_mounts += ["/Volumes", "/media"]

    for mount in potential_mounts:
        root_path = Path(mount)
        if root_path.exists():
            # Recursively search for the file
            for path in root_path.rglob(filename):
                return str(path.resolve())
    return None

def prompt_user_for_file():
    """
    Prompts the user to manually select the Results.txt file using a file dialog.
    Returns the selected file path or None.
    """
    root = tk.Tk()
    root.withdraw()
    return filedialog.askopenfilename(title="Select Results.txt file", filetypes=[("Text Files", "*.txt")])

def build_gui():
    """
    Constructs the main GUI window and its widgets.
    Returns the root Tkinter window.
    """
    root = tk.Tk()
    root.title("TDT Patient Analysis App")
    root.geometry("600x400")
    root.resizable(False, False)
    icon_path = resource_path("Calculator.ico")
    root.iconbitmap(icon_path)

    # Frame for user inputs (Patient ID)
    frm_inputs = ttk.LabelFrame(root, text="Input Patient ID", padding=12)
    frm_inputs.pack(padx=10, pady=10, fill="x")

    ttk.Label(frm_inputs, text="Patient ID:").grid(row=0, column=0, sticky="w")
    patient_id_var = tk.StringVar()
    ttk.Entry(frm_inputs, textvariable=patient_id_var).grid(row=0, column=1, padx=5, pady=5)

    # Frame for controls (buttons and checkboxes)
    frm_controls = ttk.Frame(root, padding=10)
    frm_controls.pack(padx=10, pady=(0, 10), fill="x")

    frm_controls.columnconfigure(0, weight=1)
    frm_controls.columnconfigure(1, weight=1)

    # Actions (buttons) frame
    frm_buttons = ttk.LabelFrame(frm_controls, text="Actions")
    frm_buttons.grid(row=0, column=1, sticky="nsew", padx=(0,5))

    # Checkboxes frame for trial selection
    frm_checks = ttk.LabelFrame(frm_controls, text="Select which trials to include:")
    frm_checks.grid(row=0, column=0, sticky="nsew", padx=(5,0))

    # Create checkboxes for trial options
    buttons = []
    check_vars = []
    for label in ("Left Eye", "Right Eye", "Staircase", "Random", "Include First Trial"):
        var = tk.BooleanVar()
        cb  = tk.Checkbutton(frm_checks, text=label, variable=var)
        cb.pack(anchor="w", pady=2)
        check_vars.append(var)
        buttons.append(cb)

    # Callback for "Run Analysis" button
    def on_run_analysis():
        patient_id = patient_id_var.get().strip()
        # Get checkbox values
        for i in range (5):
            flags = [var.get() for var in check_vars]
        # Ensure at least one eye and one paradigm is selected
        if (flags[0] == False and flags[1] == False):
            messagebox.showerror("Error", "Choose left eye, right eye, or both")
            on_clear()
            return
        if (flags[2] == False and flags[3] == False):
            messagebox.showerror("Error", "Choose staircase, random, or both")
            on_clear()
            return
        # Create timestamped ID string
        ID_timestamp = f"{patient_id}" +" " + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        # Try to find Results.txt automatically, otherwise prompt user
        file_path = find_file_on_usb("Results.txt")
        if not file_path:
            file_path = prompt_user_for_file()
        if not file_path:
            messagebox.showerror("Error", "No Results file found.")
            on_clear()
            return
        try:
            # Parse results from file
            test_results = extract_from_txt(file_path, patient_id, flags)
        except Exception as e:
            messagebox.showerror("File error", f"An error occurred while reading the file:\n{e}")
            return
        # Perform analysis and fitting
        avg_resps, summed_resps, resp_counter, TDT = fit.analyse_TDTs(test_results)
        mu_hat, sigma_hat = fit.Fit_to_Gaussian(test_results.all_ISIs, resp_counter, summed_resps, TDT)

        root.patient_id = patient_id
        # Bootstrap for confidence intervals
        bootstrapped_results = fit.bootstrap(mu_hat, sigma_hat, no_bootstraps, test_results)

        # Compute 95% confidence intervals using quantile method
        boot = bootstrapped_results[["PSE","JND","TDT"]].to_numpy()  # shape (B,3)
        ci_low  = np.array([np.percentile(boot[:, j], 100 * 0.025) for j in range(3)])
        ci_high = np.array([np.percentile(boot[:, j], 100 * 0.975) for j in range(3)])
        ci_array = np.column_stack((ci_low, ci_high)) # shape (3,2)

        # Build summary string for results
        type_string = ""
        for i, label in enumerate(["Left Eye", "Right Eye", "Staircase", "Random", "First Trial Included"]):
            if flags[i] == True:
                type_string += label + " "
        summary = type_string + "\n" + (
            f"TDT: {TDT:.2f} ms (95% CI: {ci_array[2,0]:.2f} – {ci_array[2,1]:.2f})\n"
            f"PSE: {mu_hat:.2f} ms (95% CI: {ci_array[0,0]:.2f} – {ci_array[0,1]:.2f})\n"
            f"JND: {sigma_hat:.2f} ms (95% CI: {ci_array[1,0]:.2f} – {ci_array[1,1]:.2f})"
        )
        root.results_var.set(summary)     # Show results in GUI
        root.update_idletasks()           # Force redraw

        # Plot results with bootstraps
        plt.figure()
        fit.plot_with_bootstraps(bootstrapped_results, mu_hat, sigma_hat, test_results.all_ISIs, avg_resps)
        plt.text(
            0, 0.8, summary,
            fontsize=9,
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray')
        )
        plt.title(ID_timestamp)
        fig = plt.gcf()
        plt.show(block=False)

        # Store analysis results in root for later export
        root.analysis_results = {
            "patient_id": patient_id,
            "no_trials": test_results.no_trials,
            "TDT": TDT,
            "PSE": mu_hat,
            "JND": sigma_hat,
            "figure": fig,
            "CI_array": ci_array,
            "ID_timestamp": ID_timestamp
        }

    # Callback for "Export Results" button
    def on_export():
        results = getattr(root, "analysis_results", None)
        if results is None:
            messagebox.showerror("Error", "No results to export. Please run analysis first.")
            return

        # Create base folder for exports
        base_dir = "TDT results"
        os.makedirs(base_dir, exist_ok=True)

        # Create timestamped subfolder for this export
        export_folder = os.path.join(base_dir, results["ID_timestamp"])
        os.makedirs(export_folder)

        # Save figure as PNG
        fig = results["figure"]
        fig_path = os.path.join(export_folder, "tdt_plot.png")
        fig.savefig(fig_path)
        location = os.path.abspath("TDT results")
        messagebox.showinfo("Export Successful", f"Results exported to:\n{location}")

    # Callback for "Clear Results" button
    def on_clear():
        # Clear results text
        root.results_var.set("Results will appear here.")
        # Clear input fields
        root.patient_id_var.set("")
        # Deselect all checkboxes
        for button in buttons:
            button.deselect()
        # Remove stored results if present
        if hasattr(root, "analysis_results"):
            del root.analysis_results

    # Add buttons to Actions frame
    ttk.Button(frm_buttons, text="Run Analysis", command=on_run_analysis).pack(fill="x", pady=5)
    ttk.Button(frm_buttons, text="Export Results", command=on_export).pack(fill="x", pady=5)
    ttk.Button(frm_buttons, text="Clear Results", command=on_clear).pack(fill="x", pady=5)

    # Frame for results display
    frm_output = ttk.LabelFrame(root, text="Results", padding=12)
    frm_output.pack(padx=10, pady=10, fill="both", expand=True)

    # Variable to hold results summary text
    results_var = tk.StringVar(value="Results will appear here.")
    ttk.Label(frm_output, textvariable=results_var, justify="left").pack(anchor="w")

    # Store state variables in root for access in callbacks
    root.results_var = results_var
    root.patient_id_var = patient_id_var

    return root

if __name__ == "__main__":
    # Run the GUI application
    app = build_gui()
    app.mainloop()

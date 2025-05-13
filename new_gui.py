"""""
GUI Main Script

for Reilly Lab TDT Analysis
Author: Spencer O'Connell

"""""


import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import tdt_fitting as fit
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scikits.bootstrap as bootci
import os
from datetime import datetime
from txt_parsing import extract_from_txt
from scipy.stats import norm
no_bootstraps = 2000 #hard-coded


def find_file_on_usb(filename="Results.txt"):#Find Results file
    from pathlib import Path
    potential_mounts = []

    if os.name == 'nt':  # Windows
        for drive in "DEFGHIJKLMNOPQRSTUVWXYZ":
            potential_mounts.append(f"{drive}:/") #potential path addresses
    else:  # macOS/Linux
        potential_mounts += ["/Volumes", "/media"]

    for mount in potential_mounts:
        root_path = Path(mount)
        if root_path.exists():
            for path in root_path.rglob(filename):
                return str(path.resolve())
    return None
#otherwise user can open file manually
def prompt_user_for_file():
    root = tk.Tk()
    root.withdraw()
    return filedialog.askopenfilename(title="Select Results.txt file", filetypes=[("Text Files", "*.txt")])


def build_gui():
    root = tk.Tk()
    root.title("TDT Patient Analysis App")
    root.geometry("600x400")
    root.resizable(False, False)

    # Frame for user inputs
    frm_inputs = ttk.LabelFrame(root, text="Input Parameters", padding=12)
    frm_inputs.pack(padx=10, pady=10, fill="x")

    ttk.Label(frm_inputs, text="Patient ID:").grid(row=0, column=0, sticky="w")
    patient_id_var = tk.StringVar()
    ttk.Entry(frm_inputs, textvariable=patient_id_var).grid(row=0, column=1, padx=5, pady=5)

    ttk.Label(frm_inputs, text="Number of Trials:").grid(row=1, column=0, sticky="w")
    trials_var = tk.StringVar()
    ttk.Entry(frm_inputs, textvariable=trials_var).grid(row=1, column=1, padx=5, pady=5)

    # Frame for buttons
    frm_buttons = ttk.Frame(root, padding=10)
    frm_buttons.pack(padx=10, pady=(0, 10), fill="x")

    def on_run_analysis():
        patient_id = patient_id_var.get()
        try:
            no_trials = int(trials_var.get())
        except ValueError:
            messagebox.showerror("Input error", "Please enter an integer number of trials.")
            return
        ID_timestamp = f"{patient_id}" +" " + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        file_path = find_file_on_usb ("Results.txt")
        if not file_path:
            file_path = prompt_user_for_file()
        if not file_path:
            messagebox.showerror("No Results file found.")
            return
        try:
            isi_list, resp_list, all_isis = extract_from_txt(file_path, no_trials)
        except Exception as e:
            messagebox.showerror("File error", f"An error occurred while reading the file:\n{e}")
            return
        avg_resps, summed_resps, resp_counter, TDT = fit.analyse_TDTs(isi_list, resp_list, all_isis, no_trials)
        mu_hat, sigma_hat = fit.Fit_to_Gaussian(all_isis, resp_counter, summed_resps, TDT)

        root.patient_id = patient_id  
        bootstrapped_results = fit.bootstrap(mu_hat, sigma_hat, all_isis, no_bootstraps, no_trials)

        """plt.figure()
        plt.hist (bootstrapped_results ["PSE"])
        plt.title("PSE")

        plt.figure()
        plt.hist (bootstrapped_results ["JND"])
        plt.title("JND")
        plt.figure()
        plt.hist (bootstrapped_results ["TDT"])
        plt.title("TDT") """

        ###95%CI calculations using scikits.bootstrap BCa method, Bias correction with acceleration
        ###Recommended by Wichmann and Hill (2001)
        z0_hats = fit.calc_z0_hats([mu_hat, sigma_hat, TDT], bootstrapped_results, no_bootstraps)
        a_hats = fit.calc_a_hats ([mu_hat, sigma_hat, TDT], isi_list, resp_list, all_isis, no_trials)
        z_low, z_high = norm.ppf(0.05/2), norm.ppf(1 - 0.05/2)
        pct_low = norm.cdf(z0_hats+ (z0_hats + z_low) / (1 - a_hats*(z0_hats + z_low)))
        pct_high = norm.cdf(z0_hats+ (z0_hats + z_high) / (1 - a_hats*(z0_hats + z_high)))
        print  ([z_low, z_high])
        print ([pct_low, pct_high])
        boot = bootstrapped_results[["PSE","JND","TDT"]].to_numpy()  # shape (B,3)

        ci_mu = bootci.ci(data=bootstrapped_results["PSE"].values, statfunction=np.mean, method='bca')
        ci_sigma = bootci.ci(data=bootstrapped_results["JND"].values, statfunction=np.mean, method='bca')
        ci_TDT = bootci.ci(data=bootstrapped_results["TDT"].values, statfunction=np.mean, method='bca')


        ci_low  = np.array([np.percentile(boot[:, j], 100 * 0.025)
                    for j in range(3)])
        ci_high = np.array([np.percentile(boot[:, j], 100 * 0.975)
                    for j in range(3)])

        # stack into a (3×2) array: rows = [PSE, JND, TDT], cols = [low, high]
        ci_array = np.column_stack((ci_low, ci_high))
        summary = (
            f"TDT: {TDT:.2f} ms (95% CI: {ci_array[2,0]:.2f} – {ci_array[2,1]:.2f})\n"
            f"PSE: {mu_hat:.2f} ms (95% CI: {ci_array[0,0]:.2f} – {ci_array[0,1]:.2f})\n"
            f"JND: {sigma_hat:.2f} ms (95% CI: {ci_array[1,0]:.2f} – {ci_array[1,1]:.2f})"
        )
        root.results_var.set(summary)     # show results immediately
        root.update_idletasks()           # force a redraw
        plt.figure()
        fit.plot_with_bootstraps (bootstrapped_results, mu_hat, sigma_hat, all_isis, avg_resps)
        plt.text(
        60, 0.7, summary,
        fontsize=9,
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray')
)
        plt.title (ID_timestamp)
        fig = plt.gcf()
        plt.show(block=False)
        root.analysis_results = { #update root attributes
            "patient_id": patient_id,
            "no_trials": no_trials,
            "TDT": TDT,
            "PSE": mu_hat,
            "JND": sigma_hat,
            "figure" : fig,
            "CI_array" : ci_array,
            "ID_timestamp" : ID_timestamp
        }


    def on_export():
        results = getattr(root, "analysis_results", None)
        if results is None:
            messagebox.showerror("Error", "No results to export. Please run analysis first.")
            return

        # Step 1: Create base folder
        base_dir = "TDT results"
        os.makedirs(base_dir, exist_ok=True)

        # Step 2: Create timestamped subfolder
        export_folder = os.path.join(base_dir, results["ID_timestamp"])
        os.makedirs(export_folder)

        # Step 3: Save figure
        fig = results["figure"]
        fig_path = os.path.join(export_folder, "tdt_plot.png")
        fig.savefig(fig_path)
        location = os.path.abspath("TDT results")
        # Step 4: Save confidence interval data
        messagebox.showinfo("Export Successful", f"Results exported to:\n{location}")


    def on_clear():
        # Clear results text
        root.results_var.set("Results will appear here.")
        
        # Clear input fields
        root.patient_id_var.set("")
        root.trials_var.set("")
        
        # Clear stored results if they exist
        if hasattr(root, "analysis_results"):
            del root.analysis_results

    ttk.Button(frm_buttons, text="Run Analysis", command=on_run_analysis).pack(fill="x", pady=5)
    ttk.Button(frm_buttons, text="Export Results", command=on_export).pack(fill="x", pady=5)
    ttk.Button(frm_buttons, text="Clear Results", command=on_clear).pack(fill="x", pady=5)

    # Frame for results display
    frm_output = ttk.LabelFrame(root, text="Results", padding=12)
    frm_output.pack(padx=10, pady=10, fill="both", expand=True)

    results_var = tk.StringVar(value="Results will appear here.")
    ttk.Label(frm_output, textvariable=results_var, justify="left").pack(anchor="w")

    # Store state
    root.results_var = results_var
    root.patient_id_var = patient_id_var
    root.trials_var = trials_var

    return root


if __name__ == "__main__":
    app = build_gui()
    app.mainloop()

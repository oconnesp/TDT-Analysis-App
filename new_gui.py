import tkinter as tk
from tkinter import ttk, messagebox
import tdt_fitting as fit
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scikits.bootstrap as bootci
import os
from datetime import datetime
csv_path = r"C:\Users\OCONNESP\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
no_bootstraps = 2000


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
        isi_list, resp_list, all_isis = fit.extract_trials(csv_path, no_trials)
        avg_resps, summed_resps, resp_counter, TDT = fit.analyse_TDTs(isi_list, resp_list, all_isis, no_trials)
        mu_hat, sigma_hat = fit.Fit_to_Gaussian(all_isis, resp_counter, summed_resps, TDT)

        root.patient_id = patient_id  
        bootstrapped_results = fit.bootstrap(mu_hat, sigma_hat, all_isis, no_bootstraps, no_trials)
        ci_mu = bootci.ci(data=bootstrapped_results["PSE"].values, statfunction=np.mean, method='bca')
        ci_sigma = bootci.ci(data=bootstrapped_results["JND"].values, statfunction=np.mean, method='bca')
        ci_TDT = bootci.ci(data=bootstrapped_results["TDT"].values, statfunction=np.mean, method='bca')
        ci_array = np.array([ci_mu, ci_sigma, ci_TDT])
        summary = (
            f"TDT: {TDT:.2f} ms (95% CI: {ci_TDT[0]:.2f} – {ci_TDT[1]:.2f})\n"
            f"PSE: {mu_hat:.2f} ms (95% CI: {ci_mu[0]:.2f} – {ci_mu[1]:.2f})\n"
            f"JND: {sigma_hat:.2f} ms (95% CI: {ci_sigma[0]:.2f} – {ci_sigma[1]:.2f})"
        )
        root.results_var.set(summary)     # show results immediately
        print("Summary set")

        root.update_idletasks()           # force a redraw if you like
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
        root.analysis_results = {
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

import tdt_fitting as fit
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import scikits.bootstrap as bootci
import tkinter as tk
# "C:\Users\OCONNESP\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
#"C:\Users\spenc\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
#Test 
csv_path = r"C:\Users\OCONNESP\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
no_trials = 4 #TODO find no_trials
no_bootstraps = 2000
# Extract files and find the average responses
isi_list, resp_list, all_isis = fit.extract_trials(csv_path, no_trials)
avg_resps, summed_resps, resp_counter, TDT = fit.analyse_TDTs(isi_list, resp_list, all_isis, no_trials)

# Print the result
print (f"Subject's TDT: {TDT:.2f}")
#print("ISI (ms):", all_isis)
#print("Mean Response:", avg_resps)
## get the fitted cumulative gaussian curve 

mu_hat, sigma_hat = fit.Fit_to_Gaussian (all_isis, resp_counter, summed_resps, TDT)

plt.figure()
fit.plot_curve(all_isis, avg_resps, sigma_hat, mu_hat)
plt.show()


bootstrapped_results = fit.bootstrap(mu_hat, sigma_hat, all_isis, no_bootstraps, no_trials)
column_means = bootstrapped_results.mean()
print (column_means)
plt.figure()
fit.plot_with_bootstraps (bootstrapped_results, mu_hat, sigma_hat, all_isis, avg_resps)
plt.show()

ci_mu = bootci.ci(data=bootstrapped_results["PSE"].values, statfunction=np.mean, method='bca')
ci_sigma = bootci.ci(data=bootstrapped_results["JND"].values, statfunction=np.mean, method='bca')
ci_TDT = bootci.ci(data=bootstrapped_results["TDT"].values, statfunction=np.mean, method='bca')

# Combine into an array
ci_array = np.array([ci_mu, ci_sigma, ci_TDT])

print("95% confidence intervals (rows: PSE, JND, TDT):\n", np.round(ci_array, 2))
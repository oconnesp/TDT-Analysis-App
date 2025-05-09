import tdt_fitting as fit
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# "C:\Users\OCONNESP\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
#"C:\Users\spenc\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
#Test 
csv_path = r"C:\Users\OCONNESP\OneDrive - Trinity College Dublin\College\College\Year 3\Reilly Lab\TDT Analysis App\EXAMPLE_24_06_2024.csv"
no_trials = 4 #TODO find no_trials
# Extract files and find the average responses
all_isis, avg_resps, summed_resps, resp_counter, TDT = fit.extract_trials(csv_path, no_trials)

# Print the result
print("ISI (ms):", all_isis)
print("Mean Response:", avg_resps)
## get the fitted cumulative gaussian curve 

mu_hat, sigma_hat = fit.Fit_to_Gaussian (all_isis, resp_counter, summed_resps, TDT)


fit.plot_curve(all_isis, avg_resps, sigma_hat, mu_hat)



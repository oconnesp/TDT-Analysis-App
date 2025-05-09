from tdt_fitting import extract_trials, Fit_to_Gaussian

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
ISIs, avg_resps = extract_trials(csv_path, no_trials)

# Print the result
print("ISI (ms):", ISIs)
print("Mean Response:", avg_resps)
## get the fitted cumulative gaussian curve 
mu_fitted, sigma_fitted = Fit_to_Gaussian(ISIs, avg_resps)





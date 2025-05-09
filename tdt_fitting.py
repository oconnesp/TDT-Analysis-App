""""
tdt_fitting.py

extracts trials from csv, fits to a cumulative gaussian
returns statistics

"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm


def extract_trials(csv_path: str, no_trials: int):
    df = pd.read_csv(csv_path)

    isi_list:  list[np.ndarray] = []
    resp_list: list[np.ndarray] = []

    for t in range(no_trials):
        isi_col  = df.columns[t * 2]       # 0, 2, 4, …
        resp_col = df.columns[t * 2 + 1]   # 1, 3, 5, …

        # Keep only numeric data, drop NaNs (headers, “TDT”, etc.)
        pair = (
            df[[isi_col, resp_col]]                     # take the two columns together
            .apply(pd.to_numeric, errors="coerce")    # make them numeric
            .dropna()                                 # drop rows where **either** is NaN
        )
        isi  = pair[isi_col].to_numpy(float) * 1000 #ms
        resp = pair[resp_col].round().astype(int).to_numpy()


        isi_list.append(isi)
        resp_list.append(resp)

    all_isis = np.unique(np.concatenate(isi_list))

#all isis with no duplciates
    length = len(all_isis)
    summed_resps = np.zeros(length) # an array of zeros for summed responses
    resp_counter = np.zeros(length) #counter for each index
    avg_resps = np.zeros (length)
    threshold_values = np.zeros (no_trials)
    #find TDT 
    for n in range (no_trials):
        for i in range (len(resp_list[n]) - 2): 
            if ( (resp_list[n][i] == 1) and (resp_list[n][i+1] == 1) and (resp_list[n][i+2] == 1)  ):
                threshold_values[n] = isi_list[n][i]
                break
    
    TDT = np.median(threshold_values)
    print (f"TDT = {TDT:.2f}")
        ### Take each response and add it into an array with indeces corresponding to each ISI
    for n in range (no_trials): 
        for i in range (len(isi_list[n])):
            idx = np.searchsorted(all_isis, isi_list[n][i])
            summed_resps [idx] += resp_list[n][i]
            resp_counter[idx] += 1
        
        ## Average them
    for j in range (length):
        avg_resps [j] = summed_resps[j]/resp_counter[j]
        

    return all_isis, avg_resps


### now have all_isis and their avg_resps


### next, fit to a cumulative Gaussian ###

#Initial Guesses
def Fit_to_Gaussian (all_isis,avg_resps):
    mu_guess = np.median (all_isis)
    sigma_guess = (np.max(all_isis) - np.min(all_isis))/4
    p0 = [mu_guess, sigma_guess]

    #bounds
    bounds = ([0, 1e-6], [np.max(all_isis), np.inf])

    popt, pcov = curve_fit (norm.cdf, all_isis, avg_resps, p0 =p0, bounds = bounds)

    mu_fitted, sigma_fitted = popt

    x_fitted = np.linspace (0, np.max(all_isis), 300)
    y_fitted = norm.cdf (x_fitted,mu_fitted,sigma_fitted)

    plt.figure()
    plt.scatter (all_isis, avg_resps, label = 'data', marker = 'o', color = 'blue')
    plt.plot (x_fitted, y_fitted, label = 'fit', color = 'red', linestyle = '--', lw = 2.0)
    plt.text (70, 0.7, f"PSE: {mu_fitted :.2f}ms\n JND: {sigma_fitted :.2f}ms")
    plt.xlabel ("Time (ms)")
    plt.ylabel ("Proportion of \"Different\" Responses")
    plt.show()
    return (popt)

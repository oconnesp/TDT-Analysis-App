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
from scipy.optimize import minimize
from scipy.stats import norm


def extract_trials(csv_path: str, no_trials: int):
    df = pd.read_csv(csv_path)

    isi_list:  list[np.ndarray] = []
    resp_list: list[np.ndarray] = []

    for t in range(no_trials):
        isi_col  = df.columns[t * 2]       # 0, 2, 4, …
        resp_col = df.columns[t * 2 + 1]   # 1, 3, 5, …

        # Keep only numeric data, drop NaNs (headers, “TDT”, etc.)
    # take the two columns, but skip row 0
        pair = (
        df.loc[1:, [isi_col, resp_col]]     # ← here!
          .apply(pd.to_numeric, errors="coerce")
        )

        # find first NaN and slice up to it
        null_mask = pair.isna().any(axis=1).to_numpy()
        if null_mask.any():
            stop_idx = int(null_mask.argmax())
            pair = pair.iloc[:stop_idx]

    # convert & scale
        isi  = pair[isi_col].to_numpy(float) * 1000
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
            summed_resps [idx] += resp_list[n][i] #k_i
            resp_counter[idx] += 1 #n_i
        
        ## Average them to find avg responses
    for j in range (length):
        avg_resps [j] = summed_resps[j]/resp_counter[j]
        
        #Create an array for ISIs, an array for 

    return all_isis, avg_resps, summed_resps, resp_counter, TDT


### now have all_isis and their avg_resps


### next, fit to a cumulative Gaussian ###


##define negative log likelihood function

def neg_log_likelihood(params, k, n, ISIs):#k is #different's per ISI, n is no trials per ISI
    mu, sigma = params
    if sigma <=0:
        return np.inf
    p = norm.cdf (ISIs, mu, sigma)
    eps = 1e-9
    p = np.clip (p, eps, 1-eps) #prevent log0

    #binomial log-likelihood
    ll = k*np.log (p) + (n-k)*np.log (1-p)
    return -np.sum(ll)


#Initial Guesses
def Fit_to_Gaussian (ISIs, n, k, TDT):
    mu_guess = TDT
    sigma_guess = 0.1*np.ptp(ISIs)

    #bounds
    bounds = ([0, np.max(ISIs)], [1e-6, np.inf])

    result = minimize(
        neg_log_likelihood,
        x0=[mu_guess, sigma_guess],
        args=(k, n, ISIs),
        bounds=bounds,
        method="L-BFGS-B"
    )
    mu_hat, sigma_hat  = result.x
    return mu_hat, sigma_hat


def plot_curve(all_isis, avg_resps, sigma_fitted, mu_fitted):
    x_fitted = np.linspace (0, np.max(all_isis), 300)
    y_fitted = norm.cdf (x_fitted,loc=mu_fitted, scale =sigma_fitted)

    plt.figure()
    plt.scatter (all_isis, avg_resps, label = 'data', marker = 'o', color = 'blue')
    plt.plot (x_fitted, y_fitted, label = 'fit', color = 'red', linestyle = '--', lw = 2.0)
    plt.text (70, 0.7, f"PSE: {mu_fitted :.2f}ms\n JND: {sigma_fitted :.2f}ms")
    plt.xlabel ("Time (ms)")
    plt.ylabel ("Proportion of \"Different\" Responses")
    plt.show()
    return



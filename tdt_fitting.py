""""
tdt_fitting.py

extracts trials from csv, fits to a cumulative gaussian
returns statistics

"""
import sys
from pathlib import Path
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm
from scikits.bootstrap import bootstrap as bootci
from datetime import datetime
from txt_parsing import TestResults

def analyse_TDTs (test_results : TestResults):
    isi_list, resp_list, all_isis, no_trials = test_results.ISIs, test_results.responses, test_results.all_ISIs, test_results.no_trials
    length = len(all_isis)
    summed_resps = np.zeros(length) # an array of zeros for summed responses
    resp_counter = np.zeros(length) #counter for each index
    avg_resps = np.zeros (length)
    left_trials = 0
    for i in range (test_results.no_trials):
        if test_results.left_eye[i] == True:#how many trials use the left eye
            left_trials += 1
    threshold_values_L = np.zeros (left_trials)
    threshold_values_R = np.zeros (test_results.no_trials - left_trials) 
    L_TDT_counter = 0
    R_TDT_counter = 0
    #find TDT 
    for n in range (no_trials):
        for i in range (len(resp_list[n]) - 2): #iterate through responses
            if ( (resp_list[n][i] == 1) and (resp_list[n][i+1] == 1) and (resp_list[n][i+2] == 1)  ): #if you get 3 "different's" in a row
                if test_results.left_eye[n] == True:
                    threshold_values_L [L_TDT_counter] = isi_list[n][i] #the first of the sequence of 3 "different's 
                    L_TDT_counter += 1
                else:
                    threshold_values_R[R_TDT_counter] = isi_list[n][i] #the first of the sequence of 3 "different's  
                    R_TDT_counter += 1
                break

    #Don't take the mean of the L and R values if you're only testing one eye
    if (len (threshold_values_L) == 0):
        TDT = np.median(threshold_values_R)

    elif (len(threshold_values_R) == 0):
        TDT = np.median(threshold_values_L)



    else:
        TDT = (np.median(threshold_values_L) + np.median(threshold_values_R))/2

        ### Take each response and add it into an array with indeces corresponding to each ISI
    for n in range (no_trials): 
        for i in range (len(isi_list[n])):
            idx = np.searchsorted(all_isis, isi_list[n][i])
            summed_resps [idx] += resp_list[n][i] #k_i
            resp_counter[idx] += 1 #n_i
        
        ## Average them to find avg responses
    for j in range (length):
        avg_resps [j] = summed_resps[j]/resp_counter[j]
        
        #Create an array for ISIs, an array for summed responses, an array for the number of responses for each ISI (not uniform in case of random) 

    return avg_resps, summed_resps, resp_counter, TDT


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
def Fit_to_Gaussian (ISIs, n, k, TDT): #n is trials per ISI, k is summed resps
    mu_guess = TDT
    sigma_guess = 0.1*np.ptp(ISIs)

    #bounds
    bounds = ([0, np.max(ISIs)], [1e-6, np.inf])# mean in [0, max ISI], stdev in (0,inf)

    result = minimize( #minimise the negative log-likelihood,
        neg_log_likelihood,
        x0=[mu_guess, sigma_guess],
        args=(k, n, ISIs),
        bounds=bounds,
        method="L-BFGS-B"
    )
    mu_hat, sigma_hat  = result.x
    return mu_hat, sigma_hat #parameters for cumulative gaussian fit


def plot_curve(all_isis, avg_resps, sigma_fitted, mu_fitted):
    x_fitted = np.linspace (0, np.max(all_isis), 300)
    y_fitted = norm.cdf (x_fitted,loc=mu_fitted, scale =sigma_fitted)
    plt.plot (x_fitted, y_fitted, label = 'fit', color = 'black', lw = 2.0)
    plt.scatter (all_isis, avg_resps, label = 'data', marker = 'o', color = 'blue', zorder = 3)
    plt.xlabel ("Time (ms)")
    plt.ylabel ("Proportion of \"Different\" Responses")
    return

def plot_one_bootstrap(mu_hat, sigma_hat, ISIs):
    x_fitted = np.linspace (0, np.max(ISIs), 300)
    y_fitted = norm.cdf (x_fitted,loc=mu_hat, scale =sigma_hat)
    plt.plot(x_fitted, y_fitted, color='lightgrey', linewidth=1)
    return

def bootstrap (mu, sigma, no_bootstraps, test_results : TestResults):
    no_trials = test_results.no_trials
    ISIs = test_results.all_ISIs
    ISI_list = [ISIs.copy() for _ in range(no_trials)] #TODO does this work for random trials?
    bootstrap_results = TestResults("Bootstrap", test_results.staircase, test_results.left_eye, [], test_results.ISIs)
    no_ISIs = len(ISIs)

    resp_list: list[np.ndarray] = []
    bootstrap_rows = []
    p = norm.cdf(ISIs, loc=mu, scale=sigma)#vectorised

    for i in range (no_bootstraps): # do this 2000 times, once per bootstrap
        resp_list: list[np.ndarray] = [] 

        for m in range (no_trials):#generate no_trials random trials, do this 8000 times (once per trial per bootstrap)
        #to find TDT, need to do each trial and assess

            resp_list.append (np.random.rand (no_ISIs) < p )#random Bernouilli trials using probabilities from the Gaussian fit


        avg_resps, summed_resps, resp_counter, TDT = analyse_TDTs (TestResults("Bootstrap", test_results.staircase, test_results.left_eye, resp_list , test_results.ISIs))#Fit a Gaussian to this bootstrapped trial
        mu_sim, sigma_sim = Fit_to_Gaussian (ISIs, resp_counter, summed_resps, TDT)#save the parameters
        bootstrap_rows.append({"PSE": mu_sim, "JND": sigma_sim, "TDT": TDT})

    return pd.DataFrame(bootstrap_rows)#return an array of each bootstrap's parameters


def plot_with_bootstraps (bootstrapped_results, mu_hat, sigma_hat, ISIs, avg_resps):#returns CI intervals
    for i in range (len(bootstrapped_results)):
        plot_one_bootstrap (bootstrapped_results.iloc[i]["PSE"], bootstrapped_results.iloc[i]["JND"], ISIs)
    plot_curve (ISIs, avg_resps, sigma_hat, mu_hat)
    return




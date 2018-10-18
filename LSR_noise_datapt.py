# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 21:16:26 2018

@author: brainfreak97
"""

import numpy as np;
import matplotlib.pyplot as plt;
import random as rand;

#%% setting up

# let's start with defining the variables
coeff_n = 4 ; # number of coefficients

timepoint = 90;
datapoint = 30;

# global variable 
Fluo_exp = np.array([0.0]*timepoint);
mRNA_exp = np.array([0.0]*timepoint);

time_list = np.array([3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90]);
#Fluo_msr = np.array([])  place observation data here, and dont run simulation

#%% simulation part
coefficients = np.array([5.5, 0.1, 8.5, 0.25]);

def positive_control_temp(coeff_temp):
    for t in range (timepoint-1):
        d_mRNA = coeff_temp[0] - coeff_temp[1]*mRNA_exp[t];
        d_Fluo = coeff_temp[2]*mRNA_exp[t] - coeff_temp[3]*Fluo_exp[t];
        mRNA_exp[t+1] = mRNA_exp[t] + d_mRNA + rand.gauss(0, 0.5);
        Fluo_exp[t+1] = Fluo_exp[t] + d_Fluo;
        if (mRNA_exp[t+1] ==0):
            mRNA_exp[t+1] = 0;
        if (Fluo_exp[t+1] ==0):
            Fluo_exp[t+1] = 0;
        
positive_control_temp(coefficients);

time_list = np.array([0]*datapoint);
Fluo_msr = np.array([0]);

for i in range (1, datapoint): #we can specify our own time_list;
    time_list[i] = 3*i;
    Fluo_msr = np.append(Fluo_msr, Fluo_exp[time_list[i]]);
    
fig, ax = plt.subplots()
ax.plot(mRNA_exp, label = 'experimental mRNA');
ax.plot(Fluo_exp, label = 'experimental Fluo');
legend = ax.legend(loc='lower right', shadow=True, fontsize='x-large');
plt.savefig('plot1.png');

# everything looks good!

#%%
# now for the regression
coefficients_guess = np.array([1.5, 0.2, 1.0, 0.15]);
#step = 0.0005;

Fluo_th = np.array([0.0]*timepoint);
mRNA_th = np.array([0.0]*timepoint);
Fluo_gs = np.array([0.0]*timepoint);

def compute_dist(fluo_th, fluo_ex):
    sum0 = 0;
    for t in range (datapoint):
        tt = time_list[t];
        sum0 += (fluo_th[tt] - fluo_ex[t])**2 / datapoint;
    return sum0;
    
def positive_control(coeff_temp):
    for t in range (timepoint-1):
        d_mRNA = coeff_temp[0] - coeff_temp[1]*mRNA_th[t];
        d_Fluo = coeff_temp[2]*mRNA_th[t] - coeff_temp[3]*Fluo_th[t];
        mRNA_th[t+1] = mRNA_th[t] + d_mRNA;
        Fluo_th[t+1] = Fluo_th[t] + d_Fluo;
        if (mRNA_exp[t+1] ==0):
            mRNA_exp[t+1] = 0;
        if (Fluo_exp[t+1] ==0):
            Fluo_exp[t+1] = 0;

# the alternative one    
# use alpha = 0.01~0.001 for best result
def run_LSR(coeff_temp, alpha, error, printq):
    flag = 1;
    cntr = 1;
    #step = 0.000001;
    while(flag == 1):
        # normal -------
        positive_control(coeff_temp);
        sum00 = compute_dist(Fluo_th, Fluo_msr);
        # +1 -----------
        positive_control(coeff_temp + np.array([alpha,0,0,0]));
        sum01 = compute_dist(Fluo_th, Fluo_msr);
        positive_control(coeff_temp + np.array([0,alpha,0,0]));
        sum02 = compute_dist(Fluo_th, Fluo_msr);
        positive_control(coeff_temp + np.array([0,0,alpha,0]));
        sum03 = compute_dist(Fluo_th, Fluo_msr);
        positive_control(coeff_temp + np.array([0,0,0,alpha]));
        sum04 = compute_dist(Fluo_th, Fluo_msr);
        # -1 -----------
        positive_control(coeff_temp - np.array([alpha,0,0,0]));
        sum05 = compute_dist(Fluo_th, Fluo_msr);
        positive_control(coeff_temp - np.array([0,alpha,0,0]));
        sum06 = compute_dist(Fluo_th, Fluo_msr);
        positive_control(coeff_temp - np.array([0,0,alpha,0]));
        sum07 = compute_dist(Fluo_th, Fluo_msr);
        positive_control(coeff_temp - np.array([0,0,0,alpha]));
        sum08 = compute_dist(Fluo_th, Fluo_msr);
        
        # now to force --------
        
        force = np.array([0.0]*4);
        force[0] = (sum01 - sum05) / alpha / 2;
        force[1] = (sum02 - sum06) / alpha / 2;
        force[2] = (sum03 - sum07) / alpha / 2;
        force[3] = (sum04 - sum08) / alpha / 2;
        
        step = 0.1*alpha/np.sqrt(np.dot(force,force));
        
        #check step to prevent negative values;
       
        coeff_temp2 = coeff_temp - force*step;
        positive_control(coeff_temp2);
        sum10 = compute_dist(Fluo_th, Fluo_msr);        
        
        #increasing step if way too small?
        
        coeff_temp = coeff_temp2;
        
        if(printq == 1):
            print("sum00 = ",sum00, " step_ord = ",int(np.log10(step)));
        
        if sum00 < error**2 and sum10 > sum00:
            flag = 0;
            #print("found!");

    return coeff_temp;
        
#%%
coefficients_guess = np.array([10, 0.05, 12.5, 0.25]);

positive_control(coefficients_guess);
for t in range(timepoint):
    Fluo_gs[t] = Fluo_th[t];

# Please use a different coefficient if it is stuck
coeff_results1 = run_LSR(coefficients_guess, 0.01, 300, 0);
coeff_results2 = run_LSR(coeff_results1, 0.005, 100, 1);
coeff_results1 = run_LSR(coeff_results2, 0.001, 40, 1);

print("Monte-Carlo complete! true value = ", coefficients, " values obtained = ", coeff_results1);

positive_control(coeff_results1);

fig, ax = plt.subplots()
ax.plot(Fluo_exp, label = 'experimental Fluo');
#ax.plot(Fluo_msr, label = 'measured Fluo');
ax.plot(Fluo_gs, label = '1st Guess Fluo');
ax.plot(Fluo_th, label = 'theoretical Fluo');
legend = ax.legend(loc='lower right', shadow=True, fontsize='x-large');
plt.savefig('plot1.png');

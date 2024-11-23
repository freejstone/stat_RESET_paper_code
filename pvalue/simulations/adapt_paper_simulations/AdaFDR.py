#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 20:45:28 2024

@author: jfre0619
"""

import adafdr.method as md
import numpy as np
import pandas as pd
import os
import random

random.seed(23052024)
np.random.seed(23052024)

path = 'data1'
alphas = [0.01, 0.05, 0.1, 0.2]
files = os.listdir(path)
files = [f for f in files if 'sim' in f]
count = 0
for file in files:
    print(file)
    df = pd.read_csv(path + '/' + file)
    pvals = np.array(df.pvals)
    x = df.loc[:, df.columns.str.contains('x')].to_numpy()
    h = np.array(df.H0)
    FDP_all = []
    power_all = []
    for alpha in alphas:    
        res = md.adafdr_test(pvals, x, alpha=alpha, fast_mode=True, random_state = 21072024 + count)
        t = res['threshold']
        D = np.sum(pvals<=t)
        FD = np.sum((pvals<=t)&(~h))
        power = np.sum((pvals<=t)&(h))/max(sum(h), 1)
        FDP = FD/max(D, 1)
        #print('# AdaFDR successfully finished!')
        #print('# D=%d, FD=%d, FDP=%0.3f'%(D, FD, FD/D))
        FDP_all.append(FDP)
        power_all.append(power)
    count += 1        
    df_out = pd.DataFrame(zip(FDP_all, power_all, alphas), columns = ['FDP', 'power', 'alpha'])
    df_out['type'] = 'adafdr'
    df_out.to_csv("result_data1/" + "adafdr_" + file, header=True, index=False, sep='\t')

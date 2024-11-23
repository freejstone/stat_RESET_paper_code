#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 23:20:12 2024

@author: jfre0619
"""

import numpy as np
import adafdr.method as md
import adafdr.data_loader as dl
import pandas as pd
import random
import time


random.seed(21072024)
np.random.seed(21072024)
alphas = [0.01, 0.05, 0.1, 0.2]


data_path = './data_files'
file_path = data_path + '/estrogen.csv'
df_data = pd.read_csv(file_path, sep=',')
p = df_data['pvals'].to_numpy()
x = df_data['ord_high'].to_numpy()
print('p:', p.shape)
print('x:', x.shape)
speed = []
FDP = []
power = []
seed = []
alpha_all = []
for i in range(1):
    for alpha in alphas:
        start_time = time.time()
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21072024 + i)
        end_time = time.time()
        n_rej = res['n_rej']
        t_rej = res['threshold']
        speed.append(end_time - start_time)
        FDP.append(alpha)
        power.append(np.sum(p<=t_rej))
        seed.append(i)
        alpha_all.append(alpha)
df_out = pd.DataFrame(zip(FDP, power, alpha_all, speed, seed), columns = ['FDP', 'power', 'alpha', 'time', 'seed'])
df_out['type'] = 'adafdr'
df_out.to_csv("results/estrogen_high_adafdr.csv", header=True, index=False, sep='\t')



data_path = './data_files'
file_path = data_path + '/estrogen.csv'
df_data = pd.read_csv(file_path, sep=',')
p = df_data['pvals'].to_numpy()
x = df_data['ord_mod'].to_numpy()
print('p:', p.shape)
print('x:', x.shape)
speed = []
FDP = []
power = []
seed = []
alpha_all = []
for i in range(1):
    for alpha in alphas:
        start_time = time.time()
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21072024 + i)
        end_time = time.time()
        n_rej = res['n_rej']
        t_rej = res['threshold']
        speed.append(end_time - start_time)
        FDP.append(alpha)
        power.append(np.sum(p<=t_rej))
        seed.append(i)
        alpha_all.append(alpha)
df_out = pd.DataFrame(zip(FDP, power, alpha_all, speed, seed), columns = ['FDP', 'power', 'alpha', 'time', 'seed'])
df_out['type'] = 'adafdr'
df_out.to_csv("results/estrogen_mod_adafdr.csv", header=True, index=False, sep='\t')




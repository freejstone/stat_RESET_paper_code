#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 14:20:47 2024

@author: jfre0619
"""

import numpy as np
import adafdr.method as md
import adafdr.data_loader as dl
import pandas as pd
import random
import time


random.seed(21062024)
np.random.seed(21062024)
alphas = [0.01, 0.05, 0.1, 0.2]

##airway
p,x = dl.data_airway()
print('p:', p.shape)
print('x:', x.shape)
speed = []
FDP = []
power = []
seed = []
alpha_all = []
mask = (x > 0)
p = p[mask.flatten()]
x = x[mask]
for i in range(1):
    print(i)
    for alpha in alphas:
        start_time = time.time()
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/airway_adafdr.csv", header=True, index=False, sep='\t')

##bottomly
p,x = dl.data_bottomly()
print('p:', p.shape)
print('x:', x.shape)
speed = []
FDP = []
power = []
seed = []
alpha_all = []
mask = (x > 0)
p = p[mask.flatten()]
x = x[mask]
for i in range(1):
    for alpha in alphas:
        start_time = time.time()
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/bottomly_adafdr.csv", header=True, index=False, sep='\t')

##pasilla
p,x = dl.data_pasilla()
print('p:', p.shape)
print('x:', x.shape)
speed = []
FDP = []
power = []
seed = []
alpha_all = []
mask = (x > 0)
p = p[mask.flatten()]
x = x[mask]
for i in range(1):
    for alpha in alphas:
        start_time = time.time()
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/pasilla_adafdr.csv", header=True, index=False, sep='\t')

#proteomics
data_path = './data_files'
file_path = data_path + '/proteomics'
df_data = pd.read_csv(file_path, sep=',')
p = df_data['p_val'].to_numpy()
x = df_data['x'].to_numpy()
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=True, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/proteomics_adafdr.csv", header=True, index=False, sep='\t')

#microbiome_enigma_al
data_path = './data_files'
file_path = data_path + '/microbiome_enigma_al'
df_data = pd.read_csv(file_path, sep=',')
p = df_data['p_val'].to_numpy()
x = df_data[['ubiquity', 'mean_abun']].to_numpy()
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=True, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/microbiome_enigma_al_adafdr.csv", header=True, index=False, sep='\t')

#microbiome_enigma_ph
data_path = './data_files'
file_path = data_path + '/microbiome_enigma_ph'
df_data = pd.read_csv(file_path, sep=',')
p = df_data['p_val'].to_numpy()
x = df_data[['ubiquity', 'mean_abun']].to_numpy()
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=True, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/microbiome_enigma_ph_adafdr.csv", header=True, index=False, sep='\t')

#fmri_auditory
data_path = './data_files'
file_path = data_path + '/fmri_auditory'
df_fmri = pd.read_csv(file_path, sep=',')
p = df_fmri['p_val'].to_numpy()
x = df_fmri['B_label'].to_numpy()
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=True, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/fmri_auditory_adafdr.csv", header=True, index=False, sep='\t')

#fmri_imagination
data_path = './data_files'
file_path = data_path + '/fmri_imagination'
df_fmri = pd.read_csv(file_path, sep=',')
p = df_fmri['p_val'].to_numpy()
x = df_fmri['B_label'].to_numpy()
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=True, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/fmri_imagination_adafdr.csv", header=True, index=False, sep='\t')

#Adipose_Subcutaneous
data_name = 'Adipose_Subcutaneous'
p,x,cate_name,cis_name = dl.data_small_gtex_chr21(opt=data_name)
print('p', p.shape)
print('x', p.shape)
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/adipose_subcutaneous_adafdr.csv", header=True, index=False, sep='\t')

#Adipose_Visceral_Omentum
data_name = 'Adipose_Visceral_Omentum'
p,x,cate_name,cis_name = dl.data_small_gtex_chr21(opt=data_name)
print('p', p.shape)
print('x', p.shape)
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
        res = md.adafdr_test(p, x, alpha=alpha, fast_mode=False, single_core=False, random_state = 21062024 + i)
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
df_out.to_csv("results/adipose_visceral_omentum_adafdr.csv", header=True, index=False, sep='\t')







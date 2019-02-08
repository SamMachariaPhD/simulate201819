#Check if there's correlation between binding motor and speed
#Sam. Prof. Nitta Lab.

import os,glob,csv,time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

current_path = os.getcwd()
filename = 'corrdata'
no = 1

dirlist = glob.glob(current_path+'/*/') #[name for name in os.listdir(".") if os.path.isdir(name)]
dirlist = sorted(dirlist, key=lambda x:x[-30:])

columns=['tym','bm','vel']
title = ['Correlation Between Binding Motors And Speed']
space = [' ']

try:
    os.remove(filename+'.csv') # if program run twice, remove the previously saved file.
except (OSError, RuntimeError, TypeError, NameError):
    pass

with open(filename+'.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow(title)
        writer.writerow(space)

for i in dirlist:
    bmsp = pd.read_csv('bmd'+str(no)+'.csv',skiprows=[0],header=None,names=columns)
    corr = bmsp['bm'].corr(bmsp['vel'])
    with open(filename+'.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow([i])
        writer.writerow([corr])
    no = no+1
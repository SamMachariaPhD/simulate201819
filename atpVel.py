"""
Go inside various complete simulation data.
Search for "v.csv".
Compute the average speed, the deviation, and return the dataframe thereof.
Regards, Sam Sirmaxford.
"""

import os, glob, csv
import numpy as np
import pandas as pd

current_path = os.getcwd()
filename = 'spd_std273ATP50'
DtFile=0.01

dirlist = glob.glob(current_path+'/*/')
dirlist = sorted(dirlist, key=lambda x:x[-18:])

for i in dirlist:
    os.chdir(i)
    os.chdir('PLOTS')
    df = pd.read_csv('SkippedTipXY_A001.csv')
    Dx_tip = np.diff(df['x_tip']); Dy_tip = np.diff(df['y_tip'])
    DD=np.sqrt((Dx_tip**2)+(Dy_tip**2))
    v=DD/(10*DtFile); Av_vel = np.mean(v)
    vSD=np.sum(((v-Av_vel)**2)/(np.size(v)-1)); vSD=np.sqrt(vSD)
    os.chdir(current_path)
    with open(filename+'.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow((Av_vel,vSD))

"""
Go inside various complete simulation data.
Search for IntParticle***.vtk
Extract motor number
Regards, Sam Sirmaxford.
"""

import os, glob, csv
import numpy as np
import pandas as pd

current_path = os.getcwd()
filename = 'bm_no'
count = 500 # smallest motor density
seed = '873'

dirlist = glob.glob(current_path+'/*/')
dirlist = sorted(dirlist, key=lambda x:x[-18:])

for i in dirlist:
    os.chdir(i)
    file_list = glob.glob('IntParticle**.vtk') # os.listdir()
    files = sorted(file_list, key=lambda x:x[-11:])
    for j in files:
        motor_no = pd.read_csv(j, delim_whitespace=True)
        extract = motor_no.iloc[3:4,1:2].values
        extract = int(np.int_(extract))
        rowElement = [str(extract)]
        os.chdir(current_path)
        with open(filename+seed+'MD'+str(count)+'.csv','a') as fd:
            writer = csv.writer(fd)
            writer.writerow(rowElement)
        os.chdir(i)
    count=count+500
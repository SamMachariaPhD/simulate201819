# prepared by Sam., Prof. Nitta Lab.
import glob, os, csv, tarfile, shutil, sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

current_path = os.getcwd() # detect the current working dir
filename = 'AnchorForce'; dtOut = 0.1

try:
    os.remove(filename+'.txt') # if program run twice, remove the previously saved file.
except (OSError, RuntimeError, TypeError, NameError):
    pass

# extract NodeData
try:
    check_dir = glob.glob('NodeData**.tar.gz')
    for item in check_dir:
        tar = tarfile.open(item, "r:gz") # must be r or w
        tar.extractall()
        tar.close()
    check_dir = glob.glob('NodeData**.tar')
    for item in check_dir:
        tar = tarfile.open(item, "r:")
        tar.extractall()
        tar.close()
except (OSError, RuntimeError, TypeError, NameError):
    sys.exit("\n=> Tar extraction failed")

file_list = glob.glob('NodeData**.txt') # os.listdir()
files = sorted(file_list, key=lambda x:x[-11:])
column_names = ['c1','c2','c3','c4','c5','c6','c7','c8']
total_files = len(files)
NumFile = 0
count = 1
NumFile_input = int(input("\n=> Out of %s files, enter the number of files to analyse: " %total_files))

if NumFile_input <= 0:
    sys.exit("\n=> Please enter a number > 0")
elif NumFile_input > total_files:
    sys.exit("\n=> Please enter a number <= %s" %total_files)
else:
    pass

for i in files: # make the required dataset
    if count < NumFile_input:
        node_data = pd.read_csv(i, names=column_names, delim_whitespace=True)
        FIx = node_data['c7']; FIy = node_data['c8']
        AnchorI = node_data['c6']
        LeftI = (node_data['c3'] < 0.0)*1 # *1 -- multiply boolean by 1
        dfl = pd.concat([FIx,(AnchorI & LeftI)],axis=1) 
        dfl.columns = ['FIx','Left']
        dfl = dfl.loc[dfl['Left'] == 1]
        TotalAnchorFxLeft = np.sum(dfl['FIx']) # <---- 2nd element
        RightI = (node_data['c3'] > 0.0)*1
        dfr = pd.concat([FIx,(AnchorI & RightI)],axis=1)
        dfr.columns = ['FIx','Right']
        dfr = dfr.loc[dfr['Right'] == 1]
        TotalAnchorFxRight = np.sum(dfr['FIx']) # <---- 3rd element
        Time = NumFile*dtOut # <---- 1st element
        rowElement = np.array([Time, TotalAnchorFxLeft, TotalAnchorFxRight])
        rowElement = rowElement.reshape(1,3)
        with open(filename+'.txt','a') as myfile:
            np.savetxt(myfile, rowElement, delimiter='\t', fmt='%.2f %.20f %.20f')
        NumFile = NumFile+1; count = count+1
    else:
        pass

# do plots
anchor_data = pd.read_csv(filename+'.txt', names=['Time','Left','Right'], delim_whitespace=True)
abs_left = anchor_data['Left'].abs(); abs_right = anchor_data['Right'].abs()
max_left = abs_left.max(); max_right = abs_right.max()
max_left_pt = anchor_data.loc[abs_left == max_left]; max_right_pt = anchor_data.loc[abs_right == max_right]
np.savetxt('max_left.txt', max_left_pt, delimiter='\t', fmt='%.2f %.20f %.20f')
np.savetxt('max_right.txt', max_right_pt, delimiter='\t', fmt='%.2f %.20f %.20f')

plt.figure(figsize=(8,6), dpi=300)
plt.plot(anchor_data['Time'], anchor_data['Left'], 'b', marker='D', markersize=1, label='F-Left')
plt.scatter(max_left_pt['Time'], max_left_pt['Left'], s=25, facecolors='green', edgecolor='blue', label='Max-F-Left')
plt.plot(anchor_data['Time'], anchor_data['Right'], 'r', marker='D', markersize=1, label='F-Right')
plt.scatter(max_right_pt['Time'], max_right_pt['Right'], s=25, facecolors='green', edgecolor='red', label='Max-F-Right')
plt.xlabel('t (s)'); plt.ylabel('F (pN)'); plt.legend() # ; plt.grid()
plt.savefig(filename+'.png', format='png'); plt.savefig(filename+'.svg', format='svg')

# go back to current dir and delete extracted files
os.chdir(current_path)

check_dir = glob.glob('NodeData**.txt')

for item in check_dir:
    os.remove(item)

try:
    shutil.rmtree('NodeData') # remove full folder
except (OSError, RuntimeError, TypeError, NameError):
    pass

print("\n=> Done.\n")
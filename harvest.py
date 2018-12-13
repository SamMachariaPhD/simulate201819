import glob, os, csv, tarfile
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

current_path = os.getcwd() # detect the current working dir

filename = 'min_max'

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
columns = ['one','two','x','y','three','four','five','six']
indexing = 0

for i in files:
    node_data = pd.read_csv(i, names=columns, delim_whitespace=True)
    x_node_data = node_data['x']
    node_data_xmin = x_node_data.min()
    node_data_xmax = x_node_data.max()
    x_minmax_diff = (node_data_xmax-node_data_xmin)
    y_node_data = node_data['y']
    node_data_ymin = y_node_data.min()
    node_data_ymax = y_node_data.max()
    y_minmax_diff = (node_data_ymax-node_data_ymin)
    rowElement = np.array([indexing, node_data_xmin, node_data_xmax, node_data_ymin, node_data_ymax, x_minmax_diff, y_minmax_diff])
    rowElement = rowElement.reshape(1,7)
    with open(filename+'.txt','a') as myfile:
        np.savetxt(myfile, rowElement, delimiter='\t', fmt='%i %.10f %.10f %.10f %.10f %.10f %.10f')
    indexing = indexing+1

os.chdir(current_path)

check_dir = glob.glob('NodeData**.txt')
for item in check_dir:
    os.remove(item)

try:
    shutil.rmtree('NodeData_A001') # remove full folder
except (OSError, RuntimeError, TypeError, NameError):
    pass

print("\n=> I'm done making %s file.\n" %filename)
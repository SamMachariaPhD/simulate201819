# prepared by Sam., Prof. Nitta Lab.
import glob, os, csv, tarfile
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

current_path = os.getcwd() # detect the current working dir

filename = 'binding_motors'

# extract IntParticle
try:
    check_dir = glob.glob('IntParticle**.tar.gz')
    for item in check_dir:
        tar = tarfile.open(item, "r:gz") # must be r or w
        tar.extractall()
        tar.close()
    check_dir = glob.glob('IntParticle**.tar')
    for item in check_dir:
        tar = tarfile.open(item, "r:")
        tar.extractall()
        tar.close()
except (OSError, RuntimeError, TypeError, NameError):
    sys.exit("\n=> Tar extraction failed")

file_list = glob.glob('IntParticle**.vtk') # os.listdir()
files = sorted(file_list, key=lambda x:x[-11:])

for i in files:
    motor_no = pd.read_csv(i, delim_whitespace=True)
    rowElement = [str(int(motor_no.iloc[3:4,1:2].values))]
    with open(filename+'.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow(rowElement)

# go back to current dir and delete extracted files
os.chdir(current_path)

ls = os.listdir(current_path)
for item in ls:
    if item.endswith(".vtk"):
        os.remove(os.path.join(current_path, item))

try:
    shutil.rmtree('IntParticle_A001') # remove full folder
except (OSError, RuntimeError, TypeError, NameError):
    pass

print("\n=> I'm done making %s.csv file.\n" %filename)
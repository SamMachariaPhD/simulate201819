import os
import glob

current_path = os.getcwd()
os.chdir(current_path+'/test/')

files = glob.glob('IntParticle**.vtk')
files = sorted(files, key=lambda x:x[-16:])
for file in files:
    try:
        os.remove(file)
    except:
        print("Error while deleting file : ", file)
        pass
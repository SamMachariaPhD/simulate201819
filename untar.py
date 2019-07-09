import glob, os, csv, tarfile, shutil, sys

current_path = os.getcwd()
os.chdir(current_path+'/test/')

try:
    check_dir = glob.glob('IntParticle**.tar')
    for item in check_dir:
        tar = tarfile.open(item, "r:")
        tar.extractall()
        tar.close()
except (OSError, RuntimeError, TypeError, NameError):
    sys.exit("\n=> Tar extraction failed")

files = glob.glob('IntParticle**.vtk')
sFiles = sorted(files, key=lambda x:x[-16:])
print(sFiles)
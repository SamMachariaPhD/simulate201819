import tarfile, glob, os, sys

current_path = os.getcwd()
os.chdir(current_path+'/test/')

grp1 = glob.glob('ContactStates**.txt')
grp1 = sorted(grp1, key=lambda x:x[-16:])
grp2 = glob.glob('Filament**.vtk')
grp2 = sorted(grp2, key=lambda x:x[-16:])
grp3 = glob.glob('IntParticle**.vtk')
grp3 = sorted(grp3, key=lambda x:x[-16:])
grp4 = glob.glob('MotorSpecie1**.vtk')
grp4 = sorted(grp4, key=lambda x:x[-16:])
grp5 = glob.glob('MotorSpecie2**.vtk')
grp5 = sorted(grp5, key=lambda x:x[-16:])
grp6 = glob.glob('MTPlusEnd**.vtk')
grp6 = sorted(grp6, key=lambda x:x[-16:])
grp7 = glob.glob('Particle**.vtk')
grp7 = sorted(grp7, key=lambda x:x[-16:])
#1=ContactStates=====================================================#
try:
    with tarfile.open('ContactStates.tar', mode='w') as tar:
        for file in grp1:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    #sys.exit("\n=> Tar creation failed")
    pass
for filePath in grp1:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('ContactStates.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)
#2=Filament=====================================================#
try:
    with tarfile.open('Filament.tar', mode='w') as tar:
        for file in grp2:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    pass
for filePath in grp2:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('Filament.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)
#3=IntParticle=====================================================#
try:
    with tarfile.open('IntParticle.tar', mode='w') as tar:
        for file in grp3:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    pass
for filePath in grp3:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('IntParticle.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)
#4=MotorSpecie1=====================================================#
try:
    with tarfile.open('MotorSpecie1.tar', mode='w') as tar:
        for file in grp4:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    pass
for filePath in grp4:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('MotorSpecie1.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)
#5=MotorSpecie2=====================================================#
try:
    with tarfile.open('MotorSpecie2.tar', mode='w') as tar:
        for file in grp5:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    pass
for filePath in grp5:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('MotorSpecie2.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)
#6=MTPlusEnd=====================================================#
try:
    with tarfile.open('MTPlusEnd.tar', mode='w') as tar:
        for file in grp6:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    pass
for filePath in grp6:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('MTPlusEnd.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)
#7=Particle=====================================================#
try:
    with tarfile.open('Particle.tar', mode='w') as tar:
        for file in grp7:
            tar.add(file)
except (OSError, RuntimeError, TypeError, NameError):
    print("\n=> Some tar was not created.")
    pass
for filePath in grp7:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
        pass
# Read the contents of the newly created archive
with tarfile.open('Particle.tar', mode='r') as t:
    for member in t.getmembers():
        print(member.name)

os.chdir(current_path)
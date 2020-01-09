#Take intParticle and MotorSpecie1 in each case
#Calculate how many motors are active out of the total binding motors

import os, glob, csv
import numpy as np
import pandas as pd

current_path = os.getcwd()
filename = 'Abm_no'
R = 0.1 
seed = '273ATP50'

dirlist = glob.glob(current_path+'/*/')
dirlist = sorted(dirlist, key=lambda x:x[-18:])

for i in dirlist:
    os.chdir(i)
    intParticleLst = glob.glob('IntParticle**.vtk') # os.listdir()
    srtdIntParticle = sorted(intParticleLst, key=lambda x:x[-11:])
    intMotorSpecieLst = glob.glob('MotorSpecie1**.vtk') # os.listdir()
    srtdMotorSpecie = sorted(intMotorSpecieLst, key=lambda x:x[-11:])

    for j in range(0, len(srtdIntParticle)-1):
        intParticle = pd.read_csv(srtdIntParticle[j], delim_whitespace=True, names=['x','y','z','nan','nan2'])
        motorno = intParticle.iloc[4:5,1:2].values
        motorNo = int( np.int_(motorno) )

        intParticle = intParticle.drop(['nan','nan2'], axis=1)
        intParticl = intParticle.drop(intParticle.index[motorNo+5:intParticle.index[-1]+1])
        intPartic = intParticl.drop(intParticl.index[0:5])
        intPartic = intPartic.values

        MotorSpecie1 = pd.read_csv(srtdMotorSpecie[j], delim_whitespace=True, names=['x','y','z','nan','nan2'])
        Sp1no = MotorSpecie1.iloc[4:5,1:2].values
        Sp1No = int( np.int_(Sp1no) )

        MotorSpecie1_ = MotorSpecie1.drop(['nan','nan2'], axis=1)
        MotorSpecie1_ = MotorSpecie1_.drop(MotorSpecie1_.index[Sp1No+5:MotorSpecie1_.index[-1]+1])
        MotorSpecie1_ = MotorSpecie1_.drop(MotorSpecie1_.index[0:5])
        MotorSpecie1_ = MotorSpecie1_.values

        activeM=0
        for k in range(0, motorNo):
            for l in range(0,MotorSpecie1_.shape[0]):
                if(np.sum(1*np.equal(MotorSpecie1_[l],intPartic[k])) == 3):
                    activeM = activeM+1

        os.chdir(current_path)
        with open(filename+seed+'R'+str(np.round(R,1))+'.csv','a') as fd:
            writer = csv.writer(fd)
            if activeM > 0:
                writer.writerow([activeM])
            else:
                writer.writerow([0])
        os.chdir(i)
    R=R+0.1
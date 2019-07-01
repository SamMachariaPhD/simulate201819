#!/usr/bin/env python3
#Check execution time
#Sam.
import fileinput, sys, shutil, os, time, socket, subprocess
#import numpy as np
import matplotlib.pyplot as plt

tic = time.time()
subprocess.call("ifort jacobi_openmp.f90 -o jcb", shell=True)
subprocess.call("./jcb >> jacobi_openmp.txt", shell=True)
toc = time.time()
ttym1 = toc-tic
print("\n No openmp time = ", ttym1)
subprocess.call("rm jcb", shell=True)

tic = time.time()
subprocess.call("ifort -openmp jacobi_openmp.f90 -o jcb", shell=True)
subprocess.call("./jcb >> jacobi_openmp.txt", shell=True)
toc = time.time()
ttym2 = toc-tic
print("\n No thread time = ", ttym2)
subprocess.call("rm jcb", shell=True)

tic = time.time()
subprocess.call("ifort -openmp jacobi_openmp.f90 -o jcb", shell=True)
subprocess.call("export OMP_NUM_THREADS=1;./jcb >> jacobi_openmp.txt", shell=True)
toc = time.time()
ttym3 = toc-tic
print("\n 1 thread time = ", ttym3)
subprocess.call("rm jcb", shell=True)

tic = time.time()
subprocess.call("ifort -openmp jacobi_openmp.f90 -o jcb", shell=True)
subprocess.call("export OMP_NUM_THREADS=2;./jcb >> jacobi_openmp.txt", shell=True)
toc = time.time()
ttym4 = toc-tic
print("\n 2 thread time = ", ttym4)
subprocess.call("rm jcb", shell=True)

tic = time.time()
subprocess.call("ifort -openmp jacobi_openmp.f90 -o jcb", shell=True)
subprocess.call("export OMP_NUM_THREADS=4;./jcb >> jacobi_openmp.txt", shell=True)
toc = time.time()
ttym5 = toc-tic
print("\n 4 thread time = ", ttym5)
subprocess.call("rm jcb", shell=True)

tic = time.time()
subprocess.call("ifort -openmp jacobi_openmp.f90 -o jcb", shell=True)
subprocess.call("export OMP_NUM_THREADS=8;./jcb >> jacobi_openmp.txt", shell=True)
toc = time.time()
ttym6 = toc-tic
print("\n 8 thread time = ", ttym6)
subprocess.call("rm jcb", shell=True)

x = [-1,0,1,2,4,8]
y = [ttym1,ttym2,ttym3,ttym4,ttym5,ttym6]

plt.figure(figsize=(10,6), dpi=500)
plt.plot(x,y) #, capsize=7, linewidth=2, markersize=7
plt.xlabel('OpenMP Threads'); plt.ylabel('Execution Time')
plt.title('CPU Threads vs. Execution Time'); plt.legend(loc='best')
plt.savefig('jacobi.svg',bbox_inches='tight', format='svg',dip=500)
plt.savefig('jacobi.png',bbox_inches='tight', format='svg',dip=500)
plt.close()
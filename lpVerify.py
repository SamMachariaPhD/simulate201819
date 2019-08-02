#!/usr/bin/env python3
# Plot persistence length graphs
# Sam

import io, csv, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import sympy as sym
plt.style.use('default')
plt.rcParams.update({'font.size': 13})

confName = 'mtpaths.txt'
v=7.19351; dt=0.1

Conf = pd.read_csv(confName, names=['mt','no','x','y','th'], delim_whitespace=True)

conf = Conf[0:600]

conf.to_csv('trConf.csv',index=False)

xmax_ = conf['x'].max(); xmin_ = conf['x'].min()
ymax_ = conf['y'].max(); ymin_ = conf['y'].min()
x1_ = (0.01*(xmax_-xmin_))+xmin_
y1_ = (0.90*(ymax_-ymin_))+ymin_
x2_ = (0.15*(xmax_-xmin_))+xmin_ 
y2_ = (0.85*(ymax_-ymin_))+ymin_

plt.figure(figsize=(10,8))
plt.plot(conf['x'],conf['y'], label='Leading tip', color='green', marker='.', linestyle='solid')
plt.text(x1_, y1_, 'avSpeed = '); plt.text(x2_, y1_, r'%.5f $\mu m/sec$'%v)
plt.xlabel('X'); plt.ylabel('Y')
plt.title('Actin Filament Leading Tip Movement'); plt.legend(loc='upper left')
plt.savefig('leadTrajectory.svg',bbox_inches='tight', format='svg', dpi=500)
plt.close()

lIndex = conf.shape[0]-2

p = 0; q = 1; df = conf

try:
    os.remove('ub.csv')
except(OSError, RuntimeError, TypeError, NameError):
    pass

while q<lIndex:
    x0 = df.loc[p,'x']
    x1 = df.loc[q,'x']
    y0 = df.loc[p,'y']
    y1 = df.loc[q,'y']
    u = np.array([(x1-x0),(y1-y0)])
    ub = np.true_divide(u,np.sqrt((x1-x0)**2 + (y1-y0)**2))
    with open('ub.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow(ub)
    p = p+1; q = q+1
ub = pd.read_csv('ub.csv', names=['dx','dy'])

no=0

while no<lIndex+2:
  try:
    os.remove('ds'+str(no)+'.csv')
    os.remove('Ds.csv')
  except(OSError, RuntimeError, TypeError, NameError):
    pass
  no=no+1

i = 0; ii = 0; il = lIndex-1; s=0; no = 1

while il>0:
    while i<il:
        ds = np.dot(ub.loc[i],ub.loc[(i+ii)])
        with open('ds'+str(no)+'.csv','a') as fd:
            writer = csv.writer(fd)
            writer.writerow([ds])
        i=i+1
        #print('i = ',i)
    ds = pd.read_csv('ds'+str(no)+'.csv',names=['ds'])
    Ds = float(ds.mean())
    s_ = v*dt*s
    rows = np.array([s_,Ds])
    with open('Ds.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow(rows)
    i=0; ii=ii+1; il=il-1; s=s+1; no=no+1
    #print('il = ',il)
    #os.remove('ds1.csv')

Ds = pd.read_csv('Ds.csv',names=['s','Ds'])

plt.figure(figsize=(10,8))
plt.scatter(Ds['s'],Ds['Ds'],facecolors='none',edgecolors='b', label=r'Correlation')
plt.title(r'$S\ vs.\ <\cos\ \Delta \Theta >$')
plt.xlabel(r'$S\ (\mu m)$'); plt.ylabel(r'$<\cos \Delta \Theta >$')
plt.legend(loc='best')
plt.savefig('allPersistence.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

tenth = int(np.around(Ds.shape[0]/10,0)); Ds=Ds[0:tenth] #10th(10%) or 20eth(5%)

plt.figure(figsize=(10,8))
plt.scatter(Ds['s'],Ds['Ds'],facecolors='none',edgecolors='b', label=r'Correlation')
plt.title(r'$S\ vs.\ <\cos \Delta \Theta >$')
plt.xlabel(r'$S (\mu m)$'); plt.ylabel(r'$<\cos \Delta \Theta >$')
plt.legend(loc='best')
plt.savefig('tenthPersistence.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

#curve fitting
plt.figure(figsize=(10,8))
plt.scatter(Ds['s'],Ds['Ds'],facecolors='none',edgecolors='b', label=r'Correlation')
xf = np.array(Ds['s'], dtype=float) #transform your data in a numpy array of floats 
yf = np.array(Ds['Ds'], dtype=float) #so the curve_fit can work
#curve function
def func(xf, Lp):
    return np.exp(-xf/(2*Lp))
#curve fit
popt, pcov = curve_fit(func, xf, yf)
plt.plot(xf, func(xf, *popt), color='green', label=r'Fitted Curve: $<\cos \Delta \Theta >=\exp \left(\frac{-S}{2Lp}\right)$')
plt.text(30, 0.8, r'$Lp = %s\ \mu m$'%np.around(popt[0],6))
plt.title(r'$S\ vs.\ <\cos \Delta \Theta >$')
plt.xlabel(r'$S (\mu m)$'); plt.ylabel(r'$<\cos \Delta \Theta >$')
plt.legend(loc='best')
plt.savefig('curveFit.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()
print("curveFitLp = %s micrometer" %np.around(popt[0],6))

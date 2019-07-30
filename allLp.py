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
#========================================================================
Ds1 = pd.read_csv('Ds1.csv',names=['s','Ds'])
Ds2 = pd.read_csv('Ds2.csv',names=['s','Ds'])
Ds3 = pd.read_csv('Ds3.csv',names=['s','Ds'])
Ds4 = pd.read_csv('Ds4.csv',names=['s','Ds'])
Ds5 = pd.read_csv('Ds5.csv',names=['s','Ds'])
#========================================================================
tenth1 = int(np.around(Ds1.shape[0]/10,0)); Ds1=Ds1[0:tenth1] #10th(10%)
tenth2 = int(np.around(Ds2.shape[0]/10,0)); Ds2=Ds2[0:tenth2] #10th(10%)
tenth3 = int(np.around(Ds3.shape[0]/10,0)); Ds3=Ds3[0:tenth3] #10th(10%)
tenth4 = int(np.around(Ds4.shape[0]/10,0)); Ds4=Ds4[0:tenth4] #10th(10%)
tenth5 = int(np.around(Ds5.shape[0]/10,0)); Ds5=Ds5[0:tenth5] #10th(10%)
#========================================================================
#curve fitting
plt.figure(figsize=(10,8))
plt.scatter(Ds1['s'],Ds1['Ds'],facecolors='none',edgecolors='r', label=None) # , label=r'AF1, R0.7'
xf = np.array(Ds1['s'], dtype=float) #transform your data in a numpy array of floats 
yf = np.array(Ds1['Ds'], dtype=float) #so the curve_fit can work
#curve function
def func(xf, Lp):
    return np.exp(-xf/(2*Lp))
#curve fit
popt, pcov = curve_fit(func, xf, yf)
plt.plot(xf, func(xf, *popt), color='red', label=r'$AF1, R0.7, Lp = %s\ \mu m$'%np.around(popt[0],6)) 
#========================================================================
plt.scatter(Ds2['s'],Ds2['Ds'],facecolors='none',edgecolors='g', label=None)
xf = np.array(Ds2['s'], dtype=float) #transform your data in a numpy array of floats 
yf = np.array(Ds2['Ds'], dtype=float) #so the curve_fit can work
#curve function
def func(xf, Lp):
    return np.exp(-xf/(2*Lp))
#curve fit
popt, pcov = curve_fit(func, xf, yf)
plt.plot(xf, func(xf, *popt), color='green', label=r'$AF1, R1.0, Lp = %s\ \mu m$'%np.around(popt[0],6))
#========================================================================
plt.scatter(Ds3['s'],Ds3['Ds'],facecolors='none',edgecolors='b', label=None)
xf = np.array(Ds3['s'], dtype=float) #transform your data in a numpy array of floats 
yf = np.array(Ds3['Ds'], dtype=float) #so the curve_fit can work
#curve function
def func(xf, Lp):
    return np.exp(-xf/(2*Lp))
#curve fit
popt, pcov = curve_fit(func, xf, yf)
plt.plot(xf, func(xf, *popt), color='blue', label=r'$AF1, R1.0, Lp = %s\ \mu m$'%np.around(popt[0],6))
#========================================================================
plt.scatter(Ds4['s'],Ds4['Ds'],facecolors='none',edgecolors='y', label=None)
xf = np.array(Ds4['s'], dtype=float) #transform your data in a numpy array of floats 
yf = np.array(Ds4['Ds'], dtype=float) #so the curve_fit can work
#curve function
def func(xf, Lp):
    return np.exp(-xf/(2*Lp))
#curve fit
popt, pcov = curve_fit(func, xf, yf)
plt.plot(xf, func(xf, *popt), color='yellow', label=r'$AF2, R1.0, Lp = %s\ \mu m$'%np.around(popt[0],6))
#========================================================================
plt.scatter(Ds5['s'],Ds5['Ds'],facecolors='none',edgecolors='k', label=None) 
xf = np.array(Ds5['s'], dtype=float) #transform your data in a numpy array of floats 
yf = np.array(Ds5['Ds'], dtype=float) #so the curve_fit can work
#curve function
def func(xf, Lp):
    return np.exp(-xf/(2*Lp))
#curve fit
popt, pcov = curve_fit(func, xf, yf)
plt.plot(xf, func(xf, *popt), color='black', label=r'$AF3, R1.0, Lp = %s\ \mu m$'%np.around(popt[0],6))
#========================================================================
#plt.text(30, 0.8, r'$Lp = %s\ \mu m$'%np.around(popt[0],6))
plt.title(r'$S\ vs.\ <\cos \Delta \Theta >$')
plt.xlabel(r'$S (\mu m)$'); plt.ylabel(r'$<\cos \Delta \Theta >$')
plt.legend(loc='best')
plt.savefig('aCurveFit.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()
#print("curveFitLp = %s micrometer" %np.around(popt[0],6))

# Plot persistence length graphs
# Sam

import io, csv, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
from scipy.interpolate import interp1d
plt.style.use('default')

confName = 'Conformation_30sec1.txt'
v=7.12534; dt=0.1

Conf = pd.read_csv(confName, names=['t','x','y','z'], delim_whitespace=True)

conf1 = Conf.iloc[0::13,:]
conf1 = conf1.reset_index(drop=True)
conf2 = conf1.iloc[0::10, :]
conf2 = conf2[['x','y']]
conf = conf2.reset_index(drop=True)

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
    rows = np.array([s,Ds])
    with open('Ds.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow(rows)
    i=0; ii=ii+1; il=il-1; s=s+1; no=no+1
    #print('il = ',il)
    #os.remove('ds1.csv')

Ds = pd.read_csv('Ds.csv',names=['s','Ds'])

plt.figure(figsize=(10,8))
plt.scatter(v*dt*Ds['s'],Ds['Ds'],facecolors='none',edgecolors='b', label='Filament 1')
#plt.scatter(v*dt*Ds2['s'],Ds2['Ds'],facecolors='none',edgecolors='b', label='Filament 2')
plt.title('S vs <cos $\Delta \Theta $>')
plt.xlabel('S'); plt.ylabel('<cos $\Delta \Theta $>')
plt.legend(loc='best')
plt.savefig('allPersistence.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

tenth = int(np.around(Ds.shape[0]/10,0)); Ds=Ds[0:tenth]

plt.figure(figsize=(10,8))
plt.scatter(v*dt*Ds['s'],Ds['Ds'],facecolors='none',edgecolors='b', label='Filament 1')
#plt.scatter(v*dt*Ds2['s'],Ds2['Ds'],facecolors='none',edgecolors='b', label='Filament 2')
plt.title('S vs <cos $\Delta \Theta $>')
plt.xlabel('S'); plt.ylabel('<cos $\Delta \Theta $>')
plt.legend(loc='best')
plt.savefig('tenthPersistence.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

x = v*dt*Ds['s']
y = np.log(Ds['Ds'])
plt.figure(figsize=(10,8))
plt.scatter(x,y,facecolors='none',edgecolors='b')
plt.title('S vs. log<cos $\Delta \Theta $>')
plt.xlabel('S')
plt.ylabel('log<cos $\Delta \Theta $>')
plt.legend(loc='best')
plt.savefig('tenthLogPersistence.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

z = np.polyfit(x, y, 1)
p = np.poly1d(z)
slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x,p(x))
x_1 = np.max(abs(x))/2 #10.0
y_ = -np.max(abs(y))/2
y_1 = slope*x_1+intercept
c_line ='y={0:.5f}x+{1:.5f}'.format(z[0],z[1])
plt.figure(figsize=(10,8))
plt.scatter(x,y, facecolors='none', edgecolors='b', label='log<cos $\Delta \Theta $>')
plt.plot(x,p(x), color='r', label='trendline', linewidth=1, marker="_")
plt.annotate('y={0:.7f}x+{1:.7f}'.format(z[0],z[1]), xy=(x_1,y_1),xycoords='data',\
             xytext=(x_1,y_),textcoords='data',arrowprops=dict(arrowstyle="->",connectionstyle='arc3,rad=-.5'))
#plt.text(10.1,7,'r value ='); plt.text(11.6,7,r_value)
#plt.text(10.1,6.6,'p value ='); plt.text(11.66,6.6,p_value)# "%.2g" %p_value for 2dp
#plt.text(10.1,6.2,'slope std error ='); plt.text(12.83,6.2,slope_std_error)
plt.title('S vs. log<cos $\Delta \Theta $>')
plt.xlabel('S')
plt.ylabel('log<cos $\Delta \Theta $>')
plt.legend(loc='best')
plt.savefig('logFitPersistence.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

Lp = -1/(2*z[0])
print("Persistence Length = ",np.around(Lp,5), "micrometer")
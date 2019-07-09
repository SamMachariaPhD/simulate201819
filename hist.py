import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
import matplotlib.mlab as mlab

df=pd.read_csv('bmotors.csv',names=['Binding Motors Distribution'])

histo = df.hist(bins=50, figsize=(12,8), grid=False)
for ax in histo.flatten():
    ax.set_xlabel("Motor Number")
    ax.set_ylabel("Frequency")
plt.savefig('bmHist.svg',bbox_inches='tight', format='svg',dip=500)
df=df.loc[(df!=0).any(axis=1)] #remove zeros from df

df=df.values

plt.figure(figsize=(12,8))
(mu,sigma)=norm.fit(df)
n, bins, patches=plt.hist(df,'auto',normed=1,facecolor='gray',alpha=0.5, histtype='bar',ec='black')
y=mlab.normpdf(bins,mu,sigma)
#plt.plot(bins,y,'r-', linewidth=1.3)
plt.xlabel('Binding Motor Number')
plt.ylabel('Probability Density')
plt.title(r'$\mathrm{Histogram\ of\ Binding\ Motors:}\ \mu=%.4f,\ \sigma=%.4f$' %(mu, sigma))
#plt.grid(True)
#plt.legend(loc='best')
plt.savefig('bmHistFitH.svg',bbox_inches='tight', format='svg',dip=500)

plt.figure(figsize=(12,8))
(mu,sigma)=norm.fit(df)
n, bins, patches=plt.hist(df,'auto',normed=1,facecolor='gray',alpha=0.5, histtype='bar',ec='black')
y=mlab.normpdf(bins,mu,sigma)
plt.plot(bins,y,'r-', linewidth=1.3, label='Gaussian P.D.F')
plt.xlabel('Binding Motor Number')
plt.ylabel('Probability Density')
plt.title(r'$\mathrm{Gaussian\ PDF\ and\ Histogram\ of\ Binding\ Motors:}\ \mu=%.4f,\ \sigma=%.4f$' %(mu, sigma))
#plt.grid(True)
plt.legend(loc='best')
plt.savefig('bmHistFitC.svg',bbox_inches='tight', format='svg',dip=500)


plt.figure(figsize=(12,8)); plt.grid()
#mu, sigma = meanVal,stdVal
s = np.random.normal(mu,sigma,6000)
count,bins,ignored = plt.hist(s,'auto',normed=True,facecolor='gray',alpha=0.5, histtype='bar',ec='black')
#plt.plot(bins, 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2)),linewidth=1.5,color='r', label='Gaussian P.D.F')
plt.title(r'$\mathrm{Gaussian\ PDF\ and\ Histogram\ of\ Binding\ Motors:}\ \mu=%.4f,\ \sigma=%.4f$' %(mu, sigma))
plt.xlabel('Binding Motor Number')
plt.ylabel('Probability Density')
plt.legend(loc='best')
plt.grid(0)
plt.savefig('bmHistFith.svg',bbox_inches='tight', format='svg',dip=500)

plt.figure(figsize=(12,8)); plt.grid()
#mu, sigma = meanVal,stdVal
s = np.random.normal(mu,sigma,6000)
count,bins,ignored = plt.hist(s,'auto',normed=True,facecolor='gray',alpha=0.5, histtype='bar',ec='black')
plt.plot(bins, 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2)),linewidth=1.5,color='r', label='Gaussian P.D.F')
plt.title(r'$\mathrm{Gaussian\ PDF\ and\ Histogram\ of\ Binding\ Motors:}\ \mu=%.4f,\ \sigma=%.4f$' %(mu, sigma))
plt.xlabel('Binding Motor Number')
plt.ylabel('Probability Density')
plt.legend(loc='best')
plt.grid(0)
plt.savefig('bmHistFitc.svg',bbox_inches='tight', format='svg',dip=500)
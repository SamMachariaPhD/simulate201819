import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt('rnd51.txt',unpack=True)
y = np.loadtxt('rnd52.txt',unpack=True)
binNo = 30
#n = y.flatten(); print(n.shape)
#n = np.fliplr(n[::-1])

plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
plt.plot(x,y,'.',color='blue')
#plt.plot(y,'.',color='green')
plt.title('randonm xy')
plt.subplot(2,2,2)
plt.hist(x,bins=binNo,color='blue')#, normed=True)
plt.title('x distribution')
plt.subplot(2,2,3)
plt.hist(y,bins=binNo,color='green')#,normed=True)
plt.title('y distribution')
plt.show()

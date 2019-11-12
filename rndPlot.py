import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('rnd.txt')

plt.figure(figsize=(10,10))
plt.subplot(1,2,1)
plt.plot(data,'.')
plt.subplot(1,2,2)
plt.hist(data, bins=200,normed=True)
plt.show()

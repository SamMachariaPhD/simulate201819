# prepared by Sam. feel free to consult (sirmaxford@gmail.com).
import fileinput, sys, shutil, os, time, socket, subprocess
import matplotlib; matplotlib.use('Agg') #set matplotlib to not use the Xwindows backend
import matplotlib.pyplot as plt; import pandas as pd; import numpy as np
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

results_dir = 'PLOTS'
plot_name = 'actin_graph'

try:  
    os.mkdir(results_dir)
except OSError:  
    print ("=> Creation of the directory: %s failed" % results_dir)
else:  
    print ("=> Successfully created %s directory." % results_dir)

columns = ['time', 'x_tip', 'y_tip']; DtFile=0.01
df_load = pd.read_csv('TipXY_A001.txt', names=columns, delim_whitespace=True)
df_nice = df_load.drop(['time'], axis=1)
df = df_nice.iloc[0::9, :] #pick every 9th row starting from the first
df.to_csv(results_dir+'/SkippedTipXY_A001.csv', encoding='utf-8', index=False)

Dx_tip = np.diff(df['x_tip']); Dy_tip = np.diff(df['y_tip'])
DD=np.sqrt((Dx_tip**2)+(Dy_tip**2))
np.savetxt(results_dir+'/DD_A001.csv', DD, delimiter=',')
v=DD/(10*DtFile); Av_vel = np.mean(v)
np.savetxt(results_dir+'/v.csv', v, delimiter=',')
vSD=np.sum(((v-Av_vel)**2)/(np.size(v)-1)); vSD=np.sqrt(vSD)

plt.figure(figsize=(10,8))
plt.text(6, 8, 'Av_Vel = ', fontdict=font); plt.text(8, 8, '%.5f'%Av_vel, fontdict=font)
plt.text(6, 7, 'Vel_SD = ', fontdict=font); plt.text(8, 7, '%.5f'%vSD, fontdict=font)
plt.plot(df['x_tip'],df['y_tip'], label='Filament', color='green', marker='o', linestyle='dashed', linewidth=2, markersize=7)
plt.xlabel('X Tip', fontdict=font); plt.ylabel('Y Tip', fontdict=font)
plt.title('Actin Filament Movement'); plt.legend(loc='upper left'); plt.grid()
plt.savefig(results_dir+'/'+plot_name+'.svg', format='svg', dpi=1200)
plt.savefig(results_dir+'/'+plot_name+'.png', format='png')

print ("=> %s successfully saved. Done!\n" %results_dir)

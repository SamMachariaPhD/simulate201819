import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }
columns = ['time', 'x_tip', 'y_tip']; DtFile=0.01
df_load = pd.read_csv('TipXY_A001.txt', names=columns, delim_whitespace=True)
df_nice = df_load.drop(['time'], axis=1)
df = df_nice.iloc[0::9, :]
df.to_csv('SkippedTipXY_A001.csv', encoding='utf-8', index=False)

Dx_tip = np.diff(df['x_tip']); Dy_tip = np.diff(df['y_tip'])
DD=np.sqrt((Dx_tip**2)+(Dy_tip**2))
np.savetxt('DD_A001.csv', DD, delimiter=',')
v=DD/(10*DtFile); Av_vel = np.mean(v)
np.savetxt('v.csv', v, delimiter=',')
vSD=np.sum(((v-Av_vel)**2)/(np.size(v)-1)); vSD=np.sqrt(vSD)

plt.figure(figsize=(10,8))
plt.text(6, 8, 'Av_Vel = ', fontdict=font); plt.text(8, 8, '%.5f'%Av_vel, fontdict=font)
plt.text(6, 7, 'Vel_SD = ', fontdict=font); plt.text(8, 7, '%.5f'%vSD, fontdict=font)
plt.plot(df['x_tip'],df['y_tip'], label='Filament', color='green', marker='o', linestyle='dashed', linewidth=2, markersize=7)
plt.xlabel('X Tip', fontdict=font); plt.ylabel('Y Tip', fontdict=font)
plt.title('Actin Filament Movement'); plt.legend(loc='upper left'); plt.grid()
plt.savefig('actin_graph.svg', format='svg', dpi=1200)

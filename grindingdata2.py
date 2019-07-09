#After simulation, get inside each folder.
#Extract the binding motor number.
#Extract speed and its deviation.
#Sort data, do instructed analysis, and later make a latex doc. using ready made template.
#Inputs => 10 complete simulations, say, R=1.0 to 1.0
#Outputs => 1. grinding_data.csv having paths for simulation, binding motors, speed and average speed.
#           2. binding motor vs. average speed plots
#           3. leading tip plots
#           4. sorted data outputs
#Sam. Prof. Nitta Lab.

import os,glob,csv,time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

current_path = os.getcwd()
filename = 'grinding_data'

dirlist = glob.glob(current_path+'/*/') #[name for name in os.listdir(".") if os.path.isdir(name)]
dirlist = sorted(dirlist, key=lambda x:x[-30:])

dt = 0.01 # change in stime between rows
conf_names = ['stime','x','y','z']
rownames = ['stime','bmotors','vel']
columnslt = ['stime', 'x_tip', 'y_tip']
avsd = ['avSpeed','std_dev']
binding_motor = 0
no = 1

try:
    os.remove(filename+'.csv') # if program run twice, remove the previously saved file.
    os.remove('bmotors.csv')
    os.remove('summarydata.csv'); os.remove('bmd1.csv'); os.remove('bmd2.csv')
    os.remove('bmd3.csv'); os.remove('bmd4.csv'); os.remove('bmd5.csv'); os.remove('bmd6.csv')
    os.remove('bmd7.csv'); os.remove('bmd8.csv'); os.remove('bmd9.csv'); os.remove('bmd10.csv')
except (OSError, RuntimeError, TypeError, NameError):
    pass

print('\n ==> Now mining data ...\n')
tic = time.time()

for i in dirlist:
    os.chdir(i)
    file_list = glob.glob('IntParticle**.vtk')
    files = sorted(file_list, key=lambda x:x[-11:])
    #print(files)
    for j in files:
        motor_no = pd.read_csv(j, delim_whitespace=True)
        binding_motor = np.asscalar(motor_no.iloc[3:4,1:2].values) #get scalar value as str
        os.chdir(current_path)
        with open('bmotors.csv','a') as fd:
            writer = csv.writer(fd)
            writer.writerow([binding_motor])
        os.chdir(i)
    os.chdir(current_path)
    bmotors = pd.read_csv('bmotors.csv', names=['bm'])
    bmotors = bmotors.iloc[0::10,:]
    bmotors = bmotors.drop(bmotors.index[30])
    bmtrs = bmotors.reset_index(drop=True)
    os.chdir(i)
    conf = pd.read_csv('Conformation_A001.txt', names=conf_names, delim_whitespace=True)
    conf = conf.iloc[0::13,:] # jump 13 starting from the first index
    conf = conf.reset_index(drop=True) # reset index to have a fresh dataframe for leading tip
    conf = conf.iloc[0::10, :] # jump every 10 -- coz recommended stime is 0.01*10 sec
    conf = conf[['x','y']]
    Dx_tip = np.diff(conf['x']); Dy_tip = np.diff(conf['y'])
    DD=np.sqrt((Dx_tip**2)+(Dy_tip**2))
    v=DD/(10*dt); Av_vel = np.mean(v)
    vSD=np.sum(((v-Av_vel)**2)/(np.size(v)-1)); vSD=np.sqrt(vSD)
    av_sd = [str(Av_vel),str(vSD)]
    v1 = v.reshape(v.shape[0],1) # v.shape[0] to generalize; rather than just 30 for 3 sec simul
    #print(av_sd)
    v1 = pd.DataFrame(v1)
    #----------------------------
    df_load = pd.read_csv('TipXY_A001.txt', names=columnslt, delim_whitespace=True)
    df_nice = df_load.drop(['stime'], axis=1)
    dflt = df_nice.iloc[0::10, :] #pick every 9th row starting from the first
    Dx_tip = np.diff(dflt['x_tip']); Dy_tip = np.diff(dflt['y_tip'])
    DDlt=np.sqrt((Dx_tip**2)+(Dy_tip**2))
    vlt=DDlt/(10*dt); Av_vellt = np.mean(vlt)
    vSDlt=np.sum(((vlt-Av_vellt)**2)/(np.size(vlt)-1)); vSDlt=np.sqrt(vSDlt)
    xmax_ = df_load['x_tip'].max(); xmin_ = df_load['x_tip'].min()
    ymax_ = df_load['y_tip'].max(); ymin_ = df_load['y_tip'].min()
    x1_ = (0.01*(xmax_-xmin_))+xmin_
    y1_ = (0.90*(ymax_-ymin_))+ymin_
    x2_ = (0.15*(xmax_-xmin_))+xmin_ 
    y2_ = (0.85*(ymax_-ymin_))+ymin_
    plt.figure(figsize=(10,6))
    plt.plot(dflt['x_tip'],dflt['y_tip'], label='Leading tip', color='green', marker='.', linestyle='solid') #, linewidth=2, markersize=7
    plt.text(x1_, y1_, 'avSpeed = '); plt.text(x2_, y1_, r'%.5f $\mu m/sec$'%Av_vellt)
    plt.text(x1_, y2_, 'sDev = '); plt.text(x2_, y2_, '%.5f'%vSDlt)
    plt.xlabel('X', fontdict=font); plt.ylabel('Y', fontdict=font)
    plt.title('Actin Filament Leading Tip Movement'); plt.legend(loc='upper left')
    os.chdir(current_path)
    plt.savefig('leadingtip'+str(no)+'.svg',bbox_inches='tight', format='svg', dpi=500)
    plt.close()
    no = round(no+1,1)
    #----------------------------
    stime = np.around(np.linspace(0,3,30,endpoint=True), decimals=1)
    tym = pd.DataFrame(stime)
    df = pd.concat([tym, bmtrs, v1], axis=1)
    df.columns=['stime','b_m','vel']
    with open(filename+'.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow([i])
        writer.writerow(rownames)
    with open(filename+'.csv','a') as myfile:
        df.to_csv(myfile,sep=",",float_format='%.14f',header=False,index=False)
    with open(filename+'.csv','a') as fd:
        writer = csv.writer(fd)
        writer.writerow(avsd)
        writer.writerow(av_sd)
    #try:
        #os.remove('bmotors.csv') # if program run twice, remove the previously saved file.
    #except (OSError, RuntimeError, TypeError, NameError):
        #pass

toc = time.time(); ttime1 = round(toc-tic,4)
print(' ... I took %s sec. for that :( ...'%ttime1)
tic = time.time()
print('\n ==> Now sorting data ...\n')
#==========================prepare data for analysis=================================
df = pd.read_csv('grinding_data.csv',skiprows=[0],header=None)
sumdata = df.iloc[32::34,:]
sumdata.columns = ['avSpeed','sdev','nan1']
sumdata = sumdata.reset_index(drop=True)
sumdata = sumdata.drop(['nan1'],axis=1)
r = pd.DataFrame(np.around(np.linspace(0.1,1.0,10,endpoint=True), decimals=1))
dfsum = pd.concat([r,sumdata],axis=1)
dfsum.columns = ['r','avSpeed','sdev']
dfsum.to_csv('summarydata.csv',index=False)

bmd1 = df.iloc[1:31]
bmd1.columns = ['tym','bm','vel']
bmd1.to_csv('bmd1.csv',index=False)

bmd2 = df.iloc[35:65] #jump 4 from where you left, +30
bmd2.columns = ['tym','bm','vel']
bmd2.to_csv('bmd2.csv',index=False)

bmd3 = df.iloc[69:99]
bmd3.columns = ['tym','bm','vel']
bmd3.to_csv('bmd3.csv',index=False)

bmd4 = df.iloc[103:133]
bmd4.columns = ['tym','bm','vel']
bmd4.to_csv('bmd4.csv',index=False)

bmd5 = df.iloc[137:167]
bmd5.columns = ['tym','bm','vel']
bmd5.to_csv('bmd5.csv',index=False)

bmd6 = df.iloc[171:201]
bmd6.columns = ['tym','bm','vel']
bmd6.to_csv('bmd6.csv',index=False)

bmd7 = df.iloc[205:235]
bmd7.columns = ['tym','bm','vel']
bmd7.to_csv('bmd7.csv',index=False)

bmd8 = df.iloc[239:269]
bmd8.columns = ['tym','bm','vel']
bmd8.to_csv('bmd8.csv',index=False)

bmd9 = df.iloc[273:303]
bmd9.columns = ['tym','bm','vel']
bmd9.to_csv('bmd9.csv',index=False)

bmd10 = df.iloc[307:337]
bmd10.columns = ['tym','bm','vel']
bmd10.to_csv('bmd10.csv',index=False)

toc = time.time(); ttime2 = round(toc-tic,4)
print(' ... I took %s sec. for that :) ...'%ttime2)
tic = time.time()
print('\n ==> Now analysing data ...\n')
#==========================do analysis=======================================
columns=['r','avSpeed','sdev']
sumdat = pd.read_csv('summarydata.csv',skiprows=[0],header=None,names=columns)
x=sumdat['r']; y=sumdat['avSpeed']; dev=sumdat['sdev']
plt.figure(figsize=(10,6), dpi=500)
plt.errorbar(x,y,yerr=dev, ecolor='b', mec='blue', color='green', marker='D') #, capsize=7, linewidth=2, markersize=7
plt.xlabel('Active motor ratio', fontdict=font); plt.ylabel('Speed ($\mu m/sec$)', fontdict=font)
plt.title('Actin Filament Grinding Movement [ATP = 2000, MD = 3000]'); plt.legend(loc='upper left')
plt.savefig('summary.svg',bbox_inches='tight', format='svg',dip=500)
plt.close()

columns=['tym','bm','vel']
bmd1 = pd.read_csv('bmd1.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd1['tym'],bmd1['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd1['tym'],bmd1['vel'], 'g', label='Change in speed')
ax2.axhline(bmd1['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd1.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd2 = pd.read_csv('bmd2.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd2['tym'],bmd2['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd2['tym'],bmd2['vel'], 'g', label='Change in speed')
ax2.axhline(bmd2['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd2.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd3 = pd.read_csv('bmd3.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd3['tym'],bmd3['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd3['tym'],bmd3['vel'], 'g', label='Change in speed')
ax2.axhline(bmd3['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd3.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd4 = pd.read_csv('bmd4.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd4['tym'],bmd4['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd4['tym'],bmd4['vel'], 'g', label='Change in speed')
ax2.axhline(bmd4['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd4.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd5 = pd.read_csv('bmd5.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd5['tym'],bmd5['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd5['tym'],bmd5['vel'], 'g', label='Change in speed')
ax2.axhline(bmd5['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd5.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd6 = pd.read_csv('bmd6.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd6['tym'],bmd6['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd6['tym'],bmd6['vel'], 'g', label='Change in speed')
ax2.axhline(bmd6['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd6.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd7 = pd.read_csv('bmd7.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd7['tym'],bmd7['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd7['tym'],bmd7['vel'], 'g', label='Change in speed')
ax2.axhline(bmd7['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd7.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd8 = pd.read_csv('bmd8.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd8['tym'],bmd8['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd8['tym'],bmd8['vel'], 'g', label='Change in speed')
ax2.axhline(bmd8['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd8.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd9 = pd.read_csv('bmd9.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd9['tym'],bmd9['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd9['tym'],bmd9['vel'], 'g', label='Change in speed')
ax2.axhline(bmd9['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd9.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd10 = pd.read_csv('bmd10.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd10['tym'],bmd10['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Speed ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd10['tym'],bmd10['vel'], 'g', label='Change in speed')
ax2.axhline(bmd10['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd10.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

toc = time.time(); ttime3 = round(toc-tic,4)
print(' ... I took %s sec. for that ;) ...'%ttime3)
ttime = round(ttime1+ttime2+ttime3,4)
print('\n ==> I think Im done. %s sec. in total!\n'%ttime)

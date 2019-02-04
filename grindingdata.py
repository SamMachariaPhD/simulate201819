#After simulation, get inside each folder.
#Extract the binding motor number.
#Extract speed and its deviation.
#Sort data and do instructed analysis
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
    os.remove('summarydata.csv'); os.remove('bmd01.csv'); os.remove('bmd02.csv')
    os.remove('bmd03.csv'); os.remove('bmd04.csv'); os.remove('bmd05.csv'); os.remove('bmd06.csv')
    os.remove('bmd07.csv'); os.remove('bmd08.csv'); os.remove('bmd09.csv'); os.remove('bmd10.csv')
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
    v1 = v.reshape(30,1)
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
    plt.figure(figsize=(10,6))
    plt.plot(dflt['x_tip'],dflt['y_tip'], label='Leading tip', color='green', marker='o', linestyle='dashed') #, linewidth=2, markersize=7
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
    try:
        os.remove('bmotors.csv') # if program run twice, remove the previously saved file.
    except (OSError, RuntimeError, TypeError, NameError):
        pass

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

bmd01 = df.iloc[1:31]
bmd01.columns = ['tym','bm','vel']
bmd01.to_csv('bmd01.csv',index=False)

bmd02 = df.iloc[35:65] #jump 4 from where you left, +30
bmd02.columns = ['tym','bm','vel']
bmd02.to_csv('bmd02.csv',index=False)

bmd03 = df.iloc[69:99]
bmd03.columns = ['tym','bm','vel']
bmd03.to_csv('bmd03.csv',index=False)

bmd04 = df.iloc[103:133]
bmd04.columns = ['tym','bm','vel']
bmd04.to_csv('bmd04.csv',index=False)

bmd05 = df.iloc[137:167]
bmd05.columns = ['tym','bm','vel']
bmd05.to_csv('bmd05.csv',index=False)

bmd06 = df.iloc[171:201]
bmd06.columns = ['tym','bm','vel']
bmd06.to_csv('bmd06.csv',index=False)

bmd07 = df.iloc[205:235]
bmd07.columns = ['tym','bm','vel']
bmd07.to_csv('bmd07.csv',index=False)

bmd08 = df.iloc[239:269]
bmd08.columns = ['tym','bm','vel']
bmd08.to_csv('bmd08.csv',index=False)

bmd09 = df.iloc[273:303]
bmd09.columns = ['tym','bm','vel']
bmd09.to_csv('bmd09.csv',index=False)

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
bmd01 = pd.read_csv('bmd01.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd01['tym'],bmd01['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd01['tym'],bmd01['vel'], 'g', label='Change in speed')
ax2.axhline(bmd01['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd01.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd02 = pd.read_csv('bmd02.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd02['tym'],bmd02['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd02['tym'],bmd02['vel'], 'g', label='Change in speed')
ax2.axhline(bmd02['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd02.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd03 = pd.read_csv('bmd03.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd03['tym'],bmd03['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd03['tym'],bmd03['vel'], 'g', label='Change in speed')
ax2.axhline(bmd03['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd03.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd04 = pd.read_csv('bmd04.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd04['tym'],bmd04['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd04['tym'],bmd04['vel'], 'g', label='Change in speed')
ax2.axhline(bmd04['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd04.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd05 = pd.read_csv('bmd05.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd05['tym'],bmd05['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd05['tym'],bmd05['vel'], 'g', label='Change in speed')
ax2.axhline(bmd05['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd05.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd06 = pd.read_csv('bmd06.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd06['tym'],bmd06['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd06['tym'],bmd06['vel'], 'g', label='Change in speed')
ax2.axhline(bmd06['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd06.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd07 = pd.read_csv('bmd07.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd07['tym'],bmd07['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd07['tym'],bmd07['vel'], 'g', label='Change in speed')
ax2.axhline(bmd07['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd07.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd08 = pd.read_csv('bmd08.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd08['tym'],bmd08['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd08['tym'],bmd08['vel'], 'g', label='Change in speed')
ax2.axhline(bmd08['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd08.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd09 = pd.read_csv('bmd09.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd09['tym'],bmd09['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
ax2.plot(bmd09['tym'],bmd09['vel'], 'g', label='Change in speed')
ax2.axhline(bmd09['vel'].mean(), linestyle='--', color='g', label='Mean speed')
ax2.legend(loc='upper right')
ax2.tick_params('y', colors='g'); ax2.set_ylim(bottom=0)
plt.title('Binding motors vs. average speed', fontdict=font)
plt.savefig('bmd09.svg',bbox_inches='tight', format='svg',dip=300)
plt.close()

bmd10 = pd.read_csv('bmd10.csv',skiprows=[0],header=None,names=columns)
fig, ax1 = plt.subplots(figsize=(10,6), dpi=500)
ax1.set_xlabel('Time (sec)', fontdict=font)
ax1.set_ylabel('Number of binding motors', fontdict=font)
ax1.plot(bmd10['tym'],bmd10['bm'], label='Change in binding motor')
ax1.legend(loc='upper left'); ax1.set_ylim(bottom=0)
ax2 = ax1.twinx()
ax2.set_xlabel('Time (sec)', fontdict=font)
ax2.set_ylabel('Velocity ($\mu m/sec$)', fontdict=font)
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
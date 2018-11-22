import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

attach_lim = 0.12 # max z for attached filament
dt = 0.01 # change in time between rows
conf_names = ['time','x','y','z']

conf1 = pd.read_csv('Conformation_A001_r0.1.txt', names=conf_names, delim_whitespace=True)
#conf01 = conf1.drop(columns=['beads'])
conf01 = conf1.iloc[0::13,:]
conf01_df = conf01.reset_index(drop=True)
print("Conf0.1 shape before filter: ", conf01.shape)
df01 = conf01_df.loc[conf01_df['z'] <= attach_lim]
#conf01_df.to_csv('testconf.csv', encoding='utf-8', index=False)
print("Conf0.1 shape after filter: ", df01.shape)

conf2 = pd.read_csv('Conformation_A001_r0.2.txt', names=conf_names, delim_whitespace=True)
conf02 = conf2.iloc[0::13,:]
conf02_df = conf02.reset_index(drop=True)
print("Conf0.2 shape before filter: ", conf02.shape)
df02 = conf02_df.loc[conf02_df['z'] <= attach_lim]
print("Conf0.2 shape after filter: ", df02.shape)

conf3 = pd.read_csv('Conformation_A001_r0.3.txt', names=conf_names, delim_whitespace=True)
conf03 = conf3.iloc[0::13,:]
conf03_df = conf03.reset_index(drop=True)
print("Conf0.3 shape before filter: ", conf03.shape)
df03 = conf03_df.loc[conf03_df['z'] <= attach_lim]
print("Conf0.3 shape after filter: ", df03.shape)

conf4 = pd.read_csv('Conformation_A001_r0.4.txt', names=conf_names, delim_whitespace=True)
conf04 = conf4.iloc[0::13,:]
conf04_df = conf04.reset_index(drop=True)
print("Conf0.4 shape before filter: ", conf04.shape)
df04 = conf04_df.loc[conf04_df['z'] <= attach_lim]
print("Conf0.4 shape after filter: ", df04.shape)

conf5 = pd.read_csv('Conformation_A001_r0.5.txt', names=conf_names, delim_whitespace=True)
conf05 = conf5.iloc[0::13,:]
conf05_df = conf05.reset_index(drop=True)
print("Conf0.5 shape before filter: ", conf05.shape)
df05 = conf05_df.loc[conf05_df['z'] <= attach_lim]
print("Conf0.5 shape after filter: ", df05.shape)

conf6 = pd.read_csv('Conformation_A001_r0.6.txt', names=conf_names, delim_whitespace=True)
conf06 = conf6.iloc[0::13,:]
conf06_df = conf06.reset_index(drop=True)
print("Conf0.6 shape before filter: ", conf06.shape)
df06 = conf06_df.loc[conf06_df['z'] <= attach_lim]
print("Conf0.6 shape after filter: ", df06.shape)

conf7 = pd.read_csv('Conformation_A001_r0.7.txt', names=conf_names, delim_whitespace=True)
conf07 = conf7.iloc[0::13,:]
conf07_df = conf07.reset_index(drop=True)
print("Conf0.7 shape before filter: ", conf07.shape)
df07 = conf07_df.loc[conf07_df['z'] <= attach_lim]
print("Conf0.7 shape after filter: ", df07.shape)

conf8 = pd.read_csv('Conformation_A001_r0.8.txt', names=conf_names, delim_whitespace=True)
conf08 = conf8.iloc[0::13,:]
conf08_df = conf08.reset_index(drop=True)
print("Conf0.8 shape before filter: ", conf08.shape)
df08 = conf08_df.loc[conf08_df['z'] <= attach_lim]
print("Conf0.8 shape after filter: ", df08.shape)

conf9 = pd.read_csv('Conformation_A001_r0.9.txt', names=conf_names, delim_whitespace=True)
conf09 = conf9.iloc[0::13,:]
conf09_df = conf09.reset_index(drop=True)
print("Conf0.9 shape before filter: ", conf09.shape)
df09 = conf09_df.loc[conf09_df['z'] <= attach_lim]
print("Conf0.9 shape after filter: ", df09.shape)

conf10 = pd.read_csv('Conformation_A001_r1.0.txt', names=conf_names, delim_whitespace=True)
conf10 = conf10.iloc[0::13,:]
conf10_df = conf10.reset_index(drop=True)
print("Conf1.0 shape before filter: ", conf10.shape)
df10 = conf10_df.loc[conf10_df['z'] <= attach_lim]
print("Conf1.0 shape after filter: ", df10.shape)
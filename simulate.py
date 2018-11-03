import fileinput, sys, shutil, os, time, socket, subprocess

pc_hostname = socket.gethostname()
date_today = time.strftime("%d-%m-%Y")
dir_name = date_today+pc_hostname

current_path = os.getcwd()  
print ("=> The current working directory is: %s" % current_path) 

try:  
    os.mkdir(dir_name)
except OSError:  
    print ("=> Creation of the directory: %s failed" % dir_name)
else:  
    print ("=> Successfully created the directory: %s " % dir_name)

#run_range = range(1,6)
#runs = int(input("Enter the number of Runs for today (1-5): "))
#if runs in run_range:
#    print ("=> Number of experimental runs: %s " % runs)
#else:
#    print ("=> Please enter a number between 1 and 5.")

#counter = 1

print("Set parameters for Run 1")
print("Enter: species_ratio, beads_number, ATP_value, MD_value: ")
species_ratio_run1, beads_number_run1, ATP_value_run1, MD_value_run1 = map(float,input().split(','))
dir_1 = 'R'+str(species_ratio_run1)+'B'+str(beads_number_run1)+'ATP'+str(ATP_value_run1)+'MD'+str(MD_value_run1)

print("Set parameters for Run 2")
print("Enter: species_ratio, beads_number, ATP_value, MD_value: ")
species_ratio_run2, beads_number_run2, ATP_value_run2, MD_value_run2 = map(float,input().split(','))
dir_2 = 'R'+str(species_ratio_run2)+'B'+str(beads_number_run2)+'ATP'+str(ATP_value_run2)+'MD'+str(MD_value_run2)

print("Set parameters for Run 3")
print("Enter: species_ratio, beads_number, ATP_value, MD_value: ")
species_ratio_run3, beads_number_run3, ATP_value_run3, MD_value_run3 = map(float,input().split(','))
dir_3 = 'R'+str(species_ratio_run3)+'B'+str(beads_number_run3)+'ATP'+str(ATP_value_run3)+'MD'+str(MD_value_run3)

print("Set parameters for Run 4")
print("Enter: species_ratio, beads_number, ATP_value, MD_value: ")
species_ratio_run4, beads_number_run4, ATP_value_run4, MD_value_run4 = map(float,input().split(','))
dir_4 = 'R'+str(species_ratio_run4)+'B'+str(beads_number_run4)+'ATP'+str(ATP_value_run4)+'MD'+str(MD_value_run4)

print("Set parameters for Run 5")
print("Enter: species_ratio, beads_number, ATP_value, MD_value: ")
species_ratio_run5, beads_number_run5, ATP_value_run5, MD_value_run5 = map(float,input().split(','))
dir_5 = 'R'+str(species_ratio_run5)+'B'+str(beads_number_run5)+'ATP'+str(ATP_value_run5)+'MD'+str(MD_value_run5)

os.chdir(dir_name)

os.mkdir(dir_1)
os.mkdir(dir_2)
os.mkdir(dir_3)
os.mkdir(dir_4)
os.mkdir(dir_5)

os.chdir(current_path)

files = ['MotilityAssayActin2MotorsParameters_v3.f90', 'MotilityAssayConfinements_v1.f90', 'MotilityAssayForceForceFunctions_v2.f90', 'MotilityAssaySubstrateDeformation_v2.f90', 'mt.f90', 'MotilityAssayActin2MotorsMain_v6.f90', 'autocomp_py.sh', 'autorun_py.sh']
for f in files:
    shutil.copy(f, dir_name+'/'+dir_1)
    shutil.copy(f, dir_name+'/'+dir_2)
    shutil.copy(f, dir_name+'/'+dir_3)
    shutil.copy(f, dir_name+'/'+dir_4)
    shutil.copy(f, dir_name+'/'+dir_5)
print ("=> Files successfully copied: %s " % files)

os.chdir(dir_name+'/'+dir_1)
param_file = "MotilityAssayActin2MotorsParameters_v3.f90"
# replace all occurrences of 'ATP = 3000.0_DP' with 'ATP = 1000.0_DP' and insert a line after the 5th
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Species1Ratio = 0.60', 'Species1Ratio = '+str(species_ratio_run1)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('NumBeads = 13', 'NumBeads = '+str(beads_number_run1)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('ATP = 2000.0', 'ATP = '+str(ATP_value_run1)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Motor_Density = 3000.0', 'Motor_Density = '+str(MD_value_run1)))
    # replace 'ATP = 3000.0_DP' and write
    #if i == 4: sys.stdout.write('\n')  # write a blank line after the 5th line
time.sleep(1)
print ("=> Files successfully updated: %s " % dir_1)
subprocess.call("./autocomp_py.sh", shell=True)
time.sleep(1)
print ("=> Program successfully compiled: %s " % param_file)
subprocess.call("./autorun_py.sh", shell=True)
print ("=> Program successfully run. Done: %s " % param_file)
os.chdir(current_path)

os.chdir(dir_name+'/'+dir_2)
param_file = "MotilityAssayActin2MotorsParameters_v3.f90"
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Species1Ratio = 0.60', 'Species1Ratio = '+str(species_ratio_run2)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('NumBeads = 13', 'NumBeads = '+str(beads_number_run2)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('ATP = 2000.0', 'ATP = '+str(ATP_value_run2)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Motor_Density = 3000.0', 'Motor_Density = '+str(MD_value_run2)))
time.sleep(1)
print ("=> Files successfully updated: %s " % dir_2)
subprocess.call("./autocomp_py.sh", shell=True)
time.sleep(1)
print ("=> Program successfully compiled: %s " % param_file)
subprocess.call("./autorun_py.sh", shell=True)
print ("=> Program successfully run. Done: %s " % param_file)
os.chdir(current_path)

os.chdir(dir_name+'/'+dir_3)
param_file = "MotilityAssayActin2MotorsParameters_v3.f90"
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Species1Ratio = 0.60', 'Species1Ratio = '+str(species_ratio_run3)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('NumBeads = 13', 'NumBeads = '+str(beads_number_run3)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('ATP = 2000.0', 'ATP = '+str(ATP_value_run3)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Motor_Density = 3000.0', 'Motor_Density = '+str(MD_value_run3)))
time.sleep(1)
print ("=> Files successfully updated: %s " % dir_3)
subprocess.call("./autocomp_py.sh", shell=True)
time.sleep(1)
print ("=> Program successfully compiled: %s " % param_file)
subprocess.call("./autorun_py.sh", shell=True)
print ("=> Program successfully run. Done: %s " % param_file)
os.chdir(current_path)

os.chdir(dir_name+'/'+dir_4)
param_file = "MotilityAssayActin2MotorsParameters_v3.f90"
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Species1Ratio = 0.60', 'Species1Ratio = '+str(species_ratio_run4)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('NumBeads = 13', 'NumBeads = '+str(beads_number_run4)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('ATP = 2000.0', 'ATP = '+str(ATP_value_run4)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Motor_Density = 3000.0', 'Motor_Density = '+str(MD_value_run4)))
time.sleep(1)
print ("=> Files successfully updated: %s " % dir_4)
subprocess.call("./autocomp_py.sh", shell=True)
time.sleep(1)
print ("=> Program successfully compiled: %s " % param_file)
subprocess.call("./autorun_py.sh", shell=True)
print ("=> Program successfully run. Done: %s " % param_file)
os.chdir(current_path)

os.chdir(dir_name+'/'+dir_5)
param_file = "MotilityAssayActin2MotorsParameters_v3.f90"
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Species1Ratio = 0.60', 'Species1Ratio = '+str(species_ratio_run5)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('NumBeads = 13', 'NumBeads = '+str(beads_number_run5)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('ATP = 2000.0', 'ATP = '+str(ATP_value_run5)))
for i, line in enumerate(fileinput.input(param_file, inplace=1)):
    sys.stdout.write(line.replace('Motor_Density = 3000.0', 'Motor_Density = '+str(MD_value_run5)))
time.sleep(1)
print ("=> Files successfully updated: %s " % dir_5)
subprocess.call("./autocomp_py.sh", shell=True)
time.sleep(1)
print ("=> Program successfully compiled: %s " % param_file)
subprocess.call("./autorun_py.sh", shell=True)
print ("=> Program successfully run. Done: %s " % param_file)
os.chdir(current_path)


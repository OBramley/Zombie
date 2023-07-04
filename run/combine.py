###################################################################################################
#
# This is the run file for the zombie states program. It checks run paramters are correct. It then
# compiles the Fortran modules; creates an execution folder and copies the relevant files to it.
# The program then begins the run requesting nodes to allow paralel execution if on the the HPC
# 
###################################################################################################
import sys

import inputs
import os
import socket
from glob import glob
import getpass
import shutil
import subprocess
import random
import math 
from pyscf import gto, scf, ao2mo
from functools import reduce
import numpy
import csv

# Checking input paramaters 
if(isinstance(inputs.run['nodes'],int)==False):
    sys.exit("Number of folders must be an integer")
elif(isinstance(inputs.run['cores'],int)==False):
    sys.exit("Number of parallel cores must be an integer")
elif(inputs.run['nodes']<1):
    sys.exit("Not enough nodes selected. Must be 1 or greater")
elif(inputs.run['nodes']>100):
    sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
elif(inputs.run['cores']>40):
    sys.exit("Too many cores selected. Maximum of 40 available")
elif(inputs.run['cores']<1):
    sys.exit("Not enough cores selected. Must be 1 or greater")
elif(inputs.zombs['norb']<1):
    sys.exit("Not enough orbitals. Must be 1 or greater")
elif(inputs.zombs['zomtyp'] not in {'ran','HF','bb','hf'}):
    sys.exit("Type of zombie state must be ran, HF or bb")
elif(inputs.run['elecs'] not in {'pyscf','mol','no'}):
    sys.exit("one and two electron integral paramter must be 'pyscf', 'mol' or 'no'")
elif(inputs.run['zomgen'] not in {'y','n'}):
    sys.exit("Generation of a new set of zombie states must be either 'y' or 'n'")
elif(inputs.run['imagprop'] not in {'y','n'}):
    sys.exit("The imaginary time propagation parameter has been inccorectly set must be either 'y' or 'n'")
elif(isinstance(inputs.run['beta'],int)==False):
    sys.exit("Beta must be an integer")
elif(isinstance(inputs.run['timesteps'],int)==False):
    sys.exit("Beta must be an integer")
elif(inputs.run['gram'] not in {'y','n'}):
    sys.exit("The Gram Schmidt pramater must be either 'y' or 'n'")
elif(inputs.run['zomhf'] not in {'y','n'}):
    sys.exit("Setting the first zombie state as a RHF determinant must be either 'y' or 'n'")
elif(isinstance(inputs.run['hfnum'],int)==False):
    sys.exit("Number of the RHF determinant must be an integer")
elif(isinstance(inputs.run['gramnum'],int)==False):
    sys.exit("Number of states to be investiated must be an integer")
elif(inputs.run['clean'] not in {'y','n','f'}):
    sys.exit("Setting cleaning on or off must be 'y' or 'n' or 'f' (to clean from a previously generated file)")


print("Arguments checked")
# Check if on HPC
Hostname=socket.gethostname()
if((Hostname==("login2.arc4.leeds.ac.uk"))or(Hostname==("login1.arc4.leeds.ac.uk"))or(Hostname==("login2.arc3.leeds.ac.uk"))or(Hostname==("login1.arc3.leeds.ac.uk"))):
    HPCFLG=1
else:
    HPCFLG=0

# Get current directory
EXDIR=os.getcwd()

if os.path.exists("rundata.csv"):
  os.remove("rundata.csv")
else:
  print("rundata.csv file does not exist")

multflg=inputs.run['multiple']
shutil.copy2(EXDIR+'/node_1/data/zombie_1001.csv',EXDIR+'/data')
cnt=2
for i in range(multflg):
    floc=EXDIR+'/node_'+str(i+1)+'/data/zombie_1'
    for j in range(2,inputs.zombs['ndet']+1):
        shutil.copy2(floc+"{:03d}".format(j)+".csv",EXDIR+'/data/zombie_1'+"{:03d}".format(cnt)+".csv")
        cnt+=1

cnt=cnt-1

zomstat='n'
hamstat='y'
cleanstat='n'

if(inputs.run['grad']=='n'):
    if(inputs.run['hamgen']=='y'):
        if os.path.exists('data/ham.csv'):
            if os.path.exists('data/ovlp.csv'):
                hamstat='n'
            else:
                os.remove("data/ham.csv")
                hamstat='y'
        else: 
            if os.path.exists('data/ovrlp.csv'):
                os.remove("data/ovrlp.csv")
                hamstat='y'
            else:
                os.remove("data/ham.csv")
                hamstat='y'

 
with open("rundata.csv",'w',newline='')as file:
    writer = csv.writer(file)
    writer.writerow([zomstat,hamstat,inputs.run['imagprop'],inputs.run['beta'],inputs.run['timesteps'],inputs.run['clean'],inputs.run['gram'],inputs.run['gramnum'],inputs.run['grad'],'n'])
    writer.writerow([inputs.zombs['norb'],inputs.zombs['nel'],inputs.zombs['spin'],cnt,inputs.zombs['zomtyp'],inputs.zombs['rhf_1'],inputs.zombs['imagflg']])
    writer.writerow([inputs.run['hamfile'],inputs.run['ovrlfile'],inputs.run['cleanham']])
   
os.environ["OMP_CANCELLATION"]="TRUE" 
if(HPCFLG==1):
    if(inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
    number=random.randint(99999,1000000)
    file1="zombie"+str(number)+".sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    if(inputs.run['cores']!=1):
        f.write("#$ -pe smp "+str(inputs.run['cores'])+" \n") #Use shared memory parallel environemnt 
    f.write("#$ -l h_rt="+inputs.run['runtime']+"\n")
    f.write("#$ -l h_vmem=5G \n")
    f.write("export OMP_CANCELLATION=true \n")
    f.write("module add mkl \n")
    # f.write('time ./d_check.exe')
    f.write('time ./ZOMBIE.exe')
    f.close()
    subprocess.call(['qsub',file1])
else:
    print(os.getcwd())
    if(inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
    subprocess.run(["./ZOMBIE.exe"])
        
        

    

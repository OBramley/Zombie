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
import getpass
import shutil
import subprocess
import random
from math import comb

# Checking input paramaters 
if(isinstance(inputs.run['nodes'],int)==False):
    sys.exit("Number of folders must be an integer")
elif(isinstance(inputs.run['cores'],int)==False):
    sys.exit("Number of parallel cores must be an integer")
elif(inputs.run['nodes']<1):
    sys.exit("Not enough nodes selected. Must be 1 or greater")
elif(inputs.run['nodes']>100):
    sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
elif(inputs.run['cores']>8):
    sys.exit("Too many cores selected. Maximum of 8 available")
elif(inputs.run['cores']<1):
    sys.exit("Not enough cores selected. Must be 1 or greater")
elif(inputs.zombs['norb']<1):
    sys.exit("Not enough orbitals. Must be 1 or greater")
elif(inputs.zombs['ndet']<2):
    sys.exit("Not enough zombie states. Must be 2 or greater")
elif((inputs.run['elecs']!='pyscf')or(inputs.run['elecs']!='mol')or(inputs.run['elecs']!='n')):
    sys.exit("one and two electron integral paramter must be 'pyscf', 'mol' or 'n'")
elif((inputs.run['zomgen']!='y')or(inputs.run['zomgen']!='n')):
    sys.exit("Generation of a new set of zombie states must be either 'y' or 'n'")
elif((inputs.run['imagprop']!='y')or(inputs.run['imagprop']!='n')):
    sys.exit("The imaginary time propagation parameter has been inccorectly set must be either 'y' or 'n'")
elif(isinstance(inputs.run['beta'],int)==False):
    sys.exit("Beta must be an integer")
elif(isinstance(inputs.run['timesteps'],int)==False):
    sys.exit("Beta must be an integer")
elif((inputs.run['gram']!='y')or(inputs.run['gram']!='n')):
    sys.exit("The Gram Schmidt pramater must be either 'y' or 'n'")
elif((inputs.run['zomhf']!='y')or(inputs.run['zomhf']!='n')):
    sys.exit("Setting the first zombie state as a RHF determinant must be either 'y' or 'n'")
elif(isinstance(inputs.run['gramnum'],int)==False):
    sys.exit("Number of states to be investiated must be an integer")

if(inputs.zombs['zomtyp']=='HF'):
    ndetcheck=0
    for i in range(inputs.zombs['norb']+1):
        ndetcheck=ndetcheck+comb(inputs.zombs['norb'],i)
    if(inputs.zombs['zomtyp']!=ndetcheck):
        sys.exit('A Hartree Fock Basis for',inputs.zombs['norb'], 'orbitals should have', ndetcheck, 'basis functions')

if(inputs.run['elecs']=='mol'):
    if not os.path.isfile('../'+ inputs.run['elecfile']):
        sys.exit("No Molpro input file")

if(inputs.run['elecs']=='n'and inputs.run['hamgen']=='y'):
    if not os.path.isfile('../'+ inputs.run['elecfile']):
        sys.exit("No electron integral input file")

if(inputs.run['elecs']=='n'and inputs.run['hamgen']=='n'):
    if not os.path.isfile('../'+ inputs.run['hamfile']):
        sys.exit("No Hamiltonian input file")

if(inputs.run['gram']=='y'):
    if(isinstance(inputs.run['gramnum'],int)==False):
        sys.exit("Number of states for Gram Schmidt must be an integer")
    elif(inputs.run['gramnum']<2):
        sys.exit("If using Gram Schmidt more than one state must investigated")

print("Arguments checked")
# Check if on HPC
Hostname=socket.gethostname()
if(Hostname==("login2.arc4.leeds.ac.uk")):
    HPCFLG=1
else:
    HPCFLG=0

# Make Execution folder
if(HPCFLG==0):
    if not os.path.exists("../EXEC"):
        os.mkdir("../EXEC")
    EXDIR="../EXEC"
else:    
    os.environ['LOGNAME']
    EXDIR="/nobackup/"+getpass.getuser()

# Check for run folder and make it
if os.path.exists(EXDIR+"/"+inputs.run['runfolder']):
    value=input("File already exists do you want to delete it? y/n\n")
    if(value=='y'):
        shutil.rmtree(EXDIR+"/"+inputs.run['runfolder'])
    else:
        sys.exit("Runfolder already exists. Change the Runfolder name or delete/move it")

os.mkdir(EXDIR+"/"+inputs.run['runfolder'])
EXDIR1=EXDIR+"/"+inputs.run['runfolder']
os.mkdir(EXDIR1="/outputs")



# Move the relevant files to the execution file

# Run the program
# 
#
# #If on a SGE machine make job submission file
if(HPCFLG==1):
    number=random.randint(99999,1000000)
    file1="MCE"+str(number)+".sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    if(inputs.run['cores']!=1):
        f.write("#$ -pe smp "+str(inputs.run['cores'])+" \n") #Use shared memory parallel environemnt 
    f.write("#$ -l h_rt=40:00:00 \n")
    f.write("#$ -l h_vmem=4G \n")
    f.write("#$ -t 1-"+str(inputs.run['nodes'])+" \n")
    f.write("date \n")
    f.write("cd "+EXDIR1+"/run-$SGE_TASK_ID/ \n")
    f.write("echo "'"Running on $HOSTNAME in folder $PWD" \n')
    f.write("time ./MCE.exe \n")
    f.write("date \n")
    f.close()
    if(inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
    subprocess.call(['qsub',file1])
else:
    if(inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
    for i in range(inputs.run['nodes']):
        SUBDIR=EXDIR1+"/run-"+str(i+1)
        subprocess.Popen('',executable=SUBDIR+"python main.py",cwd=SUBDIR) 
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
import math #import comb

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
elif(inputs.zombs['zomtyp'] not in {'ran','HF','bb'}):
    sys.exit("Type of zombie state must be ran, HF or bb")
elif(inputs.run['elecs'] not in {'pyscf','mol','n'}):
    sys.exit("one and two electron integral paramter must be 'pyscf', 'mol' or 'n'")
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


if(inputs.zombs['zomtyp']=='HF'):
    ndetcheck=0
    for i in range(inputs.zombs['norb']+1):
        ndetcheck=ndetcheck+math.comb(inputs.zombs['norb'],i)
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


# Move the relevant files to the execution file
shutil.copy2("inputs.py",EXDIR1)
shutil.copy2("../src/ham.py",EXDIR1)
shutil.copy2("../src/imgtp.py",EXDIR1)
shutil.copy2("../src/in_outputs.py",EXDIR1)
shutil.copy2("../src/main.py",EXDIR1)
shutil.copy2("../src/op.py",EXDIR1)
shutil.copy2("../src/zom.py",EXDIR1)
if(inputs.run['elecs']=='mol'):
    shutil.copy2('../'+ inputs.run['elecfile'],EXDIR1)
if(inputs.run['elecs']=='n'and inputs.run['hamgen']=='y'):
    shutil.copy2('../'+ inputs.run['elecfile'],EXDIR1)
if(inputs.run['elecs']=='n'and inputs.run['hamgen']=='n'):
    shutil.copy2('../'+ inputs.run['hamfile'],EXDIR1)

os.chdir(EXDIR1)
# Run the program
# If on a SGE machine make job submission file
if(HPCFLG==1):
    number=random.randint(99999,1000000)
    file1="zombie"+str(number)+".sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    f.write("#$ -l h_rt=00:30:00 \n")
    f.write("#$ -l h_vmem=2G \n")
    f.write('#$ -m be')#Get email at start and end of the job
    f.write('module load anaconda')
    f.write('source activate base')
    f.write('python main.py')
    f.close()
    subprocess.call(['qsub',file1])
else:
    # if(inputs.run['cores']!=1):
    #     os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
    # number=random.randint(99999,1000000)
    # file1="zombie"+str(number)+".sh"
    # f=open("../Run/"+file1,"w")
    # f.write("python " + EXDIR1+"/main.py")
    # f.close()
    # subprocess.run(['chmod', 'u+x', '../Run/zombie'+str(number)+'.sh'])
    subprocess.run(['python', 'main.py'])
    
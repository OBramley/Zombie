###################################################################################################
#
# This is the run file for the zombie states program. It checks run paramters are correct. It then
# compiles the Fortran modules; creates an execution folder and copies the relevant files to it.
# The program then begins the run requesting nodes to allow paralel execution if on the the HPC
# 
###################################################################################################
import sys

import d_check_inputs
import os
import socket
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
if(isinstance(d_check_inputs.run['nodes'],int)==False):
    sys.exit("Number of folders must be an integer")
elif(isinstance(d_check_inputs.run['cores'],int)==False):
    sys.exit("Number of parallel cores must be an integer")
elif(d_check_inputs.run['nodes']<1):
    sys.exit("Not enough nodes selected. Must be 1 or greater")
elif(d_check_inputs.run['nodes']>100):
    sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
elif(d_check_inputs.run['cores']>8):
    sys.exit("Too many cores selected. Maximum of 8 available")
elif(d_check_inputs.run['cores']<=1):
    sys.exit("Not enough cores selected. Must be 2 or greater as this program will take too long not in parallel")
elif(d_check_inputs.zombs['norb']<1):
    sys.exit("Not enough orbitals. Must be 1 or greater")
elif(d_check_inputs.zombs['ndet']<2):
    sys.exit("Not enough zombie states. Must be 2 or greater")
elif(d_check_inputs.zombs['zomtyp'] not in {'ran','HF','bb','hf'}):
    sys.exit("Type of zombie state must be ran, HF or bb")
elif(d_check_inputs.run['elecs'] not in {'pyscf','mol','no'}):
    sys.exit("one and two electron integral paramter must be 'pyscf', 'mol' or 'no'")
elif(isinstance(d_check_inputs.run['beta'],int)==False):
    sys.exit("Beta must be an integer")
elif(isinstance(d_check_inputs.run['timesteps'],int)==False):
    sys.exit("Beta must be an integer")

if(d_check_inputs.zombs['zomtyp']=='HF'):
    ndetcheck=0
    for i in range(d_check_inputs.zombs['norb']+1):
        ndetcheck=ndetcheck+math.comb(d_check_inputs.zombs['norb'],i)
    if(d_check_inputs.zombs['zomtyp']!=ndetcheck):
        sys.exit('A Hartree Fock Basis for',d_check_inputs.zombs['norb'], 'orbitals should have', ndetcheck, 'basis functions')

if(d_check_inputs.run['elecs']=='mol'):
    if not os.path.isfile('integrals/'+ d_check_inputs.run['elecfile']):
        sys.exit("No Molpro input file")

if(d_check_inputs.run['elecs']=='no'):
    if not os.path.exists("integrals"):
        sys.exit("No electron integral input file")

print("Arguments checked")
# Check if on HPC
Hostname=socket.gethostname()
if((Hostname==("login2.arc4.leeds.ac.uk"))or(Hostname==("login1.arc4.leeds.ac.uk"))):
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
if os.path.exists(EXDIR+"/"+d_check_inputs.run['runfolder']):
    value=input("File already exists do you want to delete it? y/n\n")
    if(value=='y'):
        shutil.rmtree(EXDIR+"/"+d_check_inputs.run['runfolder'])
    else:
        sys.exit("Runfolder already exists. Change the Runfolder name or delete/move it")

os.mkdir(EXDIR+"/"+d_check_inputs.run['runfolder'])
EXDIR1=EXDIR+"/"+d_check_inputs.run['runfolder']


    


shutil.copy2("d_check_inputs.py",EXDIR1)
if((d_check_inputs.run['elecs']=='no')):
    shutil.copytree('integrals',EXDIR1+"/integrals")
elif((d_check_inputs.run['elecs']=='mol')):
    os.mkdir(EXDIR1+"/integrals")
    shutil.copy2('../'+ d_check_inputs.run['elecfile'],EXDIR1+"/integrals")
if(d_check_inputs.run['elecs']=='pyscf'):
    os.mkdir(EXDIR1+"/integrals")
    mol = gto.M(
    unit = d_check_inputs.pyscf['units'],
    atom = d_check_inputs.pyscf['atoms'],
    basis = d_check_inputs.pyscf['bs'],
    verbose = d_check_inputs.pyscf['verbosity'],
    symmetry = d_check_inputs.pyscf['symmetry'],
    spin=d_check_inputs.pyscf['spin'],
    charge=d_check_inputs.pyscf['charge'],
    symmetry_subgroup = d_check_inputs.pyscf['symmetry_subgroup'], #0 is code for A1 point group
    )
    myhf=scf.RHF(mol)
    myhf.kernel()
    """Obtaining one and two electron integrals from pyscf calculation
    Code adapted from George Booth"""
    # Extract AO->MO transformation matrix
    c = myhf.mo_coeff
    # Get 1-electron integrals and convert to MO basis
    h1e = reduce(numpy.dot, (c.T, myhf.get_hcore(), c))
    # Get 2-electron integrals and transform them
    eri = ao2mo.kernel(mol, c)
    # Ignore all permutational symmetry, and write as four-index tensor, in chemical notation
    eri_full = ao2mo.restore(1, eri, c.shape[1])
    # Scalar nuclear repulsion energy
    Hnuc = myhf.energy_nuc()

    
    with open(EXDIR1+"/integrals/hnuc.csv",'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile)
        spamwriter.writerow([Hnuc,0])

    with open(EXDIR1+"/integrals/h1ea.csv",'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerows(h1e)

    for i in range(len(eri_full)):
        for j in range(len(eri_full)):
            obj=eri_full[i,j,:,:]
            with open(EXDIR1+"/integrals/h2ea_"+str(i+1)+"_"+str(j+1)+".csv",'w', newline='')as csvfile:
                spamwriter=csv.writer(csvfile, delimiter=',')
                spamwriter.writerows(obj)

os.mkdir(EXDIR1+'/data')

with open(EXDIR1+'/rundata.csv','w',newline='')as file:
    writer = csv.writer(file)
    writer.writerow(['y','y','y',d_check_inputs.run['beta'],d_check_inputs.run['timesteps'],'n','n',0])
    writer.writerow(d_check_inputs.zombs.values())


if(HPCFLG==1):
    shutil.copy2("makefile_dcheck_arc","Makefile")
    subprocess.run(["make"])
else:
    shutil.copy2("makefile_bbi","Makefile")
    # shutil.copy2("makefile_dcheck_mac","Makefile")
    subprocess.run(["make"])

# shutil.copy2("d_check.exe",EXDIR1)
# shutil.copy2("outputs.py",EXDIR1)
shutil.copy2("bbi_o.py",EXDIR1)

shutil.copy2("bbi.exe",EXDIR1)

os.chdir(EXDIR1)

if(HPCFLG==1):
    if(d_check_inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(d_check_inputs.run['cores'])
    number=random.randint(99999,1000000)
    file1="dcheck"+str(number)+".sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    if(d_check_inputs.run['cores']!=1):
        f.write("#$ -pe smp "+str(d_check_inputs.run['cores'])+" \n") #Use shared memory parallel environemnt 
    f.write("#$ -l h_rt="+d_check_inputs.run['runtime']+"\n")
    f.write("#$ -l h_vmem=12G \n")
    f.write("module add mkl \n")
    f.write('time ./d_check.exe')
    f.close()
    
    subprocess.call(['qsub',file1])
else:
    print(os.getcwd())
    if(d_check_inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(d_check_inputs.run['cores'])
    # subprocess.run(["./d_check.exe"])
    subprocess.run(["./bbi.exe"])
        
        

    

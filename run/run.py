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
import math 
from pyscf import gto, scf, ao2mo
from functools import reduce
import numpy
import csv
import gpu_code_write as gcw

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
# elif(inputs.zombs['ndet']<2):
    # sys.exit("Not enough zombie states. Must be 2 or greater")
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


if(inputs.zombs['zomtyp']=='HF'):
    ndetcheck=0
    for i in range((inputs.zombs['norb']*2)+1):
        ndetcheck=ndetcheck+math.comb(inputs.zombs['norb']*2,i)
    if(inputs.zombs['ndet']!=ndetcheck):
        sys.exit('A Hartree Fock Basis for',inputs.zombs['norb'], 'orbitals should have', ndetcheck, 'basis functions')

if(inputs.run['elecs']=='mol'):
    if not os.path.isfile('integrals/'+ inputs.run['elecfile']):
        sys.exit("No Molpro input file")

if(inputs.run['elecs']=='no'and inputs.run['hamgen']=='y'):
    if(inputs.run['language']=='python'):
        if not os.path.isfile('integrals/'+ inputs.run['elecfile']):
            sys.exit("No electron integral input file")
    elif(inputs.run['language']=='fortran'):
        if not os.path.exists("integrals"):
            sys.exit("No electron integral input file")

if(inputs.run['elecs']=='no'and inputs.run['hamgen']=='n'):
    if not os.path.isfile('data/'+ inputs.run['hamfile']):
        sys.exit("No Hamiltonian input file")

if(inputs.run['clean']=='f'):
    if not os.path.isfile('data/'+inputs.run['cleanham']):
        sys.exit("No clean hamiltonian")
    if(inputs.run['language']=='python'):
        if not os.path.isfile('../'+inputs.run['cleanzom']):
            sys.exit("No clean zombie states")

if(inputs.run['gram']=='y'):
    if(isinstance(inputs.run['gramnum'],int)==False):
        sys.exit("Number of states for Gram Schmidt must be an integer")
    elif(inputs.run['gramnum']<2):
        sys.exit("If using Gram Schmidt more than one state must investigated")


print("Arguments checked")
# Check if on HPC
Hostname=socket.gethostname()
if((Hostname==("login2.arc4.leeds.ac.uk"))or(Hostname==("login1.arc4.leeds.ac.uk"))or(Hostname==("login2.arc3.leeds.ac.uk"))or(Hostname==("login1.arc3.leeds.ac.uk"))):
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

if(inputs.run['language']=="python"):
    # Move the relevant files to the execution file
    shutil.copy2("inputs.py",EXDIR1)
    shutil.copy2("../src/ham.py",EXDIR1)
    shutil.copy2("../src/imgtp.py",EXDIR1)
    shutil.copy2("../src/in_outputs.py",EXDIR1)
    shutil.copy2("../src/main.py",EXDIR1)
    shutil.copy2("../src/op.py",EXDIR1)
    shutil.copy2("../src/zom.py",EXDIR1)
    shutil.copy2("../src/cleaning.py",EXDIR1)
    if(inputs.run['elecs']=='mol'):
        shutil.copy2('../'+ inputs.run['elecfile'],EXDIR1)
    if(inputs.run['elecs']=='n'and inputs.run['hamgen']=='y'):
        shutil.copy2('../'+ inputs.run['elecfile'],EXDIR1)
    if(inputs.run['elecs']=='n'and inputs.run['hamgen']=='n'):
        shutil.copy2('../'+ inputs.run['hamfile'],EXDIR1)
    if(inputs.run['clean']=='f'):
        shutil.copy2('../'+ inputs.run['cleanham'],EXDIR1)
        shutil.copy2('../'+ inputs.run['cleanzom'],EXDIR1)

    os.chdir(EXDIR1)
    # Run the program
    # If on a SGE machine make job submission file
    if(HPCFLG==1):
        number=random.randint(99999,1000000)
        file1="zombie"+str(number)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        f.write("#$ -l h_rt="+inputs.run['runtime']+"\n")
        f.write("#$ -l h_vmem=2G \n")
        # f.write('#$ -m be')#Get email at start and end of the job
        f.write('module load anaconda\n')
        f.write('source activate base\n')
        f.write('python main.py')
        f.close()
        subprocess.call(['qsub',file1])
    else:
        subprocess.run(['python', 'main.py'])

elif(inputs.run['language']=="fortran"):
    shutil.copy2("inputs.py",EXDIR1)
    shutil.copy2("graph.py",EXDIR1)
    shutil.copy2("restart.py",EXDIR1)
    if((inputs.run['elecs']=='no')):
        shutil.copytree('integrals',EXDIR1+"/integrals")
    elif((inputs.run['elecs']=='mol')):
        os.mkdir(EXDIR1+"/integrals")
        shutil.copy2('../'+ inputs.run['elecfile'],EXDIR1+"/integrals")
    if(inputs.run['elecs']=='pyscf'):
        os.mkdir(EXDIR1+"/integrals")
        mol = gto.M(
        unit = inputs.pyscf['units'],
        atom = inputs.pyscf['atoms'],
        basis = inputs.pyscf['bs'],
        verbose = inputs.pyscf['verbosity'],
        symmetry = inputs.pyscf['symmetry'],
        spin=inputs.pyscf['spin'],
        charge=inputs.pyscf['charge'],
        # symmetry_subgroup = inputs.pyscf['symmetry_subgroup'], #0 is code for A1 point group
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
  

    if((inputs.run['zomgen']=='n')or(inputs.run['hamgen']=='n')or(inputs.run['imagprop']=='n')or(inputs.run['clean']=='f')):
        shutil.copytree(inputs.run['datafile'],EXDIR1+'/data')
    else:
        os.mkdir(EXDIR1+'/data')
     
    gcw.ham_writer(h1e,eri_full,EXDIR1)
    # shutil.copy2("ham_2.f90",EXDIR1+'/data')
    # shutil.copy2("ham_2.f90",'../src_2')
   

    with open(EXDIR1+'/rundata.csv','w',newline='')as file:
        writer = csv.writer(file)
        writer.writerow([inputs.run['zomgen'],inputs.run['hamgen'],inputs.run['imagprop'],inputs.run['beta'],inputs.run['timesteps'],inputs.run['clean'],inputs.run['gram'],inputs.run['gramnum'],inputs.run['grad'],'n'])
        writer.writerow(inputs.zombs.values())
        writer.writerow([inputs.run['hamfile'],inputs.run['ovrlfile'],inputs.run['cleanham']])



    os.chdir("../build")
    if(inputs.run['cores']==1):
        if(HPCFLG==1):
            if(inputs.run['GPU']=='y'):
                shutil.copy2("../build/makefile_gpu","../build/Makefile")
                subprocess.run(["make"])
            else:
                shutil.copy2("../build/makefile_arc","../build/Makefile")
                subprocess.run(["make"])
        else:
            shutil.copy2("../build/makefile_mac","../build/Makefile")
            subprocess.run(["make"])
    elif(inputs.run['cores']>1):
        if(HPCFLG==1):
            if(inputs.run['GPU']=='y'):
                shutil.copy2("../build/makefile_gpu","../build/Makefile")
                subprocess.run(["make"])
            else:
                shutil.copy2("../build/makefile_arc_omp","../build/Makefile")
                subprocess.run(["make"])
        else:
            shutil.copy2("../build/makefile_mac_omp","../build/Makefile")
            subprocess.run(["make"])
    
 
    shutil.copy2("ZOMBIE.exe",EXDIR1)

    
    os.chdir(EXDIR1)
   
    if(HPCFLG==1):
        if(inputs.run['cores']!=1):
            os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
        number=random.randint(99999,1000000)
        file1="zombie"+str(number)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        if(inputs.run['cores']!=1):
            f.write("#$ -pe smp "+str(inputs.run['cores'])+" \n") #Use shared memory parallel environemnt 
        if(inputs.run['GPU']=='y'):
            f.write("#$ -l coproc_v100=1 \n")
            f.write("#$ -P feps-gpu \n")
        f.write("#$ -l h_rt="+inputs.run['runtime']+"\n")
        f.write("#$ -l h_vmem=5G \n")
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
        # subprocess.run(["./d_check.exe"])
        
        

    

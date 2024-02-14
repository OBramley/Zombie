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
from integral_write import elec_writer

# Checking input paramaters 
if(isinstance(inputs.setup['nodes'],int)==False):
    sys.exit("Number of folders must be an integer")
elif(isinstance(inputs.setup['cores'],int)==False):
    sys.exit("Number of parallel cores must be an integer")
elif(inputs.setup['nodes']<1):
    sys.exit("Not enough nodes selected. Must be 1 or greater")
elif(inputs.setup['nodes']>100):
    sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
elif(inputs.setup['cores']>40):
    sys.exit("Too many cores selected. Maximum of 40 available")
elif(inputs.setup['cores']<1):
    sys.exit("Not enough cores selected. Must be 1 or greater")
elif(inputs.zombs['norb']<1):
    sys.exit("Not enough orbitals. Must be 1 or greater")
# elif(inputs.zombs['ndet']<2):
    # sys.exit("Not enough zombie states. Must be 2 or greater")
elif(inputs.zombs['zomtyp'] not in {'ran','HF','bb','hf'}):
    sys.exit("Type of zombie state must be ran, HF or bb")
elif(inputs.setup['elecs'] not in {'pyscf','mol','no'}):
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

elif(inputs.run['clean'] not in {'y','n','f'}):
    sys.exit("Setting cleaning on or off must be 'y' or 'n' or 'f' (to clean from a previously generated file)")


if(inputs.zombs['zomtyp']=='HF'):
    ndetcheck=0
    for i in range((inputs.zombs['norb']*2)+1):
        ndetcheck=ndetcheck+math.comb(inputs.zombs['norb']*2,i)
    if(inputs.run['ndet']!=ndetcheck):
        sys.exit('A Hartree Fock Basis for',inputs.zombs['norb'], 'orbitals should have', ndetcheck, 'basis functions')

if(inputs.setup['elecs']=='mol'):
    if not os.path.isfile(inputs.setup['datafolder']+'/'+ inputs.files['elecfile']):
        sys.exit("No Molpro input file")

if(inputs.setup['elecs']=='no' and inputs.run['hamgen']=='y'):
    if not os.path.isfile(inputs.setup['datafolder']+'/integrals'):
        sys.exit("No electron integral input file")

if(inputs.setup['elecs']=='no'and inputs.run['hamgen']=='n'):
    if not os.path.isfile(inputs.setup['datafolder']+'/'+ inputs.files['hamfile']):
        sys.exit("No Hamiltonian input file")

if(inputs.run['clean']=='f'):
    if not os.path.isfile(inputs.setup['datafolder']+'/'+ inputs.files['cleanham']):
        sys.exit("No clean hamiltonian")

if(inputs.run['gram']=='y'):
    if(isinstance(inputs.gram['gramnum'],int)==False):
        sys.exit("Number of states for Gram Schmidt must be an integer")
    elif(inputs.run['gramnum']<2):
        sys.exit("If using Gram Schmidt more than one state must investigated")

if(inputs.run['grad']=='y'):
    if(inputs.run['hamgen']=='n'):
        sys.exit("If using gradient descent a new Hamiltonian must be generated")
    elif(inputs.zombs['rhf_1']=='n'):
        sys.exit("If using gradient descent the first zombie state must be the RHF state")
    elif(isinstance(inputs.grad['decay_steps'],int)==False):
        sys.exit("Decasy steps must be an integer")
    elif(inputs.grad['decay_steps']<1):
        sys.exit("Decay steps must be 1 or greater")
    elif(isinstance(inputs.grad['learning_rate_decay'],float)==False):
        sys.exit("Decay rate must be a float")
    elif(inputs.grad['learning_rate_decay']>1 or inputs.grad['learning_rate_decay']<0):
        sys.exit("Decay rate must be between 0 and 1")
    elif(isinstance(inputs.grad['epoc_max'],int)==False):
        sys.exit("Number of epocs must be an integer")
    elif(inputs.grad['epoc_max']<2):
        sys.exit("Number of epocs must be 2 or greater")
    elif(inputs.grad['clone'] not in {'y','n'}):
        sys.exit("Clone parameter must be 'y' or 'n'")
    elif(isinstance(inputs.grad['clone_max'],int)==False):
        sys.exit("Clone max must be an integer")
    elif(isinstance(inputs.grad['clone_steps'],int)==False):
        sys.exit("Clone steps must be an integer")
    elif(isinstance(inputs.grad['clone_num'],int)==False):
        sys.exit("Number of additional clones must be an integer")

if(inputs.setup['multiple']>1):
    if(isinstance(inputs.run['multiple'],int)==False):
        sys.exit("Multiple runs must be an integer")
    if(inputs.setup['subnode']>1):
        if(isinstance(inputs.setup['subnode'],int)==False):
            sys.exit("Subnodes runs must be an integer")
        if((inputs.setup['multiple']%inputs.setup['subnode'])!=0):
            sys.exit("Number of multiple basis sets must be divisible by number of subnodes so the program can split the basis evenly")


print("Arguments checked")
# Hopefully temporary solution but current one to set the size of the dual paramter is to just write a tiny fortran module to set 
# the size of the dual array
   
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
if os.path.exists(EXDIR+"/"+inputs.setup['runfolder']):
    value=input("File already exists do you want to delete it? y/n\n")
    if(value=='y'):
        shutil.rmtree(EXDIR+"/"+inputs.setup['runfolder'])
    else:
        sys.exit("Runfolder already exists. Change the Runfolder name or delete/move it")

os.mkdir(EXDIR+"/"+inputs.setup['runfolder'])
EXDIR1=EXDIR+"/"+inputs.setup['runfolder']


shutil.copy2("inputs.py",EXDIR1)
shutil.copy2("graph.py",EXDIR1)
shutil.copy2("restart.py",EXDIR1)
shutil.copy2("combine.py",EXDIR1)
if((inputs.setup['elecs']=='no')):
    shutil.copytree(inputs.setup['datafolder']+'/integrals',EXDIR1+"/integrals")
elif((inputs.setup['elecs']=='mol')):
    os.mkdir(EXDIR1+"/integrals")
    shutil.copy2(inputs.setup['datafolder']+'/'+inputs.files['elecfile'],EXDIR1+"/integrals")
if(inputs.setup['elecs']=='pyscf'):
    os.mkdir(EXDIR1+"/integrals")
    if(inputs.pyscf['units']=='atom'):
        mol = gto.M(
        atom = inputs.pyscf['atoms'],
        basis = inputs.pyscf['bs'],
        verbose = inputs.pyscf['verbosity'],
        spin=inputs.pyscf['spin'],
        charge=inputs.pyscf['charge'],
        # symmetry_subgroup = inputs.pyscf['symmetry_subgroup'], #0 is code for A1 point group
        )
        myhf=scf.RHF(mol)
        myhf.kernel()
        c = myhf.mo_coeff
        # Get 1-electron integrals and convert to MO basis
        h1e = reduce(numpy.dot, (c.T, myhf.get_hcore(), c))
        # Get 2-electron integrals and transform them
        eri = ao2mo.kernel(mol, c)
        # Ignore all permutational symmetry, and write as four-index tensor, in chemical notation
        eri_full = ao2mo.restore(1, eri, c.shape[1])
        # Scalar nuclear repulsion energy
        Hnuc = myhf.energy_nuc()
       
    else:
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


multflg=0

if(inputs.run['zomgen']=='n')or(inputs.run['clean']=='f'):
    shutil.copytree(inputs.setup['datafile'],EXDIR1+'/data')
else:
    os.mkdir(EXDIR1+'/data')

if(inputs.run['hamgen']=='n')or(inputs.run['imagprop']=='n'):
    shutil.copytree(inputs.setup['datafile']+'/'+inputs.files['hamfile'],EXDIR1+'/data/ham.csv')
    shutil.copytree(inputs.setup['datafile']+'/'+inputs.files['ovrlfile'],EXDIR1+'/data/ovlp.csv')

elec_writer(h1e,eri_full,EXDIR1)

if(HPCFLG==1)and(inputs.run['multiple']>1):
    multflg=inputs.run['multiple']
    for i in range(multflg):
        os.mkdir(EXDIR1+'/node_'+str(i+1))
        shutil.copytree((EXDIR1+'/data'),(EXDIR1+'/node_'+str(i+1)+'/data'))
        shutil.copytree((EXDIR1+'/integrals'),(EXDIR1+'/node_'+str(i+1)+'/integrals'))
        with open(EXDIR1+'/node_'+str(i+1)+'/rundata.csv','w',newline='')as file:
            writer = csv.writer(file)
            writer.writerow(inputs.run.values())
            writer.writerow(inputs.zombs.values())
            if(inputs.run['grad']=='y'):
                writer.writerow(inputs.grad.values())
            if(inputs.run['gram']=='y'):
                writer.writerow(inputs.gram.values())
else:
    with open(EXDIR1+'/rundata.csv','w',newline='')as file:
        writer = csv.writer(file)
        writer.writerow(inputs.run.values())
        writer.writerow(inputs.zombs.values())
        if(inputs.run['grad']=='y'):
            writer.writerow(inputs.grad.values())
        if(inputs.run['gram']=='y'):
            writer.writerow(inputs.gram.values())

os.chdir("../build")
if(inputs.setup['cores']==1):
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
elif(inputs.setup['cores']>1):
    if(HPCFLG==1):
        if(inputs.setup['GPU']=='y'):
            shutil.copy2("../build/makefile_gpu","../build/Makefile")
            subprocess.run(["make"])
        else:
            shutil.copy2("../build/makefile_arc_omp","../build/Makefile")
            subprocess.run(["make"])
    else:
        shutil.copy2("../build/makefile_mac_omp","../build/Makefile")
        subprocess.run(["make"])


shutil.copy2("ZOMBIE",EXDIR1)
if(multflg>1):
    for i in range(multflg):
        shutil.copy2("ZOMBIE",EXDIR1+'/node_'+str(i+1))

os.chdir(EXDIR1)

if(HPCFLG==1):
    if(inputs.run['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])
    file1="zombie_"+inputs.setup['runfolder']+".sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    if(inputs.run['cores']!=1):
        f.write("#$ -pe smp "+str(inputs.run['cores'])+" \n") #Use shared memory parallel environemnt 
    if(inputs.run['GPU']=='y'):
        f.write("#$ -l coproc_v100=1 \n")
        f.write("#$ -P feps-gpu \n")
    f.write("#$ -l h_rt="+inputs.run['runtime']+"\n")
    f.write("#$ -l h_vmem=1G \n")
    f.write("module add mkl \n")
    f.write('time ./ZOMBIE')
    f.close()
    if(multflg>1):
        for j in range(multflg):
            file2=EXDIR1+'/node_'+str(j+1)+"/zombie_"+inputs.setup['runfolder']+"_"+str(j+1)+".sh"
            shutil.copy2(file1,file2)

        for j in range(multflg):
            os.chdir(EXDIR1+'/node_'+str(j+1))
            file2="zombie_"+inputs.setup['runfolder']+"_"+str(j+1)+".sh"
            command=['qsub','-N', inputs.setup['runfolder']+"_"+str(j+1)+"_1",file2]
            for i in range(1,inputs.setup['submissions']+1):
                subprocess.call(command)
                command=['qsub','-N',inputs.setup['runfolder']+"_"+str(j+1)+"_"+str(i+1),'-hold_jid',inputs.setup['runfolder']+"_"+str(j+1)+"_"+str(i),file2]
    else:
        command=['qsub','-N',inputs.setup['runfolder']+"_1",file1]
        for i in range(1,inputs.setup['submissions']+1):
            subprocess.call(command)
            command=['qsub','-N',inputs.setup['runfolder']+"_"+str(i+1),'-hold_jid',inputs.setup['runfolder']+"_"+str(i),file1]     
else:
    print(os.getcwd())
    if(inputs.setup['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs.setup['cores'])
    subprocess.run(["./ZOMBIE"])

# if(inputs.run['cores']!=1): 
#     os.environ["OMP_NUM_THREADS"]=str(inputs.run['cores'])   
# os.environ["OMP_CANCELLATION"]="TRUE"   
# subprocess.run(["./ZOMBIE"])
 

    

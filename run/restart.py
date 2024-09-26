###################################################################################################
#
# This is the run file for the zombie states program. It checks run paramters are correct. It then
# compiles the Fortran modules; creates an execution folder and copies the relevant files to it.
# The program then begins the run requesting nodes to allow paralel execution if on the the HPC
# 
###################################################################################################
import sys
import os
import socket
import subprocess
import csv
import json 
with open('inputs.json') as f:
    inputs=json.load(f)

if(inputs['setup']['elecs']=='pyscf'):
    with open(inputs['pyscf_file']+'.json') as json_file:
        pyscf_ins=json.load(json_file)
    inputs['zombs']['norb']= pyscf_ins['norb']
    inputs['zombs']['nel']= pyscf_ins['nel']
    inputs['zombs']['spin']= pyscf_ins['spin']

# Checking input paramaters 
if(isinstance(inputs['setup']['nodes'],int)==False):
    sys.exit("Number of folders must be an integer")
elif(isinstance(inputs['setup']['cores'],int)==False):
    sys.exit("Number of parallel cores must be an integer")
elif(inputs['setup']['nodes']<1):
    sys.exit("Not enough nodes selected. Must be 1 or greater")
elif(inputs['setup']['nodes']>100):
    sys.exit("Too many nodes. Maximum of 100 simultaneous submisions")
elif(inputs['setup']['cores']>40):
    sys.exit("Too many cores selected. Maximum of 40 available")
elif(inputs['setup']['cores']<1):
    sys.exit("Not enough cores selected. Must be 1 or greater")
elif(inputs['zombs']['norb']<1):
    sys.exit("Not enough orbitals. Must be 1 or greater")
# elif(inputs['zombs']['ndet']<2):
    # sys.exit("Not enough zombie states. Must be 2 or greater")
elif(inputs['zombs']['zomtyp'] not in {'ran','HF','bb','hf'}):
    sys.exit("Type of zombie state must be ran, HF or bb")
elif(inputs['setup']['elecs'] not in {'pyscf','mol','no'}):
    sys.exit("one and two electron integral paramter must be 'pyscf', 'mol' or 'no'")
elif(inputs['run']['zomgen'] not in {'y','n'}):
    sys.exit("Generation of a new set of zombie states must be either 'y' or 'n'")
elif(inputs['run']['imagprop'] not in {'y','n'}):
    sys.exit("The imaginary time propagation parameter has been inccorectly set must be either 'y' or 'n'")
elif(isinstance(inputs['run']['beta'],int)==False):
    sys.exit("Beta must be an integer")
elif(isinstance(inputs['run']['timesteps'],int)==False):
    sys.exit("Beta must be an integer")
elif(inputs['run']['gram'] not in {'y','n'}):
    sys.exit("The Gram Schmidt pramater must be either 'y' or 'n'")
elif(inputs['run']['clean'] not in {'y','n','f'}):
    sys.exit("Setting cleaning on or off must be 'y' or 'n' or 'f' (to clean from a previously generated file)")

if(inputs['run']['gram']=='y'):
    if(isinstance(inputs['gram']['gramnum'],int)==False):
        sys.exit("Number of states for Gram Schmidt must be an integer")
    elif(inputs['gram']['gramnum']<1):
        sys.exit("If using Gram Schmidt more than one state must investigated")

if(inputs['run']['grad']=='y'):
    if(inputs['run']['hamgen']=='n'):
        sys.exit("If using gradient descent a new Hamiltonian must be generated")
    elif(inputs['zombs']['rhf_1']=='n'):
        sys.exit("If using gradient descent the first zombie state must be the RHF state")
    elif(isinstance(inputs['grad']['decay_steps'],int)==False):
        sys.exit("Decasy steps must be an integer")
    elif(inputs['grad']['decay_steps']<1):
        sys.exit("Decay steps must be 1 or greater")
    elif(isinstance(inputs['grad']['learning_rate_decay'],float)==False):
        sys.exit("Decay rate must be a float")
    elif(inputs['grad']['learning_rate_decay']>1 or inputs['grad']['learning_rate_decay']<0):
        sys.exit("Decay rate must be between 0 and 1")
    elif(isinstance(inputs['grad']['epoc_max'],int)==False):
        sys.exit("Number of epocs must be an integer")
    elif(inputs['grad']['epoc_max']<1):
        sys.exit("Number of epocs must be 1 or greater")
    elif(inputs['grad']['clone'] not in {'y','n'}):
        sys.exit("Clone parameter must be 'y' or 'n'")
    elif(isinstance(inputs['grad']['clone_max'],int)==False):
        sys.exit("Clone max must be an integer")
    elif(isinstance(inputs['grad']['clone_steps'],int)==False):
        sys.exit("Clone steps must be an integer")
    elif(isinstance(inputs['grad']['clone_num'],int)==False):
        sys.exit("Number of additional clones must be an integer")

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

with open('rundata.csv','w',newline='')as file:
    writer = csv.writer(file)
    writer.writerow(inputs['run'].values())
    writer.writerow(inputs['zombs'].values())
    if(inputs['run']['grad']=='y'):
        writer.writerow(inputs['grad'].values())
    if(inputs['run']['gram']=='y'):
        writer.writerow(inputs['gram'].values())
# HPCFLG=0
if(HPCFLG==1):
    if(inputs['setup']['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs['setup']['cores'])
    file1="zombie_"+inputs['setup']['runfolder']+".sh"
    if os.path.exists(file1):
        os.remove(file1)
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    if(inputs['setup']['cores']!=1):
        f.write("#$ -pe smp "+str(inputs['setup']['cores'])+" \n") #Use shared memory parallel environemnt 
    if(inputs['setup']['GPU']=='y'):
        f.write("#$ -l coproc_v100=1 \n")
        f.write("#$ -P feps-gpu \n")
    f.write("#$ -l h_rt="+inputs['setup']['runtime']+"\n")
    f.write("#$ -l h_vmem=1G \n")
    f.write("module add mkl \n")
    f.write('time ./ZOMBIE')
    f.close()

    command=['qsub','-N',inputs['setup']['runfolder']+"_1",file1]
    for i in range(1,inputs['setup']['submissions']+1):
        subprocess.call(command)
        command=['qsub','-N',inputs['setup']['runfolder']+"_"+str(i+1),'-hold_jid',inputs['setup']['runfolder']+"_"+str(i),file1]     
else:
    print(os.getcwd())
    if(inputs['setup']['cores']!=1):
        os.environ["OMP_NUM_THREADS"]=str(inputs['setup']['cores'])
    subprocess.run(["./ZOMBIE"])
        
        

    

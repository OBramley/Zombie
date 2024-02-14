###################################################################################################
# 
# This file contains all the paramters that can be altered to generate a system of Zombie States
# So far the run paramters are split into two types run and zomp. Run paramters control the type 
# of input file(s) if any the program needs and the type of output files it will output. Zomp are the 
# system paramters that set up the zombie state size and number of Zombie states. 
# The PyScf paramters can also be inputed here.
# Written by O.A Bramley
#
####################################################################################################

import pyscf_params

setup={
    # What is the name of the run
    'runfolder':'be_new_test',

    # Amount of time to request on HPC 
    'runtime': "48:00:00",

    # The number of nodes to request on HPC   
    'nodes':1,
    
    # The number of cores to request on HPC
    'cores':8,

    # The number of threads to request on HPC
    'submissions':25,
    # If multplie small runs are wanted to be started and then combined together
    'multiple':1,
    # If multiple small runs do they want combinign in
    'subnode':1,

    # Name of folder where any prior generated results are placed so they can 
    # be appropriately copied as required
    'datafolder':'data copy',

    #GPU flag
    'GPU':'n',

    # 1 and 2 electron integrals can be calcualted by PySCF by the program
    # or they can be inputed as a seperate file at the moment only from MOLPRO
    # as such the program takes 3 inputs 'pyscf', 'mol' or 'no' if the one and two
    # electron integrals have already been generated and placed in the run file
    'elecs':'pyscf'
}

run={
    # Number of Zombie states
    'ndet':5,

    # Set the random seed for Zombie state generation
    # Setting as zero will generate a new random seed any other number allows the 
    # same seed to be used again.
    'seed':2,
    # Do you want to generate new zombie states. The program can work using previously gerated
    # Zombie states. Takes input 'y' or 'n'.
    'zomgen':'y',

    # Do you want to generate a new Hamiltonian? Takes input 'y' or 'n'
    'hamgen':'y',

    # Do you want to clean after propagation takes 'y', 'n' or 'f' to use a previosuly generated
    # cleaning hamiltonian and zombie state files
    'clean':'n',

    # Do you want to perform imaginary time evolution? The program can be used to just generate
    # the Zombie state Hamiltonian. Takes input 'y' or 'n'.
    'imagprop':'y',
    # Length of imaginary time propagation
    'beta':500,
    # Number of timesteps to take in imaginary time propagation
    'timesteps':1000,

    # Gradient descent flag takes input 'y' or 'n'.
    'grad':'y',
    
    # Do you want to find other energy states other than the ground state. If so turn on 
    # Gram Schmidt orthogonalisation and then specify the number of states. Takes input
    # 'y' or 'n' and an integer number note in the python code gramnum=1 will not find an excited state just a 
    # single ground state. But gramnum=1 will find a single excited state in the fortran code.
    'gram':'n'
}

if(setup['elecs']=='pyscf'):
    pyscf=pyscf_params.Be_atom_ccpvdz
else:
    pyscf={
        'norb':10,
        'nel':5,
        'spin':0
    }

zombs={
    
    # Number of orbitals
    'norb':pyscf['norb'],

    # Number of electrons in the molecule
    'nel':pyscf['nel'],

    # Spin of the moleucle 

    'spin':pyscf['spin'],

    # Type of zombie states. Random (ran), Hartree Fock (HF) or biased (bb)
    'zomtyp':'bb',
    
    # Make the first Zombie state the RHF det? Takes y or n
    'rhf_1':'y',

    # Make the zombie states imaginary (y) or real (n)
    'imagflg':'n'
}

# The gradient descent paramters
grad={
    # The number of iterations to perform
    'epoc_max':1000000,
    # The learning rate
    'initial_learning_rate':2500,
    # The decay rate of the learning rate
    'learning_rate_decay':0.2,
    # The number of times the learning rate decays before returning to start
    'decay_steps':6,
    # Do you want to increase the size of the basis set as you perform Gradient Descent. If so turn on 
    # cloning and then specify the number basis set size. Takes input 'y' or 'n'
    # Also need to specify the number of steps before blind cloning occurs and the numbner of 
    #functions to increase the basis set by.
    'clone':'y',
    'clone_max':50,
    'clone_steps':450,
    'clone_num':1
} 

# The Gram Schmidt orthogonalisation paramters
gram={
    'gramnum':4,
}

# Names of files if inputting precalculated values
files={
    'hamfile':'ham.csv',
    'ovrlfile':'ovlp.csv',
    'elecfile':'integrals.pkl',
    'cleanham':'clean_ham.csv',
    'cleanzom':'exmple'
}






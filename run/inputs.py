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



run={
    # Decides if the program will run in Python or Fortran
    # Fortran program is faster and will use OpenMP for multiple threads
    # The python program is possibly easier to understand and is still effective for 
    # small runs. Takes input python or fortran
    'language':'fortran',

    # What is the name of the run
    'runfolder':'bh_76_with_overshoot',

    # Amount of time to request on HPC 
    'runtime': "00:05:00",

    'nodes':1,

    'cores':1,

    # Set the numpy random seed for Zombie state generation
    'seed':1,

    # Name of file where any prior generated results are placed so the Fortran 
    # program can access them and continue a run
    'datafile':'data copy',

    # 1 and 2 electron integrals can be calcualted by PySCF by the program
    # or they can be inputed as a seperate file at the moment only from MOLPRO
    # as such the program takes 3 inputs 'pyscf', 'mol' or 'no' if the one and two
    # electron integrals have already been generated and placed in the run file
    'elecs':'pyscf',

    'elecfile':'integrals.pkl',

    # Do you want to generate new zombie states. The program can work using previously gerated
    # Zombie states. Takes input 'y' or 'n'.
    'zomgen':'y',

    'zombiefile':'zombie_states.pkl',

    # Do you want to generate a new Hamiltonian? Takes input 'y' or 'n'
    'hamgen':'y',

    'hamfile':'ham.csv',

    'ovrlfile':'ovlp.csv',

    # Do you want to perform imaginary time evolution? The program can be used to just generate
    # the Zombie state Hamiltonian. Takes input 'y' or 'n'.
    'imagprop':'y',

    'beta':500,

    'timesteps':1000,

    # Do you want the starting energy to be the HF energy
    # Takes input 'y' or 'n' and then a number to defined the number of electrons
    'zomhf':'y',
    'hfnum': 6,

    # Do you want to clean after propagation takes 'y', 'n' or 'f' to use a previosuly generated
    # cleaning hamiltonian and zombie state files
    'clean':'n',
    'cleanham':'FCIconfigs_equilibrium.txt',
    'cleanzom':'BH_clean_zombie_states.pkl', 

    # Do you want to find other energy states other than the ground state. If so turn on 
    # Gram Schmidt orthogonalisation and then specify the number of states. Takes input
    # 'y' or 'n' and an integer number note in the python code gramnum=1 will not find an excited state just a 
    # single ground state. But gramnum=1 will find a single excited state in the fortran code.
    'gram':'n',
    'gramnum':4,

    # Gradient descent flag takes input 'y' or 'n'. Only implemented in the fortran version 
    'grad':'y',

    #GPU flag
    'GPU':'n'
}


zombs={
    # Number of orbitals
    'norb':19,

    # Number of electrons in the molecule
    'nel':6,

    # Spin of the moleucle 

    'spin':0,

    # Number of Zombie states
    'ndet':20, 

    # Type of zombie states. Random (ran), Hartree Fock (HF) or biased (bb)
    'zomtyp':'bb',
    
    # Make the first Zombie state the RHF det? Takes y or n
    'rhf_1':'y',

    # Make the zombie states imaginary (y) or real (n)
    'imagflg':'n'
}

# pyscf={
#      # The units the geometry of the molecule is set up in
#     'units':'Angstrom',
#      # The geometry of the molecule being investigated
#     'atoms': 'B 0 0 0; H 0 0 4.0',
#     # The type of basis used to generate the 1 and 2 electron integrals
#     'bs' : '6-31g**',
#     # How verbose do you want the PyScf output to be in your terminal?
#     'verbosity' : 4,
#     'symmetry' :True,
#     'spin':0,
#     'charge':0,
#     'symmetry_subgroup' : 0 #0 is code for A1 point group
# }

pyscf={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms': 'B 0 0 0; H 0 0 1.2324',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 4,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0 #0 is code for A1 point group
}

# pyscf={
#      # The units the geometry of the molecule is set up in
#     'units':'Bohr',
#      # The geometry of the molecule being investigated
#     'atoms': 'Li 0 0 0; Li 0 0 6',
#     # The type of basis used to generate the 1 and 2 electron integrals
#     'bs' : '6-31g**',
#     # How verbose do you want the PyScf output to be in your terminal?
#     'verbosity' : 4,
#     'symmetry' :True,
#     'spin':0,
#     'charge':0,
#     'symmetry_subgroup' : 0 #0 is code for A1 point group
# }

# pyscf={
#      # The units the geometry of the molecule is set up in
#     'units':'Angstrom',
#      # The geometry of the molecule being investigated
#     'atoms': 'N 0 0 0; N 0 0 1.094',
#     # The type of basis used to generate the 1 and 2 electron integrals
#     'bs' : 'cc-pVDZ',
#     # How verbose do you want the PyScf output to be in your terminal?
#     'verbosity' : 4,
#     'symmetry' :True,
#     'spin':0,
#     'charge':0,
#     'symmetry_subgroup' : 0 #0 is code for A1 point group
# }
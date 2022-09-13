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
    # What is the name of the run
    'runfolder':'d_check_ran',

    # Amount of time to request on HPC 
    'runtime': "24:00:00",

    'nodes':1,

    'cores':4,

    # Set the numpy random seed for Zombie state generation
    'seed':1,

    # 1 and 2 electron integrals can be calcualted by PySCF by the program
    # or they can be inputed as a seperate file at the moment only from MOLPRO
    # as such the program takes 3 inputs 'pyscf', 'mol' or 'no' if the one and two
    # electron integrals have already been generated and placed in the run file
    'elecs':'pyscf',

    'beta':200,

    'timesteps':2000,
}


zombs={
    # Number of orbitals
    'norb':19,

    # Number of electrons in the molecule
    'nel':6,

    # Spin of the moleucle 

    'spin':0,

    # Number of Zombie states
    'ndet':100, 

    # Type of zombie states. Random (ran), Hartree Fock (HF) or biased (bb)
    'zomtyp':'ran',
    
    # Biased basis improvement if 0 no loops to improve biased > 0 number of loops to improve the basis

    'bb_imprv':0,

    # Make the first Zombie state the RHF det? Takes y or n
    'rhf_1':'n',

    # Make the zombie states imaginary (y) or real (n)
    'imagflg':'n'
}

# Parameters for the biasing protocol
# The biasing method uses sampling of a trasformation of tanh graph
# Note no paramter chekcing has been implemented for these values yet
zom_bias={
    # Number of ocuupied spin orbitals
    # THE NUMBER OF SPIN ORBITALS IS HALF THE NUMBER OF ELECTRONS
    'alive':2,
    # Alive aplitude of 1st spin orbital
    'alive_start':0.999,
    # Alive aplitude of last "alive" orbital
    'alive_end':0.175,
    # Number of spin orbitals approximately 50/50 alive/dead aplitudes
    'mid':0,
    # Alive aplitude of 1st dead spin orbtial
    'dead_start':0.351,
    # Alive aplitude of last spin orbital
    'dead_end':0.120

}

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

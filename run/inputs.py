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
    'runfolder':'pyscf_test',

    'nodes':1,

    'cores':1,

    # Set the numpy random seed for Zombie state generation
    'seed':1,

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

    'hamfile':'bigham.csv',

    'ovrlfile':'kover.csv',

    # Do you want to perform imaginary time evolution? The program can be used to just generate
    # the Zombie state Hamiltonian. Takes input 'y' or 'n'.
    'imagprop':'y',

    'beta':200,

    'timesteps':2000,

    # Do you want the starting energy to be the HF energy
    # Takes input 'y' or 'n' and then a number to defined the number of electrons
    'zomhf':'n',
    'hfnum': 6, 

    # Do you want to find other energy states other than the ground state. If so turn on 
    # Gram Schmidt orthogonalisation and then specify the number of states. Takes input
    # 'y' or 'n' and an integer number note gramnum=1 will not find an excited state just a 
    # single ground state
    'gram':'n',
    'gramnum':4
}


zombs={
    # Number of electrons
    'norb':38,

    # Number of Zombie states
    'ndet':10, 

    # Type of zombie states. Random (ran), Hartree Fock (HF) or biased (bb)
    'zomtyp':'ran'
}

# Parameters for the biasing protocol
# The biasing method uses sampling of a trasformation of tanh graph
# Note no paramter chekcing has been implemented for these values yet
zom_bias={
    # Number of ocuupied spin orbitals
    # THE NUMBER OF SPIN ORBITALS IS HALF THE NUMBER OF ELECTRONS
    'alive':3,
    # Alive aplitude of 1st spin orbital
    'alive_start':0.99,
    # Alive aplitude of last "alive" orbital
    'alive_end':0.94,
    # Number of spin orbitals approximately 50/50 alive/dead aplitudes
    'mid':0,
    # Alive aplitude of 1st dead spin orbtial
    'dead_start':0.008,
    # Alive aplitude of last spin orbital
    'dead_end':0.00006

}

pyscf={
     # The units the geometry of the molecule is set up in
    'units':'Angstrom',
     # The geometry of the molecule being investigated
    'atoms' :'B 0 0 0; H 0 0 1.2324',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 4,
    'symmetry' :True,
    'spin':0,
    'charge':0,
    'symmetry_subgroup' : 0, #0 is code for A1 point group
}
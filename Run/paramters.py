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
    # 1 and 2 electron integrals can be calcualted by PySCF by the program
    # or they can be inputed as a seperate file at the moment only from MOLPRO
    # as such the program takes 3 inputs 'pyscf', 'mol' or 'no' if the one and two
    # electron integrals have already been generated and placed in the run file
    'elecs':'pyscf',

    # Do you want to generate new zombie states. The program can work using previously gerated
    # Zombie states. Takes input 'y' or 'n'.
    'zomgen':'y',

    # Do you want to perform imaginary time evolution? The program can be used to just generate
    # the Zombie state Hamiltonian. Takes input 'y' or 'n'.
    'imagprop':'y'
}


zomp={
    # Number of orbitals
    'norb':10,

    # Number of Zombie states
    'ndet':64 
}


pyscf={
     # The units the geometry of the molecule is set up in
    'units':'Bohr',
     # The geometry of the molecule being investigated
    'atoms' :'Li 0 0 0; Li 0 0 6',
    # The type of basis used to generate the 1 and 2 electron integrals
    'bs' : '6-31g**',
    # How verbose do you want the PyScf output to be in your terminal?
    'verbosity' : 1,
    'symmetry' :True,
    # 'spin':1,
    # 'charge':-1
    'symmetry_subgroup' : 0, #0 is code for A1 point group
}

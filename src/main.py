###################################################################################################
# 
# This is the main script for the Zombie states method. It has been written in a modular manner
# to allow it to form the basis of any molecular dynamics program. The program has a v ery simple 
# strucure: one and two electron integrals are created or loaded into the program; zombie states
# are generated; the hamiltonian matrix is generated; an initial state is set; imaginary time 
# propagation is carried out. 
# 
###################################################################################################
import inputs
import elec_integrals

# Load paramters
ndet=inputs.zombs['ndet']
norb=inputs.zombs['norb']

# Generate 1 and 2 electron integrals
if((inputs.run['elecs'])=='pyscf'):
    Hnuc, H1ei, H2ei=elec_integrals.pyscf_gen(norb)
elif((inputs.run['elecs'])=='mol'):
    Hnuc, H1ei, H2ei=elec_integrals.molpro_read(norb,inputs.run['elecfile'])
elif((inputs.run['elecs'])=='no'):
    Hnuc, H1ei, H2ei = elec_integrals.read_in(inputs.run['elecfile'])


# Generate zombie states

# Generate Hamiltonian and overlap matrix

# Imaginary time propagation
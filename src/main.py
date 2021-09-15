import elec_integrals


# Set up Parameters of simulation 
# Number of spin orbitals
norb=10

# Number of Zombie states
ndet=64

# Do you want to generate electron integrals?
elec_gen=True
# If already generated add the file name. Program currently only works for pyscf and MOLPRO
filename='example.dat'



# Generate Zombie states?



# Generate 1&2 electron integrals? or read them in?
if(elec_gen==True):
    Hnuc, H1ei, H2ei = elec_integrals.pyscf_gen(norb)
else:
    Hnuc, H1ei, H2ei = elec_integrals.molpro_read(filename,norb)
# Generate Hamiltonian 
# Initialise starting state 
# Imaginary time propagation 
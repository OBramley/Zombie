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
import ham
import zom
import in_outputs
import gc
import imgtp

# Load paramters
ndet=inputs.zombs['ndet']
norb=inputs.zombs['norb']

# Check to read in or generate Hamiltonian and overlap matrix
if((inputs.run['hamgen'])=='n'):
    Bigham=in_outputs.read_ham(inputs.run['hamfile'])
    print('Hamiltonian read in')
    Kover=in_outputs.read_ham(inputs.run['ovrlfile'])
    print('Overlap matrix read in')
elif((inputs.run['hamgen'])=='y'):
    # Generate or read in 1 and 2 electron integrals
    if((inputs.run['elecs'])=='pyscf'):
        Ham = ham.pyscf_gen(norb)
    elif((inputs.run['elecs'])=='mol'):
        Ham = ham.molpro_read(norb,inputs.run['elecfile'])
    elif((inputs.run['elecs'])=='no'):
        Ham = ham.read_in(norb, inputs.run['elecfile'])
        # elec_integrals.write(Hnuc,H1ei,H2ei,norb)

    # Generate or read in zombie states
    if((inputs.run['zomgen'])=='y'):
        zstore=zom.zom_gen(norb,ndet,inputs.zombs['zomtyp'],'zombie_states.pkl')
        print('Zombie states generated')
    elif((inputs.run['zomgen'])=='n'):
        zstore=in_outputs.read_object(inputs.run['zombiefile'])
        print('Zombie states read in')

    Bigham, Kover = ham.hamiltonian(ndet,Ham,zstore)
    # print(Bigham)


del Ham
gc.collect()

# Imaginary time propagation
if(inputs.run['imagprop']=='y'):
# Check if Gram Schmidt orthogonalisation is to be used
    if(inputs.run['gram']=='y'):
        eb=imgtp.itime_prop_gs(Bigham, Kover,inputs.run['beta'],inputs.run['timesteps'],norb,inputs.run['zomhf'],inputs.run['zomhf'],zstore,ndet, inputs.run['gramnum'])    
    elif(inputs.run['gram']=='n'):
        eb=imgtp.itime_prop(Bigham, Kover,inputs.run['beta'],inputs.run['timesteps'],norb,inputs.run['zomhf'],inputs.run['zomhf'],zstore,ndet)
elif(inputs.run['imagprop']=='n'):
    print('End of program')

if(inputs.run['gram']=='y'):
    rnum=inputs.run['gramnum']=='y'
else:
    rnum=1
in_outputs.plot(eb,rnum, inputs.run['beta'],'results.png')
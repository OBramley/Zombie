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
filenamer=inputs.run['runfolder']
ndet=inputs.zombs['ndet']
norb=inputs.zombs['norb']
beta=inputs.run['beta']
timesteps=inputs.run['timesteps']

# Check to read in or generate Hamiltonian and overlap matrix
if((inputs.run['hamgen'])=='n'):
    Bigham=in_outputs.read_ham(inputs.run['hamfile'])
    print('Hamiltonian read in')
    Kover=in_outputs.read_ham(inputs.run['ovrlfile'])
    print('Overlap matrix read in')
    if(inputs.run['zomhf']=='n'):
        zstore=[]
    else:
        if((inputs.run['zomgen'])=='y'):
            zstore=zom.zom_gen(norb,ndet,inputs.zombs['zomtyp'],filenamer,inputs.run['seed'],inputs.zom_bias)
            print('Zombie states generated')
        elif((inputs.run['zomgen'])=='n'):
            zstore=in_outputs.read_object(inputs.run['zombiefile'])
            print('Zombie states read in')
elif((inputs.run['hamgen'])=='y'):
    # Generate or read in 1 and 2 electron integrals
    if((inputs.run['elecs'])=='pyscf'):
        Ham = ham.pyscf_gen(norb,filenamer)
        print('pyscf generated integrals')
    elif((inputs.run['elecs'])=='mol'):
        Ham = ham.molpro_read(norb,inputs.run['elecfile'],filenamer)
        print('Molpro integrals generated')
    elif((inputs.run['elecs'])=='no'):
        Ham = in_outputs.read_object(inputs.run['elecfile'])
        print('Electron integrals read in')
    # Generate or read in zombie states
    if((inputs.run['zomgen'])=='y'):
        zstore=zom.zom_gen(norb,ndet,inputs.zombs['zomtyp'],filenamer,inputs.run['seed'],inputs.zom_bias)
        print('Zombie states generated')
    elif((inputs.run['zomgen'])=='n'):
        zstore=in_outputs.read_object(inputs.run['zombiefile'])
        print('Zombie states read in')

    Bigham, Kover = ham.hamiltonian(ndet,Ham,zstore,filenamer)
    del Ham
    gc.collect()


# Imaginary time propagation
if(inputs.run['imagprop']=='y'):
# Check if Gram Schmidt orthogonalisation is to be used
    if(inputs.run['gram']=='y'):
        eb=imgtp.itime_prop_gs(Bigham,Kover,beta,timesteps,norb,inputs.run['zomhf'],inputs.run['hfnum'],zstore,ndet, inputs.run['gramnum'], filenamer)    
    elif(inputs.run['gram']=='n'):
        eb=imgtp.itime_prop(Bigham, Kover,beta,timesteps,norb,inputs.run['zomhf'],inputs.run['hfnum'],zstore,ndet,filenamer)
elif(inputs.run['imagprop']=='n'):
    print('End of program')

if(inputs.run['gram']=='y'):
    rnum=inputs.run['gramnum']
else:
    rnum=1

in_outputs.plot(eb,rnum, inputs.run['beta'],'results.png')
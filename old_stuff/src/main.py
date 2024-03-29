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
import cleaning
import subprocess

# Load paramters
filenamer=inputs.run['runfolder']
ndet=inputs.zombs['ndet']
norb=inputs.zombs['norb']*2
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
            zstore=zom.zom_gen(norb,ndet,inputs.zombs,filenamer,inputs.run['seed'],inputs.zom_bias)
            print('Zombie states generated')
        elif((inputs.run['zomgen'])=='n'):
            zstore=in_outputs.read_object(inputs.run['zombiefile'])
            print('Zombie states read in')
        if(inputs.run['clean']=='y'):
            Ham = ham.pyscf_gen(norb,filenamer,inputs.pyscf)
            cleaning.clean_setup(norb,inputs.zombs['nel'],ndet,inputs.zombs['spin'], Ham, zstore, filenamer)
            print('Cleaning hamiltonian generated')
elif((inputs.run['hamgen'])=='y'):
    # Generate or read in 1 and 2 electron integrals
    if((inputs.run['elecs'])=='pyscf'):
        Ham = ham.pyscf_gen(norb,filenamer,inputs.pyscf)
        print('pyscf generated integrals')
    elif((inputs.run['elecs'])=='mol'):
        Ham = ham.molpro_read(norb,inputs.run['elecfile'],filenamer)
        print('Molpro integrals generated')
    elif((inputs.run['elecs'])=='no'):
        Ham = in_outputs.read_object(inputs.run['elecfile'])
        print('Electron integrals read in')
    # Generate or read in zombie states
    if((inputs.run['zomgen'])=='y'):
        zstore=zom.zom_gen(norb,ndet,inputs.zombs,filenamer,inputs.run['seed'],inputs.zom_bias)
        print('Zombie states generated')
    elif((inputs.run['zomgen'])=='n'):
        zstore=in_outputs.read_object(inputs.run['zombiefile'])
        print('Zombie states read in')
    if(inputs.run['clean']=='y'):
        cleaning.clean_setup(norb,inputs.zombs['nel'],ndet,inputs.zombs['spin'], Ham, zstore, filenamer)
        print('Cleaning hamiltonian generated')

    Bigham, Kover = ham.hamiltonian(ndet,Ham,zstore,filenamer)
    print('Hamiltonian generated')
    del Ham
    gc.collect()


if(inputs.zombs['bb_imprv']==0):
    # Imaginary time propagation
    if(inputs.run['imagprop']=='y'):
    # Check if Gram Schmidt orthogonalisation is to be used
        if(inputs.run['gram']=='y'):
            eb=imgtp.itime_prop_gs(Bigham,Kover,beta,timesteps,norb,inputs.run['zomhf'],inputs.run['hfnum'],zstore,ndet, inputs.run['gramnum'], filenamer)    
        elif(inputs.run['gram']=='n'):
            eb=imgtp.itime_prop(Bigham, Kover,beta,timesteps,norb,inputs.run['zomhf'],inputs.run['hfnum'],zstore,ndet,filenamer)
    elif(inputs.run['imagprop']=='n'):
        print('End of program')
else:
    for i in range(inputs.zombs['bb_imprv']):
        # Imaginary time propagation
        if(inputs.run['imagprop']=='y'):
        # Check if Gram Schmidt orthogonalisation is to be used
            if(inputs.run['gram']=='y'):
                eb=imgtp.itime_prop_gs(Bigham,Kover,beta,timesteps,norb,inputs.run['zomhf'],inputs.run['hfnum'],zstore,ndet, inputs.run['gramnum'], filenamer+'_'+str(i))    
            elif(inputs.run['gram']=='n'):
                eb=imgtp.itime_prop(Bigham, Kover,beta,timesteps,norb,inputs.run['zomhf'],inputs.run['hfnum'],zstore,ndet,filenamer+'_'+str(i))

        zom.biased_check(zstore,norb,ndet,filenamer+'_'+str(i))
        




# Cleaning Protocol
if(inputs.run['clean']=='y'):
    clean_energy, norm = cleaning.cleaner(ndet, zstore, filenamer+'_clean_hamiltonian.csv', filenamer+'_clean_zombie_states.pkl',filenamer)
    print(clean_energy)
    print(norm)
    print(clean_energy/norm)
elif(inputs.run['clean']=='f'):
    clean_energy, norm = cleaning.cleaner(ndet, zstore, inputs.run['cleanham'], inputs.run['cleanzom'],filenamer)
    print(clean_energy)
    print(norm)
    print(clean_energy/norm)


if(inputs.run['gram']=='y'):
    rnum=inputs.run['gramnum']
else:
    rnum=1

if(inputs.run['clean']=='y')or(inputs.run['clean']=='f'):
    in_outputs.clean_plot(eb,rnum, inputs.run['beta'], (clean_energy/norm),timesteps,'results.png')
else:
    in_outputs.plot(eb,rnum, inputs.run['beta'],timesteps,'results.png')



# subprocess.run(['rclone', 'copy', '/nobackup/cm14oab/'+filenamer,'onedrive:'+filenamer])



import shlex
import subprocess
from matplotlib import pyplot as plt
import matplotlib as mpl
from pylab import cm
import d_check_inputs
import numpy 

def plot(results,averages, filename):
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',4)
    fig,ax=plt.subplots(3,1,figsize=(3.37,10.11))
    # ax[0]=fig.add_axes([0,0,2,1])
    # ax[1]=fig.add_axes([0,0,2,1])
    ax[0].scatter(results[:,0],results[:,1],color=colors(0))
    ax[1].scatter(results[:,0],results[:,2],color=colors(1))
    ax[2].scatter(averages[:,0],averages[:,1],color=colors(2))
    # ax[0].set_ylim(numpy.amin(results[:,1])-0.1,numpy.amax(results[:,1])+0.1)
    # ax[1].set_ylim(numpy.amin(results[:,2])-0.1,numpy.amax(results[:,2])+0.1)
    ax[0].set_xlabel('Zombie alive coefficients',labelpad=10)
    ax[1].set_xlabel('Zombie alive coefficients',labelpad=10)
    ax[2].set_xlabel('Average Zombie alive coefficients',labelpad=10)
    ax[0].set_ylabel('energy',labelpad=10)
    ax[1].set_ylabel('d coefficient',labelpad=10)
    ax[2].set_ylabel('energy',labelpad=10)
    fig.tight_layout()
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    plt.close(fig)
    return

def plot2(results,filename):
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',4)
    fig=plt.figure(figsize=(3.37,5.055))
    ax1=fig.add_subplot(projection='3d') #plt.subplots(2,1,figsize=(3.37,11.11))
    ax1.scatter(results[:,0],results[:,1],results[:,2],color=colors(0))
    # ax1.scatter(results[0,:,0],results[0,:,1],results[0,:,2])
    ax1.set_xlabel('Zombie coefficients')
    ax1.set_ylabel('energy')
    ax1.set_zlabel('dvalue')
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    plt.close(fig)
    return

def plot3(results,average,filename):
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',38)
    fig=plt.figure(figsize=(3.37,5.055))
    ax1=fig.add_subplot(projection='3d') #plt.subplots(2,1,figsize=(3.37,11.11))
    # ax1.scatter(average[0,:,2],average[0,:,0],average[0,:,1],color=colors(0))
    ax1.scatter(average[:,:,2],average[:,:,0],average[:,:,1],color=colors(0))
    # ax1.scatter(results[:,:,2],results[:,:,0],results[:,:,1])
    ax1.set_xlabel('orbital')
    ax1.set_ylabel('Average coefficient')
    ax1.set_zlabel('energy')
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    plt.clf()
    return

results=numpy.zeros(((d_check_inputs.zombs['norb']*2),1000*d_check_inputs.zombs['ndet'],4))
averages=numpy.zeros(((d_check_inputs.zombs['norb']*2),1000,3))
for i in range(d_check_inputs.zombs['norb']*2):
    with open("electron"+f"{(i+1):04d}"+".csv")as file:
        read=numpy.loadtxt(file,delimiter=',')

    results[:,:,3]=i+1
    results[i,:,1:3]=read[:,1:]
    temp=numpy.zeros(1000*d_check_inputs.zombs['ndet'])
    numpy.arcsin(read[:,0],out=temp)
    temp=temp/(2*numpy.pi)
    results[i,:,0]=temp[:]
    # print(results[i,0,:])
    del(read)
    del(temp)
   
    for j in range(1000):
        averages[i,j,0]=numpy.mean(results[i,(j*d_check_inputs.zombs['ndet']):(j*d_check_inputs.zombs['ndet']+d_check_inputs.zombs['ndet']),0])
        averages[i,j,1]=results[i,j*d_check_inputs.zombs['ndet'],1]
        averages[i,j,2]=i+1
    
    
    plot(results[i,:,:],averages[i,:,:],"electron_result"+f"{(i+1):04d}"+".png")
    plot2(results[i,:,:],"electron_3dresult"+f"{(i+1):04d}"+".png")
   
    # plot2(results[i,:,:])
    # exit()

# plot2(results)
plot3(results,averages,"orbtials.png")

command_line="rclone copy /nobackup/cm14oab/"+d_check_inputs.run['runfolder']+" onedrive:PhD/Zombie/arc_results/"+d_check_inputs.run['runfolder']
proc=subprocess.Popen(shlex.split(command_line))
try:
    outs, errs = proc.communicate(timeout=7)
except subprocess.TimeoutExpired:
    proc.kill()
    outs, errs = proc.communicate()
exit()
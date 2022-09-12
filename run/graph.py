import shlex
import subprocess
from matplotlib import pyplot as plt
import matplotlib as mpl
from pylab import cm
import inputs
import numpy 


def plot(eb,rnum, beta, timesteps, filename):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',rnum)
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    if(rnum==1):
        ax.plot(x,eb[1,:], linewidth=2, color=colors(0),label='Converged energy = '+"{:.5f}".format(eb[1,timesteps-1]))
    else:
        for i in range(rnum):
            ax.plot(x,eb[i+1,:], linewidth=2, color=colors(i),label='State '+str(i)+' = '+"{:.5f}".format(eb[i+1,timesteps-1]))
    ax.set_xlim(0,beta)
    ax.legend()
    # ax.set_ylim(-14.8615,-14.8575)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def clean_plot(eb,rnum, beta, clean, timesteps, filename) :
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',rnum+1)
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    for i in range(rnum):
        ax.plot(x,eb[i+1,:], linewidth=2, color=colors(i), label='Converged energy = '+"{:.5f}".format(eb[i+1,timesteps-1]))
    ax.axhline(y=clean, color=colors(rnum), linewidth=2,label='Cleaned energy = '+"{:.5f}".format(clean))
    ax.legend()
    ax.set_xlim(0,beta)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return


if(inputs.run['gram']=='n'):
    with open('energy.csv','rb') as file:
        energy=numpy.loadtxt(file,delimiter=',')
    if((inputs.run['clean']=='y')or(inputs.run['clean']=='f')):
        with open('clean_energy.csv','rb') as file:
            clean= numpy.loadtxt(file,delimiter=',')
        clean_plot(energy,1,inputs.run['beta'],clean[2],inputs.run['timesteps'],'result.png')
    else:
        plot(energy,1,inputs.run['beta'],inputs.run['timesteps'],'result.png')
else:
    energy=numpy.zeros(inputs.run['gramnum']+2,inputs.run['timesteps'])
    with open('energy_state_1.csv','rb') as file:
        obj=numpy.loadtxt(file,delimiter=',')
    energy[0:2,:]=obj
    for i in range(inputs.run['gramnum']):
        with open('energy_state_'+str(i+2)+'.csv','rb') as file:
            obj=numpy.loadtxt(file,delimiter=',')
        energy[i+2,:]=obj[1,:]

    plot(energy,inputs.run['gramnum']+1,inputs.run['beta'],inputs.run['timesteps'],'result.png')


command_line="rclone copy /nobackup/cm14oab/"+inputs.run['runfolder']+" onedrive:PhD/Zombie/arc_results/"+inputs.run['runfolder']
proc=subprocess.Popen(shlex.split(command_line))
try:
    outs, errs = proc.communicate(timeout=7)
except subprocess.TimeoutExpired:
    proc.kill()
    outs, errs = proc.communicate()
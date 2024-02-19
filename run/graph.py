import shlex
import subprocess
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.colors import from_levels_and_colors
from matplotlib.collections import LineCollection
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
            ax.plot(x,eb[i+1,:], linewidth=2, color=colors(i),label='State '+str(i+1)+' = '+"{:.5f}".format(eb[i+1,timesteps-1]))
    ax.set_xlim(0,beta)
    ax.legend()
    # ax.set_ylim(-14.88,-14.830)
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

def gdplote(eb, ebf, beta, timesteps, filename):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',2)
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,eb[1,:], linewidth=2, color=colors(0),label='Energy before GD = '+"{:.5f}".format(eb[1,timesteps-1]))
    ax.plot(x,ebf[1,:], linewidth=2, color=colors(1),label='Energy after GD = '+"{:.5f}".format(ebf[1,timesteps-1]))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def epcplt(epoc,ndet):
    x=epoc[:,0]
    end=x[-1]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',1)
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,epoc[:,1], linewidth=2, color=colors(0))
    ax.set_xlim(0,end)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel('Epoc',labelpad=10)
    plt.savefig('gradient_descent_reduc.png',dpi=300, transparent=False,bbox_inches='tight')
    fig2=plt.figure(figsize=(3.37,5.055))
    ax2=fig2.add_axes([0,0,2,1])
    colors =cm.get_cmap('Set1',ndet-1)
    cmap, norm = from_levels_and_colors([0.0, 0.5, 1.5], ['red', 'black'])
    for i in range(0,int(end+1)):
        for j in range(0,ndet-1):
            print(i,j)
            ax2.scatter(x[i],epoc[i,(2+j)],cmap=cmap, norm=norm)

    ax2.set_xlim(0,x[-1])
    ax2.set_ylim(1,ndet+1)
    ax.set_ylabel('Zombie states altered',labelpad=10)
    ax.set_xlabel('Epoc',labelpad=10)
    plt.savefig('gradient_descent_zs.png',dpi=300, transparent=False,bbox_inches='tight')
    return



if(inputs.run['grad']=='y'):
    with open('energy.csv','rb') as file:
        energy=numpy.loadtxt(file,delimiter=',')
    with open('energy_final.csv','rb') as file:
        fenergy=numpy.loadtxt(file,delimiter=',')
    gdplote(energy,fenergy,inputs.run['beta'],inputs.run['timesteps']+1,'ergresult.png')
    
    with open('epoc.csv','rb')as file:
        epoc=numpy.loadtxt(file,delimiter=',')
    epcplt(epoc,inputs.zombs['ndet'])
else:
    if(inputs.run['gram']=='n'):
        with open('energy.csv','rb') as file:
            energy=numpy.loadtxt(file,delimiter=',')
        if((inputs.run['clean']=='y')or(inputs.run['clean']=='f')):
            with open('clean_energy.csv','rb') as file:
                clean= numpy.loadtxt(file,delimiter=',')
            clean_plot(energy,1,inputs.run['beta'],clean[2],inputs.run['timesteps']+1,'result.png')
        else:
            plot(energy,1,inputs.run['beta'],inputs.run['timesteps']+1,'result.png')
    else:
        energy=numpy.zeros((inputs.run['gramnum']+2,inputs.run['timesteps']+1))
        with open('energy_state_0001.csv','rb') as file:
            obj=numpy.loadtxt(file,delimiter=',')
        energy[0:2,:]=obj
        for i in range(inputs.run['gramnum']):
            with open('energy_state_'+str(i+2).zfill(4)+'.csv','rb') as file:
                obj=numpy.loadtxt(file,delimiter=',')
            energy[i+2,:]=obj[1,:]
        plot(energy,inputs.run['gramnum']+1,inputs.run['beta'],inputs.run['timesteps']+1,'result.png')



command_line="rclone copy /nobackup/cm14oab/"+inputs.run['runfolder']+" onedrive:PhD/Zombie/arc_results/"+inputs.run['runfolder']
proc=subprocess.Popen(shlex.split(command_line))
try:
    outs, errs = proc.communicate(timeout=7)
except subprocess.TimeoutExpired:
    proc.kill()
    outs, errs = proc.communicate()
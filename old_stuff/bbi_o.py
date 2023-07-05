import shlex
import subprocess
from matplotlib import pyplot as plt
import matplotlib as mpl
from pylab import cm
import d_check_inputs
import numpy 
import math

def plot(achange, iters, norb, filename):
    # x=achange[:iters+1,0]
    x=numpy.arange(1,iters+2)
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',1)
    

    for i in range(norb):
        fig=plt.figure(figsize=(3.37,5.055))
        ax=fig.add_axes([0,0,2,1])
        # y=achange[0+(i*iters):(i*iters)+iters+1,i+1]
        y=achange[:,i]
        ax.plot(x,y, linewidth=2, color=colors(1),label='orbital '+str(i+1)+' = '+"{:.5f}".format(y[-1]))
        # ax.legend()
        # plt.tight_layout()
        ax.set_ylabel('coefficient (proxy)',labelpad=10)
        ax.set_xlabel('iteration',labelpad=10)
        plt.savefig(filename+str(i+1)+".png",dpi=300, transparent=False,bbox_inches='tight')
        plt.close(fig)

    
    return

def plot2(achange,norb):
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',1)
    
    x=numpy.arange(1,norb+1)
    y=achange[-1,:]
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,y, linewidth=2, color=colors(1))
    ax.set_ylabel('coefficient (proxy)',labelpad=10)
    ax.set_xlabel('orbital',labelpad=10)
    plt.savefig("final_vals.png",dpi=300, transparent=False,bbox_inches='tight')
    plt.close(fig)
    return


def plot4(ergchange, iters, norb, filename):
    # x=achange[:iters+1,0]
    x=numpy.arange(1,iters+2)
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',1)
    

    for i in range(norb):
        fig=plt.figure(figsize=(3.37,5.055))
        ax=fig.add_axes([0,0,2,1])
        # y=achange[0+(i*iters):(i*iters)+iters+1,i+1]
        y=ergchange[:,i]
        ax.plot(x,y, linewidth=2, color=colors(1),label='orbital '+str(i+1)+' = '+"{:.5f}".format(y[-1]))
        # ax.legend()
        # plt.tight_layout()
        ax.set_ylabel('energy',labelpad=10)
        ax.set_xlabel('iteration',labelpad=10)
        plt.savefig(filename+str(i+1)+".png",dpi=300, transparent=False,bbox_inches='tight')
        plt.close(fig)

    
    return

def plot3(ergchange,iters,passes):
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',2)

    y=numpy.zeros(passes+1)
    x=numpy.arange(0,passes+1)
    for i in range(passes+1):
        print(ergchange[0,(0+(iters+1)*i)])
        y[i]=(ergchange[0,(0+(iters+1)*i)])
    
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,y, linewidth=2, color=colors(0))
    # y=-25.20649144705521*numpy.ones(1+iters*38)
    ax.plot(x,y, linewidth=2, color=colors(2))
    ax.set_ylabel('energy',labelpad=10)
    ax.set_xlabel('passes',labelpad=10)
    plt.savefig("ergoverall.png",dpi=300, transparent=False,bbox_inches='tight')
    plt.close(fig)
    return

with open("achange.csv")as file:
        achange=numpy.loadtxt(file,delimiter=',')

plot(achange,20,38,"achange")
plot2(achange,38)
del(achange)

with open("ergchange.csv")as file:
    ergchange=numpy.loadtxt(file,delimiter=',')

plot4(ergchange,21,38,"ergchange")
plot3(ergchange,20,1)
exit()
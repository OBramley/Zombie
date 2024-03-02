import matplotlib as mpl
from matplotlib import pyplot as plt
from pylab import cm
import subprocess
import numpy 
import json 
import csv
import os

with open('inputs.json') as f:
    inputs=json.load(f)
beta=inputs['run']['beta']
timesteps=inputs['run']['timesteps']

if(inputs['run']['gram'] == 'y'):
    color_num=1+inputs['run']['gramnum']
    gram_num=inputs['run']['gramnum']
    if(inputs['run']['grad']=='y'):
        color_num=color_num*2
    
elif(inputs['run']['gram'] == 'n'):
    color_num=1
    if(inputs['run']['grad']=='y'):
        color_num=color_num+1
    elif(inputs['run']['clean'] in {'y','f'}):
        color_num=color_num+1

colors =cm.get_cmap('Set1',color_num)

def plot_energy(eb, filename,c1):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,eb[1,:], linewidth=2, color=colors(c1),label='Converged energy: '+"{:.5f}".format(eb[1,timesteps-1]))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_gradient_descent(eb, ebf, filename, c1, c2):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,eb[1,:], linewidth=2, color=colors(c1),label='Initial Energy: '+"{:.5f}".format(eb[1,timesteps-1]))
    ax.plot(x,ebf[1,:], linewidth=2, color=colors(c2),label='Energy after Gradient Descent: '+"{:.5f}".format(ebf[1,timesteps-1]))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_cleaned(eb, clean, filename, c1, c2):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,eb[1,:], linewidth=2, color=colors(c1), label='Converged energy: '+"{:.5f}".format(eb[1,timesteps-1]))
    ax.axhline(y=clean, color=colors(c2), linewidth=2,label='Cleaned energy: '+"{:.5f}".format(clean))
    ax.legend()
    ax.set_xlim(0,beta)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_epoch(epoch,filename):
    x=epoch[0,:]
    end=x[-1]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,epoc[:,1], linewidth=2, color=colors(0))
    ax.set_xlim(0,end)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel('Epoc',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_gram(eb, filename,cols):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    for i in range(gram_num):
        ax.plot(x,eb[i+1,:], linewidth=2, color=colors(cols[i]),label='State '+str(i+1))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_gradient_descent_gram(eb, ebf, filename, cols1, cols2):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    for i in range(gram_num):
        ax.plot(x,eb[i+1,:], linewidth=2, color=colors(cols1[i]),label='Initial State '+str(i+1))
        ax.plot(x,ebf[i+1,:], linewidth=2, color=colors(cols2[i]),label='Final State '+str(i+1))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return


if(inputs['run']['gram'] == 'y'):
    with open('energy.csv','rb') as file:
        energy=numpy.loadtxt(file,delimiter=',')
    for i in range(gram_num):
        plot_energy(energy[i+1,:], 'energy_state_'+str(i+1)+'.png', i)
    plot_energy_gram(energy, 'energy.png', range(gram_num))

    if(inputs['run']['grad']=='y'):
        if not os.path.isfile('energy_final.csv'):
            r = csv.reader(open('rundata.csv')) 
            lines = list(r)
            lines[2][0]=2
            with open('/rundata.csv','w',newline='')as file:
                writer = csv.writer(file)
                writer.writerows(lines)
            os.environ["OMP_NUM_THREADS"]=8
            subprocess.run(["./ZOMBIE"])

        with open('energy_final.csv','rb') as file:
            fenergy=numpy.loadtxt(file,delimiter=',')
        for i in range(gram_num):
            plot_energy(fenergy[i+1,:], 'final_energy_state_'+str(i+1)+'.png', gram_num+i)

        plot_energy_gradient_descent_gram(energy, fenergy, 'initial_and_final_energy.png', range(gram_num), range(gram_num,2*gram_num))
        
        for i in range(gram_num):
            with open('epoc_0'+str(i+1)+'.csv','rb')as file:
                epoc=numpy.loadtxt(file,delimiter=',')
            plot_epoch(epoc,'epoc'+str(i+1)+'.png')

elif(inputs['run']['gram'] == 'n'):
    with open('energy.csv','rb') as file:
        energy=numpy.loadtxt(file,delimiter=',')
    plot_energy(energy, 'energy.png', 1)

    if(inputs['run']['grad']=='y'):
        if not os.path.isfile('energy_final.csv'):
            r = csv.reader(open('rundata.csv')) 
            lines = list(r)
            lines[2][0]=2
            with open('/rundata.csv','w',newline='')as file:
                writer = csv.writer(file)
                writer.writerows(lines)
            os.environ["OMP_NUM_THREADS"]=8
            subprocess.run(["./ZOMBIE"])

        with open('energy_final.csv','rb') as file:
            fenergy=numpy.loadtxt(file,delimiter=',')
        plot_energy(fenergy, 'final_energy.png', 2)

        plot_energy_gradient_descent(energy, fenergy, 'initial_and_final_energy.png', 1, 2)

        with open('epoc.csv','rb')as file:
            epoc=numpy.loadtxt(file,delimiter=',')
        plot_epoch(epoc,'epoc.png')

    elif(inputs['run']['clean'] in {'y','f'}):
        with open('clean_energy.csv','rb') as file:
            clean=numpy.loadtxt(file,delimiter=',')
        plot_energy_cleaned(energy, clean[2], 'cleaned_energy.png', 1, 2)


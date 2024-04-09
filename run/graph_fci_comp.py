import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import subprocess
import numpy 
import json 
import csv
import os

fci_erg=-14.871913783299298
fci_erg_gram=[-14.2,12,22,55]
with open('inputs.json') as f:
    inputs=json.load(f)
beta=inputs['run']['beta']
timesteps=inputs['run']['timesteps']

if(inputs['run']['gram'] == 'y'):
    color_num=1+inputs['run']['gramnum']
    gram_num=inputs['run']['gramnum']
    if(inputs['run']['grad']=='y'):
        color_num=color_num*3
    else:
        color_num=color_num*2
    
elif(inputs['run']['gram'] == 'n'):
    color_num=2
    if(inputs['run']['grad']=='y'):
        color_num=color_num+1
    # elif(inputs['run']['clean'] in {'y','f'}):
    #     color_num=color_num+1

colors =mpl.colormaps['Set1']

def plot_energy(eb, filename,c1,c2,comp):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,eb[1,:], linewidth=2, color=colors(c1),label='Converged energy: '+"{:.6f}".format(eb[1,timesteps-1]))
    ax.axhline(y=comp, color=colors(c2), linewidth=2,label='FCI energy: '+"{:.6f}".format(comp))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('% 1.5f')) 
    ax.yaxis.set_minor_formatter(FormatStrFormatter('% 1.5f'))
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_gradient_descent(eb, ebf, filename, c1, c2,c3,comp):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,eb[1,:], linewidth=2, color=colors(c1),label='Initial Energy: '+"{:6f}".format(eb[1,timesteps-1]))
    ax.plot(x,ebf[1,:], linewidth=2, color=colors(c2),label='Energy after Gradient Descent: '+"{:.6f}".format(ebf[1,timesteps-1]))
    ax.axhline(y=comp, color=colors(c3), linewidth=2,label='FCI energy: '+"{:.f}".format(comp))
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('% 1.5f')) 
    ax.yaxis.set_minor_formatter(FormatStrFormatter('% 1.5f'))
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

# def plot_energy_cleaned(eb, clean, filename, c1, c2):
#     x=eb[0,:]
#     mpl.rcParams['font.family']='DejaVu Sans'
#     plt.rcParams['font.size']=18
#     plt.rcParams['axes.linewidth']=2
#     fig=plt.figure(figsize=(3.37,5.055))
#     ax=fig.add_axes([0,0,2,1])
#     ax.plot(x,eb[1,:], linewidth=2, color=colors(c1), label='Converged energy: '+"{:.6f}".format(eb[1,timesteps-1]))
#     ax.axhline(y=clean, color=colors(c2), linewidth=2,label='Cleaned energy: '+"{:.6f}".format(clean))
#     ax.legend()
#     ax.set_xlim(0,beta)
#     ax.set_ylabel('Energy [au]',labelpad=10)
#     ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
#     plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
#     return

def plot_epoch(epoch,filename,c1,c2,comp):
    x=epoch[:,0]
    end=x[-1]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    ax.plot(x,epoc[:,1], linewidth=2, color=colors(c1))
    ax.axhline(y=comp, color=colors(c2), linewidth=2,label='FCI energy: '+"{:.6f}".format(comp))
    ax.set_xlim(0,end)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel('Epochs',labelpad=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('% 1.5f')) 
    ax.yaxis.set_minor_formatter(FormatStrFormatter('% 1.5f'))
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_gram(eb, filename,cols,comp_cols,comp):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    for i in range(gram_num):
        ax.plot(x,eb[i+1,:], linewidth=2, color=colors(cols[i]),label='State '+str(i+1))
        ax.axhline(y=comp[i], color=colors(comp_cols[i]), linewidth=2,label='State '+str(i+1)+' FCI energy')
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('% 1.6f')) 
    ax.yaxis.set_minor_formatter(FormatStrFormatter('% 1.6f'))
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return

def plot_energy_gradient_descent_gram(eb, ebf, filename, cols1, cols2,comp_cols,comp):
    x=eb[0,:]
    mpl.rcParams['font.family']='DejaVu Sans'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    for i in range(gram_num):
        ax.plot(x,eb[i+1,:], linewidth=2, color=colors(cols1[i]),label='Initial State '+str(i+1))
        ax.plot(x,ebf[i+1,:], linewidth=2, color=colors(cols2[i]),label='Final State '+str(i+1))
        ax.axhline(y=comp[i], color=colors(comp_cols[i]), linewidth=2,label='State '+str(i+1)+' FCI energy')
    ax.set_xlim(0,beta)
    ax.legend()
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    ax.yaxis.set_major_formatter(FormatStrFormatter('% 1.6f')) 
    ax.yaxis.set_minor_formatter(FormatStrFormatter('% 1.6f'))
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')
    return


if(inputs['run']['gram'] == 'y'):
    with open('energy.csv','rb') as file:
        energy=numpy.loadtxt(file,delimiter=',')
    for i in range(gram_num):
        plot_energy(energy[i+1,:], 'energy_state_'+str(i+1)+'_comp.png', i,gram_num+i,fci_erg_gram[i])
    plot_energy_gram(energy, 'energy_comp.png', range(gram_num), range(gram_num,2*gram_num), fci_erg_gram)

    if(inputs['run']['grad']=='y'):
        if not os.path.isfile('energy_final.csv'):
            value=input("File does not exists do you want to create it? y/n\n")
            if(value=='y'):
                r = csv.reader(open('rundata.csv')) 
                lines = list(r)
                lines[2][0]=2
                with open('rundata.csv','w',newline='')as file:
                    writer = csv.writer(file)
                    writer.writerows(lines)
                os.environ["OMP_NUM_THREADS"]="8"
                subprocess.run(["./ZOMBIE"])
        if os.path.isfile('energy_final.csv'):
            with open('energy_final.csv','rb') as file:
                fenergy=numpy.loadtxt(file,delimiter=',')
            for i in range(gram_num):
                plot_energy(fenergy[i+1,:], 'final_energy_state_'+str(i+1)+'_comp.png', 2*gram_num+i,gram_num+i,fci_erg_gram[i])

            plot_energy_gradient_descent_gram(energy, fenergy, 'initial_and_final_energy_comp.png', range(gram_num),range(2*gram_num,3*gram_num),range(gram_num,2*gram_num),fci_erg_gram)
        
        for i in range(gram_num):
            with open('epoc_0'+str(i+1)+'.csv','rb')as file:
                epoc=numpy.loadtxt(file,delimiter=',')
            plot_epoch(epoc,'epoc'+str(i+1)+'_comp.png',2*gram_num+i,gram_num+i,fci_erg_gram[i])

elif(inputs['run']['gram'] == 'n'):
    with open('energy.csv','rb') as file:
        energy=numpy.loadtxt(file,delimiter=',')
    plot_energy(energy, 'energy_comp.png', 1,3,fci_erg)

    if(inputs['run']['grad']=='y'):
        if not os.path.isfile('energy_final.csv'):
            value=input("File does not exists do you want to create it? y/n\n")
            if(value=='y'):
                r = csv.reader(open('rundata.csv')) 
                lines = list(r)
                lines[2][0]=2
                with open('rundata.csv','w',newline='')as file:
                    writer = csv.writer(file)
                    writer.writerows(lines)
                os.environ["OMP_NUM_THREADS"]="8"
                subprocess.run(["./ZOMBIE"])
        if os.path.isfile('energy_final.csv'):
            with open('energy_final.csv','rb') as file:
                fenergy=numpy.loadtxt(file,delimiter=',')
            plot_energy(fenergy, 'final_energy_comp.png', 2,3,fci_erg)

            plot_energy_gradient_descent(energy, fenergy, 'initial_and_final_energy_comp.png', 1, 2, 3,fci_erg)
        
        with open('epoc.csv', newline='') as csvfile:
            data = list(csv.reader(csvfile))
        epoc2=[]
        i=0
        for line in data:
            epoc2.append([i,float(line[1])])
            i=i+1
        epoc=numpy.array(epoc2)
        plot_epoch(epoc,'epoc_comp.png',2,3,fci_erg)

    # elif(inputs['run']['clean'] in {'y','f'}):
    #     with open('clean_energy.csv','rb') as file:
    #         clean=numpy.loadtxt(file,delimiter=',')
    #     plot_energy_cleaned(energy, clean[2], 'cleaned_energy.png', 1, 2)


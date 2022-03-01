from fileinput import filename
import _pickle as pickle
import numpy
import csv
from matplotlib import pyplot as plt
import matplotlib as mpl
from pylab import cm

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, -1)

def read_object(filename):
    with open(filename, 'rb') as inp:
        obj = pickle.load(inp)
    return obj

def write_ham(obj,filename):
    with open(filename,'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerows(obj)

def write_dvec(dvec,filename):
    with open(filename,'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerow(dvec)
    

def read_ham(filename):
    with open(filename, 'rb')as file:
        obj = numpy.loadtxt(file,delimiter=',')
    return obj

def plot(eb,rnum, beta, filename):
    x=eb[:,0]
    mpl.rcParams['font.family']='Avenir'
    plt.rcParams['font.size']=18
    plt.rcParams['axes.linewidth']=2
    colors =cm.get_cmap('Set1',rnum)
    fig=plt.figure(figsize=(3.37,5.055))
    ax=fig.add_axes([0,0,2,1])
    for i in range(rnum):
        ax.plot(x,eb[:,i+1], linewidth=2, color=colors(i))
    ax.set_xlim(0,beta)
    # ax.set_ylim(-14.8615,-14.8575)
    ax.set_ylabel('Energy [au]',labelpad=10)
    ax.set_xlabel(r'$\mathregular{\beta}$',labelpad=10)
    plt.savefig(filename,dpi=300, transparent=False,bbox_inches='tight')



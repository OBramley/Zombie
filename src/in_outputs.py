import _pickle as pickle
import numpy
import csv

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def read_object(filename):
    with open(filename, 'rb') as inp:
        obj = pickle.load(inp)
    return obj

def write_integrals(Hnuc, H1ei, H2ei,norb):
    with open('electron_integrals_1.csv','w',newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerow([Hnuc,0,0,0,0,0,0,0,0,0])
        spamwriter.writerows(H1ei)
        for i in range(norb):
            for j in range(norb):
                spamwriter.writerows(H2ei[i,j,:,:])

def read_in_integrals(norb, filename):
    H2ei =  numpy.zeros((norb,norb,norb,norb))
    file = open(filename)
    with open(filename,'rb')as file:
        numpy_array=numpy.loadtxt(file,delimiter=',')
        Hnuc=numpy_array[0,0]
        H1ei=numpy_array[1:norb+1,:]
        counter=norb+1
        for i in range(norb):
            for j in range(norb):
                H2ei[i,j,0:norb,0:norb]=numpy_array[counter:(counter+norb),:]
                counter=+norb

    
    return Hnuc,H1ei,H2ei

def write_ham(obj,filename):
    with open(filename,'w', newline='')as csvfile:
        spamwriter=csv.writer(csvfile, delimiter=',')
        spamwriter.writerows(obj)

def read_ham(filename):
    with open(filename, 'rb')as file:
        obj = numpy.loadtxt(file,delimiter=',')
    return obj


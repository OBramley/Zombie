import numpy
import op
import in_outputs

def hamiltonian(ndet,Ham,zstore):
    Kover = numpy.zeros((ndet,ndet))
    Bigham = numpy.zeros((ndet,ndet))
    for idet in range(ndet):
        for jdet in range(idet,ndet):
            Kover[idet,jdet] = op.overlap_f(zstore[idet].zs, zstore[jdet].zs)
            Kover[jdet,idet] = Kover[idet,jdet]
            # fover.write('{:5d}, {:5d}, {:<20.15f}\n'.format(idet,jdet,Kover[idet,jdet]))
            Bigham[idet,jdet] = Ham.HTot(zstore[idet].zs, zstore[jdet].zs)
            Bigham[jdet,idet] = Ham[idet,jdet]
            # ranham.write('{:5d}, {:5d}, {:<20.15}, {:<20.15}\n'.\
            #             format(idet,jdet,Kover[idet,jdet],Ham[idet,jdet]))

    in_outputs.write_ham(Kover,'kover.csv')
    in_outputs.write_ham(Bigham,'bigham.csv')
    return Bigham, Kover
    
#from read_iterdb_file import *
from read_iterdb_x import *
from finite_differences_x import *
#from calc_fields_from_EFIT import first_derivative
from interp import *
import matplotlib.pyplot as plt
import numpy as np
from read_EFIT import *

def dProfdR(ITERDBdict, EFITdict, profile, rhotLeft, rhotRight):
    rhot0 = EFITdict['rhotn']
    R0 = EFITdict['R']

    rhotn = ITERDBdict['rhot_'+profile]
    Prof = ITERDBdict[profile]

    leftInd = np.argmin(abs(rhotn - rhotLeft))
    rightInd = np.argmin(abs(rhotn - rhotRight))

    rhotProf = rhotn[leftInd: rightInd + 1]
    prof = Prof[leftInd: rightInd + 1]

    R = interp(rhot0, R0, rhotProf)
    uniR = np.linspace(R[0], R[-1], len(R))
    rhot_uniR = interp(R, rhotProf, uniR)
    prof_uniR = interp(R, prof, uniR)
    profPrime_uniR = - first_derivative(prof_uniR, uniR) / prof_uniR

    return rhot_uniR, prof_uniR, profPrime_uniR

    
iterdbFileName = "negOmTorDIIID98889.iterdb"
efitFileName = "g098889.04530"
EFITdict = read_EFIT(efitFileName)
#rhot0 = EFITdict['rhotn'][100:]
#R0 = EFITdict['R'][100:]

Zave = 6.
Bref_Gauss = abs(EFITdict['bcentr']) * 1.E04    #1.95451E04
Lref_m = 0.76
mref = 2.

ITERDBdict = read_iterdb_x(iterdbFileName)
#rhot, te, ti, ne, ni, nz, vrot = read_iterdb_file(file1)

#fs = 0.975
#fsInd = np.argmin(abs(rhot1 - fs))

e = 1.6E-19


if 1 == 1:
    if 1 == 1:
#        R = interp(rhot0, R0, rhot)
#        R2 = np.linspace(R[0], R[-1], len(R))
#        rhot2 = interp(R, rhot, R2)
#        te2 = interp(R, te, R2)
#        ti2 = interp(R, ti, R2)
#        ne2 = interp(R, ne, R2)
#        ni2 = interp(R, ni, R2)
#        nz2 = interp(R, nz, R2)
#        vrot2 = interp(R, vrot, R2)

#        tprime_i2 = -first_derivative(ti2, R2)/ti2
#        tprime_e2 = -first_derivative(te2, R2)/te2
#        fprime_i2 = -first_derivative(ni2, R2)/ni2
#        fprime_e2 = -first_derivative(ne2, R2)/ne2
#        fprime_z2 = -first_derivative(nz2, R2)/nz2

#def dProfdR(ITERDBdict, EFITdict, profile, rhotLeft, rhotRight):
        rhot2, te, tprime_e = dProfdR(ITERDBdict, EFITdict, 'te', 0.9, 1.0)
        rhot2, ti, tprime_i = dProfdR(ITERDBdict, EFITdict, 'ti', 0.9, 1.0)
        rhot2, ni, fprime_i = dProfdR(ITERDBdict, EFITdict, 'ni', 0.9, 1.0)
        rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(te)/Bref_Gauss
        result2 = rhoi_m * (fprime_i + tprime_i)
    if 1 == 0:
        rhot1 = np.linspace(rhot[0], rhot[-1], len(rhot))
        te1 = interp(rhot, te, rhot1)
        ti1 = interp(rhot, ti, rhot1)
        ne1 = interp(rhot, ne, rhot1)
        ni1 = interp(rhot, ni, rhot1)
        nz1 = interp(rhot, nz, rhot1)

        vrot1 = interp(rhot, vrot, rhot1)
        tprime_i1 = -first_derivative(ti1, rhot1)/ti1
        tprime_e1 = -first_derivative(te1, rhot1)/te1
        fprime_i1 = -first_derivative(ni1, rhot1)/ni1
        fprime_e1 = -first_derivative(ne1, rhot1)/ne1
        fprime_z1 = -first_derivative(nz1, rhot1)/nz1

        rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(te1)/Bref_Gauss
        result1 = rhoi_m * (fprime_i1 + tprime_i1) / Lref_m

    #ptot1 = te1 * ne1 + ti1 * (ni1 + nz1)

    #plt.plot(rhot1, fprime_i1, label='omn_i')
    #plt.plot(rhot1, fprime_e1, label='omn_e')
    #plt.xlabel('rhot')
    ##plt.axis([viewLeft,viewRight,0,100])
    #plt.legend(loc = 2)
    #plt.show()
    #plt.plot(rhot1, tprime_i1, label='omt_i')
    #plt.plot(rhot1, tprime_e1, label='omt_e')
    #plt.xlabel('rhot')
    #plt.axis([viewLeft,viewRight,0,100])
    #plt.legend(loc = 2)
    #plt.show()
    #plt.plot(rhot1, tprime_i1/fprime_i1, label='eta_i')
    #plt.plot(rhot1, tprime_e1/fprime_e1, label='eta_e')
    #plt.xlabel('rhot')
    #plt.legend(loc = 2)
    #plt.axis([viewLeft,viewRight,0,10])
    #plt.show()

    #plt.plot(rhot1, result1, label = 'radial grid: sqrt(psi_tor)')
    plt.plot(rhot2, result2, label = 'radial grid: major radius (m)')
    plt.xlabel('rhot')
    plt.legend(loc = 2)
    plt.axis([0.9,1.0,0.,0.35])
    plt.title('rho_i / L for DIIID' + efitFileName[2:7])
    plt.show()

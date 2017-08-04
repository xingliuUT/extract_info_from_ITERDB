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

    
def dProfdrhot(ITERDBdict, EFITdict, profile, rhotLeft, rhotRight):
    rhot0 = EFITdict['rhotn']
    R0 = EFITdict['R']

    rhotn = ITERDBdict['rhot_'+profile]
    Prof = ITERDBdict[profile]

    leftInd = np.argmin(abs(rhotn - rhotLeft))
    rightInd = np.argmin(abs(rhotn - rhotRight))

    rhotProf = rhotn[leftInd: rightInd + 1]
    prof = Prof[leftInd: rightInd + 1]

    uni_rhot = np.linspace(rhotProf[0], rhotProf[-1], len(rhotProf))
    prof_unirhot = interp(rhotProf, prof, uni_rhot)
    profPrime_unirhot = - first_derivative(prof_unirhot, uni_rhot) / prof_unirhot

    return uni_rhot, prof_unirhot, profPrime_unirhot

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
        rhot2, te, tprime_e = dProfdR(ITERDBdict, EFITdict, 'te', 0.9, 1.0)
        rhot2, ti, tprime_i = dProfdR(ITERDBdict, EFITdict, 'ti', 0.9, 1.0)
        rhot2, ni, fprime_i = dProfdR(ITERDBdict, EFITdict, 'ni', 0.9, 1.0)
        rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(te)/Bref_Gauss
        result2 = rhoi_m * (fprime_i + tprime_i)
    if 1 == 1:
        rhot1, te, tprime_e = dProfdrhot(ITERDBdict, EFITdict, 'te', 0.9, 1.0)
        rhot1, ti, tprime_i = dProfdrhot(ITERDBdict, EFITdict, 'ti', 0.9, 1.0)
        rhot1, ni, fprime_i = dProfdrhot(ITERDBdict, EFITdict, 'ni', 0.9, 1.0)

        rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(te)/Bref_Gauss
        result1 = rhoi_m * (fprime_i + tprime_i) / Lref_m

    plt.plot(rhot1, result1, label = 'derivative wrt sqrt(psi_tor)')
    plt.plot(rhot2, result2, label = 'derivative wrt major radius (m)')
    plt.xlabel('rhot')
    plt.legend(loc = 2)
    plt.axis([0.9,1.0,0.,0.35])
    plt.title('rho_i / L for DIIID' + efitFileName[2:7])
    plt.show()

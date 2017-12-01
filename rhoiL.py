#from read_iterdb_file import *
from read_iterdb_x import *
from finite_differences_x import *
#from calc_fields_from_EFIT import first_derivative
from interp import *
import matplotlib.pyplot as plt
import numpy as np
from read_EFIT import *
import sys

def dPressuredR(ITERDBdict, EFITdict, rhotLeft, rhotRight):
    rhot0 = EFITdict['rhotn']
    R0 = EFITdict['R']
    Pres = EFITdict['Pres']
    print(len(rhot0), len(R0), len(Pres))

    leftInd = np.argmin(abs(rhot0 - rhotLeft))
    rightInd = np.argmin(abs(rhot0 - rhotRight))

    rhotProf = rhot0[leftInd: rightInd + 1]
    prof = Pres[leftInd: rightInd + 1]
    R = R0[leftInd: rightInd + 1]

    uniR = np.linspace(R[0], R[-1], len(R))
    rhot_uniR = interp(R, rhotProf, uniR)
    prof_uniR = interp(R, prof, uniR)
    if 1 == 0:
        plt.plot(R0, rhot0, '+-', label = 'rhot0')
        plt.plot(uniR, rhot_uniR, '+-', label = 'rhot_uniR')
        plt.xlabel('R')
        plt.ylabel('rhot')
        plt.legend()
        plt.show()

        plt.plot(R0, Pres, '+-', label = 'Pres')
        plt.plot(uniR, prof_uniR, '+-', label = 'prof_uniR')
        plt.xlabel('R')
        plt.ylabel('pressure')
        plt.legend()
        plt.show()

    profPrime_uniR = - first_derivative(prof_uniR, uniR) / prof_uniR

    return rhot_uniR, uniR, prof_uniR, profPrime_uniR

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

#iterdbFileName = "negOmTorDIIID98889.iterdb"
#efitFileName = "g098889.04530"
iterdbFileName = sys.argv[1]
efitFileName = sys.argv[2]
EFITdict = read_EFIT(efitFileName)
#rhot0 = EFITdict['rhotn'][100:]
#R0 = EFITdict['R'][100:]

Zave = 6.
Bref_Gauss = abs(EFITdict['bcentr']) * 1.E04    #1.95451E04
Lref_m = 0.764  #off7 = 0.7567 #Callen = 0.764, Diallo = 0.756
mref = 2.

ITERDBdict = read_iterdb_x(iterdbFileName)
#rhot, te, ti, ne, ni, nz, vrot = read_iterdb_file(file1)

#fs = 0.975
#fsInd = np.argmin(abs(rhot1 - fs))

e = 1.6E-19

if 1 == 1:
    uni_rhot, shat, Ls = magneticShear(EFITdict)
    rhot_uniR, uni_R, pres, dpresdR = dPressuredR(ITERDBdict, EFITdict, 0.9, 1.0)
    Lped = 1. / dpresdR


    if 1 == 1:
        # convert Ls onto uniform R grid
        Ls_uniR = interp(uni_rhot, Ls, rhot_uniR)
        # convert Lped onto uniform rhot grid
        Lped_unirhot = interp(rhot_uniR, Lped, uni_rhot)
        if 1 == 0:
            plt.plot(rhot_uniR, Lped, '+-')
            plt.plot(uni_rhot, Lped_unirhot, '+-')
            plt.ylabel('Lped')
            plt.xlabel('rhot')
            plt.xlim([0.9, 1.])
            plt.show()

        fig, ax = plt.subplots(figsize = (6, 4.5))
        y1 = mref * 1836 * (Lped_unirhot / Ls)**2
        ax.semilogy(uni_rhot, y1, 'o-', linewidth = 2., ms = 4., \
        label = r'$\frac{m_i}{m_e} (\frac{L_{ped}}{L_s})^2$')
        ax.set_xlabel(r'$\rho_{t}$', fontsize = 12)
        ax.set_xlim([0.9, 1.])
        ax.set_ylim([1E-5, 20.])
        #ax.set_xlim([0.94, 1.])
        #ax.set_ylim([5E-6, 10.])
        ax.semilogy(rhot_uniR, Lped / uni_R, 'o-', linewidth = 2., \
        ms = 4., label = r'$\frac{L_{ped}}{R}$')
#        ax.set_title(efitFileName)
        ax.legend(fontsize = 'large')
        ax.tick_params(axis='both',which='major',labelsize=12,\
        length=3,width=2,direction='out')
        plt.savefig(filename = 'Callen_log_plot1.pdf', format = 'pdf')
        #plt.savefig(filename = 'Diallo_off7_log_plot1.pdf', format = 'pdf')
        #plt.show()

if 1 == 0:
    if 1 == 1:
        rhot2, te, tprime_e = dProfdR(ITERDBdict, EFITdict, 'te', 0.9, 1.0)
        rhot2, ti, tprime_i = dProfdR(ITERDBdict, EFITdict, 'ti', 0.9, 1.0)
        rhot2, ni, fprime_i = dProfdR(ITERDBdict, EFITdict, 'ni', 0.9, 1.0)
        #rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(te)/Bref_Gauss
        rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(ti)/Bref_Gauss
        result2 = rhoi_m * (fprime_i + tprime_i)
    if 1 == 1:
        rhot1, te, tprime_e = dProfdrhot(ITERDBdict, EFITdict, 'te', 0.9, 1.0)
        rhot1, ti, tprime_i = dProfdrhot(ITERDBdict, EFITdict, 'ti', 0.9, 1.0)
        rhot1, ni, fprime_i = dProfdrhot(ITERDBdict, EFITdict, 'ni', 0.9, 1.0)
        #rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(te)/Bref_Gauss
        rhoi_m = 1.02*np.sqrt(mref)*np.sqrt(ti)/Bref_Gauss
        result1 = rhoi_m * (fprime_i + tprime_i) / Lref_m

    plt.plot(rhot1, result1, label = 'derivative wrt sqrt(psi_tor)')
    plt.plot(rhot2, result2, label = 'derivative wrt major radius (m)')
    plt.xlabel('rhot')
    plt.legend(loc = 2)
    plt.axis([0.9,1.0,0.,0.5])
    plt.title('rho_i / L_i for DIIID' + efitFileName[1:7])
#    plt.title('rho_i / L_i for DIIID off7')
    plt.show()
#    plt.savefig('Callen_rhoiLi.pdf', format = 'pdf')

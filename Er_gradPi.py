# python Er_gradPi.py newDIIID153764.iterdb g153764.01359 Diallo_exp_er.txt
# python Er_gradPi.py fullDIIID98889.iterdb g098889.04530 Callen_exp_er.txt
# python Er_gradPi.py CMOD1120815027May08.iterdb g1120815027.01075_001 CMOD_Hmode_exp_er.txt
# python Er_gradPi.py CMOD1120907032May12.iterdb g1120907032.01012 CMOD_Imode_exp_er.txt

from read_iterdb_x import *
from read_EFIT import *
from finite_differences_x import *
from interp import *
import matplotlib.pyplot as plt
import numpy as np
import sys

iterdbFileName = sys.argv[1]
efitFileName = sys.argv[2]
er_data = np.genfromtxt(sys.argv[3])
use_psip = False

EFITdict = read_EFIT(efitFileName)
ITERDBdict = read_iterdb_x(iterdbFileName)

Rgrid = EFITdict['R']

uniR = np.linspace(Rgrid[0], Rgrid[-1], len(Rgrid) * 8)
rhot_uniR = interp(Rgrid, EFITdict['rhotn'], uniR)
psip_uniR = interp(Rgrid, EFITdict['psipn'], uniR)

if 1 == 0:
    plt.plot(ITERDBdict['rhot_ni'], ITERDBdict['ni'] * 1E-19, label = 'ni (1E19)')
    plt.plot(ITERDBdict['rhot_ti'], ITERDBdict['ti'] * 1E-03, label = 'ti (1E3)')
    plt.plot(ITERDBdict['rhot_ne'], ITERDBdict['ne'] * 1E-19, label = 'ne (1E19)')
    plt.plot(ITERDBdict['rhot_te'], ITERDBdict['te'] * 1E-03, label = 'te (1E3)')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(Rgrid, EFITdict['rhotn'], label = 'before')
    plt.plot(uniR, rhot_uniR, '.', label = 'after')
    plt.xlabel('R')
    plt.ylabel('rhot')
    plt.legend()
    plt.show()

if 1 == 0:
    plt.plot(Rgrid, EFITdict['psipn'], label = 'before')
    plt.plot(uniR, psip_uniR, '.', label = 'after')
    plt.xlabel('R')
    plt.ylabel('psip')
    plt.legend()
    plt.show()

# ti in ITERDB file is in the unit of eV
ti_uniR = interp(ITERDBdict['rhot_ti'], ITERDBdict['ti'], rhot_uniR)
ni_uniR = interp(ITERDBdict['rhot_ni'], ITERDBdict['ni'], rhot_uniR)
pi_uniR = ti_uniR * ni_uniR

te_uniR = interp(ITERDBdict['rhot_te'], ITERDBdict['te'], rhot_uniR)
ne_uniR = interp(ITERDBdict['rhot_ne'], ITERDBdict['ne'], rhot_uniR)
pe_uniR = te_uniR * ne_uniR

if 1 == 0:
    plt.plot(ITERDBdict['rhot_ti'], ITERDBdict['ti'], label = 'before')
    plt.plot(rhot_uniR, ti_uniR, '.', label = 'after')
    plt.xlabel('rhot')
    plt.ylabel('ti')
    plt.legend()
    plt.show()

if 1 == 0:
    plt.plot(ITERDBdict['rhot_ni'], ITERDBdict['ni'], label = 'before')
    plt.plot(rhot_uniR, ni_uniR, '.', label = 'after')
    plt.xlabel('rhot')
    plt.ylabel('ni')
    plt.legend()
    plt.show()

gradPioverNie = first_derivative(pi_uniR, uniR) / ni_uniR
gradPeoverNee = first_derivative(pe_uniR, uniR) / ne_uniR

if use_psip:
    er_uniR = interp(er_data[:, 0], er_data[:, 1], psip_uniR)
else:
    er_uniR = interp(er_data[:, 0], er_data[:, 1], rhot_uniR)

if 1 == 1:
    plt.plot(rhot_uniR, gradPioverNie/1000., label = r'$\frac{\nabla P_i}{e n_i}$')
    plt.plot(rhot_uniR, gradPeoverNee/1000., label = r'$\frac{\nabla P_e}{e n_e}$')
    plt.plot(rhot_uniR, er_uniR, label = r'$E_r$')
    plt.xlabel('rhot')
    plt.ylabel('kV/m')
    plt.legend(fontsize = 'large')
    plt.xlim([0.9, 1.])
    plt.ylim([-300., 50.])
    plt.show()
if 1 == 0:
    plt.plot(uniR, gradPioverNie/1000., label = r'$\frac{\nabla P_i}{e n_i}$')
    plt.plot(uniR, gradPeoverNee/1000., label = r'$\frac{\nabla P_e}{e n_e}$')
    plt.plot(uniR, er_uniR, label = r'$E_r$')
    plt.xlabel('R')
    plt.ylabel('kV/m')
    plt.legend(fontsize = 'large')
    plt.show()

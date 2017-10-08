import sys
import numpy as np
import matplotlib.pyplot as plt
from interp import *
from finite_differences_x import *
from read_iterdb_x import *
from read_EFIT import *
from calc_volume_from_EFIT import *
from calc_fields_from_EFIT import *

EFITdict = read_EFIT(sys.argv[1])
ITERDBdict = read_iterdb_x(sys.argv[2])

rhotn = EFITdict['rhotn']
psipn = EFITdict['psipn']

e_charge = 1.6E-19
zave = 6.
zeff = 1.5

te0 = ITERDBdict['te']
rhot_te0 = ITERDBdict['rhot_te']
ti0 = ITERDBdict['ti']
rhot_ti0 = ITERDBdict['rhot_ti']
ne0 = ITERDBdict['ne']
rhot_ne0 = ITERDBdict['rhot_ne']
ni0 = ITERDBdict['ni']
rhot_ni0 = ITERDBdict['rhot_ni']
nz0 = ITERDBdict['nz']
rhot_nz0 = ITERDBdict['rhot_nz']

te = interp(rhot_te0, te0, rhotn)
ti = interp(rhot_ti0, ti0, rhotn)
ne = interp(rhot_ne0, ne0, rhotn)
ni1 = interp(rhot_ni0, ni0, rhotn)
nz1 = interp(rhot_nz0, nz0, rhotn)

ni = ne * (zave - zeff) / (zave - 1.)
nz = ne * (zeff - 1.) / (zave - 1.) / zave

pressure_MKS = (ne * te + (ni + nz) * ti) * e_charge

if 1 == 0:
    plt.plot(psipn, EFITdict['Pres'], label = 'Pressure EFIT')
    plt.plot(psipn, pressure_MKS, label = 'Pressure ITERDB')
    plt.xlabel('psip')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(psipn, ne * 1E-20, label = 'ne (10^20 m^-3)')
    plt.plot(psipn, ti * 1E-3, label = 'ti (KeV)')
    plt.plot(psipn, te * 1E-3, label = 'te (KeV)')
    plt.xlabel('psip')
    plt.legend()
    plt.show()

if 1 == 0:
    f = open('inputs.txt','w')
    f.write('#1.psi_pol 2.pressure(MKS) 3.ne(10^20 m^-3) 4.ti(KeV) 5.te(KeV)\n')
    np.savetxt(f,np.column_stack((psipn, pressure_MKS, ne * 1E-20, ti * 1E-3, te * 1E-3)))
    f.close()

psip_fs = 0.995
ntheta = 512
R_fs, Z_fs, B_pol, B_tor, B_tot = BfieldsFS(EFITdict, psip_fs, ntheta, True)
if 1 == 1:
    f = open('rz_formatted_0995.txt','w')
    f.write('#1.R 2.Z\n')
    np.savetxt(f,np.column_stack((R_fs, Z_fs)), fmt = '%.18e', delimiter = ' ')

from read_EFIT import *
from read_iterdb_x import *
from finite_differences_x import *
import matplotlib.pyplot as plt
from interp import *
import sys


iterdb_file_name = sys.argv[1]
efitFileName = sys.argv[2]
EFITdict = read_EFIT(efitFileName)

# Callen
# aGENE_m = 0.76378
# n0_global = 16
# kymin = 0.1511
# x0 = 0.9725

# Diallo off7
aGENE_m = 0.7567
n0_global = 13
kymin = 0.1032
x0 = 0.97

e = 1.6*10**(-19)
mref = 2.
M_kg = 3.3*10**(-27)

#rhot0, te0, ti0, ne0, ni0, nz0, vrot0 = read_iterdb_file(iterdb_file_name)
ITERDBdict = read_iterdb_x(iterdb_file_name)
print(list(ITERDBdict.keys()))
#['rhot_te', 'te', 'rhot_ti', 'ti', 'rhot_ne', 'ne', \
#'rhot_ni', 'ni', 'rhot_nz', 'nz', 'rhot_vrot', 'vrot']

uni_rhot = np.linspace(min(ITERDBdict['rhot_te']), \
max(ITERDBdict['rhot_te']), len(ITERDBdict['rhot_te']) * 10)

te_u = interp(ITERDBdict['rhot_te'], ITERDBdict['te'], uni_rhot)
ne_u = interp(ITERDBdict['rhot_ne'], ITERDBdict['ne'], uni_rhot)
vrot_u = interp(ITERDBdict['rhot_vrot'], ITERDBdict['vrot'], uni_rhot)

tprime_e = -first_derivative(te_u,uni_rhot)/te_u
nprime_e = -first_derivative(ne_u,uni_rhot)/ne_u

#kyGENE = 102*np.sqrt(mref)*np.sqrt(te_u)/Bref_Gauss*k_theta_cm*k2Factor
x0Ind = np.argmin(abs(uni_rhot - x0))
te_mid = te_u[x0Ind]
kyGENE = kymin * np.sqrt(te_u / te_mid)
omMTM = kyGENE * (tprime_e + nprime_e)
vref = 9.79E3 / np.sqrt(mref) * np.sqrt(te_u)
freqref = vref /aGENE_m
mtmFreq = omMTM * freqref / 2. / np.pi / 1000.

# vrot_u is in the unit of rad/s
omegaDoppler = vrot_u * n0_global / 2. / np.pi / 1E3

if 1 == 1:
    x0scan = np.arange(0.96, 0.985, 0.005)
    for xi in x0scan:
        index = np.argmin(abs(uni_rhot - xi))
        print("x0 =", uni_rhot[index])
        print("omega_ExB_kHz =", omegaDoppler[index])
    plt.plot(uni_rhot, omegaDoppler, label = 'omega_ExB (kHz)')
    plt.xlabel('rhot')
    plt.axis([0.9, 1., -100, 100])
    plt.title(iterdb_file_name + ', n0_global = ' + str(n0_global))
    plt.legend()
    plt.show()

#pedtopInd = np.argmin(abs(EFITdict['rhotn'] - 0.93)) + 1
pedtopInd = np.argmin(abs(EFITdict['rhotn'] - 0.89)) + 1
if 1 == 0:
    fig, ax = plt.subplots(figsize = (6,4.5))
    ax.plot(uni_rhot, omegaDoppler / mtmFreq, \
    linewidth = 2., label=r'$\frac{\omega \,\, in \,\, plasma \,\, frame}{\omega^*_e}$')
    ax.plot(EFITdict['rhotn'], \
    EFITdict['Pres'] / EFITdict['Pres'][pedtopInd], \
    linewidth = 2., label = 'Pressure (a.u.)')
    ax.axhline(y = 0.,\
        color='k',linewidth = 2.5)
    ax.axis([0.9, 1., -2., 1.])
    #ax.axis([0.94, 1., -1.5,1.])
    ax.set_xlabel(r'$\rho_t$', fontsize = 14)
    ax.legend(fontsize = 'x-large')
    #plt.title(iterdb_file_name)
    ax.tick_params(axis='both',which='major',labelsize=12,\
        length=3,width=2,direction='out')
    #plt.savefig('Diallo_off7_linear1.pdf', format='pdf')
    #plt.savefig('Callen_linear1.pdf', format='pdf')
    plt.show()
    #file_name = 'DIIID_Diallo_freq'


if 1 == 0:
    plt.plot(uni_rhot,omegaDoppler,label='Doppler Shift')
    plt.plot(uni_rhot,mtmFreq,label='Electron Diamagnetic (MTM in plasma frame)')
    plt.plot(uni_rhot,mtmFreq + omegaDoppler,label='Diamagnetic plus Doppler (MTM in lab frame)')
    plt.axhline(y = 100.,\
        color='red',linestyle='dashed',linewidth = 4.5)
    plt.axhline(y = 135.,\
        color='red',linestyle='dashed',linewidth = 4.5)
    plt.axis([0.96,1.,-50.,500.])
    plt.xlabel('rhot')
    plt.ylabel('frequency (kHz)')
    plt.legend(loc = 2, prop = {'size':12})
    plt.grid(True)
    plt.title(iterdb_file_name)
    #plt.tick_params(axis='both',which='major',labelsize=18,\
    #    length=5,width=2,direction='out')
    plt.show()
    #file_name = 'DIIID_Diallo_freq'
    #plt.savefig(file_name+'.pdf', format='pdf')

#from read_iterdb_file import *
from read_iterdb_x import *
import sys
from interp import *
from finite_differences_x import *
from read_EFIT import *
from calc_volume_from_EFIT import *

#rhot_te1, te1, ti1, ne1, ni1, nb1, vrot1 = read_iterdb_file(sys.argv[1])
ITERDBdict = read_iterdb_x(sys.argv[1])
EFITdict = read_EFIT(sys.argv[2])

# compute d Te / d R
uni_R = np.linspace(EFITdict['R'][0], EFITdict['R'][-1], len(EFITdict['R'])*10)
rhot_uniR = interp(EFITdict['R'], EFITdict['rhotn'], uni_R)
te_uniR = interp(ITERDBdict['rhot_te'], ITERDBdict['te'], rhot_uniR)
d_te_d_R = - first_derivative(te_uniR, uni_R)
ne_uniR = interp(ITERDBdict['rhot_ne'], ITERDBdict['ne'], rhot_uniR)

# compute chi
Qheating_MW = 1.
ntheta = 512
area_m2 = surfaceArea(EFITdict, ntheta)
print('surface area: {}'.format(area_m2))
# convert eV into J: multiply by e
charge_e = 1.6E-19
chi_uniR = Qheating_MW * 1.E6 / area_m2 / (d_te_d_R / 2.) / charge_e / ne_uniR

# compute diffusivity from B_tilda / B
p_kg = 1.673E-27
me_kg = 9.11E-31
mD_kg = 2. * p_kg
vth_uniR = np.sqrt(te_uniR * charge_e / me_kg)

rhot_hires, shat_hires, Ls_hires = magneticShear(EFITdict, True)
vth_hires = interp(rhot_uniR, vth_uniR, rhot_hires)
# Bdot [T/s] and f [kHz] is given from experiment
# angular velocity omega = 2 pi f
# total magnitude of the magnetic field is B ~ 5 T
# the decay factor is about 0.25
# B_tilda_n is B_tilda / B in the plasma
B_tilda_n = 3.E-5
D_Btilda = np.pi * vth_hires * Ls_hires * B_tilda_n**2

if 1 == 1:
    plt.plot(rhot_uniR, vth_uniR, label = 'v_th')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.semilogy(rhot_hires, D_Btilda, label = 'D')
    plt.semilogy(rhot_uniR, chi_uniR, label = 'chi')
    plt.axis([0.85, 1., 0.01, 10.])
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

if 1 == 1:
    plt.plot(rhot_hires, D_Btilda, label = 'D')
    plt.plot(rhot_uniR, chi_uniR, label = 'chi')
    plt.axis([0.85, 1., 0.01, 10.])
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

if 1 == 0:
    #plt.plot(rhot_uniR, chi_uniR, '.', label = 'chi')
    plt.semilogy(rhot_uniR, chi_uniR, '.', label = 'chi')
    plt.axis([0.85, 1., 0.01, 1.])
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(rhot_uniR, d_te_d_R, '.', label = 'd Te / d R')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(rhot_uniR, ne_uniR, '.', label = 'ne')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 0:
    #plt.plot(rhot_te1, te1, '.', label = 'old')
    plt.plot(rhot_uniR, te_uniR, '.', label = 'uniR')
    plt.plot(ITERDBdict['rhot_te'], ITERDBdict['te'], label = 'new')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(rhot_te1, ti1, '.', label = 'old')
    plt.plot(ITERDBdict['rhot_ti'], ITERDBdict['ti'], label = 'new')
    plt.legend()
    plt.show()
    plt.plot(rhot_te1, ne1, '.', label = 'old')
    plt.plot(ITERDBdict['rhot_ne'], ITERDBdict['ne'], label = 'new')
    plt.legend()
    plt.show()
    plt.plot(rhot_te1, ni1, '.', label = 'old')
    plt.plot(ITERDBdict['rhot_ni'], ITERDBdict['ni'], label = 'new')
    plt.legend()
    plt.show()
    plt.plot(rhot_te1, vrot1, '.', label = 'old')
    plt.plot(ITERDBdict['rhot_vrot'], ITERDBdict['vrot'], label = 'new')
    plt.legend()
    plt.show()

#from read_iterdb_file import *
from read_iterdb_x import *
import sys
from interp import *
from finite_differences_x import *
from read_EFIT import *
from calc_volume_from_EFIT import *
import matplotlib.pyplot as plt

#rhot_te1, te1, ti1, ne1, ni1, nb1, vrot1 = read_iterdb_file(sys.argv[1])
ITERDBdict = read_iterdb_x(sys.argv[1])
EFITdict = read_EFIT(sys.argv[2])

# compute d ne / d R
uni_R = np.linspace(EFITdict['R'][0], EFITdict['R'][-1], len(EFITdict['R'])*10)
rhot_uniR = interp(EFITdict['R'], EFITdict['rhotn'], uni_R)
ne_uniR = interp(ITERDBdict['rhot_ne'], ITERDBdict['ne'], rhot_uniR)
d_ne_d_R = abs(first_derivative(ne_uniR, uni_R))
te_uniR = interp(ITERDBdict['rhot_te'], ITERDBdict['te'], rhot_uniR)
d_te_d_R = - first_derivative(te_uniR, uni_R)

# compute chi
Gamma_LCFS = 1.E20
Qheating_MW = 4.
ntheta = 512
area_m2 = surfaceArea(EFITdict, ntheta)
print('surface area: {}'.format(area_m2))
D_uniR = Gamma_LCFS / area_m2 / (d_ne_d_R / 2.)
# convert eV into J: multiply by e
charge_e = 1.6E-19
chi_uniR = Qheating_MW * 1.E6 / area_m2 / (d_te_d_R / 2.) / charge_e / ne_uniR

print('LCFS: ', uni_R[-1], ne_uniR[-1])
x1 = uni_R[-1] - 0.007
ix1 = np.argmin(abs(uni_R - x1))
print('top', uni_R[ix1], ne_uniR[ix1])
print('avg gradient = ', (ne_uniR[ix1] - ne_uniR[-1]) / (uni_R[-1] - uni_R[ix1]))

x2 = uni_R[-1] - 0.003
ix2 = np.argmin(abs(uni_R - x2))
print('LCFS: ', uni_R[-1], te_uniR[-1])
print('top', uni_R[ix2], te_uniR[ix2])
print('avg gradient = ', (te_uniR[ix2]*1E-3 - te_uniR[-1]*1E-3) / (uni_R[-1] - uni_R[ix2]))

if 1 == 1:
#    plt.plot(uni_R, ne_uniR, '.', label = 'ne')
    #plt.scatter([uni_R[-1], uni_R[ix1]], [ne_uniR[-1], ne_uniR[ix1]], color = 'red')
    plt.plot(uni_R, te_uniR * 1E-3, '.', label = 'te (keV)')
    plt.scatter([uni_R[-1], uni_R[ix2]], [te_uniR[-1]*1E-3, te_uniR[ix2]*1E-3], color = 'red')
    #plt.xlabel('rhot')
#    plt.xlabel('R (m)')
#    plt.ylabel('ne (m^-3)')
    plt.axis([0.84, 0.89, 0., 1.5])
#    plt.axis([0.84, 0.89, 0.5E20, 1.8E20])
    plt.tick_params(axis='both', which='major',labelsize=18,\
    length=5,width=2,direction='out')
#    plt.legend()
#    plt.show()
    plt.savefig('Imode_te02.png')
if 1 == 0:
    plt.plot(rhot_uniR, d_ne_d_R, '.', label = 'd ne / d R')
    #plt.plot(rhot_uniR, d_te_d_R * 1E-3, '.', label = 'd te / d R')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(rhot_uniR, D_uniR, label = 'D')
    plt.axis([0.85, 1., 0., .05])
    plt.ylabel('m^2/s')
    plt.xlabel('rhot')
    plt.title('I mode')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.plot(rhot_uniR, D_uniR, label = 'D')
    plt.plot(rhot_uniR, chi_uniR, label = 'chi')
    plt.axis([0.85, 1., 0., 5.])
    plt.ylabel('m^2/s')
    plt.xlabel('rhot')
    plt.title('I mode')
    plt.legend()
    plt.show()
if 1 == 0:
    plt.semilogy(rhot_uniR, D_uniR, label = 'D')
    plt.semilogy(rhot_uniR, chi_uniR, label = 'chi')
    plt.axis([0.85, 1., 0., 5.])
    plt.ylabel('m^2/s')
    plt.xlabel('rhot')
    plt.title('I mode')
    plt.legend()
    plt.show()



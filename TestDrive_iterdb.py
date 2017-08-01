#from read_iterdb_file import *
from read_iterdb_x import *
import sys
from interp import *
from finite_differences_x import *
from read_EFIT import *

#rhot_te1, te1, ti1, ne1, ni1, nb1, vrot1 = read_iterdb_file(sys.argv[1])
ITERDBdict = read_iterdb_x(sys.argv[1])
EFITdict = read_EFIT(sys.argv[2])

uni_R = np.linspace(EFITdict['R'][0], EFITdict['R'][-1], len(EFITdict['R'])*10)
rhot_uniR = interp(EFITdict['R'], EFITdict['rhotn'], uni_R)
te_uniR = interp(ITERDBdict['rhot_te'], ITERDBdict['te'], rhot_uniR)
d_te_d_R = - first_derivative(te_uniR, uni_R)
ne_uniR = interp(ITERDBdict['rhot_ne'], ITERDBdict['ne'], rhot_uniR)

if 1 == 1:
    plt.plot(rhot_uniR, d_te_d_R, '.', label = 'd Te / d R')
    plt.legend()
    plt.show()
if 1 == 1:
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

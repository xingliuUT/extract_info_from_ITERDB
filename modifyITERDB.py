#from read_iterdb_file import *
from read_iterdb_x import *
import sys
from interp import *
from finite_differences_x import *
from read_EFIT import *
from calc_volume_from_EFIT import *
from write_iterdb import *
import numpy as np

ITERDBdict = read_iterdb_x(sys.argv[1])
EFITdict = read_EFIT(sys.argv[2])
vrot0 = ITERDBdict['vrot']
rhot_vrot0 = ITERDBdict['rhot_vrot']
shift = 0.005
vrot_shifted = interp(rhot_vrot0 - shift, vrot0, rhot_vrot0)
plt.plot(rhot_vrot0, vrot0, label = 'before')
plt.plot(ITERDBdict['rhot_te'], vrot_shifted, label = 'new')
plt.legend()
plt.show()
file_out_base = 'efit_DIIID153764_' + 'off7Zp2bsx11_120_exp_vrot_left'
base_number = '153764'
time_str = '0.0'
psipn = interp(EFITdict['rhotn'], EFITdict['psipn'], ITERDBdict['rhot_te'])
rhop = np.sqrt(abs(psipn))
output_iterdb(ITERDBdict['rhot_te'], \
              rhop, \
              ITERDBdict['ne']*1.E-19, \
              ITERDBdict['te']*1.E-3, \
              ITERDBdict['ni']*1.E-19, \
              ITERDBdict['ti']*1.E-3, \
              file_out_base+base_number, \
              base_number, \
              time_str, \
              vrot=vrot_shifted, \
              nimp=ITERDBdict['nz']*1.E-19)

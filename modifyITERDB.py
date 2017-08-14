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
file_out_base = 'new' + 'DIIID' 
base_number = '98889'
time_str = '4530'
psipn = interp(EFITdict['rhotn'], EFITdict['psipn'], ITERDBdict['rhot_te'])
plt.plot(EFITdict['rhotn'], EFITdict['psipn'], label = 'old')
plt.plot(ITERDBdict['rhot_te'], psipn, label = 'new')
plt.legend()
plt.show()
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
              vrot=ITERDBdict['vrot'], \
              nimp=ITERDBdict['nz']*1.E-19)

#from read_iterdb_file import *
from read_iterdb_x import *
import sys
from interp import *
from finite_differences_x import *
from read_EFIT import *
from calc_volume_from_EFIT import *

ITERDBdict = read_iterdb_x(sys.argv[1])
EFITdict = read_EFIT(sys.argv[2])

ntheta = 128
V, N = totalN(EFITdict, ITERDBdict, sys.argv[3], ntheta)
print('# Total volume = {}'.format(V))
print('# Total particle = {}'.format(N))

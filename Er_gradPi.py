from read_iterdb_x import *
from read_EFIT import *
from finite_differences_x import *
from interp import *
import matplotlib.pyplot as plt
import numpy as np
import sys

iterdbFileName = sys.argv[1]
efitFileName = sys.argv[2]
EFITdict = read_EFIT(efitFileName)
ITERDBdict = read_iterdb_x(iterdbFileName)

Rgrid = EFITdict['R']

uniR = np.linspace(Rgrid[0], Rgrid[-1], len(Rgrid))
rhot_uniR = interp(Rgrid, EFITdict['rhotn'], uniR)

# ti in ITERDB file is in the unit of eV
ti_uniR = interp(ITERDBdict['rhot_ti'], ITERDBdict['ti'], rhot_uniR)
ni_uniR = interp(ITERDBdict['rhot_ni'], ITERDBdict['ni'], rhot_uniR)
pi_uniR = ti_uniR * ni_uniR

gradPioverNie = first_derivative(pi_uniR, uniR) / ni_uniR



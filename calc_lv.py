#! /usr/bin/python

from pylab import *
from sys import argv,exit,stdout
import optparse as op
import matplotlib.pyplot as plt
import numpy as np
from interp import *
from read_iterdb_file import *
import sys

sys.stdout = open('lvlw.txt', 'w')

parser = op.OptionParser()
options,args = parser.parse_args()
c_buffer_x = float(args[0])
buffer_size = float(args[1])
physW = float(args[2])
file1 = args[3]

l_buffer_x = float(c_buffer_x) - physW*float(buffer_size)
r_buffer_x = float(c_buffer_x) + physW*float(buffer_size)

#pdata=np.genfromtxt('p_info.dat')
#pdata=np.genfromtxt(prof_file_name)
#rhot=pdata[:,0]
#te=pdata[:,2]

rhot, te, ti1, ne1, ni1, nb1, vrot1 = read_iterdb_file(file1)
e = 1.6*10**(-19)

rhot_fine = linspace(rhot[0],rhot[-1],10*len(rhot))
te_fine = interp(rhot,te,rhot_fine)

l_ind = np.argmin(abs(rhot_fine - float(l_buffer_x)))
c_ind = np.argmin(abs(rhot_fine - float(c_buffer_x)))
r_ind = np.argmin(abs(rhot_fine - float(r_buffer_x)))

te_l = te_fine[l_ind]
te_c = te_fine[c_ind]
te_r = te_fine[r_ind]

lv = 3*np.sqrt(te_l/te_c)
lw = lv**2

print 'lv = ', lv
print 'te_l = ', te_l*0.001/e
print 'te_c = ', te_c*0.001/e
print 'te_r = ', te_r*0.001/e
print 'lw = ', lw
print 'buffer', l_buffer_x, r_buffer_x

nv = 48*np.sqrt(te_l/te_r)
nw = 16*te_l/te_r

print 'nv = ', nv
print 'nw = ', nw


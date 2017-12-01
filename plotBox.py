from read_iterdb_x import *
import sys
from interp import *
from finite_differences_x import *
from read_EFIT import *
from calc_volume_from_EFIT import *
from write_iterdb import *
import numpy as np

ITERDBdict = read_iterdb_x('./CMOD1120907032May12.iterdb')
EFITdict = read_EFIT('../extract_info_from_EFIT/g1120907032.01012')

#ITERDBdict = read_iterdb_x(sys.argv[1])
#EFITdict = read_EFIT(sys.argv[2])

if 1 == 1:

    plt.figure(figsize = (9, 3.5))

    mid_ped = 0.97
    mid_ped_ind = np.argmin(abs(ITERDBdict['rhot_te'] - mid_ped))

    plt.plot(ITERDBdict['rhot_te'], ITERDBdict['te'] * 1E-3, linewidth = 4.5, color = 'blue')
    plt.plot(ITERDBdict['rhot_ti'], ITERDBdict['ti'] * 1E-3, linewidth = 4.5, color = 'red')
    plt.plot(ITERDBdict['rhot_ne'], ITERDBdict['ne'] * 1E-20, linewidth = 4.5, color = 'green')
    plt.plot(ITERDBdict['rhot_ni'], ITERDBdict['ni'] * 1E-20, linewidth = 4.5, color = 'magenta')
    #plt.plot(ITERDBdict['rhot_vrot'], ITERDBdict['vrot'] * 1E-5, linewidth = 4.5, color = 'black')

#    plt.plot(ITERDBdict['rhot_te'], ITERDBdict['te'] / ITERDBdict['te'][mid_ped_ind], linewidth = 4.5, color = 'blue')
#    plt.plot(ITERDBdict['rhot_ti'], ITERDBdict['ti'] / ITERDBdict['ti'][mid_ped_ind], linewidth = 4.5, color = 'green')
#    plt.plot(ITERDBdict['rhot_ne'], ITERDBdict['ne'] / ITERDBdict['ne'][mid_ped_ind], linewidth = 4.5, color = 'red')

    x_ls = 0.951; x_ld = 0.9605; x_rd = 0.9795; x_rs = 0.989
#    x_ls = 0.9455; x_ld = 0.9514; x_rs = 0.9945; x_rd = 0.9886
#    x_ls = 0.9464; x_ld = 0.96; x_rs = 0.9986; x_rd = 0.985
    plt.axvline(x = x_ls, linewidth = 2., color = 'black')
    plt.axvline(x = x_ld, linewidth = 2., color = 'black', linestyle = '--')
    plt.axvline(x = x_rd, linewidth = 2., color = 'black', linestyle = '--')
    plt.axvline(x = x_rs, linewidth = 2., color = 'black')


    plt.axis([x_ls, x_rs, 0., 1.4])
#    plt.axis([0.94, 1., 0., 0.7])
#    plt.grid()
#    plt.xlabel('x')

    plt.tick_params(axis='both',which='both',labelsize=20,\
        length=5,width=2,direction='out')
    #plt.show()
    #plt.savefig('Diallo_global_box.pdf', format = 'pdf')
    plt.savefig('Imode_global_box_2.pdf', format = 'pdf')

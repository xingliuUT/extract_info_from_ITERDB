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
#GENE_shat_data = 'shat_summary_Callen_old_EFIT.txt'
#EFIT_shat_data = 'shat_Callen_old.txt'

GENE_shat_data = 'shat_summary_Diallo_EFIT.txt'
EFIT_shat_data = 'Diallo_EFIT_shat.txt'

#GENE_shat_data = 'shat_summary_Diallo_vmec04.txt'
#EFIT_shat_data = 'shat_VMEC04.txt'

shat_gene = np.genfromtxt(GENE_shat_data)
shat_efit = np.genfromtxt(EFIT_shat_data)

if 1 == 1:
    fig, ax = plt.subplots()
    ax.plot(ITERDBdict['rhot_ne'], ITERDBdict['ne']*1E-20*10, linewidth = 2., label = 'ne')
    ax.plot(ITERDBdict['rhot_te'], ITERDBdict['te']*1E-3*10, linewidth = 2., label = 'te')
    ax.plot(shat_gene[:,0], shat_gene[:,1], marker = 'o', markersize = 10, linewidth = 4., label = 'gene shat')
    ax.plot(shat_efit[:,0], shat_efit[:,2], marker = 'o', markersize = 10, linewidth = 4., label = 'efit shat')
    #ax.plot(shat_efit[:,0], shat_efit[:,1], marker = 'X', markersize = 4, linewidth = 2., label = 'efit shat')
    #ax.scatter(shat_efit[:,0], shat_efit[:,2], marker = 'o', s = 100, label = 'efit shat', color = 'red')
    ax.axhline(y = 0, linewidth = 2.5, color = 'black')
    ax.axis([0.94,1.,-2., 9.])
    #ax.xticks(np.arange(0.94, 1., 0.1))
    #plt.legend()

    ax.tick_params(axis='both',which='major',labelsize=22,\
        length=5,width=2,direction='out')
    plt.locator_params(nbins = 7)
    #plt.show()
    #plt.savefig('Diallo_vmec04_shat_comparison.pdf', format = 'pdf')
    plt.savefig('Diallo_shat_comparison.pdf', format = 'pdf')
    #plt.savefig('Callen_shat_comparison.pdf', format = 'pdf')

if 1 == 0:
    plt.plot(EFITdict['rhotn'], EFITdict['Pres'])
    plt.show()
if 1 == 0:
    np.savetxt('vrot.txt', np.column_stack((ITERDBdict['rhot_vrot'], ITERDBdict['vrot'])))
if 1 == 0:
    psipn = interp(EFITdict['rhotn'], EFITdict['psipn'], ITERDBdict['rhot_te'])
    plt.plot(psipn, ITERDBdict['te'])
    plt.xlabel('psi_pol_n')
    plt.ylabel('te')
    plt.axis([0.99, 1., 0., 300.])
    plt.grid()
    plt.show()
    psipn = interp(EFITdict['rhotn'], EFITdict['psipn'], ITERDBdict['rhot_ne'])
    plt.plot(psipn, ITERDBdict['ne'])
    plt.xlabel('psi_pol_n')
    plt.ylabel('ne')
    plt.axis([0.99, 1., 0., 3.E19])
    plt.grid()
    plt.show()
if 1 == 0:
    file_out_base = 'efit_DIIID153764_' + 'vmec04_exp_vrot' 
    base_number = '153764'
    time_str = '0.0'
    psipn = interp(EFITdict['rhotn'], EFITdict['psipn'], ITERDBdict['rhot_te'])
    plt.plot(EFITdict['rhotn'], EFITdict['psipn'], label = 'old')
    plt.plot(ITERDBdict['rhot_te'], psipn, label = 'new')
    plt.legend()
    plt.show()
    rhop = np.sqrt(abs(psipn))
    if 1 == 1:
        exp_data = np.genfromtxt('exp_vrot.txt')
        exp_vrot = interp(exp_data[:,0], exp_data[:,1], ITERDBdict['rhot_vrot'])
    output_iterdb(ITERDBdict['rhot_te'], \
              rhop, \
              ITERDBdict['ne']*1.E-19, \
              ITERDBdict['te']*1.E-3, \
              ITERDBdict['ni']*1.E-19, \
              ITERDBdict['ti']*1.E-3, \
              file_out_base, \
              base_number, \
              time_str, \
              vrot=exp_vrot, \
              nimp=ITERDBdict['nz']*1.E-19)

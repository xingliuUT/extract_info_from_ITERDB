from read_pfile import *
from read_EFIT import *
from finite_differences_x import *
import matplotlib.pyplot as plt
from write_iterdb import *
import sys

def comparePlot(x1, y1, x2, y2, label1, label2, xl, yl, lc = 1):
    plt.plot(x1, y1, label = label1)
    plt.plot(x2, y2, label = label2)
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.legend(loc = lc)
    plt.show()
def singlePlot(x, y, labl, xl, yl, lc = 1):
    plt.plot(x, y, label = labl)
    plt.xlabel(xl)
    plt.ylabel(yl)
    plt.legend(loc = int(lc))
    plt.show()

#fs = 0.975
#aGENE_m = 0.76
#k_theta_cm = 0.19
#Bref_Gauss = 19545.
mref = 2.
#k2Factor = 3.3
e_charge = 1.6E-19
M_kg = 3.3*10**(-27)
impurityCharge = 6.
efit_file_name = sys.argv[1]
p_file_name = sys.argv[2]

rbsProfs = np.genfromtxt('rbsProfs', skip_header = 3, skip_footer = 0)
rhot0 = rbsProfs[:, 0]
psi0 = rbsProfs[:, 1]
Pres0 = rbsProfs[:, 2]
ne0 = rbsProfs[:,3]
ti0 = rbsProfs[:, 4]
te0 = rbsProfs[:, 5]
Bpol0 = rbsProfs[:, 25]
Btor0 = rbsProfs[:, 26]

if 1 == 0:
    # check input profiles are consistent with output profiles in rbsprofs
    input_profs = np.genfromtxt(sys.argv[3])
    plt.plot(input_profs[:, 0], input_profs[:, 2], label = 'ne input')
    plt.plot(psi0, ne0, label = 'ne rbsProfs')
    plt.legend()
    plt.show()

zeff_prof = np.genfromtxt('Zeff_out')

EFITdict = read_EFIT(efit_file_name)
print(list(EFITdict.keys()))
Zeff_out = interp(zeff_prof[:, 0], zeff_prof[:, 1], psi0)
ni0 = ne0 * (impurityCharge - Zeff_out) / (impurityCharge - 1.)
nz0 = ne0 * (Zeff_out - 1.) / (impurityCharge - 1.) / impurityCharge

psipne,ne,psipte,te,psipni,ni,psipti,ti,psiper,er = read_pfile_raw(p_file_name)
er0 = interp(psiper, er, psi0)
if 1 == 1:
#    singlePlot(EFITdict['rhotn'], EFITdict['Bpol'], '', 'rhot', 'Bpol')
#    singlePlot(EFITdict['rhotn'], EFITdict['Btor'], '', 'rhot', 'Btor')
    # rhotn in EFIT doesn't start exactly from psip = 0
    comparePlot(EFITdict['psipn'], EFITdict['rhotn'], psi0, \
                rhot0, 'rhot EFIT', 'rhot rbsprofs', 'psip', 'rhot')
    # slight difference between rbsProfs Bpol and my calculation from EFIT
    comparePlot(EFITdict['rhotn'], EFITdict['Bpol'], rhot0, \
                Bpol0, 'Bpol EFIT', 'Bpol rbsprofs', 'rhot', 'Bpol')
#    comparePlot(EFITdict['rhotn'], EFITdict['Btor'], rhot0, \
#                -Btor0, 'Btor EFIT', 'Btor rbsprofs', 'rhot', 'Btor')
#    comparePlot(psipne, ne, psi0, \
#                ne0, 'ne p-file', 'ne rbsProf', 'rhot', 'ne')
#    comparePlot(psipni, ni, psi0, \
#                ni0, 'ni p-file', 'ni rbsProf', 'rhot', 'ni')


sepInd = np.argmin(abs(EFITdict['psipn'] - 1.))
# print('index at psipn = 1 is ', sepInd)
Rsep = EFITdict['R'][sepInd]
# print('major R at psipn = 1 is ', Rsep)
# print('major R at index = 1 is ', EFITdict['R'][0])
# print('R grid length is ', len(EFITdict['R']))

uni_R = np.linspace(EFITdict['R'][0],Rsep,EFITdict['nw']*10)
psip_uniR = interp(EFITdict['R'], EFITdict['psipn'], uni_R)
rhot_uniR = interp(EFITdict['R'], EFITdict['rhotn'], uni_R)

if 1 == 0:
    comparePlot(EFITdict['rhotn'], EFITdict['R'], rhot_uniR, \
                uni_R, 'R raw', 'uniform R grid', 'rhot', 'R(m)')

# rhot0 = interp(EFITdict['psipn'], EFITdict['rhotn'], psi0)
if 1 == 0:
    comparePlot(EFITdict['psipn'], EFITdict['rhotn'], psi0, \
                rhot0, 'rhot EFIT', 'rhot p-file', 'psip', 'rhot')
pi0 = ni0 * ti0
pi_uniR = interp(rhot0,pi0,rhot_uniR)
ni_uniR = interp(rhot0,ni0,rhot_uniR)
ti_uniR = interp(rhot0,ti0,rhot_uniR)
pe0 = ne0 * te0
pe_uniR = interp(rhot0,pe0,rhot_uniR)
ne_uniR = interp(rhot0,ne0,rhot_uniR)
te_uniR = interp(rhot0,te0,rhot_uniR)

if 1 == 0:
    comparePlot(rhot0, pi0, rhot_uniR, pi_uniR, 'p_i from p-file', \
                'p_i on uniform R_grid', 'rhot', 'N/m^2')
    comparePlot(rhot0, ni0, rhot_uniR, ni_uniR, 'n_i from p-file', \
                'n_i on uniform R_grid', 'rhot', '10^20 m^-3')
    comparePlot(rhot0, ti0, rhot_uniR, ti_uniR, 't_i from p-file', \
                't_i on uniform R_grid', 'rhot', 'keV')

gradPioverNe = first_derivative(pi_uniR,uni_R)/ni_uniR 
gradPeoverNe = first_derivative(pe_uniR,uni_R)/ne_uniR 
gradPoverNe2 = first_derivative(ni_uniR,uni_R)*ti_uniR/ni_uniR +\
               first_derivative(ti_uniR,uni_R)

#fsInd = np.argmin(abs(rhot0-fs))
#print('tau = zeff * Te / Ti =', zeff[fsInd]*te0[fsInd]/ti0[fsInd])

if 1 == 0:
    plt.plot(rhot0,nz0*6./ne0,label='nz * 6 / ne')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(rhot0,ne0,label='ne (10^20 m^-3)')
    plt.plot(rhot0,te0,label='te (KeV)')
    plt.plot(rhot0,ni0,label='ni (10^20 m^-3)')
    plt.plot(rhot0,nz0*6.,label='nz * 6')
    plt.plot(rhot0,ti0,label='ti (KeV)')
    #plt.axvline(x = fs,color='black',label='mid te ped')
    plt.axis([0.85,1.,0.,3.])
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

#kyGENE = 102*np.sqrt(mref)*np.sqrt(te0*1E3)/Bref_Gauss*k_theta_cm*k2Factor
if 1 == 0:
    singlePlot(rhot0, kyGENE, '', 'rhot', 'ky*rho', 3)

uni_rhot = np.linspace(min(rhot0),max(rhot0),len(rhot0))
ti_u = interp(rhot0,ti0,uni_rhot)
te_u = interp(rhot0,te0,uni_rhot)
ne_u = interp(rhot0,ne0,uni_rhot)
ni_u = interp(rhot0,ni0,uni_rhot)
nz_u = interp(rhot0,nz0,uni_rhot)
p_u = (ni_u + nz_u) * ti_u + ne_u * te_u
if 1 == 1:
    comparePlot(uni_rhot, p_u * 1.E03 * 1.6E-19 * 1.E20, EFITdict['rhotn'], EFITdict['Pres'], 'Pressure ITERDB', \
                'Pressure EFIT', 'rhot', '', 1)
if 1 == 0:
    comparePlot(uni_rhot, ne_u, rhot0, ne0, 'ne ITERDB', \
                'ne pfile', 'rhot', '', 1)

tprime_i = -first_derivative(ti_u,uni_rhot)/ti_u
tprime_e = -first_derivative(te_u,uni_rhot)/te_u
nprime_e = -first_derivative(ne_u,uni_rhot)/ne_u
nprime_i = -first_derivative(ni_u,uni_rhot)/ni_u
nprime_z = -first_derivative(nz_u,uni_rhot)/nz_u
eta_i = tprime_i / nprime_i
eta_e = tprime_e / nprime_e
eta_z = tprime_i / nprime_z

if 1 == 1:
    plt.plot(uni_rhot,nprime_e,label='nprime_e')
    plt.plot(uni_rhot,nprime_i,label='nprime_i')
    plt.plot(uni_rhot,nprime_z,label='nprime_z')
    plt.xlabel('rhot')
    plt.legend()
    plt.show()

#jtor = np.pi * 4.E-7 * (rmag * pprime + ffprime / rmag)
#jtor_u = interp(rhot_n,jtor,uni_rhot)
#jtorprime = -first_derivative(jtor_u,uni_rhot)
if 1 == 0:
    singlePlot(rhot_n, jtor, '', 'rhot', 'jtor')
    singlePlot(uni_rhot, jtorprime, '', 'rhot', 'jtor prime')

#fsInd = np.argmin(abs(uni_rhot-fs))
#print 'At rhot =', fs
#print 'tprime_i =', tprime_i[fsInd]
#print 'tprime_e =', tprime_e[fsInd]
#print 'nprime_e =', nprime_e[fsInd]
#print 'nprime_i =', nprime_i[fsInd]

#kyGENE_u = interp(rhot0,kyGENE,uni_rhot)
#omMTM = kyGENE_u*(tprime_e+nprime_e)
#omMTM = kyGENE_u*nprime_e
#gyroFreq = v_th / a
#gyroFreq = 9.79E3/np.sqrt(mref)*np.sqrt(te_u*1E3)/aGENE_m
#convert from radians/s to kHz
#mtmFreq = omMTM*gyroFreq/2./np.pi/1000.

if 1 == 0:
    singlePlot(uni_rhot, mtmFreq, 'mtm frequency', 'rhot', 'kHz', 2)

if 1 == 0:
    comparePlot(uni_rhot, tprime_i, uni_rhot, tprime_e, 'tprime_i', \
                'tprime_e', 'rhot', '', 2)
    comparePlot(uni_rhot, nprime_e, uni_rhot, nprime_i, 'nprime_e', \
                'nprime_i', 'rhot', '', 2)
    singlePlot(rhot0, ti0/te0, 'ti/te', 'rhot', '', 2)
    comparePlot(uni_rhot, eta_e, uni_rhot, eta_i, 'eta_e', \
                'eta_i', 'rhot', '', 1)
    
# convert from kV/m to V/m
Er_Vm = interp(rhot0,er0,uni_rhot)*1E3

if 1 == 0:
    comparePlot(rhot0, er0, uni_rhot, Er_Vm*1E-3, 'Er old', 'Er new', \
                'rhot', 'kV/m',3)
if 1 == 0:
    comparePlot(rhot_uniR,gradPioverNe,rhot_uniR,gradPeoverNe,'grad Pe / ne /e',\
                'grad Pi / ni / e', 'rhot', 'kV/m', 3)
    comparePlot(uni_rhot,Er_Vm/1E3, rhot_uniR, gradPioverNe, 'experimental Er',\
                'grad Pi / ni / e', 'rhot', 'kV/m', 3)

# vExB m/s
#vExB = Er_Vm/Bref_Gauss*1E4
# negative Er shift in electron diamagnetic direction
#omegaDoppler = -vExB*k_theta_cm*1E2/2./np.pi/1E3

if 1 == 0:
    plt.plot(uni_rhot,omegaDoppler,label='Doppler Shift')
    plt.plot(uni_rhot,mtmFreq,label='Electron Diamagnetic (MTM in plasma frame)')
    plt.plot(uni_rhot,mtmFreq + omegaDoppler,label='Diamagnetic plus Doppler (MTM in lab frame)')
    #plt.axvline(x = fs,color='black')
    plt.axis([0.92,1.,-50.,300.])
    plt.xlabel('rhot')
    plt.ylabel('frequency (kHz)')
    plt.legend(loc = 2, prop = {'size':12})
    plt.show()
    #file_name = 'DIIID_Diallo_freq'
    #plt.savefig(file_name+'.pdf', format='pdf')

R_u = interp(EFITdict['rhotn'],EFITdict['R'],uni_rhot)
if 1 == 0:
    plt.plot(R_u,omegaDoppler,label='Doppler Shift',linewidth = 2)
    plt.plot(R_u,mtmFreq,label='Electron Diamagnetic (MTM in plasma frame)',linewidth = 2)
    plt.plot(R_u,mtmFreq + omegaDoppler,label='Diamagnetic plus Doppler (MTM in lab frame)',linewidth = 2)
    #plt.axvline(x = fs,color='black')
    plt.axis([2.24,2.3,-50.,300.])
    plt.xlabel('R(m)')
    plt.ylabel('frequency (kHz)')
    plt.legend(loc = 2, prop = {'size':12})
    plt.show()
    #file_name = 'DIIID_Diallo_freq'
    #plt.savefig(file_name+'.pdf', format='pdf')
    np.savetxt('mtmFreq.dat',np.column_stack((uni_rhot,R_u,omegaDoppler,mtmFreq)))


Bpol_u = interp(EFITdict['rhotn'],EFITdict['Bpol'],uni_rhot)
if 1 == 0:
    comparePlot(EFITdict['rhotn'], EFITdict['Bpol'], uni_rhot, Bpol_u, 'Bpol EFIT', 'Bpol uniform rhot', 'rhot', '')


#te_u_J = te_u*e*1E3
#q_u = interp(rhot_n,qpsi,uni_rhot)

# minus sign for consistency
omega_tor = - Er_Vm / (R_u * Bpol_u)
if 1 == 1:
    singlePlot(uni_rhot, omega_tor, 'omega_tor', 'rhot', '')

if 1 == 0:
    np.savetxt('omegaTor.txt',np.column_stack((uni_rhot,omega_tor)))

#cs = np.sqrt(te_u_J/M_kg)

#gammaE = - (uni_rhot / q_u)*first_derivative(omega_tor,uni_rhot) / (cs / aGENE_m)
#print 'ExBrate =', gammaE[fsInd]

if 1 == 0:
    singlePlot(uni_rhot, abs(gammaE), 'abs gammaE', 'rhot', '')

if 1 == 0:
    #rhop=np.sqrt(psi0)
    psi_u = interp(rhot0,psi0,uni_rhot)
    rhop_u = np.sqrt(psi_u)
    f = open('profiles_e','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Te(kev) 4.ne(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((uni_rhot,rhop_u,te_u,ne_u*10.)))
    f.close()

    f = open('profiles_i','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((uni_rhot,rhop_u,ti_u,ni_u*10.)))
    f.close()

    f = open('profiles_z','w')
    f.write('# 1.rho_tor 2.rho_pol 3.Ti(kev) 4.ni(10^19m^-3)\n#\n')
    np.savetxt(f,np.column_stack((uni_rhot,rhop_u,ti_u,nz_u*10.)))
    f.close()

if 1 == 1:
    rhop=np.sqrt(psi0)
    file_out_base = sys.argv[3]
    base_number = '_new'
    time_str = '01380'
    psi_u = interp(rhot0,psi0,uni_rhot)
    rhop_u = np.sqrt(psi_u)
    output_iterdb(uni_rhot,rhop_u,ne_u*10.,te_u,ni_u*10.,ti_u,file_out_base+base_number,base_number,time_str,vrot=omega_tor,nimp=nz_u*10.)


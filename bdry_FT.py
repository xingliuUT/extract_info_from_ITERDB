import numpy as np
import matplotlib.pyplot as plt
import math
import sys

ncoeff = 49

#bndy = np.genfromtxt(sys.argv[1],skip_header=1,dtype = 'float')
bndy = np.loadtxt(sys.argv[1], dtype = 'float')
r_bndy = bndy[:,0]
z_bndy = bndy[:,1]
#plt.plot(r_bndy,z_bndy)
#plt.show()

n=len(r_bndy)
r_m = np.fft.fft(r_bndy)/n
z_m = np.fft.fft(z_bndy)/n
freq = np.fft.fftfreq(n)
#plt.plot(freq, r_m.imag,'r.')
#plt.plot(freq, r_m.real)
#plt.show()
#plt.plot(freq, z_m.imag,'r.')
#plt.plot(freq, z_m.real)
#plt.show()

rbc = np.empty(0)
rbs = np.empty(0)
zbc = np.empty(0)
zbs = np.empty(0)

rbc = np.append(rbc,r_m[0])
rbs = np.append(rbs,0)
zbc = np.append(zbc,z_m[0])
zbs = np.append(zbs,0)

for i in np.arange(1,n/2):
	rbc_i = (r_m[i].real+r_m[-i].real)
        rbs_i = -(r_m[i].imag-r_m[-i].imag)
	zbc_i = (z_m[i].real+z_m[-i].real)
        zbs_i = -(z_m[i].imag-z_m[-i].imag)
	rbc = np.append(rbc,rbc_i)
	rbs = np.append(rbs,rbs_i)
	zbc = np.append(zbc,zbc_i)
	zbs = np.append(zbs,zbs_i)

if 1 == 0:
    plt.plot(rbc, label = 'rbc')
    plt.plot(rbs, label = 'rbs')
    plt.plot(zbc, label = 'zbc')
    plt.plot(zbs, label = 'zbs')
    plt.legend()
    plt.show()

f = open('coeff.dat','w')
f.write('# rbc rbs zbc zbs\n')
np.savetxt(f,np.column_stack((rbc,rbs,zbc,zbs)))
f.close()

#f = open('','w')
def format_D(fl):
    a = '%.12E' % fl
    return a.split('E')[0]+'D'+a.split('E')[1]
if 1 == 1:
    f = open('coeff_' + sys.argv[1],'w')
    for i in range(ncoeff):
        f.write('{0}{1}{2}= {3:.24}    {4}{5}{6}= {7:.24}'.
        format('RBC(0,', i, ')', format_D(rbc[i]), 'ZBS(0,', i, ')', format_D(zbs[i])))
        f.write('\n') 
    for i in range(ncoeff):
        f.write('{0}{1}{2}= {3:.24}    {4}{5}{6}= {7:.24}'.
        format('RBS(0,', i, ')', format_D(rbs[i]), 'ZBC(0,', i, ')', format_D(zbc[i])))
        f.write('\n') 
    f.close()

r = np.empty(0)
z = np.empty(0)
theta = np.empty(0)

for i in range(2000):
    this_theta = i*2*math.pi/2000
    this_r = 0.0
    this_z = 0.0
    for j in range(len(rbc)):
        arg = j*this_theta
        this_r = this_r + rbc[j]*math.cos(arg) + rbs[j]*math.sin(arg)
        this_z = this_z + zbc[j]*math.cos(arg) + zbs[j]*math.sin(arg)
    r = np.append(r, this_r)
    z = np.append(z, this_z)
    theta = np.append(theta, this_theta)

plt.plot(r,z,'.')
plt.plot(r_bndy,z_bndy,'r.')
plt.xlabel('r')
plt.ylabel('z')
plt.axis('equal')
plt.show()

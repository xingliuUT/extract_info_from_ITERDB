from read_pfile import *
import matplotlib.pyplot as plt
import sys

p_file_name = sys.argv[1]
psipne,ne,psipte,te,psipni,ni,psipti,ti,psiper,er = read_pfile_raw(p_file_name)
if 1 == 1:
    plt.plot(psipne, ne, label = 'ne')
    plt.plot(psipni, ni, label = 'ni')
    plt.plot(psipte, te, label = 'te')
    plt.plot(psipti, ti, label = 'ti')
    plt.xlabel('psip')
    plt.legend()
    plt.show()
if 1 == 1:
    plt.plot(psiper, er, label = 'er')
    plt.xlabel('psip')
    plt.legend()
    plt.show()

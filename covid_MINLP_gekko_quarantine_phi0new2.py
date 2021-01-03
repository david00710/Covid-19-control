import numpy as np
import matplotlib.pyplot as plt

#import pyscipopt as scip

import numpy as np

N=8.9

sigma=0.001

gamma2=0.05

# gamma1=0.063
# gamma2=0.05

beta=0.2967

Pmax=0.1
Pmax0=8.9
Pmin=0
Pmaxper=0.001
gamma1 = 0.063

Ts=1

n = 200

S = [0] * n
U = [0] * n
Q = [0] * n
C = [0] * n
Cper = [0] * n

S[0]= 8.9
U[0]= 0.001
Q[0]= 0
C[0]= 0

s=0

for j in range(n-1):
    S[j+1]= S[j]- Ts*beta*S[j]*U[j]/N
    U[j+1]= U[j]+ Ts*beta*S[j]*U[j]/N - Ts*gamma1*U[j]
    Q[j+1]= Q[j] + Ts*gamma1*U[j] - Ts*(gamma2+(1-gamma2)*sigma)*Q[j]
    C[j+1]= C[j]+ Ts*(gamma2+(1-gamma2)*sigma)*Q[j]
    Cper[j+1]=C[j + 1] - C[j]


n=200


time = range(n)
time2 = range(n-1)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

#plt.plot(time,S,'r-',label='susceptible')
plt.plot(time,U,'m-',label='un-quarantined infected',linewidth=3.0)
plt.plot(time,Q,'b-',label='quarantined infected',linewidth=3.0)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('number of individuals (M)', fontsize=18)
plt.xlim(0, n)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig('case_quarantine_phi0a.png')
plt.show()

plt.plot(time,C,'k-',label='confirmed infected',linewidth=3.0)
plt.xlim(0, n)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('confirmed infected (M)', fontsize=18)
plt.tight_layout()
plt.savefig('case_quarantine_phi0c.png')
plt.show()

plt.plot(time,Cper,'k-',label='confirmed infected per day',linewidth=3.0)
plt.xlim(0, n)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('confirmed infected (M) per day', fontsize=18)
plt.tight_layout()
plt.savefig('case_quarantine_phi0d.png')
plt.show()

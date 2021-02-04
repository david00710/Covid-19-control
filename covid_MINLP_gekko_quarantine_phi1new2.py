import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
import pyomo.environ as pyo
#import pyscipopt as scip

import numpy as np
from gekko import GEKKO

# SUQC


m = GEKKO()


N=8.9

sigma=0


gamma2=0.05

# gamma1=0.063
# gamma2=0.05

beta=0.2967

Pmax=0.1
Pmax0=8.9
Pmin=0
Pmaxper=0.001

Ts=1

n = 200

S = m.Array(m.Var,(n))
U = m.Array(m.Var,(n))
Q = m.Array(m.Var,(n))
C = m.Array(m.Var,(n))
Cper = m.Array(m.Var,(n-1))
u = m.Array(m.Var,(n-1))

m.Equation(S[0] == 8.9)
m.Equation(U[0] == 0.001)
m.Equation(Q[0] == 0)
m.Equation(C[0] == 0)

s=0

for j in range(n-1):
    m.Equation(S[j+1] == S[j]- Ts*beta*S[j]*U[j]/N)
    m.Equation(U[j+1] == U[j]+ Ts*beta*S[j]*U[j]/N - Ts*u[j]*U[j])
    m.Equation(Q[j+1] == Q[j]+ Ts*u[j]*U[j] - Ts*(gamma2+(1-gamma2)*sigma)*Q[j])
    m.Equation(C[j+1] == C[j]+ Ts*(gamma2+(1-gamma2)*sigma)*Q[j])
    m.Equation(Cper[j] == C[j + 1] - C[j])
    m.Equation(0 <= u[j])
    m.Equation(u[j] <= 1)
    m.Equation(Pmin <= S[j])
    m.Equation(S[j] <= Pmax0)
    m.Equation(Pmin <= U[j])
    m.Equation(U[j] <= Pmax)
    m.Equation(Pmin <= Q[j])
    m.Equation(Q[j] <= Pmax)
    m.Equation(Pmin <= C[j])
    m.Equation(C[j] <= Pmax)
    m.Equation(Cper[j] <= Pmaxper)
    s = s+u[j]*u[j]

j=n-1
m.Equation(Pmin <= S[j])
m.Equation(S[j] <= Pmax0)
m.Equation(Pmin <= U[j])
m.Equation(U[j] <= Pmax)
m.Equation(Pmin <= Q[j])
m.Equation(Q[j] <= Pmax)
m.Equation(Pmin <= C[j])
m.Equation(C[j] <= Pmax)


m.Obj(s)
m.solve()


Svalue=[]
for j in range(n):
    Svalue.append(S[j].value)

Uvalue=[]
for j in range(n):
    Uvalue.append(U[j].value)

Qvalue = []
for j in range(n):
    Qvalue.append(Q[j].value)

Cvalue = []
for j in range(n):
    Cvalue.append(C[j].value)

Cpervalue=[]
for j in range(n-1):
    Cpervalue.append(Cper[j].value)

uvalue=[]
for j in range(n-1):
    uvalue.append(u[j].value)
print(uvalue)

time = range(n)
time2 = range(n-1)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

plt.plot(time,Uvalue,'m-',label='un-quarantined infected',linewidth=3.0)
plt.plot(time,Qvalue,'b-',label='quarantined infected',linewidth=3.0)
plt.ylim(-0.01, 0.5)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('number of individuals (M)', fontsize=18)
plt.xlim(0, n)
plt.legend(fontsize=16)
plt.tight_layout()
plt.savefig('case_quarantine_phi1a.png')
plt.show()

plt.plot(time2,uvalue,'k-',label='quarantine',linewidth=3.0)
plt.xlim(0, n)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('quarantine rate', fontsize=18)
plt.tight_layout()
plt.savefig('case_quarantine_phi1b.png')
plt.show()

plt.plot(time,Cvalue,'k-',label='confirmed infected',linewidth=3.0)
plt.xlim(0, n)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('confirmed infected (M)', fontsize=18)
plt.tight_layout()
plt.axhline(0.1, color='r', linestyle=':')
plt.savefig('case_quarantine_phi1c.png')
plt.show()

plt.plot(time2,Cpervalue,'k-',label='confirmed infected per day',linewidth=3.0)
plt.xlim(0, n)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('confirmed infected (M) per day', fontsize=18)
plt.tight_layout()
plt.axhline(0.001, color='r', linestyle=':')
plt.savefig('case_quarantine_phi1d.png')
plt.show()

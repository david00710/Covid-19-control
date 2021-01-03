import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
import pyomo.environ as pyo
#import pyscipopt as scip

import numpy as np
from gekko import GEKKO

# shield immunity

m = GEKKO()

N=10
beta=0.75
epsilon=1/5
mu=1/30295
gamma=1/5
alpha = 0.006

Ts=1

DeltaDmax = 0.002
Dmax = 0.07

lambda1 = mu*N

n = 100
I = m.Array(m.Var,(n))
E = m.Array(m.Var,(n))
S = m.Array(m.Var,(n))
R = m.Array(m.Var,(n))
D = m.Array(m.Var,(n))
Dper = m.Array(m.Var,(n-1))
u = m.Array(m.Var,(n-1))

m.Equation(I[0] == 0.001)
m.Equation(E[0] == 0.02)
m.Equation(S[0] == 9.979)
m.Equation(R[0] == 0)
m.Equation(D[0] == 0)

s=0

for j in range(n-1):
    m.Equation(I[j+1] == I[j]+ Ts*epsilon*E[j] - Ts*(gamma+mu+alpha)*I[j])
    m.Equation(E[j+1] == E[j]+ Ts*beta*S[j]*I[j]/(N + u[j]*R[j]) - Ts*(mu+epsilon)*E[j])
    m.Equation(S[j+1] == S[j]+ Ts*mu*(I[j] + E[j] + S[j] + R[j]) - Ts*mu*S[j] - Ts*beta*S[j]*I[j]/(N + u[j]*R[j]) )
    m.Equation(R[j+1] == R[j]+ Ts*gamma*I[j] - Ts*mu*R[j] )
    m.Equation(D[j] == 10 - I[j] - E[j] - S[j] - R[j])
    m.Equation(Dper[j] == D[j + 1] - D[j])
    m.Equation(0 <= u[j])
    m.Equation(u[j] <= 200)
    m.Equation(0 <= I[j])
    m.Equation(I[j] <= 10)
    m.Equation(0 <= E[j])
    m.Equation(E[j] <= 10)
    m.Equation(0 <= S[j])
    m.Equation(S[j] <= 10)
    m.Equation(0 <= R[j])
    m.Equation(R[j] <= 10)
    m.Equation(0 <= D[j])
    m.Equation(D[j] <= Dmax)
    m.Equation(Dper[j] <= DeltaDmax)
    s = s+u[j]*u[j]

j=n-1
m.Equation(D[j] == 10 - I[j] - E[j] - S[j] - R[j])
m.Equation(0 <= I[j])
m.Equation(I[j] <= 10)
m.Equation(0 <= E[j])
m.Equation(E[j] <= 10)
m.Equation(0 <= S[j])
m.Equation(S[j] <= 10)
m.Equation(0 <= R[j])
m.Equation(R[j] <= 10)
m.Equation(0 <= D[j])
m.Equation(D[j] <= Dmax)

#Objective

m.Obj(s)
m.solve()

Ivalue=[]
for j in range(n):
    Ivalue.append(I[j].value)

Evalue=[]
for j in range(n):
    Evalue.append(E[j].value)

Svalue=[]
for j in range(n):
    Svalue.append(S[j].value)

Rvalue=[]
for j in range(n):
    Rvalue.append(R[j].value)

Dvalue=[]
for j in range(n):
    Dvalue.append(D[j].value)

Dpervalue=[]
for j in range(n-1):
    Dpervalue.append(Dper[j].value)

uvalue=[]
for j in range(n-1):
    uvalue.append(u[j].value)

time = range(n)
time2 = range(n-1)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

plt.plot(time,Ivalue,'r-',label='infectious',linewidth=3.0)
plt.plot(time,Evalue,'m-',label='exposed',linewidth=3.0)
plt.plot(time,Svalue,'b-',label='susceptible',linewidth=3.0)
plt.plot(time,Rvalue,'g-',label='immune',linewidth=3.0)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('number of individuals (M)', fontsize=18)
plt.legend(fontsize=18)
plt.axhline(1, color='r', linestyle=':')
plt.xlim(0, 100)
plt.tight_layout()
plt.savefig('case_shield_phi2a.png')
plt.show()

plt.plot(time2,uvalue,'k-',label='vaccinated individuals per day',linewidth=3.0)
plt.xlim(0, 100)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('shield strength', fontsize=18)
plt.tight_layout()
plt.savefig('case_shield_phi2b.png')
plt.show()

plt.plot(time,Dvalue,'k-',label='dead individuals',linewidth=3.0)
plt.xlim(0, 100)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('deaths (M)', fontsize=18)
plt.tight_layout()
plt.axhline(0.07, color='r', linestyle=':')
plt.savefig('case_shield_phi2c.png')
plt.show()

plt.plot(time2,Dpervalue,'k-',label='dead individuals per day',linewidth=3.0)
plt.xlim(0, 100)
plt.ylim(-0.001, 0.0025)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('deaths (M) per day', fontsize=18)
plt.tight_layout()
plt.axhline(0.002, color='r', linestyle=':')
plt.savefig('case_shield_phi2d.png')
plt.show()
import numpy as np
import matplotlib.pyplot as plt


N=10
beta=0.75
epsilon=1/5
mu=1/30295
gamma=1/5
alpha = 0.006


Ts=1

lambda1 = mu*N
n = 100

S = [0] * n
I = [0] * n
E = [0] * n
R = [0] * n
D = [0] * n
Dper = [0] * n

S[0]= 9.979
I[0]= 0.001
E[0]= 0.02
R[0]= 0
D[0]= 0

print(I[0] + E[0] + S[0] + R[0])

s=0

for j in range(n-1):
    S[j+1]= S[j]+ Ts*mu*(I[j] + E[j] + S[j] + R[j]) - Ts*mu*S[j]- Ts*beta*S[j]*I[j]/N
    I[j+1]= I[j]+ Ts*epsilon*E[j] - Ts*(gamma+mu+alpha)*I[j]
    E[j+1]= E[j]+ Ts*beta*S[j]*I[j]/N - Ts*(mu+epsilon)*E[j]
    R[j+1]= R[j]+ Ts*gamma*I[j] - Ts*mu*R[j]
    D[j+1]= 10 - I[j] - E[j] - S[j] - R[j]
    Dper[j+1]=D[j + 1] - D[j]

# I[j] + E[j] + S[j] + R[j]

time = range(n)
time2 = range(n-1)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)

plt.plot(time,I,'r-',label='infectious',linewidth=3.0)
plt.plot(time,E,'m-',label='exposed',linewidth=3.0)
plt.plot(time,S,'b-',label='susceptible',linewidth=3.0)
plt.plot(time,R,'g-',label='immune',linewidth=3.0)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('number of individuals (M)', fontsize=18)
plt.legend(fontsize=18)
plt.xlim(0, 100)
plt.tight_layout()
plt.savefig('case0_a.png')
plt.show()


plt.plot(time,D,'k-',label='dead individuals',linewidth=3.0)
plt.xlim(0, 100)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('deaths (M)', fontsize=18)
plt.tight_layout()
plt.savefig('case0_c.png')
plt.show()

plt.plot(time,Dper,'k-',label='dead individuals per day',linewidth=3.0)
plt.xlim(0, 100)
plt.xlabel('time (days)', fontsize=18)
plt.ylabel('deaths (M) per day', fontsize=18)
plt.tight_layout()
plt.savefig('case0_d.png')
plt.show()


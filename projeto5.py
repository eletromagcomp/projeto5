#%% Importando
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#%% Variáveis

def F():#Fonte
    return 10**7 

def Factor(caso):
    if caso == 1:
        return 0.75
    if caso == 2:
        return 1
    if caso == 3:
        return 1.25

def omega():
    return 10**8

def alpha():
    return omega()*Factor(1)

'''
Escolhendo:
    
    Underdamped
    alpha > omega
    
    Critically damped
    alpha = omega
    
    Overdamped
    alpha < omega
'''
#%% 

def RLC(Vx, t):
    V,x = Vx
    res = np.array([x,(-2*alpha()*x - omega()**2*V + F())])
    return res
    
t = np.linspace(0.0,0.6e-6,1001)    
v,x = odeint(RLC,[0.0,0.0],t).T

#%% Plot Charge

fig, ax = plt.subplots(figsize=(8, 4))

ax.set(xlabel='Tempo', ylabel='Carga',
       title=r'Carga $\alpha$ = '+str(Factor(1))+'$\omega_0$')

ax.plot(t,x)
#ax.legend(loc = 'upper left')
ax.grid(which='both')

plt.savefig('Cargas_RLC_Simples.png',dpi=800,bbox_inches='tight')
plt.show()

#%% Plot Current
fig, ax = plt.subplots(figsize=(8, 4))

ax.set(xlabel='Tempo', ylabel='Carga',
       title=r'Carga $\alpha$ = '+str(Factor(1))+'$\omega_0$')

ax.plot(t,v)
#ax.legend(loc = 'upper left')
ax.grid(which='both')

plt.savefig('Correntes_RLC_Simples.png',dpi=800,bbox_inches='tight')
plt.show()

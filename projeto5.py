#%% Importando
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#%% Variáveis

def F():#Fonte
    return 10**7 

def Factor():
    return 0.25

def omega():
    return 10**8

def alpha():
    return omega()*Factor()

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

fontsize = 16
fig, ax = plt.subplots(figsize=(8, 4))
ax.set_xlabel('Tempo', fontsize = fontsize)
ax.plot(t,x)
#ax.legend(loc = 'upper left')
ax.grid(which='both')

plt.savefig('RLC_Simples.png',dpi=800,bbox_inches='tight')
plt.show()

#%% Plot Current
fig, ax = plt.subplots(figsize=(8, 4))
ax.set_xlabel('Tempo', fontsize = fontsize)
ax.plot(t,v)
#ax.legend(loc = 'upper left')
ax.grid(which='both')

plt.savefig('RLC_Simples.png',dpi=800,bbox_inches='tight')
plt.show()

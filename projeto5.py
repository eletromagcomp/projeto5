#%% Importando
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#%% Variáveis

def F():#Fonte/L com auto-indutância L = 1
    return 10**7
    
def alpha(): #R/2L
    return 2.5*10**10

def omega(): #(LC)^(-1/2) com capacitância C = 10^(-16)
    return 10**8

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
    res = np.array([x,(F() - ((omega())**2)*V - 2*alpha()*x)])
    return res
    
t = np.linspace(0.0,0.6e-6,1001)    
v,x = odeint(RLC,[0.0,0.0],t).T

plt.plot(t,v)
plt.show()

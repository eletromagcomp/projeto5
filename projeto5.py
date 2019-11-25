print("Hello World")

# Para um sistema RLC simples em série:

#%% Importando
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


#%% Variáveis

def F():#Fonte
    return int(1) 

def RC():
    R = 5.0
    C = 1.0e-9
    return R*C

def LC():
    L = 100.0e-9
    C = 1.0e-9
    return L*C

#%% 

def RLC(Vx, t):
    V,x = Vx
    res = np.array([x,(F() - V - RC()*x)/(LC())])
    return res
    
t = np.linspace(0.0,0.6e-6,1001)    
v,x = odeint(RLC,[0.0,0.0],t).T

plt.plot(t,v)
plt.show()

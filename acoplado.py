#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 17:01:46 2019

@author: hiro
"""

#%% Importando
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

#%% Variáveis para o RLC acoplado
def fontes():
    F1 = 0
    F2 = 0
    return [F1, F2]

def indutancia():
    L1 = 1e-2
    L2 = 1e-2
    return [L1, L2]

def capacitancia():
    C1 = 1e-5
    C2 = 1e-3
    return [C1, C2]

def resistencia():
    R1 = 500
    R2 = 0.5
    return [R1, R2]

def acoplamento_mutua():
    return 0.5

def mutua():
    k = acoplamento_mutua()
    L1, L2 = indutancia()
    M = k*np.sqrt(L1*L2)
    return M

def RLC_coupled_inicial():
    Q1_inicial = 1e-4
    Q2_inicial = 0
    I1_inicial = 0
    I2_inicial = 0
    return [Q1_inicial, Q2_inicial, I1_inicial, I2_inicial]

#%%
def RLC_coupled(QI_coupled, t):
    Q1, Q2, I1, I2 = QI_coupled
    
    F1, F2 = fontes()
    R1, R2 = resistencia()
    L1, L2 = indutancia()
    C1, C2 = capacitancia()
    M = mutua()
    
    dI1 = -L2*R1 * I1 + M*R2 *I2 - L2/C1 * Q1 + M/C2 * Q2 + L2*F1 - M*F2
    dI1 = dI1/(L1*L2 - M**2)
    dI2 = -M*R1 * I1 + L1*R2 * I2 - M/C1 * Q1 + L2/C2 * Q2 + M*F1 - L1*F2
    dI2 = dI2/(M**2 - L1*L2)
    
    res = [I1, I2, dI1, dI2]
    return res


#%%
def simulate_RLC_Coupled():
    #Parâmetros
    F1, F2 = fontes()
    R1, R2 = resistencia()
    L1, L2 = indutancia()
    C1, C2 = capacitancia()
    M = mutua()
    
    alpha1 = R1/(2*L1)
    alpha2 = R2/(2*L2)
    
    omega1 = np.sqrt(1/(L1*C1))
    omega2 = np.sqrt(1/(L2*C2))
    
    print(alpha1/omega1)
    print(alpha2/omega2)

    
    #Modos normais
    A = (C1*L1 + C2*L2)/(2*C1*C2*(L1*L2 - M**2))
    B = 1/(C1*C2 * (L1*L2 - M**2)) 
    
    modo_normal1 = np.sqrt(A + np.sqrt(A**2 - B))
    modo_normal2 = np.sqrt(A - np.sqrt(A**2 - B))
    print('omega_+: ' + str(modo_normal1))
    print('omega_-: ' + str(modo_normal2))
    
    #Condição inicial
    QI_coupled_inicial = RLC_coupled_inicial()
    
    #Simulação
    t = np.linspace(0, 0.1, 100001)
    Q1, Q2, I1, I2 = odeint(RLC_coupled, QI_coupled_inicial, t).T 
    
    
    #Plot Carga
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.set(xlabel='Tempo', ylabel='Carga')

    ax.plot(t,Q1, label='Circuito 1')
    ax.plot(t,Q2, label='Circuito 2')
    ax.grid(which='both')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.legend()

    plt.savefig('Acop_Carga.png',dpi=800,bbox_inches='tight')
    plt.show()
    
    
    
    #Plot Corrente
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set(xlabel='Tempo', ylabel='Corrente')

    ax.plot(t,I1, label='Circuito 1')
    ax.plot(t,I2, label='Circuito 2')
    ax.grid(which='both')
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax.legend()

    plt.savefig('Acop_Corrente.png',dpi=800,bbox_inches='tight')
    plt.show()
    
    
    
    
    #Modos normais
    normal1 = Q2 + (L1/M - 1/(modo_normal1**2 * M * C1 )) * Q1
    normal2 = Q2 + (L1/M - 1/(modo_normal2**2 * M * C1 )) * Q1
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(8, 4))
    
    ax1.set(xlabel='Tempo', ylabel='Carga', title = '$Q_+$')
    ax1.plot(t,normal2, label =r'$\omega \approx $' + str(int(modo_normal1)) + '$ s^{-1}$')
    ax1.grid(which='both')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax1.legend()
    
    ax2.set(xlabel='Tempo', ylabel='Carga', title = '$Q_-$')
    ax2.plot(t,normal1, label=r'$\omega \approx $' + str(int(modo_normal2)) + '$ s^{-1}$')
    ax2.grid(which='both')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax2.legend()
    
    fig.tight_layout()
    plt.savefig('Acop_Carga_Normal.png',dpi=800,bbox_inches='tight')
    plt.show()
    
#%%
    
simulate_RLC_Coupled()

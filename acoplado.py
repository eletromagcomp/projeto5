#%% Variáveis para o RLC acoplado
def fontes():
    F1 = 1
    F2 = 1
    return [F1, F2]

def indutancia():
    L1 = 1
    L2 = 1
    return [L1, L2]

def capacitancia():
    C1 = 1
    C2 = 1
    return [C1, C2]

def resistencia():
    R1 = 5
    R2 = 0.2
    return [R1, R2]

def acoplamento_mutua():
    return 0.5

def mutua():
    k = acoplamento_mutua()
    L1, L2 = indutancia()
    M = k*np.sqrt(L1*L2)
    return M

def RLC_coupled_inicial():
    Q1_inicial = 1
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
    
    dI1 = -L2*R1 * I1 + M*R2 *I2 - L2/C1 * Q1 + M/C2 * Q2 + L2*F2 - M*F1
    dI1 = dI1/(L1*L2 - M**2)
    dI2 = -M*R1 * I1 + L1*R2 * I2 - M/C1 * Q1 + L2/C2 * Q2 + M*F1 - L1*F2
    dI2 = dI2/(M**2 - L1*L2)
    
    res = [I1, I2, dI1, dI2]
    return res
    
#%%
def simulate_RLC_Coupled():
    #Condição inicial
    QI_coupled_inicial = RLC_coupled_inicial()
    R1, R2 = resistencia()
    #Simulação
    t = np.linspace(0,80,100001)
    Q1, Q2, I1, I2 = odeint(RLC_coupled, QI_coupled_inicial, t).T 
    
    #Plot Carga
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.set(xlabel='Tempo', ylabel='Carga')

    ax.plot(t,Q1)
    ax.plot(t,Q2)
    ax.grid(which='both')

    plt.savefig('Acop_Corrente_R1_' + str(R1) + '_R2_' + str(R2) + '.png',dpi=800,bbox_inches='tight')
    plt.show()
    
    #Plot Corrente
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.set(xlabel='Tempo', ylabel='Corrente')

    ax.plot(t,I1)
    ax.plot(t,I2)
    ax.grid(which='both')

    plt.savefig('Acop_Carga_R1_' + str(R1) + '_R2_' + str(R2) + '.png',dpi=800,bbox_inches='tight')
    plt.show()
    
#%%
simulate_RLC_Coupled()

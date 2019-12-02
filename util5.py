import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from celluloid import Camera
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)


class RLC_Simples(object):

    def __init__(self):
        self.Fator = 0.1
        self.Omega = 10**4
        self.F = 10**2 #Fonte/L
        self.t = np.linspace(0.0,0.4e-2,1001)
        self.Intial_Conditions = [0.0, 0.0]
        self.Flag = 1
        
    def Alpha(self):     
        '''
        Escolhendo:
            
            Underdamped
            alpha > omega
            
            Critically damped
            alpha = omega
            
            Overdamped
            alpha < omega
        '''
        return self.Omega*self.Fator
    
    def RLC(self, Vx, t):
        Q,I = Vx
        res = np.array([I,(-2*self.Alpha()*I - self.Omega**2*Q + self.F)])
        return res
    
def Create_GIF(RLC, factors):
    fig_Carga, ax_Carga = plt.subplots()
    fig_Carga.set_size_inches((16,9))
    ax_Carga.grid(which='both')
    camera_Carga = Camera(fig_Carga)
    
    fig_Corrente, ax_Corrente = plt.subplots()
    fig_Corrente.set_size_inches((16,9))
    ax_Corrente.grid(which='both')
    camera_Corrente = Camera(fig_Corrente)
    
    fsize = 26

    def plot(Q,I,RLC):
        if RLC.Fator >= 1 and RLC.Flag != 0:
            RLC.Flag = 0
            xy = [2e-3, 1.2e-6]#usar 2e-3 para Carga
            arr_img = plt.imread('/home/lordemomo/EletroComp/projeto5/comic_boom_explosion.png', format='png')
            imagebox = OffsetImage(arr_img, zoom=0.15)
            imagebox.image.axes = ax_Corrente
            ab = AnnotationBbox(imagebox, xy,
                                xybox=(120., -80.),
                                xycoords='data',
                                boxcoords="offset points",
                                pad=0.5,
                                )
            ax_Corrente.add_artist(ab)
            
            
        
        ax_Carga.plot(RLC.t, Q, color='b', lw=4)
        ax_Carga.set_xlabel('Tempo', fontsize=fsize)
        ax_Carga.set_ylabel('Carga', fontsize=fsize)
        ax_Carga.legend([r'$\alpha$ = {0:.2f}$\omega_0$'.format(RLC.Fator)], fontsize=fsize)
        ax_Carga.set_ylim((-8.0e-3,10.0e-3))
        ax_Carga.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax_Carga.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax_Carga.tick_params(labelsize=20)
        fig_Carga.tight_layout()
        camera_Carga.snap()
        
        ax_Corrente.plot(RLC.t, I, color='g', lw=4)
        ax_Corrente.set_xlabel('Tempo', fontsize=fsize)
        ax_Corrente.set_ylabel('Corrente', fontsize=fsize)
        ax_Corrente.legend([r'$\alpha$ = {0:.2f}$\omega_0$'.format(RLC.Fator)], fontsize=fsize)
        ax_Corrente.set_ylim((0.0,2.0e-6))
        ax_Corrente.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax_Corrente.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax_Corrente.tick_params(labelsize=20)
        fig_Corrente.tight_layout()
        camera_Corrente.snap()

    for i in range(len(factors)):
        RLC.Fator = factors[i]
        I,Q = odeint(RLC.RLC, RLC.Intial_Conditions, RLC.t).T
        print(r'Solved for alpha = {0:.2f} * omega_0'.format(RLC.Fator))
        plot(Q,I,RLC)
    

    animation_Carga = camera_Carga.animate()
    animation_Corrente = camera_Corrente.animate()
    
    animation_Carga.save('RLC_Simples_Carga.gif' , writer = 'imagemagick')
    animation_Corrente.save('RLC_Simples_Corrente.gif' , writer = 'imagemagick')
    
def plot_simples(RLC):
    
    fig_Carga, ax_Carga = plt.subplots()
    fig_Carga.set_size_inches((10,6))
    ax_Carga.grid(which='both')
    
    fig_Corrente, ax_Corrente = plt.subplots()
    fig_Corrente.set_size_inches((10,6))
    ax_Corrente.grid(which='both')
    fsize = 20
    
    I,Q = odeint(RLC.RLC, RLC.Intial_Conditions, RLC.t).T
    
    ax_Carga.plot(RLC.t, Q, color='b', lw=3)
    ax_Carga.set_xlabel('Tempo', fontsize=fsize)
    ax_Carga.set_ylabel('Carga', fontsize=fsize)
    ax_Carga.legend([r'$\alpha$ = {0:.2f}$\omega_0$'.format(RLC.Fator)], fontsize=fsize)
    ax_Carga.set_ylim((-8.0e-3,10.0e-3))
    ax_Carga.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax_Carga.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax_Carga.tick_params(labelsize=fsize-8)
    fig_Carga.tight_layout()
    plt.savefig('Simples_Carga.png',dpi=800,bbox_inches='tight')
    
    ax_Corrente.plot(RLC.t, I, color='g', lw=3)
    ax_Corrente.set_xlabel('Tempo', fontsize=fsize)
    ax_Corrente.set_ylabel('Corrente', fontsize=fsize)
    ax_Corrente.legend([r'$\alpha$ = {0:.2f}$\omega_0$'.format(RLC.Fator)], fontsize=fsize)
    ax_Corrente.set_ylim((0.0,2.0e-6))
    ax_Corrente.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax_Corrente.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax_Corrente.tick_params(labelsize=fsize-8)
    fig_Corrente.tight_layout()
    plt.savefig('Simples_Corrente.png',dpi=800,bbox_inches='tight')
    
#%%
def fontes():
    F1 = 100
    F2 = 0
    return [F1, F2]

def indutancia():
    L1 = 1e-2
    L2 = 1e-2
    return [L1, L2]

def capacitancia():
    C1 = 1e-6
    C2 = 1e-6
    return [C1, C2]

def resistencia():
    R1 = 1e1
    R2 = 1e1
    return [R1, R2]

def acoplamento_mutua():
    return 0.5

def mutua():
    k = acoplamento_mutua()
    L1, L2 = indutancia()
    M = k*np.sqrt(L1*L2)
    return M

def RLC_coupled_inicial():
    Q1_inicial = 0
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
    t = np.linspace(0,0.01,100001)
    Q1, Q2, I1, I2 = odeint(RLC_coupled, QI_coupled_inicial, t).T 
    #Plot Carga
    fsize = 20
    
    fig_Carga, ax_Carga = plt.subplots()
    fig_Carga.set_size_inches((10,6))
    ax_Carga.grid(which='both')
    ax_Carga.plot(t,Q1, lw=3)
    ax_Carga.plot(t,Q2, lw=3)
    ax_Carga.set_xlabel('Tempo', fontsize=fsize)
    ax_Carga.set_ylabel('Carga', fontsize=fsize)
    ax_Carga.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
    ax_Carga.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax_Carga.tick_params(labelsize=fsize-8)
    fig_Carga.tight_layout()
    plt.savefig('Acoplado_Carga.png',dpi=800,bbox_inches='tight')
    plt.show()
    
    fig_Corrente, ax_Corrente = plt.subplots()
    fig_Corrente.set_size_inches((10,6))
    ax_Corrente.grid(which='both')    
    ax_Corrente.plot(t,I1, lw=3)
    ax_Corrente.plot(t,I2, lw=3)
    ax_Corrente.set_xlabel('Tempo', fontsize=fsize)
    ax_Corrente.set_ylabel('Corrente', fontsize=fsize)
    ax_Corrente.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax_Corrente.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax_Corrente.tick_params(labelsize=fsize-8)
    fig_Corrente.tight_layout()
    plt.savefig('Acoplado_Corrente.png',dpi=800,bbox_inches='tight')
    plt.show()

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from celluloid import Camera


class RLC_Simples(object):

    def __init__(self):
        self.Fator = 1
        self.Omega = 10**4
        self.F = 10**2 #Fonte/L
        self.t = np.linspace(0.0,0.6e-2,1001)
        self.Intial_Conditions = [0.0, 0.0]
        
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
        ax_Carga.plot(RLC.t, Q, color='b', lw=4)
        ax_Carga.set_xlabel('Tempo', fontsize=fsize)
        ax_Carga.set_ylabel('Carga', fontsize=fsize)
        ax_Carga.legend([r'$\alpha$ = {0:.2f}$\omega_0$'.format(RLC.Fator)], fontsize=fsize)
        ax_Carga.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax_Carga.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
        ax_Carga.tick_params(labelsize=20)
        fig_Carga.tight_layout()
        camera_Carga.snap()
        
        ax_Corrente.plot(RLC.t, I, color='b', lw=4)
        ax_Corrente.set_xlabel('Tempo', fontsize=fsize)
        ax_Corrente.set_ylabel('Corrente', fontsize=fsize)
        ax_Corrente.legend([r'$\alpha$ = {0:.2f}$\omega_0$'.format(RLC.Fator)], fontsize=fsize)
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

import util5 as u5
import cProfile

#%%
#cProfile.run('RLC = u5.RLC_Simples()')

#cProfile.run('u5.plot_simples(RLC)')

#cProfile.run('u5.simulate_RLC_Coupled()')

#%%
RLC = u5.RLC_Simples()

u5.plot_simples(RLC)
                
u5.simulate_RLC_Coupled()

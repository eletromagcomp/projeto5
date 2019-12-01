import numpy as np
import util5 as u5
import cProfile

cProfile.run('RLC = u5.RLC_Simples()')

cProfile.run('factors = np.linspace(0.25, 1.75, 20)')

cProfile.run('u5.Create_GIF(RLC, factors)')

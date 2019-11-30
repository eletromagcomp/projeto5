import numpy as np
import util5 as u5


RLC = u5.RLC_Simples()

factors = np.linspace(0.25, 1.75, 10)

u5.Create_GIF(RLC, factors)

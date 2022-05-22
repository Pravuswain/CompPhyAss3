import monte_carlo as mc
import numpy as np

f2 = lambda x : np.exp(-x**2)

#setting a seed
mc.multiplicative_LCG.current = 2000

out2 = mc.montecarlo(f2,0,1,1000,16381,572)

print(f'Without using Importance Sampling = {out2}')

def func3(x):
    a = np.exp(1)/(np.exp(1)-1)
    y = -np.log(1-(x/a))
    p = a*np.exp(-y)
    f = np.exp(-y**2)

    return f/p




#setting a seed
mc.multiplicative_LCG.current = 2000

out2 = mc.montecarlo(func3,0,1,1000,16381,572)

print(f'Using importance Sampling = {out2}')

'''
OUTPUT
Without using Importance Sampling = 0.7495792821873427
Using importance Sampling = 0.7451042287724318as
'''
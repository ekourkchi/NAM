import numpy as np
import pandas as pd

##########################################

df = pd.read_csv('VX.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('VX.npy', A)
#C = np.load('test.npy')
print 'VX'

df = pd.read_csv('VY.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('VY.npy', A)
print 'VY'

df = pd.read_csv('VZ.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('VZ.npy', A)
print 'VZ'


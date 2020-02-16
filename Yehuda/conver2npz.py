import numpy as np
import pandas as pd

##########################################

df = pd.read_csv('VX_600_256.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('VX_600_256.npy', A)
#C = np.load('test.npy')
print('VX')

df = pd.read_csv('VY_600_256.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('VY_600_256.npy', A)
print('VY')

df = pd.read_csv('VZ_600_256.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('VZ_600_256.npy', A)
print('VZ')


df = pd.read_csv('X_600_256.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('X_600_256.npy', A)
print('X')

df = pd.read_csv('Y_600_256.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('Y_600_256.npy', A)
print('Y')

df = pd.read_csv('Z_600_256.txt',  dtype=float, header=None)
A = df.values
A = A.flatten()
np.save('Z_600_256.npy', A)
print('Z')

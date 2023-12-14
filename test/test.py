import numpy as np
import scipy.signal

A=np.array([[1,2,3,4],
[1,2,3,4],
[1,2,3,4],
[1,2,3,4]])

B=np.array([[1,2],[1,2]])

print(scipy.signal.convolve2d(A,B, mode='same', boundary='wrap'))
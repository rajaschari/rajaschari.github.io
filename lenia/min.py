import numpy as np
import scipy.signal
import matplotlib.pylab as plt

import sys
np.set_printoptions(threshold=sys.maxsize)

# orbium configuration
orbium = {
    "R": 13,
    "T": 10,
    "m": 0.15,
    "s": 0.015,
    "b": [1],
    "cells": np.array([
      [0,0,0,0,0,0,0.1,0.14,0.1,0,0,0.03,0.03,0,0,0.3,0,0,0,0],
      [0,0,0,0,0,0.08,0.24,0.3,0.3,0.18,0.14,0.15,0.16,0.15,0.09,0.2,0,0,0,0],
      [0,0,0,0,0,0.15,0.34,0.44,0.46,0.38,0.18,0.14,0.11,0.13,0.19,0.18,0.45,0,0,0],
      [0,0,0,0,0.06,0.13,0.39,0.5,0.5,0.37,0.06,0,0,0,0.02,0.16,0.68,0,0,0],
      [0,0,0,0.11,0.17,0.17,0.33,0.4,0.38,0.28,0.14,0,0,0,0,0,0.18,0.42,0,0],
      [0,0,0.09,0.18,0.13,0.06,0.08,0.26,0.32,0.32,0.27,0,0,0,0,0,0,0.82,0,0],
      [0.27,0,0.16,0.12,0,0,0,0.25,0.38,0.44,0.45,0.34,0,0,0,0,0,0.22,0.17,0],
      [0,0.07,0.2,0.02,0,0,0,0.31,0.48,0.57,0.6,0.57,0,0,0,0,0,0,0.49,0],
      [0,0.59,0.19,0,0,0,0,0.2,0.57,0.69,0.76,0.76,0.49,0,0,0,0,0,0.36,0],
      [0,0.58,0.19,0,0,0,0,0,0.67,0.83,0.9,0.92,0.87,0.12,0,0,0,0,0.22,0.07],
      [0,0,0.46,0,0,0,0,0,0.7,0.93,1,1,1,0.61,0,0,0,0,0.18,0.11],
      [0,0,0.82,0,0,0,0,0,0.47,1,1,0.98,1,0.96,0.27,0,0,0,0.19,0.1],
      [0,0,0.46,0,0,0,0,0,0.25,1,1,0.84,0.92,0.97,0.54,0.14,0.04,0.1,0.21,0.05],
      [0,0,0,0.4,0,0,0,0,0.09,0.8,1,0.82,0.8,0.85,0.63,0.31,0.18,0.19,0.2,0.01],
      [0,0,0,0.36,0.1,0,0,0,0.05,0.54,0.86,0.79,0.74,0.72,0.6,0.39,0.28,0.24,0.13,0],
      [0,0,0,0.01,0.3,0.07,0,0,0.08,0.36,0.64,0.7,0.64,0.6,0.51,0.39,0.29,0.19,0.04,0],
      [0,0,0,0,0.1,0.24,0.14,0.1,0.15,0.29,0.45,0.53,0.52,0.46,0.4,0.31,0.21,0.08,0,0],
      [0,0,0,0,0,0.08,0.21,0.21,0.22,0.29,0.36,0.39,0.37,0.33,0.26,0.18,0.09,0,0,0],
      [0,0,0,0,0,0,0.03,0.13,0.19,0.22,0.24,0.24,0.23,0.18,0.13,0.05,0,0,0,0],
      [0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.09,0.07,0.05,0.01,0,0,0,0,0]])
}

def bell(x, m, s):
    return np.exp(-((x - m) / s) ** 2 / 2)

# world size
size = 128

# scaling factor and initial position
scale = 1
cx, cy = 20, 20
cx_1, cy_1 = 50, 50
# load orbium configuration (including parameters R, T, m, s, and cells)
globals().update(orbium)
C = np.asarray(cells)

# prepare empty world
A = np.zeros([size, size])

C = scipy.ndimage.zoom(C, scale, order=0)
R *= scale
A[cx:cx + C.shape[0], cy:cy + C.shape[1]] = C
A[2*cx:2*cx + C.shape[0], 2*cy:2*cy + C.shape[1]] = np.fliplr(C)
A[-2*cx:-2*cx + C.shape[0], -2*cy:-2*cy + C.shape[1]] = np.flipud(np.fliplr(C))
A[10:10 + C.shape[0], 10:10 + C.shape[1]] = np.flipud(C)
D = np.linalg.norm(np.asarray(np.ogrid[-R+1:R+1, -R+1:R+1], dtype="object")) / R
# print(D)
# print(np.size(D)): o/p: 26x26 grid
# print(D[25][25]): test-> 1.4142135 (Correct)

K = (D < 1) * bell(D, 0.5, 0.15)
K = K / np.sum(K)
# print(np.size(K)): o/p: 26x26 grid
# print(K[12][12]): test-> 1.9388359465263956e-05 (Correct)

def growth(U):
    return bell(U, m, s) * 2 - 1

# render video
fig, ax = plt.subplots()
img = ax.imshow(A, cmap="jet", interpolation="nearest", vmin=0)
ax.set_title('World A')


for i in range(50):
    U = scipy.signal.convolve2d(A, K, mode='same', boundary='wrap')

    A = np.clip(A + 1 / T * growth(U), 0, 1)
    img.set_array(A)
    plt.pause(0.003)  # Pause to make the animation visible

plt.show()

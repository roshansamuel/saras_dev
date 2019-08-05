import numpy as np
import h5py


#Reduces the size of the array to a lower level, 2^(n-1)+1.
def restrict(function):
    global sInd, sLst

    [rx, ry, rz] = [sLst[sInd[0]], sLst[sInd[1]], sLst[sInd[2]]]
    
    restricted = np.zeros([rx + 1, ry + 1, rz + 1])

    for i in range(2, rx-1):
        for j in range(2, ry-1):
            for k in range(2, rz-1):
                restricted[i, j, k] = function[2*i - 1, 2*j - 1, 2*k - 1]

    return restricted


#Increases the size of the array to a higher level, 2^(n+1)+1.
def prolong(function):
    global sInd, sLst

    [rx, ry, rz] = [sLst[sInd[0]], sLst[sInd[1]], sLst[sInd[2]]]
    prolonged = np.zeros([rx + 1, ry + 1, rz + 1])

    [lx, ly, lz] = np.shape(function)
    for i in range(2, lx-2):
        for j in range(2, ly-2):
            for k in range(2, lz-2):
                prolonged[i*2 - 1, j*2 - 1, k*2 - 1] = function[i, j, k]
    
    for i in range(2, rx-1, 2):
        for j in range(2, ry-1, 2):
            for k in range(3, rz-1, 2):
                prolonged[i, j, k] = (prolonged[i, j, k-1] + prolonged[i, j, k+1])/2

    for i in range(2, rx-1, 2):
        for j in range(3, ry-1, 2):
            for k in range(2, rz-1):
                prolonged[i, j, k] = (prolonged[i, j-1, k] + prolonged[i, j+1, k])/2

    for i in range(3, rx-1, 2):
        for j in range(2, ry-1):
            for k in range(2, rz-1):
                prolonged[i, j, k] = (prolonged[i-1, j, k] + prolonged[i+1, j, k])/2

    return prolonged



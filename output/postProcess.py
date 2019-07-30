#!/usr/bin/python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import polyfit
import pyfftw as pf
import numpy as np
import h5py as hp
import yaml as yl
import sys

# Pyplot-specific directives
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams["font.weight"] = "medium"

Nx = 0
Ny = 0
Nz = 0
probType = 1

def parseYAML():
    global Nx, Ny
    global probType

    yamlFile = open("../parameters.yaml", 'r')
    yamlData = yl.load(yamlFile)

    Nx = 2**yamlData["Mesh"]["X Index"] + 1
    Ny = 2**yamlData["Mesh"]["Z Index"] + 1

    probType = yamlData["Program"]["Problem Type"]


def loadData(timeVal):
    global Nx, Ny
    global probType

    fileName = "Soln_{0:09.4f}.h5".format(timeVal)
    print("Processing file " + fileName + "\n")

    try:
        f = hp.File(fileName, 'r')
    except:
        print("Could not open file " + fileName + "\n")
        exit()

    if (probType == 5):
        F_data = np.zeros([Nx, Ny, 4])
        for i in range(4):
            F_key = list(f.keys())[i]
            F_data[:,:,i] = np.matrix(f[F_key]).transpose()
    elif (probType == 1):
        F_data = np.zeros([Nx, Ny, 3])
        for i in range(3):
            F_key = list(f.keys())[i]
            F_data[:,:,i] = np.matrix(f[F_key]).transpose()

    return F_data


def plotSnap(plotData):
    global Nx, Ny
    global probType

    # Plot a data frame
    fig, axes = plt.subplots(1, 1)

    qSkip = 5
    qxCount, qyCount = plotData[::qSkip,::qSkip,0].shape
    x = np.linspace(0.0, 1.0, qxCount)
    y = np.linspace(0.0, 1.0, qyCount)
    X, Y = np.meshgrid(x, y)

    if probType == 5:
        im = axes.imshow(plotData[:,:,1], cmap='jet', extent=[0,1,0,1], origin='lower')
        axes.quiver(X, Y, plotData[::qSkip,::qSkip,2], plotData[::qSkip,::qSkip,3], color='black')

    axes.tick_params(labelsize=15)

    cbar = plt.colorbar(im, ax=axes)
    cbar.ax.tick_params(labelsize=20)

    plt.show()


def plotEvolution(coeffMat, modeNumber):
    global fStr, fEnd, nSamples

    timeAxis = np.linspace(fStr, fEnd, nSamples)

    fig, axes = plt.subplots(1, 1, figsize=(12, 8))

    axes.plot(timeAxis, np.array(coeffMat[modeNumber,:])[0], linewidth=4)
    axes.axhline(y=0, color='black', linewidth=2)
    axes.axvline(x=0, color='black', linewidth=2)

    titleString = "Evolution of mode {0:1d}".format(modeNumber+1)
    axes.set_title(titleString, fontsize=20)
    axes.tick_params(labelsize=20)

    axes.set_xlabel("t", fontsize=25)
    axes.set_ylabel("Amplitude", fontsize=20, labelpad=20)
    axes.yaxis.offsetText.set_fontsize(20)

    plt.show()


def plotEnergy(modeEnergy):
    fig, axes = plt.subplots(1, 1, figsize=(10, 8))
    axes.plot(modeEnergy, linewidth=5)
    axes.set_xlabel("N", fontsize=25)
    axes.set_ylabel(r"$\frac{\lambda_m}{\sum_m \lambda_m}$", fontsize=25, rotation=0, labelpad=20)
    axes.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axes.set_title("$E_m$", fontsize=20)
    axes.tick_params(labelsize=20)
    axes.yaxis.offsetText.set_fontsize(20)

    plt.show()


def plotSpectrum(fftData):
    k = getFreqRange()

    #m, A0 = computeFit(k, fftData, pl2)

    fig = plt.figure(figsize=(12,10))  #(12,7) for laptop screen
    plt.clf()

    #ax = plt.subplot(1, 1, 1, projection='3d')
    #ax.plot_surface(k, k, np.absolute(fftData))

    nFTModes = 20
    ax = plt.subplot(1, 1, 1)

    #fftPlot = ax.imshow(np.absolute(fftData)[:nFTModes,:nFTModes], vmin=-1, vmax=1, origin='lower', cmap='jet')
    fftPlot = ax.imshow(np.absolute(fftData)[:nFTModes,:nFTModes], origin='lower', cmap='jet')

    ax.set_xticks(np.arange(0, nFTModes, nFTModes/10))
    ax.set_yticks(np.arange(0, nFTModes, nFTModes/10))
    ax.set_xticklabels(np.arange(1, nFTModes, nFTModes/10))
    ax.set_yticklabels(np.arange(1, nFTModes, nFTModes/10))
    #plt.colorbar(fftPlot)

    ax.tick_params(labelsize=20)

    cbar = plt.colorbar(fftPlot, ax=ax)
    cbar.ax.tick_params(labelsize=20)

    plt.show()


def main():
    global Nx, Ny

    print("\nPython postprocessor for RBC runs using Saras\n")

    solTime = 1.0
    try:
        solTime = float(sys.argv[1])
    except:
        print("Time value for solution file to be processed not specifed\n")

    parseYAML()

    # Matrix of snapshots - argument to loadData specifies the solution time
    A = loadData(solTime)

    plotSnap(A)

    # Perform POD analysis on the input
    #podAnalysis(A)

    # Perform Fourier analysis on the input
    #fourierAnalysis(A)


main()

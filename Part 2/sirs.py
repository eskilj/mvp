import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as pl
import re
import scipy as S
from scipy import interpolate
import os
import sys
import subprocess

#Array for import data
p1 = []
p3 = []
contourData = []
waveContourData = []

#Call "free simulator"
subprocess.call("javac Ising.java", shell = True)
subprocess.call("java Ising", shell = True)

#Call contour plot data collection:
subprocess.call("javac Sirs.java", shell = True)
subprocess.call("java Sirs", shell = True)

#read plot data from files
p1File = open("p1.txt", "r")
for line in p1File:
    p1.append(float(line))
p1File.close()

p3File = open("p3.txt", "r")
for line in p3File:
        p3.append(float(line))
p3File.close()

contourDataFile = open("data.txt", "r")
for line in contourDataFile:
    
    if re.search( r" ", line):
        words = line.split()
        contourData.append(words)
contourDataFile.close()

waveContourDataFile = open("waveData.txt", "r")
for line in waveContourDataFile:
    
    if re.search( r" ", line):
        words = line.split()
        waveContourData.append(words)
waveContourDataFile.close()

p1 = np.array(p1)
p3 = np.array(p3)
waveContourData = np.array(waveContourData)
contourData = np.array(contourData)

sp = interpolate.RectBivariateSpline(p3, p1, contourData, kx=1 , ky=1 , s = 0)
spWave = interpolate.RectBivariateSpline(p3, p1, waveContourData, kx=1 , ky=1 , s = 0)

#Make contour plot: p3 vs p1

n = 20
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
X,Y = np.meshgrid(x, y)

ax = pl.axes([0.0, 0.0, 1.0, 1.0])

pl.contourf(x, y, sp([x],[y]), 30, alpha=.75, cmap=pl.cm.hot)
C = pl.contour(x, y, sp([x],[y]), 5, colors='black', linewidth=.5)
pl.clabel(C, inline=1, fontsize=10)


pl.xticks((np.linspace(0, 1, 9)))
pl.yticks(())
pl.show()

#Make contour plot: Wave-like behaviour: p3 vs p1

pl.contourf(x, y, spWave([x],[y]), 30, alpha=.75, cmap=pl.cm.hot)
C = pl.contour(x, y, spWave([x],[y]), 5, colors='black', linewidth=.5)
pl.clabel(C, inline=1, fontsize=10)


pl.xticks((np.linspace(0, 1, 9)))
pl.yticks(())
pl.show()

#Call wave-like I/N data collection:
subprocess.call("javac Waves.java", shell = True)
subprocess.call("java Waves", shell = True)

sweeps = []
waveFunc = []

waveFuncFile = open("waveFunc.txt", "r")
for line in waveFuncFile:
    if re.search(r" ", line):
        words = line.split()
        sweeps.append(float(words[0]))
        waveFunc.append(float(words[1]))
waveFuncFile.close()

waveFunc = np.array(waveFunc)
sweeps = np.array(sweeps)

x = sweeps
y = waveFunc
pl.plot(x,y)
pl.show()

#Call wave-like I/N data collection:
subprocess.call("javac Immune.java", shell = True)
subprocess.call("java Immune", shell = True)

p11 = []
immuneList = []

waveFuncFile = open("immune.txt", "r")
for line in waveFuncFile:
    if re.search(r" ", line):
        words = line.split()
        p11.append(float(words[0]))
        immuneList.append(float(words[1]))
waveFuncFile.close()

p11 = np.array(p11)
immuneList = np.array(immuneList)

xx = p11
yy = immuneList
pl.plot(xx,yy)
pl.show()


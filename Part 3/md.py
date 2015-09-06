import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as pl
import re
import scipy as S
import os
import sys

aveE = []
aveK = []
aveU = []
count = []

waveFuncFile = open("energyData.txt", "r")
for line in waveFuncFile:
    if re.search(r" ", line):
        words = line.split()
        aveE.append(float(words[0]))
        aveK.append(float(words[1]))
        aveU.append(float(words[2]))
        count.append(len(aveE)-1)
waveFuncFile.close()

aveE = np.array(aveE)
aveK = np.array(aveK)
aveU = np.array(aveU)
count = np.array(count)



x = count
pl.plot(x,aveE, 'b', label='Total E' )
pl.plot(x,aveK, 'g',label='Kinetic E')
pl.plot(x,aveU, 'r', label='Potential E')
pl.legend()
pl.show()




tempArray = []
tempCount = []


tempFile = open("tempData.txt", "r")
for line in tempFile:
        tempArray.append(float(line))
        tempCount.append(len(tempArray) - 1)

tempFile.close()

tempArray = np.array(tempArray)
tempCount = np.array(tempCount)

pl.plot(tempCount,tempArray, 'r', label='Temperature')
pl.show()










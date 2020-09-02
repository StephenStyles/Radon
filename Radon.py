# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 09:57:47 2020

@author: sjsty
"""

import numpy as np
import numpy.random as rnd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression



def expcount(n,time,step, *args):
    countlist = np.array([0]*int((time/step)))
    
    for _ in range(n):
        templist=[int(x)//step for x in np.cumsum([rnd.exponential(a) for a in args])]
        
        for j in templist:
            if j <(time/step):
                countlist[j] = countlist[j]+1
    return countlist

def get_inputs(time,step, *args):
    
    #times = [0,15,30,45,60,75,90,105,120,135,150,165,180]
    times = list(range(0,int(time+step),int(step)))
    inputs = [0]*int((time/step))
    for i in range(int(time/step)):
        inputs[i] = inputs[i]+(np.exp(-args[0]*times[i])-np.exp(-args[0]*times[i+1]))
        progeny_sum = 0
        for r in range(len(args)):
            tmp=1
            for q in range(len(args)):
                if q==r:
                    continue
                tmp *= args[q]/(args[q]-args[r])
            progeny_sum = progeny_sum + tmp*(np.exp(-args[r]*times[i])-np.exp(-args[r]*times[i+1]))
        inputs[i] = inputs[i]+progeny_sum
    return inputs

LRn222 = 3.825*24*60*60 # the half-life for Rn-222 (in seconds)
LPo218 = 3.05*60        # the half-life for Po-218
DC1HL = np.array([LRn222,LPo218])
DC1Lambda = np.log(2)/DC1HL # the decay constants for this first decay chain in units of 1/min
DC1Lambda

LRn220 = 54.5 # the half-life for Rn-220 (in seconds)
LPo216 = 0.158        # the half-life for Po-216
DC2HL = np.array([LRn220,LPo216])
DC2Lambda = np.log(2)/DC2HL # the decay constants for this first decay chain in units of 1/min
DC2Lambda

n=8000
m=8000

step = 60
time = 3600
counts = expcount(n,time,step,1/DC1Lambda) + expcount(m,time,step,1/DC2Lambda)

X = np.array([get_inputs(time,step,*(DC1Lambda)), get_inputs(time,step,*(DC2Lambda))]).transpose()


#indices=[0,6,len(counts)-1]
#fitX = [X[index] for index in indices]
#fitY = [counts[index] for index in indices]


Rn222 = []
Rn220 = []

for _ in range(10):
    counts = expcount(n,time,step,1/DC1Lambda) + expcount(m,time,step,1/DC2Lambda)
    reg = LinearRegression(fit_intercept = False).fit(X, counts)
    Rn222.append(reg.coef_[0])
    Rn220.append(reg.coef_[1])
    
print("The mean for Rn222 is {0} with standard deviation {1}".format(np.mean(Rn222), np.std(Rn222)))
print("The mean for Rn220 is {0} with standard deviation {1}".format(np.mean(Rn220), np.std(Rn220)))


time = list(range(300,1500,60))
step = [10,15,20,30,60]


stdmatrix = np.zeros((len(time),len(step)))

for i in range(len(time)):
    for j in range(len(step)):
        X = np.array([get_inputs(time[i],step[j],*(DC1Lambda)), get_inputs(time[i],step[j],*(DC2Lambda))]).transpose()

        Rn222 = []

        for _ in range(5):
            counts = expcount(n,time[i],step[j],1/DC1Lambda) + expcount(m,time[i],step[j],1/DC2Lambda)
            reg = LinearRegression(fit_intercept = False).fit(X, counts)
            Rn222.append(reg.coef_[0])
        
        tmp = np.std(Rn222)
        print(tmp)   
        stdmatrix[i,j] = tmp
stdmatrix
np.savetxt("std.txt",stdmatrix)

from mpl_toolkits.mplot3d import Axes3D

X2, Y2 = np.meshgrid(step, time)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X2, Y2, stdmatrix, cmap="coolwarm")

plt.show()



n = 10
step = 60
time = 600
inputs = []
times = list(range(0,int(time+step),int(step)))
inputs = [0]*int((time/step))
for i in range(int(time/step)):
    inputs[i] = inputs[i]+(np.exp(-DC1Lambda[1]*times[i])-np.exp(-DC1Lambda[1]*times[i+1]))


counts = expcount(n,time,step,1/DC1Lambda[1])
X = np.array([inputs]).transpose()
reg = LinearRegression(fit_intercept = False).fit(X, counts)
reg.coef_[0]

Po218 = []

for _ in range(100):
   counts = expcount(n,time,step,1/DC1Lambda[1])
   reg = LinearRegression(fit_intercept = False).fit(X, counts)
   Po218.append(reg.coef_[0])
    
print("The mean for Po218 is {0} with standard deviation {1}".format(np.mean(Po218), np.std(Po218)))

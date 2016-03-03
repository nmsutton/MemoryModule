'''
Fitting change in rate of spikes to 
synapse weight between NN layers

Weight: 0.0
Layer 2 Spikes: 600
_
0.5
1750
_
1.0
2900
________
Found that curve is linear!
y = 0 x2 + 2300 x + 600
'''

import numpy as np
from scipy.optimize import curve_fit

def func(x, a, b, c):
	return a * np.exp(-b * x) + c

xdata = np.linspace(0, 1, 3)#[0.0,0.5,1.0]#
print(xdata)
y = func(xdata, 2.5, 1.0, 0.0)
ydata = y + [600.0,1750.0,2900.0]#y + 0.2 * np.random.normal(size=len(xdata))

popt, pcov = curve_fit(func, xdata, ydata)

print(popt)
print(pcov)
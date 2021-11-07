# -*- coding: utf-8 -*-
"""
This program imports data from the file specified by the string filename
The first line of the file is ignored (assuming it's the name of the variables)
After that the data file needs to be formatted: 
number space number space number space number newline
Do NOT put commas in your data file!!
The data file should be in the same directory as this python file
The data should be in the order:
x_data y_data x_uncertainty y_uncertainty

Then this program tries to fit a function to the data points
It plots the data as dots with errorbars and the best fit line
It then prints the best fit information
After that, it plots the "residuals": ydata - fittedfunction(xdata)
That is it subtracts your fitted ideal values from your actual data to see 
what is "left over" from the fit
Ideally your residuals graph should look like noise, otherwise there is more
information you could extract from your data (or you used the wrong fitting
function)

If you want to change the file name, that's the next line below this comment
"""

segmented = True

filename = ""
if segmented:
    filename = "lab_2_segmented_peaks.txt"
else:
    filename = "lab_2_processed_peaks.txt"
# change this if your filename is different

import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt
from pylab import loadtxt
import math

data=loadtxt(filename, usecols=(0,1,2,3), skiprows=0, unpack=True)
# load filename, take columns 0 & 1 & 2 & 3, skip 0 rows, unpack=transpose x&y

xdata=data[0]
ydata=data[1]
xerror=data[2]
yerror=data[3]
# finished importing data, naming it sensibly

def my_func(p,a,b,c):
    return a+b*p+c*p**2
# this is the function we want to fit. the first variable must be the
# x-data (period), the rest are the unknown constants we want to determine

init_guess=(1.0,0.0,2.0)
# your initial guess of (b,c)

popt, pcov = optimize.curve_fit(my_func, xdata, ydata, p0=init_guess)
# we have the best fit values in popt[], while pcov[] tells us the uncertainties

a=popt[0]
b=popt[1]
c=popt[2]
# best fit values are named nicely
u_a=pcov[0,0]**(0.5)
u_b=pcov[1,1]**(0.5)
u_c=pcov[2,2]**(0.5)
# uncertainties of fit are named nicely

def fitfunction(p):
    return a+b*p+c*p**2
#fitfunction(p) gives you your ideal fitted function, i.e. the line of best fit

start=min(xdata)
stop=max(xdata)
xs=np.arange(start,stop,(stop-start)/(20 if segmented else 9671)) # fit line has 9671 points
curve=fitfunction(xs)
# (xs,curve) is the line of best fit for the data in (xdata,ydata)

curve_up_error = curve + 0.02
curve_down_error = curve - 0.02

fig, (ax1,ax2) = plt.subplots(2, 1)
fig.subplots_adjust(hspace=0.6)
#hspace is horizontal space between the graphs

t_0liney=[0.848255,0.848255]
t_0linex=[start,stop]

ax1.errorbar(xdata,ydata,yerr=yerror,xerr=xerror,fmt=".",zorder=1)
# plot the data, fmt makes it data points not a line
ax1.plot(xs,curve,zorder=2)
ax1.plot(xs,curve_up_error,zorder=2, color="tab:gray")
ax1.plot(xs,curve_down_error,zorder=2, color="tab:gray")
ax1.plot(t_0linex,t_0liney,zorder=2, color="tab:green")
# plot the best fit curve on top of the data points as a line

#ax1.set_ylim(0, 1.2)
ax1.set_xlim(-1.6, 1.6)

ax1.set_xlabel("Starting Angle (rad)")
ax1.set_ylabel("Period (s)")
ax1.set_title("Best + Constant Fit of Starting Angle (rad) vs. Period (s)")
# HERE is where you change how your graph is labelled

print("A:", a, "+/-", u_a)
print("B:", b, "+/-", u_b)
print("C:", c, "+/-", u_c)
# prints the various values with uncertainties

residual=ydata-fitfunction(xdata)

# find the residuals
zeroliney=[0,0]
zerolinex=[start,stop]
# create the line y=0

RSS = 0
TSS = 0

mean = sum(ydata) / len(ydata)

for r in residual:
    RSS += r**2

for y in ydata:
    TSS += (y - mean)**2

R2 = 1 - RSS / TSS
print("R2:", R2)

chi2_terms = []

for i in range(len(residual)):
    chi2_terms.append((residual[i] ** 2) / fitfunction(xdata[i]))

chi2_value = sum(chi2_terms)
print("Chi2:", chi2_value)

ax2.errorbar(xdata,residual,yerr=yerror,xerr=xerror,fmt=".",zorder=1)
# plot the residuals with error bars
ax2.plot(zerolinex,zeroliney,zorder=2)
# plotnthe y=0 line on top

ax2.set_xlabel("Starting Angle (rad)")
ax2.set_ylabel("Residuals of Period (s)")
ax2.set_title("Residuals of Starting Angle (rad) vs. Period (s)")
# HERE is where you change how your graph is labelled

plt.show()
# show the graph
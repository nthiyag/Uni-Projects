import numpy as np
import scipy
import matplotlib.pyplot as plt
import pandas as pd

'''
Remember that the j_0 function is kinda like e^ix
It's a bit different bc instead of i in the argument it's more like (-i)^0.5
But then you have a similar euler approximation type vibe going on for ber0 and bei0

NORMALIZING THE OUTPUT OF BESSEL results in a complex number
This is multiplied by e^iwt
So it acts like a constant phase shift - that's the whole point of the R(r) term
'''

ID = 8      #radius in mm
OD = 16     #radius in mm

def load_file(filename):
    raw_data = pd.read_csv(filename, sep='\t', header=0, index_col=0)
    raw_data = raw_data.drop(0)

    data = {}
    
    for col in list(raw_data):
        data[col] = (raw_data[col].copy()).to_frame()

        data[col] = data[col].rename(columns={col: "inside"})

        out_data = [(100 if i % (int(col) * 2) > (int(col) - 1) else 0) if i < data[col].count()[0] else None for i in range(0, data[col].size)]

        data[col]["outside"] = out_data
    
    return data

def find_direct_phaseshifts(data):
    allshifts = {}
    
    for col in data:
        inside = [e for e in list(data[col]["inside"]) if not np.isnan(e)]
        outside = [e for e in list(data[col]["outside"]) if not np.isnan(e)]

        rawshifts = []

        for i in range(len(outside)-1):
            if outside[i] != outside[i+1]:
                wait_for = np.sign(outside[i+1] - outside[i])
            else:
                continue
            
            j = 0
            early_break = False
            while np.sign(inside[i+j+1] - inside[i+j]) != wait_for:
                j += 1

                if i+j+1 > len(inside) - 1:
                    early_break = True
                    break
            
            if not early_break:
                rawshifts.append((j / int(col)) * np.pi)
        
        rawshift = sum(rawshifts) / len(rawshifts)
        print(rawshifts)
        allshifts[col] = rawshift % (2 * np.pi) if rawshift % (2 * np.pi) <= np.pi else rawshift % (2 * np.pi) - (2 * np.pi)

    print(allshifts)
        
    return allshifts

def find_sinusoid_phaseshifts(data):
    phaseshifts = {}
    errors = {}

    for col in data:
        inside = [e for e in list(data[col]["inside"]) if not np.isnan(e)]
        outside = [e for e in list(data[col]["outside"]) if not np.isnan(e)]
        time = [list(data[col].index)[i] for i in range(len(list(data[col]["inside"]))) if not np.isnan(list(data[col]["inside"])[i])]

        rawshift, covar = scipy.optimize.curve_fit(lambda t, a: sinusoid_approx(t, a, col), time, inside)
        phaseshift = rawshift[0] % (2 * np.pi) if rawshift[0] % (2 * np.pi) <= np.pi else rawshift[0] % (2 * np.pi) - (2 * np.pi)
        phaseshifts[col] = phaseshift

        errors[col] = np.sqrt(np.diag(covar))[0]

        sin_in_data = ([sinusoid_approx(t, phaseshift, col) for t in time]) + ([None for i in range(len(list(data[col].index))-len(time))])

        data[col]["sin_approx_inside"] = sin_in_data

    return data, phaseshifts, errors

def sinusoid_approx(t, a, col):
    return 50 * np.sin(-(np.pi/(int(col) * 10)) * t + a) + 50

def visualize(data):
    for col in data:
        plt.figure()
        plt.plot(data[col])

    plt.show()

def ber0(z, terms=200):
    result = 1
    for i in range(1, terms-1):
        denom = 1
        for j in range(1, 2*i + 1):
            denom *= (2*j)**2
        try:
            term = (z ** (4*i) / denom)
        except OverflowError:
            term = 0
        result = result + term if i % 2 == 0 else result - term
    return result

def bei0(z, terms=200):
    result = 0
    for i in range(0, terms):
        denom = 1
        for j in range(1, 2*i + 2):
            denom *= (2*j)**2
        try:
            term = (z ** (4*i + 2) / denom)
        except OverflowError:
            term = 0
        result = result + term if i % 2 == 0 else result - term
    return result

def j0(z, terms=200):
    result = 1
    for i in range(1, terms-1):
        denom = 1
        for j in range(1, i + 1):
            denom *= (2*j)**2
        try:
            term = (z ** (2*i) / denom)
        except OverflowError:
            term = 0
        result = result + term if i % 2 == 0 else result - term
    return result

def dber0(z, terms=200):
    result = 1
    for i in range(1, terms-1):
        denom = 1
        for j in range(1, 2*i + 1):
            denom *= (2*j)**2
        try:
            term = (4*i) * (z ** (4*i - 1) / denom)
        except OverflowError:
            term = 0
        result = result + term if i % 2 == 0 else result - term
    return result

def dbei0(z, terms=200):
    result = 0
    for i in range(0, terms):
        denom = 1
        for j in range(1, 2*i + 2):
            denom *= (2*j)**2
        try:
            term = (4*i + 2) * (z ** (4*i + 1) / denom)
        except OverflowError:
            term = 0
        result = result + term if i % 2 == 0 else result - term
    return result

def phase(z, terms=200):
    return np.arctan2(bei0(z, terms), ber0(z, terms))

def vphase(z, terms=200):
    vphase_sig = np.vectorize(phase)
    return vphase_sig(z, terms)

def inv_phase(phi):
    return scipy.optimize.root_scalar(lambda z: phase(z) - phi, x0=1, x1=0.9).root

def solve_m_old(phi, omega, r):
    z = abs(inv_phase(phi))

    return (omega * (r ** 2)) / (z ** 2)

def delta_phase(x, r_in_mm, r_out_mm):
    rawdshift = phase(x * r_out_mm) - phase(x * r_in_mm)
    return rawdshift % (2 * np.pi) if rawdshift % (2 * np.pi) <= np.pi else rawdshift % (2 * np.pi) - (2 * np.pi)

def vdelta_phase(x, r_in_mm, r_out_mm):
    vdelta_phase_sig = np.vectorize(delta_phase)
    return vdelta_phase_sig(x, r_in_mm, r_out_mm)

def solve_m(dphi, omega, r_in_mm, r_out_mm):
    guess = 0
    while guess < 100:
        x = scipy.optimize.root_scalar(lambda x: (delta_phase(x, r_in_mm, r_out_mm) - dphi), x0=guess, x1=guess+0.1).root
        if not np.isnan(x):
            break
        else:
            guess += 0.1
    
    return ((omega) / (x ** 2)) * (10 ** -6)



# print("48", phase(1.791))
# print("24", phase(2.533))
# print("12", phase(3.582))
# print("6", 2*np.pi+phase(5.065))

# print("48", (phase(1.791) / (2*np.pi)) * 48 * 2)
# print("24", (phase(2.533) / (2*np.pi)) * 24 * 2)
# print("12", (phase(3.582) / (2*np.pi)) * 12 * 2)
# print("6", (1+(phase(5.065) / (2*np.pi))) * 6 * 2)

# x = np.linspace(0, 5, num=100)
# plt.plot(x, vdelta_phase(x, ID, OD))
# plt.show()

data = load_file("data.csv")
# phaseshifts = find_avg_phaseshifts(data)
# print(phaseshifts)

data, phaseshifts, errors = find_sinusoid_phaseshifts(data)
print(phaseshifts, errors)

for col in phaseshifts:
    omega = (2 * np.pi) / (int(col) * 10 * 2)
    dphi = phaseshifts[col]

    print(str(col) + ":", solve_m(dphi, omega, ID, OD))

visualize(data)
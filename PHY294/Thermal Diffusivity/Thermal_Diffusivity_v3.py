import math
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


#data functions
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

def find_sinusoid_phaseshifts(data):
    phaseshifts = {}
    errors = {}
    chi_squares = {}

    for col in data:
        inside = [e for e in list(data[col]["inside"]) if not np.isnan(e)]
        outside = [e for e in list(data[col]["outside"]) if not np.isnan(e)]
        time = [list(data[col].index)[i] for i in range(len(list(data[col]["inside"]))) if not np.isnan(list(data[col]["inside"])[i])]

        y_unc = 3.5

        rawshift, covar = scipy.optimize.curve_fit(lambda t, x_shift: sinusoid_approx(t, x_shift, 50, 50, 0, col), time, inside, sigma=[y_unc for i in range(len(time))], absolute_sigma=True)
        phaseshift = rawshift[0] % (2 * np.pi) if rawshift[0] % (2 * np.pi) <= np.pi else rawshift[0] % (2 * np.pi) - (2 * np.pi)
        phaseshifts[col] = phaseshift

        overrawshift, overcovar = scipy.optimize.curve_fit(lambda t, x_shift, y_shift, y_scale, a1: sinusoid_approx(t, x_shift, y_shift, y_scale, a1, col), time, inside, sigma=[y_unc for i in range(len(time))], absolute_sigma=True)
        overrawshift[0] = overrawshift[0] + np.pi if overrawshift[2] < 0 else overrawshift[0]
        overrawshift[2] = abs(overrawshift[2])
        overphaseshift = overrawshift[0] % (2 * np.pi) if overrawshift[0] % (2 * np.pi) <= np.pi else overrawshift[0] % (2 * np.pi) - (2 * np.pi)
        
        std_dev = np.sqrt(np.diag(covar))[0]
        model_dev = abs(phaseshift-overphaseshift)
        print(col, std_dev, model_dev)

        #print(std_dev, model_dev)

        errors[col] = std_dev + model_dev

        sin_in_data = ([sinusoid_approx(t, phaseshift, 50, 50, 0, col) for t in time]) + ([None for i in range(len(list(data[col].index))-len(time))])
        over_sin_in_data = ([sinusoid_approx(t, overphaseshift, overrawshift[1], overrawshift[2], overrawshift[3], col) for t in time]) + ([None for i in range(len(list(data[col].index))-len(time))])

        chi_squares[col] = sum([((inside[i] - sinusoid_approx(time[i], phaseshift, 50, 50, 0, col)) ** 2) / (y_unc ** 2) for i in range(len(time))]) / (len(time) - 1)

        data[col]["residuals"] = [inside[i] - sin_in_data[i] for i in range(len(time))] + ([None for i in range(len(list(data[col].index))-len(time))])
        #data[col]["sin_approx_inside"] = sin_in_data
        #data[col]["over_sin_approx_inside"] = over_sin_in_data

    return data, phaseshifts, errors, chi_squares

def sinusoid_approx(t, x_shift, y_shift, y_scale, a1, col):
    return a1* t + y_scale * np.sin(-(np.pi/(int(col) * 10)) * t + x_shift) + y_shift

def visualize(data):
    for col in data:
        plt.figure()
        plt.plot(data[col]["residuals"])
        plt.title("Temp. vs. Time (sinusoidal fit residuals) | Period: " + str(int(col) * 20) + "s")
        plt.xlabel("Time (s)")
        plt.ylabel("Temperature (°C)")
    plt.show()


#math functions
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
    result = 0
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

def delta_phase(x, r_in, r_out):
    rawdshift = phase(x * r_out) - phase(x * r_in)
    return rawdshift % (2 * np.pi) if rawdshift % (2 * np.pi) <= np.pi else rawdshift % (2 * np.pi) - (2 * np.pi)

def vdelta_phase(x, r_in, r_out):
    vdelta_phase_sig = np.vectorize(delta_phase)
    return vdelta_phase_sig(x, r_in, r_out)

def dp_dz(omega, m, r):
    x = (omega / m) ** (1/2)
    return  (1/(1 + (bei0(x * r) / ber0(x * r)) ** 2)) * \
            ((dbei0(x * r) / ber0(x * r)) - ((bei0(x * r) * dber0(x * r)) / (ber0(x * r) ** 2)))

def dp_dm(omega, m, r_in, r_out):
    x = (omega / m) ** (1/2)
    return  dp_dz(omega, m, r_out) * (-(1/2) * (omega ** (1/2)) * (m ** (-3/2)) * r_out) - \
            dp_dz(omega, m, r_out) * (-(1/2) * (omega ** (1/2)) * (m ** (-3/2)) * r_in)


#numerical derivatives
def dm_domega(delta_phi, omega, r_in, r_out, d=0.00001):
    m1 = solve_m(delta_phi, omega, r_in, r_out)
    m2 = solve_m(delta_phi, omega+d, r_in, r_out)
    return (m2 - m1) / d

def dm_dr_in(delta_phi, omega, r_in, r_out, d=0.00001):
    m1 = solve_m(delta_phi, omega, r_in, r_out)
    m2 = solve_m(delta_phi, omega, r_in+d, r_out)
    return (m2 - m1) / d

def dm_dr_out(delta_phi, omega, r_in, r_out, d=0.00001):
    m1 = solve_m(delta_phi, omega, r_in, r_out)
    m2 = solve_m(delta_phi, omega, r_in, r_out+d)
    return (m2 - m1) / d

def solve_m(delta_phi, omega, r_in, r_out):
    guess = 0
    while guess < 10000:
        x = scipy.optimize.root_scalar(lambda x: (delta_phase(x, r_in, r_out) - delta_phi), x0=guess, x1=guess+100).root
        if not np.isnan(x):
            break
        else:
            guess += 100
    
    return ((omega) / (x ** 2))

def solve_m_unc(omega, m, r_in, r_out, delta_phi_unc, omega_unc, r_in_unc, r_out_unc):
    return (((1/dp_dm(omega, m, r_in, r_out)) * delta_phi_unc) ** 2 + (dm_domega(delta_phi, omega, r_in, r_out) * omega_unc) ** 2 + (dm_dr_in(delta_phi, omega, r_in, r_out) * r_in_unc) ** 2 + (dm_dr_out(delta_phi, omega, r_in, r_out) * r_out_unc) ** 2) ** (1/2)


#main

# print("48", phase(1.791))
# print("24", phase(2.533))
# print("12", phase(3.582))
# print("6", 2*np.pi+phase(5.065))

# print("48", (phase(1.791) / (2*np.pi)) * 48 * 2)
# print("24", (phase(2.533) / (2*np.pi)) * 24 * 2)
# print("12", (phase(3.582) / (2*np.pi)) * 12 * 2)
# print("6", (1+(phase(5.065) / (2*np.pi))) * 6 * 2)

# x = np.linspace(0, 5, num=100)
# plt.plot(x, vdelta_phase(x, ID * 0.001, OD * 0.001))
# plt.plot(x, bei0(x))
# plt.plot(x, ber0(x))
# plt.figure()
# plt.plot(x, dbei0(x))
# plt.plot(x, dber0(x))
# plt.show()

data = load_file("data.csv")
# phaseshifts = find_avg_phaseshifts(data)
# print(phaseshifts)

data, phaseshifts, phaseshift_errors, chi_squares = find_sinusoid_phaseshifts(data)
print(phaseshifts, phaseshift_errors)

trials = {'6': 4, '12': 2, '24': 2, '48': 1}

for col in phaseshifts:
    period = (int(col) * 10 * 2)
    omega = (2 * np.pi) / period
    delta_phi = phaseshifts[col]
    r_in = ID * 0.001
    r_out = OD * 0.001
    
    period_unc = 5
    omega_unc = (((2 * np.pi) / ((int(col) * 10 * 2) ** 2)) * period_unc)
    delta_phi_unc = phaseshift_errors[col]
    r_in_unc = 0.0005
    r_out_unc = 0.0005

    m = solve_m(delta_phi, omega, r_in, r_out)
    m_unc = solve_m_unc(omega, m, r_in, r_out, delta_phi_unc, omega_unc, r_in_unc, r_out_unc)
    print(str(col) + " | m:", m, "m_unc:", m_unc, "| Chi-Squared:", chi_squares[col])

visualize(data)
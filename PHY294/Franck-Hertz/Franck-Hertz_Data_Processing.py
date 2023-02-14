import matplotlib.pyplot as plt
import numpy as np
import scipy

def find_peak_idx(data, window_len=8, err_patience=2):
    peaks = []
    for i in range(window_len, len(data)-window_len):
        window_up = data[i-window_len:i+1]
        window_down = data[i:i+window_len+1]

        if sum([window_up[j] < window_up[j+1] for j in range(len(window_up)-1)]) >= len(window_up)-1 - err_patience and sum([window_down[j] > window_down[j+1] for j in range(len(window_up)-1)]) >= len(window_up)-1 - err_patience:
            peaks.append(i)
    
    return peaks

def find_contiguous(lst):
    blocks = []
    tmp_block = []
    for i in range(len(lst)):
        if len(tmp_block) == 0:
            tmp_block.append(lst[i])
        elif lst[i] - 1 == lst[i-1]:
            tmp_block.append(lst[i])
        else:
            blocks.append(tmp_block)
            tmp_block = [lst[i]]
    
    if len(tmp_block) > 0:
        blocks.append(tmp_block)
    
    return blocks

def moving_average(a, window_len=2):
    out = []
    max_error = 0
    for i in range(window_len, len(a)-window_len):
        out.append(sum(a[i-window_len:i+window_len+1]) / (window_len*2 + 1))
        error = abs(a[i] - out[-1])
        max_error = error if error > max_error else max_error
    
    return out, max_error

def process(filename):
    file = open(filename)

    for i in range(3):
        line = file.readline()

    x, y = [], []
    while line:
        line_data = line.split()
        x.append(float(line_data[0]))
        y.append(float(line_data[1]))
        line = file.readline()

    x, x_ma_error = moving_average(x, 2)
    y, y_ma_error = moving_average(y, 2)

    raw_peaks = find_peak_idx(y)

    x_raw_peaks = [x[j] for j in raw_peaks]
    y_raw_peaks = [y[j] for j in raw_peaks]

    x_peaks = []
    y_peaks = []
    for block in find_contiguous(raw_peaks):
        x_peaks.append(x[block[len(block)//2]])
        y_peaks.append(y[block[len(block)//2]])
    
    d_volts = [x_peaks[j] - x_peaks[j-1] for j in range(1, len(x_peaks))]

    x_diffs = [round(x[j] - x[j-1], 5) for j in range(1, len(x))]

    plt.figure()
    plt.errorbar(x, y, xerr=0.08, fmt='c', ecolor='gray', capsize=5, zorder=0)
    plt.scatter(x, y, s=[25 for i in range(len(x))], c='c', marker='x', zorder=1)
    plt.plot(x_raw_peaks, y_raw_peaks, 'x', color='yellow', zorder=2)
    plt.plot(x_peaks, y_peaks, 'o', color='orange', zorder=3)
    plt.title("Franck-Hertz Current vs. Voltage: Trial 1")
    plt.xlabel("Voltage (V)")
    plt.ylabel("Current (10^-8 A)")

    # plt.errorbar(("1-2", "2-3", "3-4", "4-5"), [d_volt-5.0031210937 for d_volt in d_volts], fmt='o', yerr=0.08, ecolor='gray', capsize=5)
    # plt.title("Residuals of Constant ΔE Fit")
    # plt.axhline(0)
    # plt.xlabel("Peak #")
    # plt.ylabel("Energy Separation (eV)")

    # plt.errorbar(("1-2", "2-3", "3-4", "4-5"), [d_volts[i] - (0.14053710936000036 * i + 4.7923154296599995) for i in range(len(d_volts))], fmt='o', yerr=0.08, ecolor='gray', capsize=5)
    # plt.title("Residuals of Linear ΔE Fit")
    # plt.axhline(0)
    # plt.xlabel("Peak #")
    # plt.ylabel("Energy Separation (eV)")
    
    return d_volts, x_diffs, x_ma_error

d_volts = []
x_diffs = set()
x_ma_errors = []

for i in range(3):
    result = process("fh_data_" +str(i+1) +".txt")
    d_volts.extend(result[0])
    for j in result[1]:
        x_diffs.add(j)
    x_ma_errors.append(result[2])

x_diffs = sorted(list(x_diffs))
x_ma_max_error = max(x_ma_errors)

print("X Jumps:", x_diffs)

avg_volt_sep = sum(d_volts) / len(d_volts)
print("X Chi-Squared:", scipy.stats.chisquare(d_volts, [avg_volt_sep for i in range(len(d_volts))], 0))
print("X MA Max Error:", x_ma_max_error)

print("Avg Volt Sep:", avg_volt_sep)
print("Wavelength:", 4.135*10**-15 * 2.99*10**8 / avg_volt_sep)
print("Wavelength Unc:", abs((-4.135*10**-15 * 2.99*10**8 / (avg_volt_sep**2)) * 0.08))

x_linear_fit = scipy.stats.linregress((0,1,2,3,0,1,2,3,0,1,2,3), d_volts)
print("X Linear Fit:", x_linear_fit)
#plt.plot([i for i in range(4)], [x_linear_fit[0] * (i%4) + x_linear_fit[1] for i in range(4)])

observed = d_volts
expected = [x_linear_fit[0] * (i%4) + x_linear_fit[1] for i in range(len(d_volts))]

chi_squared = sum([(observed[i] - expected[i])**2 / expected[i] for i in range(len(d_volts))])
print(chi_squared)
print("X Linear Chi-Squared:", scipy.stats.chisquare(d_volts, [x_linear_fit[0] * (i%4) + x_linear_fit[1]  for i in range(len(d_volts))], 0))

plt.show()
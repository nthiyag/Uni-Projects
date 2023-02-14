import json
import math

import matplotlib.pyplot as plt
import better_regression_optimizer as reg

segmented_data = open('segmented_data.json')
samples = json.load(segmented_data)
segmented_data.close()

Q1, Q2, rs = [], [], []
for (number, sample) in samples.items():
    data = sample["data"]
    u_range = sample["upslope_range"]
    d_range = sample["downslope_range"]
    data_u = data[u_range[0]:u_range[1]]
    data_d = data[d_range[0]:d_range[1]]

    u_opt_lin = reg.opt_linfit(list(range(len(data_u))), data_u)
    d_opt_lin = reg.opt_linfit(list(range(len(data_d))), data_d)

    p_o = 875.3
    p_a = 1.204
    g = 9.81
    eta = 1.827*10**-5
    d = 6*10**-3

    C = 9*math.pi*d*(2*eta**3/(g*(p_o-p_a)))**(1/2)

    fps = 5
    m_per_px = 1/(540 * 10 ** 3)
    v_t = u_opt_lin[1][0] * fps * m_per_px
    v_2 = -d_opt_lin[1][0] * fps * m_per_px
    V_stop = sample["stp_volt"]
    V_up = 699

    q1 = C * (v_t ** 1.5) / V_stop
    q2 = C * (v_t + v_2) * (v_t ** 0.5) / V_up

    r = 3*(eta*v_t/(2*g*(p_o-p_a)))**(1/2)
    rs.append(r)

    Q1.append(q1)
    Q2.append(q2)

print("Average drop size:", sum(rs)/len(rs))

bounds = [(1.602*(0.675))*10**-19+i*1.602*10**-19 for i in range(14)]
Q2_flattened = []

for Q in Q2:
    for i in range(len(bounds)-1):
        if Q > bounds[i] and Q < bounds[i+1]:
            Q2_flattened.append(Q / (i+1))
            break

print(Q2_flattened)

print(sum(Q2_flattened)/len(Q2_flattened))

# plt.scatter(Q1, Q2)
# plt.xlabel("Q1: Method 1")
# plt.ylabel("Q2: Method 2")

plt.scatter(rs, Q2)
plt.xlabel("Drop Radius")
plt.ylabel("Q2")

#plt.hlines(bounds, xmin=0, xmax=max(Q1), color='tab:orange')

plt.show()
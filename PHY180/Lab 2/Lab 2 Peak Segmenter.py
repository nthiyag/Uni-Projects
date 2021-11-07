import math

peaks = open("lab_2_processed_peaks.txt", "r")

segmented = open("lab_2_segmented_peaks.txt", "w")

x_uncertainty = "0.06"
y_uncertainty = "0.02"

discard_lines = 0
for i in range(discard_lines):
    peaks.readline()

#degrees
targets_array = [-80, 80, -75, 75, -65, 65, -60, 60, -55, 55, -50, 50, -45, 45, -40, 40, -35, 35, -30, 30, -25, 25, -20, 20, -15, 15, -10, 10, -5, 5]
targets = {}
threshold = 3

for i in targets_array:
    targets[math.radians(i)] = []

threshold = math.radians(threshold)

avg_periods = []

while True:
    line = peaks.readline().split()
    if not line:
        break

    amplitude = float(line[0])
    period = float(line[1])

    print("amplitude (deg):", math.degrees(amplitude))

    for bucket in targets:
        if abs(amplitude - bucket) < threshold:
            targets[bucket].append(period)

for bucket in targets:
    if len(targets[bucket]) != 0:
        avg_periods.append(str(math.degrees(bucket)) +": " +f"{sum(targets[bucket]) / len(targets[bucket]): 4f}")
        segmented.write(str(bucket) +" " +str(sum(targets[bucket]) / len(targets[bucket])) +" " +x_uncertainty +" " +y_uncertainty +"\n")

for i in targets:
    print(len(targets[i]))

peaks.close()
segmented.close()

print("done")
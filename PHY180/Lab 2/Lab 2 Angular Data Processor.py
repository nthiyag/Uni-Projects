import math

data = open("lab_2_raw_data.txt", "r")

data.readline()
data.readline()

processed = open("lab_2_processed_data.txt", "w")
peaks = open("lab_2_processed_peaks.txt", "w")

t0 = 3.175911111111111

x_uncertainty = "0.06"
y_uncertainty = "0.02"

discard_lines = 0
for i in range(discard_lines):
    data.readline()

last_angle = 0
last_rising = False
last_anti_rising = False

last_peak_time = t0
last_peak = None

last_anti_peak_time = t0
last_anti_peak = None

while True:
    line = data.readline().split()
    if not line:
        break

    time = float(line[0]) - t0
    
    origin = (11.09, 349.6)

    rest = (24.39, -153.3)
    vertical = (rest[0]-origin[0], rest[1]-origin[1])

    point = (float(line[1]), float(line[2]))
    string = (point[0]-origin[0], point[1]-origin[1])

    vertical_mag = math.sqrt(vertical[0]**2 + vertical[1]**2)
    string_mag = math.sqrt(string[0]**2 + string[1]**2)

    vertical_dot_string = vertical[0]*string[0] + vertical[1]*string[1]

    angle = math.acos(vertical_dot_string / (vertical_mag * string_mag))
    if math.isnan(angle):
        print(string, point, point[0], origin[0], line)
    angle = -angle if point[0] < origin[0]+(point[1]-origin[1])/(vertical[1]/vertical[0]) else angle

    processed.write(str(time) +" " +str(angle) +" " +str(x_uncertainty) +" " +str(y_uncertainty) +"\n")

    rising = angle > last_angle
    anti_rising = angle < last_angle

    if angle == 0 or last_angle == 0:
        peak_reached = False
    elif (angle / abs(angle)) != (last_angle / abs(last_angle)):
        peak_reached = False

    if (last_rising and not rising and angle > 0) and not peak_reached:
        peak_reached = True
        if last_peak != None:
            peak_elapsed_time = time - last_peak_time
            peaks.write(str(last_peak) +"  " +str(peak_elapsed_time)  +" " +str(x_uncertainty) +" " +str(y_uncertainty) +"\n")
        last_peak_time = time
        last_peak = angle
    
    if (last_anti_rising and not anti_rising and angle < 0) and not peak_reached:
        peak_reached = True
        if last_anti_peak != None:
            anti_peak_elapsed_time = time - last_anti_peak_time
            peaks.write(str(last_anti_peak) +"  " +str(anti_peak_elapsed_time)  +" " +str(x_uncertainty) +" " +str(y_uncertainty) +"\n")
        last_anti_peak_time = time
        last_anti_peak = angle

    last_angle = angle
    last_rising = rising
    last_anti_rising = anti_rising

data.close()
processed.close()
peaks.close()

print("done")
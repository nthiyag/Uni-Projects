import math

data = open("raw_data.txt", "r")

data.readline()
data.readline()

processed = open("processed_data.txt", "w")
peaks = open("processed_peaks.txt", "w")

t0 = 10.49330599

x_uncertainty = "0.02"
y_uncertainty = "0.06"

discard_lines = 16
for i in range(discard_lines):
    data.readline()

last_amplitude = 0
last_rising = False

while True:
    line = data.readline().split()
    if not line:
        break

    time = float(line[0]) - t0
    
    origin = (156.4, 273.0)

    rest = (168.7, -371.1)
    vertical = (rest[0]-origin[0], rest[1]-origin[1])

    point = (float(line[1]), float(line[2]))
    string = (point[0]-origin[0], point[1]-origin[1])

    vertical_mag = math.sqrt(vertical[0]**2 + vertical[1]**2)
    string_mag = math.sqrt(string[0]**2 + string[1]**2)

    vertical_dot_string = vertical[0]*string[0] + vertical[1]*string[1]

    amplitude = math.acos(vertical_dot_string / (vertical_mag * string_mag))
    if math.isnan(amplitude):
        print(string, point, point[0], origin[0], line)
    amplitude = -amplitude if point[0] < origin[0]+(point[1]-origin[1])/(vertical[1]/vertical[0]) else amplitude

    processed.write(str(time) +" " +str(amplitude) +" " +str(x_uncertainty) +" " +str(y_uncertainty) +"\n")

    rising = amplitude > last_amplitude

    if last_rising and not rising and amplitude > 0:
        peaks.write(str(time) +"  " +str(amplitude)  +" " +str(x_uncertainty) +" " +str(y_uncertainty) +"\n")

    last_amplitude = amplitude
    last_rising = rising

data.close()
processed.close()
peaks.close()
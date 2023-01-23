def open_file(filename):
    file= open(filename)
    line = file.readline()
    lines = []
    while line:
        lines.append(float(line))
        line = file.readline()
    
    return lines

elements = {"argon": None, "helium": None, "krypton": None, "mercury": None, "neon": None, "nitrogen": None, "xenon": None}

measured = open_file("measured.txt")

for ele in elements:
    baseline = open_file(ele + ".txt")

    error = 0
    for m in measured:
        lowest_sep = (m-baseline[0])**2
        for b in baseline:
            lowest_sep = (m-b)**2 if (m-b)**2 < lowest_sep else lowest_sep
        lowest_sep = lowest_sep if lowest_sep < 100 else 0
        error += lowest_sep
    
    elements[ele] = error

print(elements)
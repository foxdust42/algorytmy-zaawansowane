import numpy as np
import matplotlib.pyplot as plt
import csv

x = []
y = []
z1 = []
z2 = []

w = 1

with open('res2.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')
    
    for row in plots:
        x.append(row[0])
        y.append(int(row[1]))
        z1.append(2 * w * w * w)
        z2.append(5 * w * w * w)
        w = w + 1
    

w, xt = 0.3, np.arange(len(x))

xt = np.arange(len(x))

cl_w = 0.7
br_w = cl_w/3

fig, ax = plt.subplots()

ax.plot(x, z1, label = "2 * n^3")
ax.plot(x, y , label = "czas obliczeń", marker="o")
ax.plot(x, z2, label = "5 * n^3")

# ax.bar(xt + (br_w * 0) - cl_w/2, z1, width=br_w, label="2 * n^3")
# ax.bar(xt + (br_w * 1) - cl_w/2, y , width=br_w, label="czas obliczeń")
# ax.bar(xt + (br_w * 2) - cl_w/2, z2, width=br_w, label="5 * n^3")

ax.set_xticks(xt)
ax.set_xticklabels(x)
ax.set_ylabel('milisekundy')
ax.set_xlabel('ilość studni')
ax.set_title('czas obliczenia 1000 iteracji w porównaniu do n^3')
ax.legend()
plt.show()
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


file = 'thresh_scan'
file = 'tt'
np.set_printoptions(suppress=True, precision=6, linewidth=1500)

data_path = file+'.csv'
with open(data_path, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader)
    #print(list(reader))
    #data = np.array(list(reader)).astype(float)
    data = list(reader)
    #print(data)



energy_var = {}
energy_pt2 = {}
num_thresh = 6
for i in range(len(data)):
    if  i > 0:
        energy_var[i] = np.array([float(a) for a in data[i][1:6]])
        energy_pt2[i] = np.array([float(a) for a in data[i][6:12]])


print(energy_var)
print(energy_pt2)

cb1 = 'k'
cb2 = 'grey'
cb3 = 'tab:blue'
cb4 = 'tab:red'
cb5 = 'tab:orange'
cb6 = 'tab:orange'


cb = ['0','k','tab:blue','tab:purple','tab:red','tab:orange','tab:green','tab:cyan','k']
label = ['0','$S_0$','$T_1$','$T_2$','$T_3$','$T_4$','$T_5$','$T_6$']

#full spectra
cb = ['0', "black", "darkgrey", "rosybrown", "indianred", "lightsalmon", "salmon", "darkorange", "tan", "darkkhaki", "olive", "darkseagreen", "mediumseagreen", "green", "mediumaquamarine", "darkcyan", "lightblue", "cadetblue", "navy", "royalblue", "blue", "mediumslateblue", "cornflowerblue", "blueviolet", "mediumslateblue", "plum", "purple", "orchid", "hotpink", "palevioletred", "crimson", "lightpink"]

#greens, blue, purple
cb = ['black', 'k', 'forestgreen', 'limegreen', 'green', 'springgreen', 'mediumspringgreen', 'mediumaquamarine', 'aquamarine', 'turquoise', 'lightseagreen', 'darkslategrey', 'teal', 'c', 'aqua', 'cadetblue', 'powderblue', 'lightblue', 'deepskyblue', 'skyblue', 'lightskyblue', 'steelblue', 'dodgerblue', 'slategrey', 'lightsteelblue', 'cornflowerblue', 'midnightblue', 'slateblue', 'darkslateblue', 'mediumslateblue', 'mediumpurple', 'rebeccapurple']

fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
extrap = []


for key in energy_var:
    #print(energy_var[key],energy_pt2[key])

    x = energy_var[key] - energy_pt2[key]
    z = energy_pt2[key]
    y = energy_var[key]
    m, b = np.polyfit(x, y, 1)
    print(x,z)
    m2,b = np.polyfit(x, z, 1)
    #m, c, b = np.polyfit(x, energy_var[key], 2)
    print(b)
    extrap.append(b)


    plt.rcParams.update({'font.size': 10})


    #ax.plot(x, y, marker='.',linestyle='-' ,markersize=10)
    ax.plot(x, y, marker='.',linestyle='-' ,markersize=10,color = cb[key]     )
    #ax.plot(x, y, marker='.',linestyle='-' ,label=label[key] ,markersize=10,color = cb[key]     )
    ax.plot(x, z, marker='.',linestyle=' ' ,markersize=10,color = cb[key]     )
    #ax.plot(x, z, marker='.',linestyle=' ' ,markersize=10)
    #if key ==3:
    #    ax.plot(x, y, marker='x',linestyle=' ' ,label=headers[1] ,markersize=8,color = cb3     )
    #    ax.plot(x, z, marker='x',linestyle=' ' ,label=headers[1] ,markersize=8,color = cb4     )
    #else:
    #    ax.plot(x, y, marker='.',linestyle=' ' ,label=headers[1] ,markersize=10,color = cb3     )
    #    ax.plot(x, z, marker='.',linestyle=' ' ,label=headers[1] ,markersize=10,color = cb4     )

    x2 = np.array([0,1])*0.02
    line = m*x2+b
    #ax.plot(x2, line, 'r',alpha=1.0 ,linestyle='-')
    ax.plot(x2, line, 'r',alpha=1.0 ,color =cb[key],linestyle='-')
    line = m2*x2+b
    #ax.plot(x2, line, 'r',alpha=0.5 ,linestyle='--')
    ax.plot(x2, line, 'r',alpha=0.5 ,color =cb[key],linestyle='--')
    print("Extrap  %14.8f"%b)
    print("Var root",key,y)
    print("PT  root",key,z)
    print("DIFFF",x)

#full
ax.set_ylim(-2754.775, -2754.61)
#middle
#ax.set_ylim(-2754.705, -2754.696)
#top
#ax.set_ylim(-2754.640, -2754.62)
#ax.set_ylim(-90.445,-90.425)
ax.set_xlim(0,0.005)
#ax.set_ylim(-90.51,-90.50)
#ax.set_xlim(0,0.005)

ax.set_xlabel('$\Delta$E$_2$ (E$_h$) ')
ax.set_ylabel('Energy (E$_h$) ')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)

ax.legend()
ax.set_title("Tetracence Tetramere \n (40o, 40e) \n ε = {0.0005, 0.0007, 0.001, 0.005, 0.01} \n 31 roots")
if 'shci' in file:
    plt.gca().invert_xaxis()
else:
    ax.yaxis.tick_left()
    #ax.yaxis.tick_right()

fig = plt.gcf()
fig.set_size_inches(9.,9.5)
#fig.set_size_inches(2.,4.5)
fig.savefig(file+'full_ext.pdf', dpi=100)
plt.show()

#print(extrap)

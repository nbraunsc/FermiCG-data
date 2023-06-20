import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib import colors as mcolors


#file = 'thresh_scan'
#file = 'tt'
file = 'tt_root_tracking'

#file = 'tt_tucker'
#file = 'tucker_data_clip/tucker_clip'
file = 'tucker_data_clip/tucker_clip_ev'
#file = 'notuck_clip'
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
        energy_var[i] = np.array([float(a) for a in data[i][1:4]])
        energy_pt2[i] = np.array([float(a) for a in data[i][4:7]])


print(energy_var)
print(energy_pt2)

blue = '#3E6D9C'
orange = '#FD841F'
red_orange = '#E14D2A'
dark_blue = '#001253'

cb = ['0', dark_blue, blue, blue, blue, blue, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, red_orange, orange, orange, blue, blue, orange, blue, blue, orange]

fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
extrap = []


for key in energy_var:
    x = energy_var[key] - energy_pt2[key]
    z = energy_pt2[key]
    y = energy_var[key]
    # use the ones below if you want to shift by extrapolated value
    #z = energy_pt2[key] + 74961.16852588826
    #y = energy_var[key] + 74961.16852588826
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
    ax.plot(x, z, marker='.',linestyle=' ' , markersize=10,color = cb[key]     )
    #ax.plot(x, z, marker='.',linestyle=' ' ,markersize=10)
    #if key ==3:
    #    ax.plot(x, y, marker='x',linestyle=' ' ,label=headers[1] ,markersize=8,color = cb3     )
    #    ax.plot(x, z, marker='x',linestyle=' ' ,label=headers[1] ,markersize=8,color = cb4     )
    #else:
    #    ax.plot(x, y, marker='.',linestyle=' ' ,label=headers[1] ,markersize=10,color = cb3     )
    #    ax.plot(x, z, marker='.',linestyle=' ' ,label=headers[1] ,markersize=10,color = cb4     )

    x2 = np.array([0,1])*0.9
    line = m*x2+b
    #ax.plot(x2, line, 'r',alpha=1.0 ,linestyle='-')
    ax.plot(x2, line, 'r',alpha=1.0 ,color =cb[key],linestyle='-', linewidth=1.5)
    line = m2*x2+b
    #ax.plot(x2, line, 'r',alpha=0.5 ,linestyle='--')
    ax.plot(x2, line, 'r',alpha=0.5 ,color =cb[key],linestyle='--', linewidth=1.5)
    print("Extrap  %14.8f"%b)
    print("Var root",key,y)
    print("PT  root",key,z)
    print("DIFFF",x)

#full
#ax.set_ylim(-2754.775, -2754.61)
#middle
#ax.set_ylim(-2754.702, -2754.685)
#ax.set_ylim(-2754.709, -2754.693)

#new middle
#ax.set_ylim(-2754.706, -2754.698)
#biexcitons
#ax.set_ylim(-2754.640, -2754.632)


#top
#ax.set_ylim(-2754.640, -2754.62)
#ax.set_ylim(-2754.639, -2754.623)
#new top
#ax.set_ylim(-2754.632, -2754.624)

#ax.set_ylim(-0.5,4.5)
ax.set_ylim(1.7,2.1)
#ax.set_ylim(1.75,2)
#ax.set_ylim(3.6,4)
#ax.set_ylim(-90.51,-90.50)
ax.set_xlim(0,0.04)

#ax.set_xlabel('$\Delta$E$_2$ (E$_h$) ')
#ax.set_ylabel('Energy (E$_h$) ')
ax.set_xlabel('$\Delta$E$_{PT2}$ (eV) ')
ax.set_ylabel('E – E$_{\infty}$ (eV) ')
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
ax.set_yticklabels([])

#ax.legend()
#ax.set_title("Tetracence Tetramere \n (40o, 40e) \n ε = {0.0005, 0.0007, 0.001, 0.005, 0.01} \n 31 roots")
#ax.set_title("(a) All 31 States")
ax.set_title("(b) 4 Lowest Triplets")
#ax.set_title("TPSCI with HOSVD")
#ax.set_title("TPSCI")
#ax.set_title("(b) TPSCI with HOSVD")
#ax.set_title("(b) TPSCI with HOSVD")
#ax.set_title("(c) TPSCI")
#ax.set_title("(c) TPSCI")
#ax.set_title("(c) 18 Biexcitons &\n Higher Excited States")
#ax.set_title("(d) 4 Bright States &\n 4 Triplets")
#ax.set_title("(a) NO HOSVD, Full")
#ax.set_title("(b) NO HOSVD, Middle")
#ax.set_title("(c) NO HOSVD, Top")
#ax.set_title("(b) Middle Triplet States")
#ax.set_title("(a) All 31 States")

black_patch = mpatches.Patch(color=dark_blue, label='Ground')
blue_patch = mpatches.Patch(color=blue, label='Triplets')
green_patch = mpatches.Patch(color=orange, label='Singlets')
red_patch = mpatches.Patch(color=red_orange, label='Biexcitons')
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='center right')
#ax.legend(handles=[black_patch, blue_patch, green_patch, red_patch], loc='upper center')
#ax.legend()

if 'shci' in file:
    plt.gca().invert_xaxis()
else:
    #plt.gca().invert_xaxis()
    ax.yaxis.tick_left()
    #ax.yaxis.tick_right()
#ax.set_facecolor('gainsboro')
fig = plt.gcf()
#fig.set_size_inches(9.,9.5)
fig.set_size_inches(2.,4.5)
#fig.savefig(file+'nohosvd_top_ext.pdf', dpi=100)
#fig.savefig(file+'nohosvd_middle_ext.pdf', dpi=100)
#fig.savefig(file+'nohosvd_full_ext.pdf', dpi=100)
#fig.savefig(file+'_colors_large_top_ext.pdf', dpi=100)
#fig.savefig(file+'_colors_legend_ext.pdf', dpi=100)

#fig.savefig(file+'_colors_biexc_ext.pdf', dpi=100)
#fig.savefig(file+'_colors_large_middle_ext.pdf', dpi=100)
#fig.savefig(file+'_colors_full_ext.pdf', dpi=100)
#fig.savefig(file+'_full_ext.pdf', dpi=100)
fig.savefig(file+'_yaxis_ext.pdf', dpi=100)
#fig.savefig(file+'_middle_ext.pdf', dpi=100)
#fig.savefig(file+'_top_ext.pdf', dpi=100)

plt.show()

print(extrap)

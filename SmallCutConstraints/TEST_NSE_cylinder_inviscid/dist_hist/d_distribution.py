
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import cm


## get data
file_name = "distances_r-1e-10"
h_average = 0.001
distances = []


with open(file_name+".csv", 'r', newline ='') as file:
    reader = csv.DictReader(file)
    for row in reader:
        distances.append(float(row['DISTANCE']))
n_cut_nodes = len(distances)
distances = [d/h_average for d in distances]

## distance modification
#distmod = 1e-7
#distances = [d if abs(d) > 1e-7 else -distmod for d in distances]

## set color
viridis = cm.get_cmap('viridis', 8)
c = viridis(0.25)

## use LaTex
plt.rc('text', usetex=True)
latex_params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(latex_params)

## format plots
math_label_fontsize = 18
ticks_fontsize = 14

## plot distribution
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
# An "interface" to matplotlib.axes.Axes.hist() method
n, bins, patches = ax1.hist(x=distances, bins=40, weights=[100.0/n_cut_nodes for d in distances], color=c, alpha=1.0, rwidth=0.9)

plt.xlabel(r'$d_\text{cut} / h_\text{avg}$ [-]', fontsize=math_label_fontsize)
plt.ylabel(r'Frequency [$\%$]', fontsize=math_label_fontsize)
# Set a clean upper y-axis limit.
#maxfreq = n.max()
#plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.ylim(ymax=np.ceil(n.max()))
ax1.tick_params(axis="both", which="both", direction="in")
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

# save plot and clear
#plt.savefig("d_distribution.png",format='png', bbox_inches='tight', dpi=300)
plt.tight_layout()
plt.savefig(file_name+"_hist.pdf")
#plt.clf()

## plot distribution using log scale (positive and negative distances)
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
min_log = 11
logbins_pos = [1e-20,5*10**(-min_log)]
logbins_pos.extend((5*np.logspace(-min_log+1,0,min_log+1)).tolist())
logbins = ((-1)*np.array(logbins_pos)).tolist()
logbins.reverse()
logbins.extend(logbins_pos)
n, bins, patches = ax2.hist(x=distances, bins=logbins, weights=[100.0/n_cut_nodes for d in distances], log=False, color=c, alpha=1.0, rwidth=0.9)

plt.xlabel(r'$d_\text{cut} / h_\text{avg}$ [-]', fontsize=math_label_fontsize)
plt.ylabel(r'Frequency [$\%$]', fontsize=math_label_fontsize)
plt.xscale('symlog', linthreshx=1e-12)
plt.xlim(-9.0,9.0)
x_ticks_pos = np.logspace(-min_log,0,min_log+1)
x_ticks = ((-1)*np.array(x_ticks_pos)).tolist()
x_ticks.reverse()
x_ticks.extend(x_ticks_pos)
ax2.set_xticks(x_ticks)
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 5) * 5 if maxfreq % 5 else maxfreq + 5)
ax2.tick_params(axis="both", which="both", direction="in")
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

# save plot
fig2.set_figwidth(20)
#plt.tight_layout()
plt.savefig(file_name+"_hist_log.pdf")

## plot distribution using log scale (only positive distances)
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
min_log = 11
distances = np.array(distances)
distances = np.delete(distances, np.where(distances<0.0))
logbins_pos = [1e-20,5*10**(-min_log)]
logbins_pos.extend((5*np.logspace(-min_log+1,0,min_log+1)).tolist())
n, bins, patches = ax3.hist(x=distances, bins=logbins_pos, weights=[100.0/n_cut_nodes for d in distances], log=False, color=c, alpha=1.0, rwidth=0.9)

plt.xlabel(r'$d_\text{cut} / h_\text{avg}$ [-]', fontsize=math_label_fontsize)
plt.ylabel(r'Frequency [$\%$]', fontsize=math_label_fontsize)
plt.xscale('symlog', linthreshx=1e-12)
plt.xlim(2*10**(-min_log-1),9.0)
x_ticks_pos = np.logspace(-min_log,0,min_log+1)
ax3.set_xticks(x_ticks_pos)
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 5) * 5 if maxfreq % 5 else maxfreq + 5)
ax3.tick_params(axis="both", which="both", direction="in")
plt.xticks(fontsize=ticks_fontsize)
plt.yticks(fontsize=ticks_fontsize)

# save plot
fig3.set_figwidth(10)
#plt.tight_layout()
plt.savefig(file_name+"_hist_log_pos.pdf")
plt.clf()

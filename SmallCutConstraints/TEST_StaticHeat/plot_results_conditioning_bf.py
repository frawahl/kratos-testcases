## Results post-process
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


# CONDITIONING without ref5
h_vect = [0.2, 0.1, 0.05, 0.025, 0.0125]
cond_bf = [7.2665768599018055, 30.98363686598401, 125.8778758192609, 505.46120447576783, 2023.7961045026027]


## use LaTex
matplotlib.rc('text', usetex=True)
latex_params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(latex_params)

## format plots
colors_line=['-k','-k','s--','x--','v--','+--','o--'] # Mixed lines
msize=[1,1,5,7,5,7,5]
locations_list = ['lower right','lower right','upper right','upper right', 'best']
linewidth = 1.0
legend_fontsize = 11
math_label_fontsize = 18

# get color map
viridis = cm.get_cmap('viridis', 8)
# Indices to step through colormap
cmap_x = np.linspace(0.0, 1.0, 5)

# Plot data
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(h_vect, cond_bf, colors_line[6], c='black', linewidth=linewidth, markersize=msize[6], alpha=0.9, label=r'Bodyfitted')

# Set legend
plt.legend(numpoints=1, markerscale=1, loc='right', fontsize=legend_fontsize, frameon=False)

# Set plot limits
#plt.xlim(-0.11, 0.11)
plt.ylim(5, 5e17)

# Set labels
x_label = r'$h$'
y_label = r'$(\lambda_{max}/\lambda_{min})_{\Tilde{\Omega}}$'
plt.xlabel(x_label, fontsize = math_label_fontsize)
plt.ylabel(y_label, fontsize = math_label_fontsize)
ax.set_xticks(h_vect)
ax.set_xticklabels(h_vect)
ax.tick_params(which='minor', length=0)

# invert x axis
plt.gca().invert_xaxis()

# Set axis ticks
# x_ticks = []
# for i in range(0,5):
#     x_ticks.append(i*0.5*math.pi)
# x_labels = [r"$0$",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"]
# ax.set_xticks(x_ticks)
# ax.set_xticklabels(x_labels,fontsize=14)

plt.yticks(fontsize = 14)

# Set brackground grid
# plt.grid()

# Save plot and clear
plt.savefig("conditioning_bf_limits.png",format='png', bbox_inches='tight', dpi=300)
plt.clf()

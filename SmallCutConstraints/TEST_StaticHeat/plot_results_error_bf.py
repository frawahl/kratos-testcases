## Results post-process
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

## Convergence results
h_vect = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625] #(old mesh 0: for y=1.5 h=0.1875)
#y_total = [1.6, 1.5, 1.5, 1.5, 1.5, 1.5]
#divisions = [8, 15, 30, 60, 120, 240]

# ERROR NORM AT GAUSS POINTS  
err_norm_bf = [6.088327493078292e-17, 3.297454721554084e-16, 1.0176142735678096e-15, 5.569505164642452e-10, 7.4360811551859516e-09, 8.350682188287363e-09] 
err_norm_bf_relativeError = [1.24907902183471e-16, 9.750161817133875e-16, 2.7195718574455487e-15, 2.5707380857327813e-09, 2.980494948537671e-08, 2.811349016293396e-08]

posCut = False
negCut = False

## Calculate convergence Parameters
err_linear = []
err_quadratic = []
if posCut:
    err_0_linear = 0.085
    err_0_quadratic = 0.16
elif negCut:
    err_0_linear = 1.0
    err_0_quadratic = 0.2	
else:
    err_0_linear = 0.5
    err_0_quadratic = 0.13
h_vect_err = [0.025, 0.0125, 0.00625]	
for h in h_vect_err:
    err_linear.append(err_0_linear*h)
    err_quadratic.append(err_0_quadratic*h**2)

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
#ax.loglog(h_vect, flux_err_norm, colors_line[3], linewidth=linewidth, markersize=4, alpha=0.9, label='')
ax.loglog(h_vect_err, err_linear, colors_line[0], linewidth=linewidth*0.5, markersize=msize[0], alpha=0.5, label='')
ax.loglog(h_vect_err, err_quadratic, colors_line[1], linewidth=linewidth*0.5, markersize=msize[1], alpha=0.5, label='')
ax.loglog(h_vect, err_norm_bf, colors_line[6], c='black', linewidth=linewidth, markersize=msize[6], alpha=0.9, label=r'Bodyfitted')

# Set description of convergence rate
textstr_h = r'$h$'
textstr_h2 = r'$h^2$'
if posCut:
    ax.text(0.92, 0.65, textstr_h, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')
    ax.text(0.92, 0.185, textstr_h2, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')
elif negCut:
    ax.text(0.92, 0.75, textstr_h, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')
    ax.text(0.92, 0.18, textstr_h2, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')
else:
    ax.text(0.92, 0.96, textstr_h, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')
    ax.text(0.92, 0.78, textstr_h2, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')

# Set legend
plt.legend(numpoints=1, markerscale=1, loc='center left',fontsize=legend_fontsize,frameon=False)

# Set plot limits
#plt.xlim(-0.11, 0.11)
#plt.ylim(-4.5, 1.5)

# Set labels
x_label = r'$h$'
y_label = r'$\left\vert T - T_{ex} \right\vert_{\Tilde{\Omega}}$'
plt.xlabel(x_label, fontsize = math_label_fontsize)
plt.ylabel(y_label, fontsize = math_label_fontsize)
ax.set_xticks(h_vect)
ax.set_xticklabels([r'${}$'.format(h) for h in h_vect])
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

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
ax.tick_params(axis="both", which="both", direction="in")

# Set brackground grid
# plt.grid()

# Save plot and clear
plt.savefig("convergence_bf.png",format='png', bbox_inches='tight', dpi=300)
plt.clf()

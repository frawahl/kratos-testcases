## Results post-process
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

## Convergence results
h_vect = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625] #(old mesh 0: for y=1.5 h=0.1875)
#y_total = [1.6, 1.5, 1.5, 1.5, 1.5, 1.5]
#divisions = [8, 15, 30, 60, 120, 240]

# ERROR NORM AT GAUSS POINTS  
#### cut at y=0.75+1e-12, penalty=1e1, avoid zeros: 1e-12, DistMod: 1e-3 ####
posCut = True
negCut = True
if posCut:
    err_norm_none = [0.004022941583044258, 0.0009829970413894, 0.00023394505198513288, 5.7177714590514655e-05, 1.4112615439091186e-05, 3.5142318355672077e-06]
    err_norm_distmod = [0.00394520384113419, 0.0009267004547239497, 0.00022145717479965692, 5.3900170759006546e-05, 1.3012391157034536e-05, 3.0813031719710002e-06] 
    err_norm_mls = [0.004336860127170321, 0.0009029008769875675, 0.000217950136703944, 5.4828632016566934e-05, 1.3823943273743341e-05, 3.4697029300197285e-06]
    err_norm_local = [0.014835729666603039, 0.006277936688845564, 0.00288962420365917, 0.0013885711752331734, 0.0006812003348473048, 0.000337447470817805]
    err_norm_localabc = [0.013353755704182704, 0.005578327615275564, 0.0025497428218194157, 0.0012209282969972727, 0.0005978775185315119, 0.00029590145600283505]
#### cut at y=0.75-1e-12, penalty=1e1, avoid zeros: 1e-12, DistMod: 1e-3 ####
elif negCut:
    err_norm_none = [0.0014467533422534123, 0.0005735190153192788, 0.00018301097722400934, 5.197987726898036e-05, 1.3915522479411525e-05, 3.6162073239275877e-06]
    err_norm_distmod = [0.0014397366522313044, 0.0005693670114212557, 0.00018050379464443843, 5.0637412617332015e-05, 1.3267704876957333e-05, 3.3304009124515578e-06]
    err_norm_mls = [0.011199990191221676, 0.0019445495986157586, 0.0002464139447079262, 5.5496869441589774e-05, 1.7108540192916653e-05, 4.8986435630016125e-06]
    err_norm_local = [0.11614143774225856, 0.06424324374580052, 0.03357560695174651, 0.017168118886039162, 0.008684417541067148, 0.004368337469372038]
    err_norm_localabc = [0.0017404905452673963, 0.0006253515402847296, 0.00018680847422996376, 5.10964840819385e-05, 1.3370736719641071e-05, 3.4216644264392867e-06]
#### cut at y=0.8+h/2, penalty=1e1, avoid zeros: 1e-12, DistMod: 1e-3 ####
else:
    err_norm_none = [0.0012989611574178073, 0.000499020071959081, 0.00015119302477026062, 4.167789516958652e-05, 1.0985417103784238e-05, 2.8355537564438224e-06]
    err_norm_distmod = [0.0012989611574178073, 0.000499020071959081, 0.00015119302477026062, 4.167789516949135e-05, 1.0985417103796193e-05, 2.8355537564438224e-06]
    err_norm_mls = [0.005666442438199272, 0.0009783644059135972, 0.00017575214736955053, 4.206315585013436e-05, 1.1262832292843172e-05, 2.9723242615108283e-06]
    err_norm_local = [0.06370794899048574, 0.032057140347804804, 0.016104262742918415, 0.00808073138481119, 0.004049601963101194, 0.002027488027646042]
    err_norm_localabc = [0.0012764703774185216, 0.0004837579869205104, 0.00014533959513748194, 3.983994943571971e-05, 1.0455461050044096e-05, 2.685484999859093e-06]    

# CONDITIONING without ref5
#h_vect = [0.2, 0.1, 0.05, 0.025, 0.0125]
#### cut at y=0.75+1e-12, penalty=1e1, avoid zeros: 1e-12, DistMod: 1e-3 ####
# None: [1.424735820419105e+17, 1.270129405053078e+16, 5283687478499783.0, 5558906861393745.0, 1220878890859023.8]
# DistMod: [33248783.844061542, 4335846.385450855, 523180.9391737967, 57592.06069444402, 86753.16846462486]
# MLSConstraint: [14.218414023153452, 56.666212058646835, 225.0564577802069, 898.8068235765604, 3593.8862224395507]
# LocalConstraint: [14.842012680609422, 60.30058107670217, 241.4715267986745, 964.8448857684218, 3855.736630875529]
# LocalConstraintABC: [17.06540933201767, 68.82838584316848, 274.3920020213558, 1093.7378273840961, 4365.455325159425]
#### cut at y=0.75-1e-12, penalty=1e1, avoid zeros: 1e-12, DistMod: 1e-3 ####
# None: [12.885408891117036, 51.71224097466204, 207.52887946651916, 831.0573015624786, 3325.307804547582]
# DistMod: [12.739713151428775, 50.59650744150634, 198.98914919757, 766.4317857915853, 2852.1190185604974]
# MLSConstraint: [24.190755373315845, 121.39752747212853, 541.7682315040784, 2237.1343942292556, 9026.148633648156]
# LocalConstraint: [13.551188144656173, 50.56305450832397, 222.35381972106808, 927.2817714310014, 3781.4336670776956]
# LocalConstraintABC: [7.266576859893723, 30.98363686594063, 125.87787581907047, 505.46120447551647, 2023.796104522145]


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
ax.loglog(h_vect, err_norm_none, colors_line[2], c=viridis(cmap_x[0]), linewidth=linewidth, markersize=msize[2], alpha=0.9, label=r'None')
ax.loglog(h_vect, err_norm_distmod, colors_line[3], c=viridis(cmap_x[1]), linewidth=linewidth, markersize=msize[3], alpha=0.9, label=r'DistMod')
ax.loglog(h_vect, err_norm_mls, colors_line[4], c=viridis(cmap_x[2]), linewidth=linewidth, markersize=msize[4], alpha=0.9, label=r'MLSConstraints')
ax.loglog(h_vect, err_norm_local, colors_line[5], c=viridis(cmap_x[3]), linewidth=linewidth, markersize=msize[5], fillstyle='none', alpha=0.9, label=r'LocalConstraints')
ax.loglog(h_vect, err_norm_localabc, colors_line[6], c=viridis(cmap_x[4]), linewidth=linewidth, markersize=msize[6], fillstyle='none', alpha=0.9, label=r'LocalConstraintsBC')

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
    ax.text(0.92, 0.73, textstr_h, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')
    ax.text(0.92, 0.163, textstr_h2, transform=ax.transAxes, fontsize=11, color="gray", verticalalignment='top')

# Set legend
plt.legend(numpoints=1, markerscale=1, loc='lower left',fontsize=legend_fontsize,frameon=False)

# Set plot limits
#plt.xlim(-0.11, 0.11)
#plt.ylim(-4.5, 1.5)

# Set labels
x_label = r'$h$'
y_label = r'$\left\vert T - T_{ex} \right\vert_{\Tilde{\Omega}}$'
plt.xlabel(x_label, fontsize = math_label_fontsize)
plt.ylabel(y_label, fontsize = math_label_fontsize)

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
plt.savefig("rectangle_1_15_str_convergence_conforming_linear.png",format='png', bbox_inches='tight', dpi=300)
plt.clf()

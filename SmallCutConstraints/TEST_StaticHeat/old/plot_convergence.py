## Results post-process
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

## Convergence results
h_vect = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625]
err_norm = [0.004810336395373692, 0.0012726819896440428, 0.00019766740674027414, 5.233736690663701e-05, 1.2204588593605655e-05, 5.3280043415437274e-06]
#[0.00473648180088828, 0.0009055590782249239, 0.00023209566467160707, 4.467408943539313e-05, 1.4084216467379718e-05, 1.692989402714388e-05]
#[0.004296387844618439, 0.0015271910390128802, 0.00022361260895186293, 0.00035327014165438753, 1.2673877043561406e-05, 1.536674456807305e-05]  #[0.005660581419761063, 0.0017246469999579058, 0.00039595009755584795, 6.989874112912418e-05, 1.177054254974085e-05, 2.8782768950048156e-06]
flux_err_norm = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

## Calculate convergence Parameters
err_linear = []
err_quadratic = []
err_0_linear = 1.2
err_0_quadratic = 0.8
for h in h_vect:
    err_linear.append(err_0_linear*h)
    err_quadratic.append(err_0_quadratic*h**2)

## Plot the computed pressure coefficients solutions
matplotlib.rc('text', usetex=True) ## Makes matplotlib use LaTeX instead of the internal mathtext engine
latex_params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(latex_params)

colors_line=['x-k','--b','--r','x-g'] # Mixed lines
locations_list = ['lower right','lower right','upper right','upper right', 'best']
linewidth = 1.0
legend_fontsize = 14
math_label_fontsize = 18

# Plot data
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(h_vect, err_norm, colors_line[0], linewidth=linewidth, markersize=4, alpha=0.9, label='')
ax.loglog(h_vect, flux_err_norm, colors_line[3], linewidth=linewidth, markersize=4, alpha=0.9, label='')
ax.loglog(h_vect, err_linear, colors_line[1], linewidth=linewidth, markersize=4, alpha=0.9, label='Linear')
ax.loglog(h_vect, err_quadratic, colors_line[2], linewidth=linewidth, markersize=4, alpha=0.9, label='Quadratic')

# Set legend
plt.legend(numpoints=1, markerscale=2, loc='best',fontsize=legend_fontsize)

# Set plot limits
# plt.xlim(-0.11, 0.11)
# plt.ylim(-4.5, 1.5)

# Set labels
x_label = r'$h$'
y_label = r'$\left\lVert u - u_{ex} \right\rVert_{\Tilde{\Omega}}$'
plt.xlabel(x_label, fontsize = math_label_fontsize)
plt.ylabel(y_label, fontsize = math_label_fontsize)

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
plt.savefig("rectangle_1_15_str_convergence_conforming_linear.png",format='png', dpi=300)
plt.clf()

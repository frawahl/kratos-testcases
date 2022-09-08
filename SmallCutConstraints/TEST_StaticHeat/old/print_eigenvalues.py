from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
from scipy.linalg import eigvals

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

path = "comparison_smallCut_e-12_ref1-0.05/"
folder_names = ["split", "split+distMod", "split+mpc"]

h_vect = [0.2, 0.1, 0.05, 0.025, 0.0125]  #, 0.00625]

compare_ratios = []


for fn in folder_names:

    ratios = [] 

    for count, h in enumerate(h_vect):

        lhs_matrix = mmread(path + fn + "_ref" + str(count) + "/A.mm")
        rhs_matrix = mmread(path + fn + "_ref" + str(count) + "/b.mm.rhs")
        
        lhs_matrix = lhs_matrix.tocsr().toarray()
        
        vals = eigvals(lhs_matrix, b=None, overwrite_a=False, check_finite=True, homogeneous_eigvals=False)
        
        eig_max = max(abs(vals))
        eig_min = min(abs(vals)) 
        ratios.append(eig_max/eig_min)
    
    print("\nEigenvalueratios of LHS matrix of " + fn + ":")
    print(ratios)
    
    compare_ratios.append(ratios)
    
    
## Plot the computed pressure coefficients solutions
matplotlib.rc('text', usetex=True) ## Makes matplotlib use LaTeX instead of the internal mathtext engine
latex_params = {'text.latex.preamble' : [r'\usepackage{amsmath}']}
plt.rcParams.update(latex_params)

colors_line=['x-k','x-b','x-r','x-g'] # Mixed lines
locations_list = ['lower right','lower right','upper right','upper right', 'best']
linewidth = 1.0
legend_fontsize = 14
math_label_fontsize = 18

# Plot data
fig = plt.figure()
ax = fig.add_subplot(111)
for count, fn in enumerate(folder_names):
    ax.loglog(h_vect, compare_ratios[count], colors_line[count], linewidth=linewidth, markersize=4, alpha=0.9, label=fn)

# Set legend
plt.legend(numpoints=1, markerscale=2, loc='best',fontsize=legend_fontsize)

# Set labels
x_label = r'$h$'
y_label = r'$\lambda_{max} / \lambda_{min}$'
plt.xlabel(x_label, fontsize = math_label_fontsize)
plt.ylabel(y_label, fontsize = math_label_fontsize)

plt.yticks(fontsize = 14)

# Set brackground grid
# plt.grid()

# Save plot and clear
plt.savefig("rectangle_e-12_matrixConditioning.png",format='png', dpi=300)
plt.clf()

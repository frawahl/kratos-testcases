from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs
from scipy.linalg import eigvals

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

path = ""  #"comparison_smallCut_e-12_ref1-0.05/"



lhs_matrix = mmread(path + "A.mm")
rhs_matrix = mmread(path + "b.mm.rhs")

lhs_matrix = lhs_matrix.tocsr().toarray()

vals = eigvals(lhs_matrix, b=None, overwrite_a=False, check_finite=True, homogeneous_eigvals=False)

eig_max = max(abs(vals))
eig_min = min(abs(vals)) 

print("\nEigenvalueratios of LHS matrix :" + str(eig_max/eig_min))


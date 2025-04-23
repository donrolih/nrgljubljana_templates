from nrgljubljana_interface import Solver, Flat
from h5 import *
from triqs.utility import mpi
import numpy as np
# Parameters
D, V, U = 1.0, np.sqrt(2*0.1/np.pi), 0.3
e_f = -0.1
T = 1e-8

# Set up the Solver
templatedir = "/Users/drolih/sources/nrgljubljana_templates/new_templates"
S = Solver(templatedir = templatedir, model = "SIAM", symtype = "QS", mesh_max = 2.0, mesh_min = 1e-6, mesh_ratio = 1.01)
# Solve Parameters
sp = { "T": T, "Lambda": 2., "Nz": 4, "Tmin": 1e-9, "keep": 1000, "bandrescale": 1.0, "alpha": 0.6}

# Model Parameters
mp = { "U1": U, "eps1": e_f }
sp["model_parameters"] = mp

# Initialize hybridization function
S.Delta_w['imp'] << V**2 * Flat(D)

# nrgp = {"removefiles": False}
#nrgp["bandrescale"] = 5.0
# S.set_nrg_params(**nrgp)

# Solve the impurity model
S.solve(**sp)

# Store the Result
if mpi.is_master_node():
    with HDFArchive("aim_solution.h5", 'w') as arch:
        arch["A_w"] = S.A_w
        arch["G_w"] = S.G_w
        arch["F_l_w"] = S.F_l_w
        arch["F_r_w"] = S.F_r_w
        arch["I_w"] = S.I_w
        arch["SigmaHartree_w"] = S.SigmaHartree_w
        arch["Sigma_w"] = S.Sigma_w
        arch["expv"] = S.expv
        arch["V"] = V
        arch["sp"] = sp
        arch["mp"] = mp
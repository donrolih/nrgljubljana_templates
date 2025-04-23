from nrgljubljana_interface import Solver, SemiCircular
from h5 import *
from triqs.utility import mpi

# Construct Parameters
cp = {}
cp["templatedir"] = "/Users/drolih/sources/nrgljubljana_templates/new_templates"
cp["model"] = "2orb-UJ"
cp["symtype"] = "QS"
cp["mesh_max"] = 2.0
cp["mesh_min"] = 1e-2
cp["mesh_ratio"] = 1.01

# Set up the Solver
S = Solver(**cp)

# Solve Parameters
sp = {}
sp["T"] = 1e-1
sp["Lambda"] = 4.0
sp["Nz"] = 2
sp["Tmin"] = 0.5
sp["keep"] = 2000

# Model Parameters
mp = {}
mp["U1"] = 1.0
mp["U2"] = 0.9
mp["eps1"] = -0.5
mp["eps2"] = -0.4
mp["U12"] = 0.1
mp["J12"] = 0.05
sp["model_parameters"] = mp

# Low-energy NRG parameters
np = {}
np["bins"] = 50
S.set_nrg_params(**np)

# # Initialize hybridization function
S.Delta_w['imp'][0,0] << 0.5 * SemiCircular(1.0)
S.Delta_w['imp'][1,1] << 0.4 * SemiCircular(1.0)
# Out-of-diagonal Delta is zero

# Solve the impurity model
S.solve(**sp)

# # Store the Result
if mpi.is_master_node():
    with HDFArchive("2orb-UJ_QS.h5", 'w') as arch:
        arch["A_w"] = S.A_w
        arch["G_w"] = S.G_w
        arch["F_l_w"] = S.F_l_w
        arch["F_r_w"] = S.F_r_w
        arch["I_w"] = S.I_w
        arch["Sigma_w"] = S.Sigma_w
        arch["SigmaHartree_w"] = S.SigmaHartree_w
        arch["expv"] = S.expv
        arch["sp"] = sp

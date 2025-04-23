# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
from triqs.gf import *
from nrgljubljana_interface import MeshReFreqPts
plt.style.use('don-custom')
%config InlineBackend.figure_format = 'svg'

import numpy as np
from h5 import HDFArchive

with HDFArchive('aim_solution.h5', 'r') as A:
    A_w = A['A_w']['imp']
    G_w = A['G_w']['imp']
    F_l_w = A['F_l_w']['imp']
    F_r_w = A['F_r_w']['imp']
    I_w = A['I_w']['imp']
    SigmaHartree_w = A['SigmaHartree_w']['imp']
    Sigma_w = A['Sigma_w']['imp']
    expv = A['expv']
    sp = A['sp']
    mp = A['mp']
    Gamma = (np.pi * 0.5 * A['V'] ** 2)

w_mesh = np.array(list(A_w.mesh.values()))

Sigma_FG_w = F_l_w * inverse(G_w)

# %%
fig, ax = plt.subplots()

ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

ax.plot(w_mesh, -Sigma_FG_w.data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{FG}$')
ax.plot(w_mesh, -Sigma_w.data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{IFG}$')
ax.set_ylim(-1.5e-3, 1e-3)
ax.set_xlim(-0.001, 0.001)

ax.set_xlabel(r'$\omega / D$')
ax.set_ylabel(r'$-\mathrm{Im}\,\Sigma(\omega) / \Gamma$')

ax.set_title(f'$N_\mathrm{{ kp }} = { sp["keep"] }$, $U = { mp["U1"] }$, $\Gamma = { np.format_float_positional(Gamma, 2) }$, $\epsilon_d = { mp["eps1"] }$, $D = 1.$')

ax.legend()
fig.savefig("sigma_im_compare_low_freq.pdf")
# %%
fig, ax = plt.subplots()

ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
# ax.plot(w_mesh, Sigma_w.data[:, 0, 0].imag, label=r'old')
ax.plot(w_mesh, -Sigma_FG_w.data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{FG}$')
ax.plot(w_mesh, -Sigma_w.data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{IFG}$')
# ax.set_ylim(-1.5e-3, 1e-3)
ax.set_xlim(-1, 1)

ax.set_xlabel(r'$\omega / D$')
ax.set_ylabel(r'$-\mathrm{Im}\,\Sigma(\omega) / \Gamma$')

ax.set_title(f'$N_\mathrm{{ kp }} = { sp["keep"] }$, $U = { mp["U1"] }$, $\Gamma = { np.format_float_positional(Gamma, 2) }$, $\epsilon_d = { mp["eps1"] }$, $D = 1.$')

ax.legend()
fig.savefig("sigma_im_compare.pdf")
# %%
fig, ax = plt.subplots()

SigmaHartree = mp['U1'] * expv['n_d'] * 0.5
ax.plot(w_mesh, (Sigma_FG_w.data[:, 0].real - SigmaHartree_w.data[:, 0].real), label=r'$\Sigma^\mathrm{FG}$')
ax.plot(w_mesh, (Sigma_w.data[:, 0].real - SigmaHartree_w.data[:, 0].real), label=r'$\Sigma^\mathrm{IFG}$')

ax.set_xlim(-1, 1)

ax.legend()

ax.set_xlabel(r'$\omega / D$')
ax.set_ylabel(r'$(\mathrm{Re}\,\Sigma(\omega) - \Sigma^\mathrm{H}) / \Gamma$')
ax.set_title(f'$N_\mathrm{{ kp }} = { sp["keep"] }$, $U = { mp["U1"] }$, $\Gamma = { np.format_float_positional(Gamma, 2) }$, $\epsilon_d = { mp["eps1"] }$, $D = 1.$')

fig.savefig("sigma_re_compare.pdf")
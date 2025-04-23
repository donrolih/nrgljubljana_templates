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
    A_w = A['A_w']
    G_w = A['G_w']
    F_l_w = A['F_l_w']
    F_r_w = A['F_r_w']
    I_w = A['I_w']
    Sigma_w = A['Sigma_w']
    SigmaHartree_w = A['SigmaHartree_w']
    expv = A['expv']
    sp = A['sp']
    mp = A['mp']
    Gamma = (np.pi * 0.5 * A['V'] ** 2)

w_mesh = np.array(list(A_w.mesh.values()))

Sigma_FG_w = F_l_w * inverse(G_w)
# %%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax1.plot(w_mesh, -Sigma_FG_w['up'].data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{FG}_\mathrm{\uparrow}$')
ax1.plot(w_mesh, -Sigma_w['up'].data[:, 0].imag / Gamma, "-", label=r'$\Sigma^\mathrm{IFG}_\mathrm{\uparrow}$')
ax1.set_ylim(-1e-4, 1e-4)
ax1.set_xlim(-3e-3, 3e-3)
# ax1.set_xlabel(r'$\omega / D$')
ax1.set_ylabel(r'$-\mathrm{Im}\,\Sigma_\mathrm{\uparrow}(\omega) / \Gamma$')
ax1.legend()

ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax2.plot(w_mesh, -Sigma_FG_w['dn'].data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{FG}_\mathrm{\downarrow}$')
ax2.plot(w_mesh, -Sigma_w['dn'].data[:, 0].imag / Gamma, "-", label=r'$\Sigma^\mathrm{IFG}_\mathrm{\downarrow}$')
ax2.set_ylim(-1e-3, 2.5e-4)
ax2.set_xlim(-1e-3, 1e-3)
ax2.set_xlabel(r'$\omega / D$')
ax2.set_ylabel(r'$-\mathrm{Im}\,\Sigma_\mathrm{\downarrow}(\omega) / \Gamma$')
ax2.legend()

fig.tight_layout()

fig.suptitle(f'$N_\mathrm{{ kp }} = { sp["keep"] }$, $U = { mp["U1"] }$, $\Gamma = { np.format_float_positional(Gamma, 2) }$, $\epsilon_d = { mp["eps1"] }$, $B = {mp["B1"]}$, $D = 1.$', y=1.01)

fig.savefig("sigma_im_compare_low_freq.pdf")
# %%
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax1.plot(w_mesh, -Sigma_FG_w['up'].data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{FG}_\mathrm{\uparrow}$')
ax1.plot(w_mesh, -Sigma_w['up'].data[:, 0].imag/ Gamma, "-", label=r'$\Sigma^\mathrm{IFG}_\mathrm{\uparrow}$')
# ax1.set_ylim(-1e-4, 1e-4)
# ax1.set_xlim(-3e-3, 3e-3)
# ax1.set_xlabel(r'$\omega / D$')
ax1.set_ylabel(r'$-\mathrm{Im}\,\Sigma_\mathrm{\uparrow}(\omega) / \Gamma$')
ax1.legend()

ax2.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax2.plot(w_mesh, -Sigma_FG_w['dn'].data[:, 0].imag / Gamma, label=r'$\Sigma^\mathrm{FG}_\mathrm{\downarrow}$')
ax2.plot(w_mesh, -Sigma_w['dn'].data[:, 0].imag / Gamma, "-", label=r'$\Sigma^\mathrm{IFG}_\mathrm{\downarrow}$')
# ax2.set_ylim(-1e-3, 2.5e-4)
ax2.set_xlim(-1, 1)
ax2.set_xlabel(r'$\omega / D$')
ax2.set_ylabel(r'$-\mathrm{Im}\,\Sigma_\mathrm{\downarrow}(\omega) / \Gamma$')
ax2.legend()

fig.tight_layout()

fig.suptitle(f'$N_\mathrm{{ kp }} = { sp["keep"] }$, $U = { mp["U1"] }$, $\Gamma = { np.format_float_positional(Gamma, 2) }$, $\epsilon_d = { mp["eps1"] }$, $B = {mp["B1"]}$, $D = 1.$', y=1.01)

fig.savefig("sigma_im_compare.pdf")
# %%
fig, ax = plt.subplots()

ax.plot(w_mesh, (Sigma_FG_w['up'].data[:, 0].real - SigmaHartree_w['up'].data[:, 0].real), label=r'$\Sigma^\mathrm{FG}_\mathrm{\uparrow}$')
ax.plot(w_mesh, (Sigma_w['up'].data[:, 0].real - SigmaHartree_w['up'].data[:, 0].real), label=r'$\Sigma^\mathrm{IFG}_\mathrm{\uparrow}$')

ax.plot(w_mesh, (Sigma_FG_w['dn'].data[:, 0].real - SigmaHartree_w['dn'].data[:, 0].real), "--", label=r'$\Sigma^\mathrm{FG}_\mathrm{\downarrow}$')
ax.plot(w_mesh, (Sigma_w['dn'].data[:, 0].real - SigmaHartree_w['dn'].data[:, 0].real), "--", label=r'$\Sigma^\mathrm{IFG}_\mathrm{\downarrow}$')

ax.set_xlim(-1, 1)

ax.legend()

ax.set_xlabel(r'$\omega / D$')
ax.set_ylabel(r'$(\mathrm{Re}\,\Sigma_\sigma(\omega) - \Sigma^\mathrm{H}_\sigma) / \Gamma$')
ax.set_title(f'$N_\mathrm{{ kp }} = { sp["keep"] }$, $U = { mp["U1"] }$, $\Gamma = { np.format_float_positional(Gamma, 2) }$, $\epsilon_d = { mp["eps1"] }$, $D = 1.$')

fig.savefig("sigma_re_compare.pdf")
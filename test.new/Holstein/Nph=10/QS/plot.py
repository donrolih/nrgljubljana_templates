# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
from triqs.gf import *
from nrgljubljana_interface import MeshReFreqPts
plt.style.use('don-custom')
%config InlineBackend.figure_format = 'svg'

import numpy as np
from h5 import HDFArchive

with HDFArchive('holstein_solution.h5', 'r') as A:
    A_w = A['A_w']['imp']
    G_w = A['G_w']['imp']
    F_l_w = A['F_l_w']['imp']
    F_r_w = A['F_r_w']['imp']
    I_w = A['I_w']['imp']
    Sigma_w = A['Sigma_w']['imp']
    SigmaHartree_w = A['SigmaHartree_w']['imp']
    expv = A['expv']

w_mesh = np.array(list(A_w.mesh.values()))

Sigma_FG_w = F_l_w * inverse(G_w)
# %%
fig, ax = plt.subplots(nrows=2, sharex=True)

ax[0].plot(w_mesh, Sigma_FG_w.data[:, 0, 0].real, label=r'$\Sigma^\mathrm{FG}$')
ax[0].plot(w_mesh, Sigma_w.data[:, 0, 0].real, label=r'$\Sigma^\mathrm{IFG}$')
ax[0].set_ylabel(r'$\mathrm{Re}\Sigma(\omega)$')
ax[0].set_xlim(-0.5, 0.5)


ax[1].plot(w_mesh, -Sigma_FG_w.data[:, 0, 0].imag, label=r'$\Sigma^\mathrm{FG}$')
ax[1].plot(w_mesh, -Sigma_w.data[:, 0, 0].imag, label=r'$\Sigma^\mathrm{IFG}$')
ax[1].set_ylabel(r'$-\mathrm{Im}\Sigma(\omega)$')
ax[1].set_xlabel(r'$\omega$')
ax[1].legend()
fig.suptitle(r'$N_\mathrm{keep} = 1000$, $N_\mathrm{ph}=10$, $T=10^{-5}$, $\omega=0.2$, $g=0.2$, $n=1$, $U = 0.$', y=0.95)
fig.tight_layout()
fig.savefig("sigma_re_im_compare.pdf")
# %%
fig, ax = plt.subplots()

ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

ax.plot(w_mesh, -Sigma_FG_w.data[:, 0, 0].imag, label=r'$\Sigma^\mathrm{FG}$')
ax.plot(w_mesh, -Sigma_w.data[:, 0, 0].imag, label=r'$\Sigma^\mathrm{IFG}$')

ax.set_xlim(-1e-4, 1e-4)
ax.set_ylim(-1e-3, 1e-3)

ax.legend()

ax.set_ylabel(r'$-\mathrm{Im}\Sigma(\omega)$')

ax.set_xlabel(r'$\omega$')

ax.set_title(r'$N_\mathrm{keep} = 1000$, $N_\mathrm{ph}=10$, $T=10^{-5}$, $\omega=0.2$, $g=0.2$, $n=1$, $U = 0.$', y=1.05)

fig.savefig("sigma_im_compare_low_freq.pdf")
# %%

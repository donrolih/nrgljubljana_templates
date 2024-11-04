# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
from triqs.gf import *
from nrgljubljana_interface import MeshReFreqPts
plt.style.use('don-custom')
%config InlineBackend.figure_format = 'svg'

import numpy as np
from h5 import HDFArchive

with HDFArchive('2orb-UJ_QS.h5', 'r') as A:
    A_w = A['A_w']['imp']
    G_w = A['G_w']['imp']
    F_l_w = A['F_l_w']['imp']
    F_r_w = A['F_r_w']['imp']
    I_w = A['I_w']['imp']
    SigmaHartree_w = A['SigmaHartree_w']['imp']
    Sigma_w = A['Sigma_w']['imp']
    expv = A['expv']

w_mesh = np.array(list(A_w.mesh.values()))

Sigma_FG_w = F_l_w * inverse(G_w)
# %%
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(w_mesh, -Sigma_FG_w.data[:, 0, 0].imag, "-", label=r'FG')
ax1.plot(w_mesh, -Sigma_w.data[:, 0, 0].imag, "-", label=r'IFG')
ax1.legend()

ax1.set_ylabel(r'$-\mathrm{Im}\,\Sigma_{11}(\omega)$')

ax2.plot(w_mesh, -Sigma_FG_w.data[:, 1, 1].imag, "-", label=r'FG')
ax2.plot(w_mesh, -Sigma_w.data[:, 1, 1].imag, "-", label=r'IFG')
ax2.legend()

ax2.set_ylabel(r'$-\mathrm{Im}\,\Sigma_{22}(\omega)$')
ax2.set_xlabel(r'$\omega$')

fig.suptitle(r'$N_\mathrm{keep} = 2000$', y=0.95)
fig.tight_layout()

fig.savefig("sigma_im_compare.pdf")
# %%
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(w_mesh, Sigma_FG_w.data[:, 0, 0].real, "-", label=r'FG')
ax1.plot(w_mesh, Sigma_w.data[:, 0, 0].real, "-", label=r'IFG')
ax1.legend()

ax1.set_ylabel(r'$\mathrm{Re}\,\Sigma_{11}(\omega)$')

ax2.plot(w_mesh, Sigma_FG_w.data[:, 1, 1].real, "-", label=r'FG')
ax2.plot(w_mesh, Sigma_w.data[:, 1, 1].real, "-", label=r'IFG')
ax2.legend()

ax2.set_ylabel(r'$\mathrm{Re}\,\Sigma_{22}(\omega)$')
ax2.set_xlabel(r'$\omega$')

fig.suptitle(r'$N_\mathrm{keep} = 2000$', y=0.95)
fig.tight_layout()

fig.savefig("sigma_re_compare.pdf")
# %%
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

ax1.plot(w_mesh, -Sigma_FG_w.data[:, 0, 0].imag, "-", label=r'FG')
ax1.plot(w_mesh, -Sigma_w.data[:, 0, 0].imag, "-", label=r'IFG')
ax1.legend()

ax1.set_ylabel(r'$-\mathrm{Im}\,\Sigma_{11}(\omega)$')

ax2.plot(w_mesh, -Sigma_FG_w.data[:, 1, 1].imag, "-", label=r'FG')
ax2.plot(w_mesh, -Sigma_w.data[:, 1, 1].imag, "-", label=r'IFG')
ax2.legend()

ax2.set_ylabel(r'$-\mathrm{Im}\,\Sigma_{22}(\omega)$')
ax2.set_xlabel(r'$\omega$')

fig.suptitle(r'$N_\mathrm{keep} = 2000$', y=0.95)
fig.tight_layout()

ax2.set_xlim(-0.3, 0.3)
fig.savefig("sigma_im_compare_low_freq.pdf")
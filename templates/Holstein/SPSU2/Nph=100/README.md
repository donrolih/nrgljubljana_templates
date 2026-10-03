# Anderson-Holstein impurity, SPSU2 symmetry, Nph=100

A single impurity level coupled to one phonon mode and to a superconducting
bath. Only the spin SU(2) symmetry is used, so the bath may have a pairing
term (`scdelta`, `sckappa`).

    H_imp = delta n_d + U/2 (n_d-1)^2 + omega b^dag b + g (n_d-1)(b + b^dag)

The phonon cutoff is fixed when the template is generated (`Nph=100` in
`src/param`, 101 phonon states), hence the directory name. The three invariant
subspaces have dimensions 505, 404 and 101.

Based on the template in `/project/teodor/anderson-holstein/template`. The
Hamiltonian, the operators and the parameter names are the same; the numerical
matrices (`opf0_*`, `op.*`) are identical file by file.

## Generating

    sbatch generate_template_files.sh

from this directory. Mathematica needs a little under two hours, 75 minutes
of which go into the Hamiltonian.

## Using

Copy `data.in`, `ham_*`, `opf0_*`, `op.*_*_*` and `instantiate` to the run
directory, add the files listed below and run `./instantiate`. It writes
`data`; no Mathematica is involved. For Nph=100 this takes about 15 seconds.

`param` must contain

    [extra]
    U=
    delta=
    omega=
    g=

    [param]
    symtype=SPSU2
    data_has_rescaled_energies=false
    polarized=false
    Nmax=

followed by the usual `nrg` settings. `instantiate` stops if `Nmax` or
`polarized` is missing.

The Wilson chain goes in `xi1.dat`, `zeta1.dat`, `scdelta1.dat` and
`sckappa1.dat`, each with `Nmax+1` values. `instantiate` deletes these four
files when it is done.

**The hybridisation enters through `theta1.dat`, not through `[extra]`.** The
model uses the standard `Hc`, that is `gammaPolCh[1] hop[f[0], d[]]`, and the
`matrix` tool evaluates `gammaPolCh[1]` as `sqrt(theta/pi)` with `theta` read
from `theta1.dat`. This is the integral of the hybridisation function: `2
Gamma` for a flat band of half-width 1. A `Gamma=` line in `[extra]` is
ignored.

## What differs from Teodor's template

`model.m` has two more lines. `Conjugate[coefdelta[i__]] ^= coefdelta[i]`
removes `Conjugate[...]` from `ham_1` and `ham_2`, which `matrix` cannot
evaluate. `snegrealconstants[delta, U, omega, g]` declares the parameters;
without it, and with the `[extra]` entries commented out in `src/param`, the
phonon terms are left unevaluated.

`data_has_rescaled_energies=false` gives `SCALE 1` in `data.in`, so the
template does not depend on `Lambda`.

`data.in` lists all sixteen operators.

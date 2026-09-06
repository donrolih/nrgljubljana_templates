# SIAM, U(1) charge symmetry

Only charge is conserved. The impurity field may point in any direction in the
x-z plane (`Bz1`, `Bx1`), and with `pol2x2=true` the Wilson chain has a full
2x2 structure in the spin space: the up and down chains are coupled by cross
hopping coefficients.

What differs from the other SIAM templates:

**Four coefficient sets** instead of one or two: `1`=up-up, `2`=down-down,
`3`=up-down, `4`=down-up, for `coefxi` and `coefzeta`. So `xi{1..4}.dat` and
`zeta{1..4}.dat` are needed at run time, and `instantiate` calls `matrix -c4`.
Set `pol2x2=true`, not `polarized` — they are mutually exclusive.

`zeta3.dat` and `zeta4.dat` must be equal, otherwise the seed Hamiltonian is
not Hermitian.

**The hybridisation is a 2x2 matrix `coefV`, not `theta`.** `Hc` from `HC[]`
is built from `gammaPolCh[ch] = sqrt(theta_ch/Pi) >= 0` and so cannot carry a
negative or vanishing off-diagonal entry — which the channel-mixing
discretization routinely produces, since `V` is fixed only up to a unitary on
the bath index. `SIAM.m` therefore writes `H_hyb` out term by term using
`coefV[i,j]`, the symbol `matrix -s` fills in from `V{i}{j}{ch}.dat`
(`tools/matrix/parser.cc`, `load_discretization_sc`).

Two consequences. `coefV` is a `double`, so a genuinely complex hybridisation
is out of reach without an upstream change. And `coefV[i,j]` is used raw,
unlike `hybV[i,j] = Sqrt[1/Pi] V[i,j]` in `nrginit`, so these files hold the
physical amplitude directly — √π smaller than the `V{i}{j}.dat` of the
`band=manual_V` route.

`matrix -s` also insists on `V{i}{j}{ch}.dat` for every `ch=1..4` plus
`scdelta{ch}.dat` and `sckappa{ch}.dat`, even though `coefV` reads only
`ch=1` and nothing here references `coefdelta`/`coefkappa`. Supply the
`ch=2..4` copies and empty `scdelta`/`sckappa` files.

**`instantiate` is not the QSZ one.** Subspaces carry a single quantum number,
there are four coefficient sets, and U(1) emits *two* `f` blocks per channel
(`f 0 0`, `f 0 1`) — spin-up and spin-down matrix elements are not related by
the Wigner-Eckart theorem. The C++ `instantiate` tool cannot be used: it
rejects spin-polarized coefficient tables.

**Doublet operators are spin-resolved** for the same reason: `A_d_u`/`A_d_d`
and `self_d_u`/`self_d_d` replace `A_d`/`self_d`. Triplet operators
(`sigma_d`) are not available.

**G and Sigma are 2x2 matrices.** The off-diagonal components come from
cross-spectra between the operators already present, so no extra `ops=`
entries are needed:

    specd=A_d_u-A_d_u A_d_u-A_d_d A_d_d-A_d_u A_d_d-A_d_d
    specd=self_d_u-A_d_u self_d_u-A_d_d self_d_d-A_d_u self_d_d-A_d_d

The Hartree term is likewise a matrix, hence `SigmaHartree-{uu,dd,ud,du}`. The
off-diagonal part is proportional to `U <d^dag_up d_do>` and is non-zero
whenever `Bx1 != 0` or the bath mixes spin.

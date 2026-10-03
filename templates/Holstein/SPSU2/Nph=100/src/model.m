def1ch[1];

snegrealconstants[delta, U, omega, g];

(* Currently the parser does not support conjugated coefficients *)
Conjugate[coefdelta[i__]] ^= coefdelta[i];

nph = ToExpression @ optionvalue["Nph"];

Himp = delta number[d[]] + U/2 pow[number[d[]]-1, 2];
Himp = Himp + omega phononnumber[nph] + g nc[number[d[]]-1, phononplus[nph] + phononminus[nph]];
MAKEPHONON = 1; (* One phonon mode *)

PR[i_] := PR[i] = nc[ket[i], bra[i]];
PRlast = PR[nph];

(* phonon sign *)
nphsign = Sum[ (-1)^i PR[i], {i, 0, nph}];

(* sigma_x and sigma_z in the standard notation of quantum Rabi model, translated into ABS subspace *)
sigmaz = nc[d[CR,UP], d[CR,DO]] + nc[d[AN,DO], d[AN,UP]]; (* 2*Ix_d *)
sigmax = number[d[]] - 1;

(* parity operator *)
par = nc[ sigmaz, nphsign ];

(* Hc is the hybridization, gammaPolCh[1] hop[f[0], d[]]; H0 is the first site of Wilson chain *)
H = Himp + Hc + H0;

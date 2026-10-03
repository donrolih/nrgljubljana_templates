MyPrint["modelops begin"];
Module[{t},
  t = {};
  MPVCFAST = False;
  t = Join[t, mtSingletOp["PR0",     PR[0] ]];
  t = Join[t, mtSingletOp["PR1",     PR[1] ]];
  t = Join[t, mtSingletOp["PR2",     PR[2] ]];
  t = Join[t, mtSingletOp["PR3",     PR[3] ]];
  t = Join[t, mtSingletOp["PR4",     PR[4] ]];
  t = Join[t, mtSingletOp["PRlast",  PRlast ]];
  t = Join[t, mtSingletOp["nphsign", nphsign ]];
  t = Join[t, mtSingletOp["sigmax",  sigmax ]];
  t = Join[t, mtSingletOp["sigmaz",  sigmaz ]];
  t = Join[t, mtSingletOp["par",     par ]];
  MPVCFAST = True;
  texportable = t;
];
MyPrint["modelops end"];
texportable

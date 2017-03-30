(* definitions *)
ClearAll[
  SpaceDim,
  C0,
  CalcMSE,
  NormToNearest,
  STH,
  HTH,
  WaveletF2D,
  WDECoefficientF2D,
  WaveletSupport,
  WaveletZSupport1D,
  WaveletPartialSum,
  WaveletEstimator2D
  ];

(* dimension of space, c_0 *)
SpaceDim = 2;
C0 = Pi^(SpaceDim/2)/Gamma[SpaceDim/2 + 1];

(* distance to k-th neightbor in data *)

NormToNearest[data_, k_] :=
    With[
        {
            f = Nearest[DeleteDuplicates[data]]
        },
        Function[x, Norm[(x - Last[f[x, k + 1]])]]
    ];

(* thresholding *)
STH[th_, v_] := If[Abs[v] < th, 0., v - Sign[v] th];
HTH[cc_, n_, j0_, j_] := cc Sqrt[(j - j0)/n];

(* R^d wavelet function for j,q,z - tensor product *)

WaveletF2D[wave_, j_Integer, qs_, zs_] :=
    Module[
        {waveTs},
        waveTs = (If[# == 0, WaveletPhi, WaveletPsi] @ wave) & /@ qs;
        Function[
            {xs}, 
            With[
                {zip = {waveTs, zs, xs}, twoJ = 2.^ j },
                2.^(SpaceDim j/2) (Times @@ MapThread[#1[twoJ #3 - #2] &, zip])
            ]
        ]
    ];

(* coefficient calculator w/ memoization *)

WDECoefficientF2D[wave_, data_List, k_] :=
    Module[
        {
            nearF,
            numData,
            factor,
            f
        },
        nearF = NormToNearest[data, 1];
        numData = Length[data];
        factor = (Gamma[k]/Gamma[k + 1/2]) (1/Sqrt[numData]) Sqrt[C0];
        Print["Def WDECoeffs"];
        f[j_Integer, qs_, zs_] :=
            f[j, qs, zs] =
            With[
                {
                    wf = WaveletF2D[wave, j, qs, zs],
                    dHalf = (SpaceDim/2)
                },
                With[
                    {v = factor Sum[wf[p] nearF[p]^dHalf, {p, data}]}, 
                    v
                ]
            ];
        f
    ];

(* support in 1D for wavelet families; it uses filter coefficients *)

WaveletSupport[wave_, q_] :=
    Module[
        {
            kind,
            offset,
            coeffs,
            minX,
            maxX
        },
        kind = If[q == 0, "PrimalLowpass", "PrimalHighpass"];
        offset = If[q == 0, 0, 1];
        coeffs = WaveletFilterCoefficients[wave, kind];
        minX = Min[First[#] & /@ coeffs] + offset;
        maxX = Max[First[#] & /@ coeffs] + offset;
        {minX, maxX}
    ];

WaveletZSupport1D[j_, min_, max_, supMin_, supMax_] :=
    {Floor[(2^j) min - supMax], Ceiling[(2^j) max - supMin]} ;

WaveletPartialSum[wave_, coeffs_, mins_, maxs_, j_, qs_] :=
    Module[
        {
            supps,
            rangeZ,
            zsAll
        },
        supps = WaveletSupport[wave, #] & /@ qs;
        rangeZ = MapThread[
            WaveletZSupport1D[j, #1, #2, First[#3], Last[#3]] &,
            {mins, maxs, supps}
            ];
        zsAll = Tuples @ ((Range @@ # &) /@ rangeZ );
        (*Print[{"WaveletPartialSum",j,qs,supps,mins,maxs,{"rangeZ",rangeZ},{"zsAll",zsAll}}];*)
        Function[
            {xs},
            Sum[
                With[
                    {
                        c1 = coeffs[j, qs, zs],
                        c2 = WaveletF2D[wave, j, qs, zs][xs]
                    },
                    c1 c2
                ],
                {zs, zsAll}
            ]
        ]
    ];

(* estimator *)

WaveletEstimator2D[wave_, data_, k_, j0_, j1_] :=
    Module[
        {
            coeffsF,
            minX0,
            minXs,
            maxXs,
            q0s,
            qss,
            psum
        },
        coeffsF = WDECoefficientF2D[wave, data, k];
        minXs = Min /@ Transpose[data];
        maxXs = Max /@ Transpose[data];
        q0s = ConstantArray[0, SpaceDim];
        qss = Rest [Tuples @ ConstantArray[{0, 1}, SpaceDim]];
        psum[j_, qs_] := 
            psum[j, qs] = 
            (
                Print[{"psum", j, qs, ToString[TimeObject[]]}]; 
                WaveletPartialSum[wave, coeffsF, minXs, maxXs, j, qs]
            );
        Function[
            {xs},
            psum[j0, q0s][xs] +
            Sum[
                psum[j, qs][xs],
                {j, j0, j1},
                {qs, qss}
            ]
        ]
   ];

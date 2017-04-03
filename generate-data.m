SetOptions[$Output,FormatType->OutputForm];

<< wde.m;
<< distprocess.m;
<< mydist_bumps.m;
<< utils.m;

If[Length[$ScriptCommandLine] < 4, Print["Usage: generate.m <wave> <j0> <J1>"]; Quit[]];
{WAVE, J0, J1} = ToExpression /@ $ScriptCommandLine[[-3 ;; -1]];
If[Head[WaveletPhi[WAVE]] == WaveletPhi, Quit[]];


RunOnce[dname_, i_, j0_, j1_, n_] :=
    Module[
        {
            data,
            estimator,
            truePdf,
            table,
            rval,
            ise
        },
        data = DistData[MyDist, n];
        estimator = WaveletEstimator2D[WAVE, data, 1, j0, j1];
        truePdf = DistPDF[MyDist];
        table = Flatten[
            Table[
                With[
                    {v = estimator[{x0,x1}]^2},
                    {N[x0], N[x1], v, (v - truePdf[{x0,x1}])^2 }
                ],
                {x0,0,1,1/32},
                {x1,0,1,1/32}
            ],
            1
        ];
        rval = IntegerString[IntegerPart[RandomVariate[UniformDistribution[{0,1000000}]]],10,7];
        Export[FileNameJoin[{dname, "/data-"<>IntegerString[i, 10, 4]<>"-"<>rval<>".csv"}], table];
        Export[FileNameJoin[{dname, "/table-"<>IntegerString[i, 10, 4]<>"-"<>rval<>".csv"}], table];
        ise = N[Total[table[[All,4]]]/1024];
        SaveISE[j0, j1, n, i, ise, rval];
    ];

DoLevel[j0_,j1_,n_] := Module[
    {
        dname
    },
    dname = MkResultDir[WAVE, j0, j1, n];
    ParallelDo[
        RunOnce[dname, i, j0, j1, n],
        {i,1,100}
    ]
];

(* main *)
LaunchKernels[4];
Do[
    Do[
        DoLevel[j0, j1, n],
        {j1, j0-1, Min[j0+2,J1]}
    ],
    {j0, J0, J1},
    {n, {50, 75, 100, 150, 200, 225, 300, 400}}
]

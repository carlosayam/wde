SetOptions[$Output,FormatType->OutputForm];

<< wde.m;
<< distprocess.m;
<< mydist_bumps.m;
<< utils.m;

If[Length[$ScriptCommandLine] < 3, Print["Usage: generate.m <wave> <j0> <J1> <sample_size>"]; Quit[]];
{WAVE, J0, J1, SAMPLESIZE} = ToExpression /@ $ScriptCommandLine[[-4 ;; -1]];


RunOnce[dname_, i_] :=
    Module[
        {
            data,
            estimator,
            truePdf,
            table,
            ise
        }
        data = DistData[MyDist, SAMPLESIZE];
        estimator = WaveletEstimator2D[WAVE, data, 1, J0, J1];
        truePdf = DistPDF[MyDist];
        table = Flatten[
            Table[
                With[
                    {v = gHat[{x0,x1}]^2},
                    {x0, x1 , v, (v - truePdf[{x0,x1}])^2 }
                ],
                {x0,0,1,1/32},
                {x1,0,1,1/32}
            ],
            1
        ];
        Export[FileNameJoin[{dname, "/sample-"<>IntegerString[i, 10, 4]<>".csv"}], table];
        ise = Total[table[[All,4]]]/1024;
        SaveISE[i, ise];
    ];

Main[] := Module[
    {
        dname
    },
    dname = MkResultDir[WAVE, J0, J1, SAMPLESIZE];
    Do[
        RunOnce[dname, i],
        {i,1,500}
    ]
];

Main[];

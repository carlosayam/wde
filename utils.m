(* utils for main .m scripts *)

MkDirP[{}] := "";
MkDirP[{x_, xs_}] := f[xs_] := StringJoin[ Table[
    (
     Print[x];
     x <> "/"
     ), {x, xs}
    ]
   ];
    (
        If[!DirectoryQ[x],CreateDirectory[x]];
        x <> "/" <> MkDirP[xs];
    );

WaveFName[wave_] := StringReplace[ToString[wave], {"[" -> "_", "]" -> ""}];

MkResultDir[wave_, j0_, j1_, sampleSize_] :=
    With[
        {
            dname = FileNameJoin[
                {
                    "data",
                    WaveFName[wave],
                    "opts-j0_" <> IntegerString[j0, 10, 3] <> "j1_" <> IntegerString[j1, 10, 3] <> "n_" <> IntegerString[sampleSize, 10, 5]
                }
            ]
        },
        CreateDirectory[dname, CreateIntermediateDirectories -> True]
    ];

PrintISE[i_,ise_] :=
    Module[
        {
            file
        }
        file = OpenAppend["data/ise.csv"];
        Export[file, {WaveFName[WAVE],J0,J1,SAMPLESIZE,i,ise}, "CSV"];
        WriteString[file, "\n"];
        Close[file];
    ];
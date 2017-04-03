(* utils for main .m scripts *)

WaveFName[wave_] := StringReplace[ToString[wave], {"[" -> "_", "]" -> ""}];

MkResultDir[wave_, j0_, j1_, sampleSize_] :=
    With[
        {
            dname = FileNameJoin[
                {
                    "data",
                    WaveFName[wave],
                    "opts-j0_" <> IntegerString[j0, 10, 3] <> ",j1_" <> IntegerString[j1, 10, 3] <> ",n_" <> IntegerString[sampleSize, 10, 5]
                }
            ]
        },
        Quiet[CreateDirectory[dname, CreateIntermediateDirectories -> True],{CreateDirectory::filex}]
    ];

SaveISE[i_,ise_, rval_] :=
    Module[
        {
            string
        },
        string = "\"" <> WaveFName[WAVE] <> "\","
          <> ToString[J0] <> ","
          <> ToString[J1] <> ","
          <> ToString[SAMPLESIZE] <> ","
          <> "\"" <> ToString[i] <> "-" <> rval <> "\","
          <> ToString[ise] <> "\n";

        s = OpenAppend["data/ise.csv"];
        WriteString[s, string];
        Close[s];
    ];
(* utils for main .m scripts *)

WaveFName[wave_] := StringReplace[ToString[wave], {"[" -> "_", "]" -> ""}];

MkResultDir[wave_, j0_, j1_, sampleSize_] :=
    With[
        {
            dname = FileNameJoin[
                {
                    "data",
                    WaveFName[wave],
                    "opts-j0_" <> IntegerString[j0, 10, 3]
                        <> ",j1_" <> If[j1 >= j0, IntegerString[j1, 10, 3], "___"]
                        <> ",n_"
                        <> IntegerString[sampleSize, 10, 5]
                }
            ]
        },
        Quiet[CreateDirectory[dname, CreateIntermediateDirectories -> True],{CreateDirectory::filex}]
    ];

SaveISE[j0_, j1_, n_, i_, ise_, rval_] :=
    Module[
        {
            string
        },
        string = "\"" <> WaveFName[WAVE] <> "\","
          <> ToString[j0] <> ","
          <> If[j1 >= j0, ToString[j1],""] <> ","
          <> ToString[n] <> ","
          <> "\"" <> ToString[i] <> "-" <> rval <> "\","
          <> ToString[ise] <> "\n";

        s = OpenAppend["data/ise.csv"];
        WriteString[s, string];
        Close[s];
    ];
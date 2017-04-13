(* utils for main .m scripts *)

MyIntegerString[i_, ds_] := If[i < 0, "-", "0"] <> IntegerString[i, 10, ds];

WaveFName[wave_] := StringReplace[ToString[wave], {"[" -> "_", "]" -> ""}];

MkResultDir[wave_, j0_, j1_, sampleSize_] :=
    With[
        {
            dname = FileNameJoin[
                {
                    "data",
                    WaveFName[wave],
                    "opts-j0_" <> MyIntegerString[j0, 3]
                        <> ",j1_" <> If[j1 >= j0, MyIntegerString[j1, 3], "____"]
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
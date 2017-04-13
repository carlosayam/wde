(* ::Package:: *)

(* dist processing *)

ClearAll[DistPDF, DistData];

DistPDF[dist_] :=
    Module[
        {
            total
        },
        total = Total[Flatten[Table[PDF[dist,{x0,x1}],{x0,0,1,1/32},{x1,0,1,1/32}]]]/1024.0;
        Function[{p}, PDF[dist,p]/total]
    ];

DistData[dist_,n_] :=
    Module[
        {
            makeC,
            dataC
        },
        makeC[d_] := 
            Module[
                {r = {-1,-1}},
                (
                    While[!RegionMember[Rectangle[{0,0}]][r], r = RandomVariate[d]];
                    r
                ) &
            ];
        dataC := makeC[dist];
        Table[dataC[],{ni,1,n}]
    ];

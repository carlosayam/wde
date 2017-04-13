(* dist easy *)

ClearAll[MyDist];

MyDist :=
    ProductDistribution[
        TransformedDistribution[(x - 0.125) / 0.5, x \[Distributed] BetaDistribution[2,4]],
        TransformedDistribution[(x - 0.125) / 0.5, x \[Distributed] BetaDistribution[2,4]]
    ];

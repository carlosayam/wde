(* ::Package:: *)

SetOptions[$Output,FormatType->OutputForm];

ClearAll[Main];

(* 2-d *)
Main[] := Module[
   {
    \[Rho],
    \[Sigma],
    dist,
    makeC,
    dataC,
    data
    },
   \[Rho] = 0.5; \[Sigma] = 0.05;
   
   dist = MixtureDistribution[{1, 
       8}, {MultinormalDistribution[{0.2, 
         0.3}, {{\[Sigma]/6, 0}, {0, \[Sigma]/6}}], 
       MultinormalDistribution[{0.7, 
         0.7}, {{0.1, \[Sigma]/8}, {\[Sigma]/8, 0.1}}]}];
   data = Flatten[Table[{x,y,PDF[dist,{x,y}]},{x,0,1,0.1},{y,0,1,0.1}],1];
   ListPlot3D[data]
];

Main[];

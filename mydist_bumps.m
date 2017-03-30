(* dist difficult *)
    
MyDist :=
    With[
        {
            \[Rho],
            \[Sigma],
        },
        \[Rho] = 0.5;
        \[Sigma] = 0.05;
        dist = MixtureDistribution[
            {1,8}, 
            {
                MultinormalDistribution[{0.2, 0.3}, {{\[Sigma]/6, 0}, {0, \[Sigma]/6}}], 
                MultinormalDistribution[{0.7, 0.7}, {{0.1, \[Sigma]/8}, {\[Sigma]/8, 0.1}}]
            }
        ]
    ]

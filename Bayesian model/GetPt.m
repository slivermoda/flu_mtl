function res = GetPt( z, t, T0, T1, T2, Alpha, Sigma, Sigma_far )
res = z * normpdf(t, T0, Sigma) + ( 1 - z ) * ...
            ( Alpha * ( t >= T1 ) * normpdf(t, T1, Sigma_far) + ( 1-Alpha ) * ( t <= T2 ) * normpdf(t, T2, Sigma_far) );
end
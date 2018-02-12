function res = GetPd( z, d, D1, D2, Delta )
res = z * normpdf(d, D1, Delta) + ( 1 - z ) * normpdf(d, D2, Delta);
end
clear; clc;
load Data/SpatialUse_2009_2016.mat
%% Get samples within one drift
% cluster order: 1994-2009 [1 5 6 2 4], 2009-2016 [2 5 1 3]
new_id = 1;
c_old = mean( X0(ClusterIDX == 5, :) );
c_new = mean( X0(ClusterIDX == 1, :) );
%%
BoundaryNode = 6091;
for i = 1:length(VNode)
    if VNode(i).dn >= VNode(BoundaryNode).dn - 30*6 % half year
        LeftPoint = i;
        break;
    end
end
for i = 1:length(VNode)
    if VNode(i).dn >= VNode(BoundaryNode).dn + 30*6 % half year
        RightPoint = i;
        break;
    end
end
Selection = (LeftPoint:RightPoint);
%%
V = VNode(Selection);
coor = X0(Selection, :);
n = length(V);
t = zeros(n, 1);
d = zeros(n, 1);
for i = 1:n
    t(i) = V(i).dn;
    d(i) = norm(coor(i, :) - c_new, 2) - norm(coor(i, :) - c_old, 2);
end
N = 1e4;
z = zeros(n, N);
posterior = zeros(n, N);
%% Sigma (20) and Sigma_far (50) are defined for half year
T0 = VNode(BoundaryNode).dn;
T1 = V(1).dn;
T2 = V(end).dn;
D1 = -norm(c_old-c_new, 2);
D2 = norm(c_old-c_new, 2);
Sigma = 20;
Sigma_far = 50;
Delta = 1;
p0 = 0.5;
%%
k = 0;
while true
    for i = 1:n
        %% z^(k+1)
        a = rand;
        if a <= p0
            z(i, k+1) = 0;
        else
            z(i, k+1) = 1;
        end
        %% Alpha
        if t(i) <= T0
            Alpha = 1;
        else
            Alpha = 0;
        end
        %%
        Pt_new = GetPt( z(i, k+1), t(i) - T0, 0, T1 - T0, T2 - T0, Alpha, Sigma, Sigma_far );
        Pd_new = GetPd( z(i, k+1), d(i), D1, D2, Delta );
        posterior(i, k+1) = Pt_new*Pd_new;
        
        z_opp = 1 - z(i, k+1);
        Pt_opp = GetPt( z_opp, t(i) - T0, 0, T1 - T0, T2 - T0, Alpha, Sigma, Sigma_far );
        Pd_opp = GetPd( z_opp, d(i), D1, D2, Delta );
        Rho = min( (Pt_new*Pd_new) / (Pt_opp*Pd_opp), 1 );
        a = rand;
        if ( a > Rho )
            z(i, k+1) = z_opp;
            posterior(i, k+1) = Pt_opp*Pd_opp;
        end
        if V(i).ClusterID ~= new_id
            z(i, k+1) = 0;
        end
    end
    k = k + 1;
    if mod(k, 1000) == 0
        fprintf('iter = %d\n', k);
    end
    if k > N
        break;
    end
end
%%
z_star = mean( z, 2 );
V_star = V(z_star > 0);
posterior_star = posterior(z_star > 0);
z_star = z_star(z_star > 0);
for i = 1:length(V_star)
    V_star(i).z = z_star(i);
    V_star(i).posterior = posterior_star(i);
end
%%
continents = {'Africa', 'Asia', 'Australia', 'Europe', 'NorthAmerica', 'SouthAmerica'};
p = containers.Map();
p_sum = 0;
for i = 1:length(continents)
    p(continents{i}) = 0;
end
for i = 1:length(V_star)
    if isempty( V_star(i).Continent )
        continue;
    end
    p(V_star(i).Continent) = p(V_star(i).Continent) + V_star(i).posterior;
    p_sum = p_sum + V_star(i).posterior;
end
for i = 1:length(continents)
    p(continents{i}) = p(continents{i}) / p_sum;
    fprintf('%s: %.4f\n', continents{i}, p(continents{i}) );
end
%% time
p_time = containers.Map();
for i = 1:length(V_star)
    y = V_star(i).Y;
    m = V_star(i).M;
    key = strcat( num2str(y), '-', num2str(m) );
    if ~p_time.isKey(key)
        p_time( key ) = 0;
    else
        p_time( key ) = p_time( key ) + V_star(i).posterior;
    end
end
keyList = p_time.keys;
p_time_sum = 0;
for i = 1:length( keyList )
    p_time_sum = p_time_sum + p_time( keyList{i} );
end
for i = 1:length( keyList )
    p_time( keyList{i} ) = p_time( keyList{i} ) / p_time_sum;
end
for i = 1:length( keyList )
    fprintf('%s: %.4f\n', keyList{i}, p_time( keyList{i} ) );
end








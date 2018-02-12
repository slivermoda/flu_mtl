function [X, pos] = Remove_Rare_Feature(X, Rate)
X_bool = (X > 0);
n = size(X, 1);
vec = sum(X_bool);
vec = vec / n;
pos = find( vec < Rate );
X(:, pos) = [];
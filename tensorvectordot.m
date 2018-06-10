function [M] = tensorvectordot(T,v)
%
% serves for weighted matrix summation
% 
% Input
%	T: x by y by z tensor. Each frontal slice of T is a x by y matrix.
%	v: z by 1 vector. Each element of v(k) is the weight of the corresponding frontal slice matrix of T(:,:,k).
% Output
%	M: x by y matrix, the weighted summation of all frontal slices of T, each weighted by an element of v.

[x,y,z] = size(T);

% [z,1] = size(v);

M = zeros(x,y);

for i = 1:x
for j = 1:y
temp = 0;
for k = 1:z
temp = temp + T(i,j,k)*v(k,1);
end
M(i,j) = temp;
end
end

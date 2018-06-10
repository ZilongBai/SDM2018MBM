function [M] = tensorvectordot(T,v)
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

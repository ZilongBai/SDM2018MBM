function [F,M,alpha] = InitFMAlpha(N,k,p,r)
%
% randomly initializes F, M, and alpha given N, k, p, and r
%
% Input
%       N: the number of vertices of each single graph
%       k: the number of blocks in each block structure
%       p: the number of different block structures to find
%       r: the number of single graphs in the multigraph X
% Output
%	F: randomly initialized N x k x p nonnegative tensor. 
%%	With column vectors of each frontal slice matrix normalized.
% 	M: randomly initialized k x k x p x r nonnegative tensor. 
%%	Conditioned to be positive semi-definitive. And then absorbed original norms of column vectors of each frontal slice matrix of F.
%	alpha: randomly generated p x r nonnegative matrix with normalized column vectors. NOTE: randomly initialized alpha does not participate in the later part of the algorithm implementation.

F = abs(rand(N,k,p));
M = abs(rand(k,k,p,r));
alpha = abs(rand(p,r));

for j = 1:r
for i = 1:p
temp = M(:,:,i,j)*M(:,:,i,j)';
temp = temp./(max(max(temp)));
M(:,:,i,j) = temp;
end
sumalphaj = sum(alpha(:,j));
alpha(:,j) = alpha(:,j)./sumalphaj;
end

for i = 1:p
	for b = 1:k
		F(:,b,i) = F(:,b,i)./max(F(:,b,i));
	end
end

%------ Begin Rescaling without changing obj
for j = 1:r

sumalphaj = sum(alpha(:,j));
alpha(:,j) = alpha(:,j)./sumalphaj;

        for a = 1:p

                M(:,:,a,j) = M(:,:,a,j).*(sumalphaj);

        end

end

normF = zeros(k,p);
for a = 1:p
        for b = 1:k
                tmpnormF(b,a) = norm(F(:,b,a));
                F(:,b,a) = F(:,b,a)./norm(F(:,b,a));
        end
end

for j = 1:r
        for a = 1:p
                for c = 1:k
                        for d = 1:k
                                M(c,d,a,j) = M(c,d,a,j)*tmpnormF(c,a)*tmpnormF(d,a);
                        end
                end
        end
end


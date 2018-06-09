function [alpha,F,M] = AlphaExtractor_normM(F,M,N,k,p,r)
%
% extracts alphas from weighted mixing matrices M (Equation 4.17 of the paper Mixtures of Block Models for Brain Networks, SDM 2018)
% (By Zilong Bai, KDD Lab @ University of California, Davis)
% 
% Input
% 	F: N x k x p tensor of set of block structures learnt with the alternative solving framework
% 	M: k x k x p x r tensor of sets of weighted mixing matrices learnt with the alternative solving framework
%	N: the number of vertices of each single graph
%	k: the number of blocks in each block structure
%	p: the number of different block structures to find
%	r: the number of single graphs in the multigraph X
%
% Output
%	alpha: p x r matrix. The set of block model weights. alpha(:,j) are the block model weights for modeling graph X(:,:,j). alpha (i,j) is corresponded to block structure F(:,:,i) and mixing matrix M(:,:,i,j).
%	F: N x k x p tensor. With normalized block indicator vectors, i.e., column vectors in each frontal slice of F.
%	M: k x k x p x r tensor. Sets of rescaled mixing matrices.

%% Firstly rescale F and M without changing optimality to normalize columns of F.
Fnorm = zeros(k,p);

for i = 1:p

	for a = 1:k

		Fnorm(a,i) = norm(F(:,a,i));
	end

end

for i = 1:p

	for a = 1:k

		F(:,a,i) = F(:,a,i)./Fnorm(a,i);

	end

end

for j = 1:r

	for i = 1:p

		for a = 1:k

			for b = 1:k

				M(a,b,i,j) = M(a,b,i,j)*Fnorm(a,i)*Fnorm(b,i);

			end

		end

	end

end

%% Extracting alpha without changing objetive function value. alpha(:,j) is Frobenius norm 1.

alpha = zeros(p,r);

for j = 1:r
	
	for i = 1:p
	
		alpha(i,j) = norm(M(:,:,i,j),'fro');
	
	end

end

for j = 1:r

	for i = 1:p

		M(:,:,i,j) = M(:,:,i,j)./alpha(i,j);

	end

end

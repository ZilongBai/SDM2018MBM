function [F,R] = updateF_continuity(X,Xbar,M,F, N,k,p,r,pc,theta)
% 
% serves to use multiplicative update rules to update nonnegative tensor F (Equation 4.15 in paper Mixtures of Block Models for Brain Networks, SDM 2018)
% (By Zilong Bai, KDD Lab @ University of California, Davis)
% Input:
%	F: N x k x p tensor. N is the number of vertices from the original graphs, k is the number of blocks in each F(:,:,i), p is the number of latent block models.
%	M: k x k x p x r tensor. r is the number of graphs in the multigraph X.
% 	alpha: p x r matrix. alpha(:,j) is the weights for the block models in the linear combination that approximates certain graph from the stack.
%	N: the number of vertices of each single graph
%	k: the number of blocks in each block structure
%       p: the number of different block structures to find
%       r: the number of single graphs in the multigraph X
%	pc: nonnegative scalar, Weight of the continuity regularization term, i.e., beta in the paper.
%       Theta: N x N spatial closeness penalty matrix. In particular for the applciation of our method to networks constructed with brain imaging data in this paper, the penalty matrix Theta is calculated according to the spatial coordinates of the centroids of the anatomical regions in the brain. 
%%	The Spatial Continuity Regularization (Bai et al. KDD 2017) is introduced to faciliate interpretable discovery of nodes.
% Output
%	F: N x k x p tensor. Updated F.
% 	R: N x N x p x r tensor. Intermediate variables that record the residual between given graphs and the reconstructed graphs based on current F and M.

epsilon = 1e-6; % ignorable positive value to avoid issues in multiplicative update rules caused by intermediate zero elements.
R = zeros(N,N,p,r);

for j = 1:r
	for i = 1:p
		R(:,:,i,j) = X(:,:,j);
		for a = 1:p
			if a~= i
				R(:,:,i,j) = R(:,:,i,j) - F(:,:,i)*M(:,:,i,j)*F(:,:,i)';
			end
		end
	end
end


for i = 1:p
tmp1 = zeros(N,k,r);
tmp2 = tmp1;
for j = 1:r
	tmp1(:,:,j) = R(:,:,i,j)*F(:,:,i)*M(:,:,i,j);
	tmp2(:,:,j) = F(:,:,i)*F(:,:,i)'*R(:,:,i,j)*F(:,:,i)*M(:,:,i,j);
end
numerator = tensorvectordot(tmp1,ones(r,1));
denominator = tensorvectordot(tmp2,ones(r,1));

denominator = denominator + theta*F(:,:,i).*pc - F(:,:,i)*F(:,:,i)'*theta*F(:,:,i).*pc;

denominator = denominator + ones(size(denominator)).*epsilon;

SQ = numerator./denominator;

SQ(find(SQ<=0)) = epsilon;
F(:,:,i) = F(:,:,i).*sqrt(SQ);

end
